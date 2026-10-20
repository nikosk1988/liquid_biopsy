#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <htslib/sam.h>
#include <sqlite3.h>

const char* nt_map = "=ACMGRSVTWYHKDBN";

int main(int argc, char *argv[]) {
    if (argc < 2) return 1;
    int n_threads = (argc > 2) ? atoi(argv[2]) : 4;
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    samFile *fp = sam_open(argv[1], "r");
    bam_hdr_t *header = sam_hdr_read(fp);
    if (n_threads > 1) hts_set_threads(fp, n_threads);
    bam1_t *aln = bam_init1();

    sqlite3 *db;
    sqlite3_open("liquid_biopsy.db", &db);

    // ΑΝΑΓΚΑΣΤΙΚΟ REFRESH SCHEMA (ΠΡΟΣΘΗΚΗ ΑΥΤΩΝ ΤΩΝ ΓΡΑΜΜΩΝ)
    sqlite3_exec(db, "DROP TABLE IF EXISTS Mutations;", 0,0,0);
    sqlite3_exec(db, "CREATE TABLE Mutations (WT_Count INT, Mut_Count INT, Avg_Qual REAL, Avg_GC REAL);", 0,0,0);
    sqlite3_exec(db, "CREATE TABLE IF NOT EXISTS FragmentSizes (Size INT, Count INT);", 0,0,0);
    sqlite3_exec(db, "CREATE TABLE IF NOT EXISTS Motifs (Motif TEXT PRIMARY KEY, Count INT);", 0,0,0);
    sqlite3_exec(db, "CREATE TABLE IF NOT EXISTS ProcessingMetrics (TotalReads INT, Runtime REAL, Efficiency REAL, Threads INT);", 0,0,0);

    sqlite3_exec(db, "DELETE FROM FragmentSizes; DELETE FROM Motifs; DELETE FROM ProcessingMetrics;", 0,0,0);
    sqlite3_exec(db, "BEGIN TRANSACTION;", 0,0,0);

    long size_dist[1000] = {0}, total_reads = 0, wt = 0, mut = 0, qual_sum_cnt = 0;
    double sum_qual = 0, sum_gc = 0;

    while (sam_read1(fp, header, aln) >= 0) {
        total_reads++;
        uint8_t *seq = bam_get_seq(aln);
        int len = aln->core.l_qseq;

        // GC-Content
        int gc_chars = 0;
        for (int i = 0; i < len; i++) {
            char b = nt_map[bam_seqi(seq, i)];
            if (b == 'G' || b == 'C') gc_chars++;
        }
        sum_gc += ((double)gc_chars / len) * 100.0;

        // Targeted SNV (Position 10) - QC Filtering (Q30)
        uint8_t *q_score = bam_get_qual(aln);
        if (len >= 15 && aln->core.qual >= 30) {
            if (q_score[10] >= 30) { // Μόνο αν η βάση είναι υψηλής ποιότητας
                char b = nt_map[bam_seqi(seq, 10)];
                if (b == 'T') mut++;
                else if (b == 'G') wt++;
                sum_qual += q_score[10]; qual_sum_cnt++;
            }
        }

        // Fragments & Motifs
        int32_t isize = abs(aln->core.isize);
        if (isize > 0 && isize < 1000) size_dist[isize]++;
        if (len >= 4) {
            char m[5] = {0}; for(int j=0; j<4; j++) m[j]=nt_map[bam_seqi(seq, j)];
            char sql_m[128]; sprintf(sql_m, "INSERT INTO Motifs VALUES ('%s', 1) ON CONFLICT(Motif) DO UPDATE SET Count=Count+1;", m);
            sqlite3_exec(db, sql_m, 0,0,0);
        }
    }

    for (int i=0; i<1000; i++) if (size_dist[i]>0) {
        char sql_s[128]; sprintf(sql_s, "INSERT INTO FragmentSizes VALUES (%d, %ld);", i, size_dist[i]);
        sqlite3_exec(db, sql_s, 0,0,0);
    }

    char sql_mut[256];
    sprintf(sql_mut, "INSERT INTO Mutations VALUES (%ld, %ld, %.2f, %.2f);", wt, mut, (qual_sum_cnt?sum_qual/qual_sum_cnt:0), (total_reads?sum_gc/total_reads:0));
    sqlite3_exec(db, sql_mut, 0,0,0);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double runtime = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
    char sql_p[256]; sprintf(sql_p, "INSERT INTO ProcessingMetrics VALUES (%ld, %f, %f, %d);", total_reads, runtime, (total_reads?runtime*1000.0/total_reads:0), n_threads);
    sqlite3_exec(db, sql_p, 0,0,0);

    sqlite3_exec(db, "COMMIT;", 0,0,0);
    sqlite3_close(db); bam_destroy1(aln); bam_hdr_destroy(header); sam_close(fp);
    return 0;
}
