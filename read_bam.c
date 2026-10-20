/*
 * read_bam.c  (v2)
 * ctDNA-Discoverer Pro - C-Engine
 *
 * Χρήση:  ./read_bam <αρχείο.bam> [threads]
 *
 * Αλλαγές σε σχέση με την v1:
 *  1. Φιλτράρισμα με βάση τα BAM flags: απορρίπτονται unmapped, secondary,
 *     supplementary, QC-fail και PCR/optical duplicates (πίνακας QCMetrics).
 *  2. Κάθε θραύσμα μετριέται ΜΙΑ φορά στην κατανομή μεγεθών
 *     (μόνο η ανάγνωση με θετικό TLEN ενός proper pair).
 *  3. Τα 4-mer end-motifs μετρώνται σε πίνακα μνήμης [256] και γράφονται
 *     στη βάση μία φορά στο τέλος (πριν: 1 εντολή SQL ανά ανάγνωση).
 *  4. Το end-motif λαμβάνεται από το 5' άκρο του θραύσματος σε κάθε κλώνο
 *     (reverse complement για αναγνώσεις στον αρνητικό κλώνο).
 *  5. Έλεγχοι σφαλμάτων (αρχείο, header, βάση) και προστασία από len = 0.
 *
 * Ο πίνακας Mutations εξακολουθεί να μετράει τη βάση στη θέση 10 ΤΗΣ ΑΝΑΓΝΩΣΗΣ
 * (όχι σε γονιδιωματική συντεταγμένη). Αντικαθίσταται στο Βήμα 5, μαζί με τη
 * νέα γεννήτρια δεδομένων.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <htslib/sam.h>
#include <sqlite3.h>

#define MAX_ISIZE 1000
#define TARGET_READ_POS 10
#define MIN_MAPQ 30
#define MIN_BASEQ 30

static const char *nt_map = "=ACMGRSVTWYHKDBN";

/* Κωδικοποίηση 2-bit: A=0 C=1 G=2 T=3, -1 για οτιδήποτε άλλο (π.χ. N) */
static int base2bit(char b) {
    switch (b) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; }
    return -1;
}

/* Δείκτης 0..255 του 5' end-motif του θραύσματος, ή -1 αν περιέχει μη-ACGT βάση */
static int end_motif_index(const bam1_t *aln) {
    const uint8_t *seq = bam_get_seq(aln);
    int len = aln->core.l_qseq, idx = 0;
    if (len < 4) return -1;
    int reverse = (aln->core.flag & BAM_FREVERSE) != 0;
    for (int j = 0; j < 4; j++) {
        int code;
        if (!reverse) {
            code = base2bit(nt_map[bam_seqi(seq, j)]);
        } else {
            /* 5' άκρο στον αρνητικό κλώνο = τελευταίες βάσεις, συμπληρωματικές */
            code = base2bit(nt_map[bam_seqi(seq, len - 1 - j)]);
            if (code >= 0) code = 3 - code;   /* A<->T, C<->G */
        }
        if (code < 0) return -1;
        idx = idx * 4 + code;
    }
    return idx;
}

static void motif_string(int idx, char out[5]) {
    const char *b = "ACGT";
    for (int j = 3; j >= 0; j--) { out[j] = b[idx % 4]; idx /= 4; }
    out[4] = '\0';
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <file.bam> [threads]\n", argv[0]);
        return 1;
    }
    int n_threads = (argc > 2) ? atoi(argv[2]) : 4;
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    samFile *fp = sam_open(argv[1], "r");
    if (!fp) { fprintf(stderr, "Cannot open %s\n", argv[1]); return 1; }
    /* Τα threads χρησιμοποιούνται από το htslib για παράλληλη αποσυμπίεση BGZF */
    if (n_threads > 1) hts_set_threads(fp, n_threads);
    sam_hdr_t *header = sam_hdr_read(fp);
    if (!header) { fprintf(stderr, "Cannot read BAM header\n"); sam_close(fp); return 1; }
    bam1_t *aln = bam_init1();

    long size_dist[MAX_ISIZE] = {0}, motif_counts[256] = {0};
    long total_reads = 0, passed_reads = 0, fragments = 0;
    long f_unmapped = 0, f_secondary = 0, f_supplementary = 0, f_qcfail = 0, f_dup = 0;
    long wt = 0, mut = 0, qual_sum_cnt = 0, gc_reads = 0;
    double sum_qual = 0, sum_gc = 0;

    int ret;
    while ((ret = sam_read1(fp, header, aln)) >= 0) {
        total_reads++;
        uint16_t flag = aln->core.flag;

        /* ---- 1. QC με βάση τα flags ---- */
        if (flag & BAM_FUNMAP)         { f_unmapped++;      continue; }
        if (flag & BAM_FSECONDARY)     { f_secondary++;     continue; }
        if (flag & BAM_FSUPPLEMENTARY) { f_supplementary++; continue; }
        if (flag & BAM_FQCFAIL)        { f_qcfail++;        continue; }
        if (flag & BAM_FDUP)           { f_dup++;           continue; }
        passed_reads++;

        const uint8_t *seq = bam_get_seq(aln);
        const uint8_t *q_score = bam_get_qual(aln);
        int len = aln->core.l_qseq;

        /* ---- 2. GC-content ανά ανάγνωση ---- */
        if (len > 0) {
            int gc = 0;
            for (int i = 0; i < len; i++) {
                char b = nt_map[bam_seqi(seq, i)];
                if (b == 'G' || b == 'C') gc++;
            }
            sum_gc += 100.0 * gc / len;
            gc_reads++;
        }

        /* ---- 3. SNV στη θέση 10 της ανάγνωσης (MAPQ >= 30, BaseQ >= 30) ---- */
        if (len > TARGET_READ_POS && aln->core.qual >= MIN_MAPQ
            && q_score[0] != 0xff && q_score[TARGET_READ_POS] >= MIN_BASEQ) {
            char b = nt_map[bam_seqi(seq, TARGET_READ_POS)];
            if (b == 'T') mut++;
            else if (b == 'G') wt++;
            sum_qual += q_score[TARGET_READ_POS];
            qual_sum_cnt++;
        }

        /* ---- 4. Μέγεθος θραύσματος: μία φορά ανά proper pair (TLEN > 0) ---- */
        if ((flag & BAM_FPAIRED) && (flag & BAM_FPROPER_PAIR)) {
            int64_t isize = aln->core.isize;
            if (isize > 0) {
                fragments++;
                if (isize < MAX_ISIZE) size_dist[isize]++;
            }
        }

        /* ---- 5. End-motif (5' άκρο, ανά ανάγνωση = ανά άκρο θραύσματος) ---- */
        int mi = end_motif_index(aln);
        if (mi >= 0) motif_counts[mi]++;
    }
    if (ret < -1) fprintf(stderr, "Warning: BAM read error (%d); results may be truncated\n", ret);

    /* ---- Αποθήκευση στη SQLite (μία συναλλαγή, prepared statements) ---- */
    sqlite3 *db;
    if (sqlite3_open("liquid_biopsy.db", &db) != SQLITE_OK) {
        fprintf(stderr, "Cannot open database: %s\n", sqlite3_errmsg(db));
        return 1;
    }
    sqlite3_exec(db,
        "DROP TABLE IF EXISTS Mutations; DROP TABLE IF EXISTS FragmentSizes;"
        "DROP TABLE IF EXISTS Motifs; DROP TABLE IF EXISTS ProcessingMetrics;"
        "DROP TABLE IF EXISTS QCMetrics;"
        "CREATE TABLE Mutations (WT_Count INT, Mut_Count INT, Avg_Qual REAL, Avg_GC REAL);"
        "CREATE TABLE FragmentSizes (Size INT, Count INT);"
        "CREATE TABLE Motifs (Motif TEXT PRIMARY KEY, Count INT);"
        "CREATE TABLE ProcessingMetrics (TotalReads INT, Runtime REAL, Efficiency REAL, Threads INT);"
        "CREATE TABLE QCMetrics (TotalReads INT, PassedReads INT, Fragments INT, Unmapped INT,"
        " Secondary INT, Supplementary INT, QCFail INT, Duplicates INT);"
        "BEGIN TRANSACTION;", 0, 0, 0);

    sqlite3_stmt *st;
    sqlite3_prepare_v2(db, "INSERT INTO FragmentSizes VALUES (?, ?);", -1, &st, 0);
    for (int i = 0; i < MAX_ISIZE; i++) if (size_dist[i] > 0) {
        sqlite3_bind_int(st, 1, i); sqlite3_bind_int64(st, 2, size_dist[i]);
        sqlite3_step(st); sqlite3_reset(st);
    }
    sqlite3_finalize(st);

    sqlite3_prepare_v2(db, "INSERT INTO Motifs VALUES (?, ?);", -1, &st, 0);
    for (int i = 0; i < 256; i++) if (motif_counts[i] > 0) {
        char m[5]; motif_string(i, m);
        sqlite3_bind_text(st, 1, m, -1, SQLITE_TRANSIENT); sqlite3_bind_int64(st, 2, motif_counts[i]);
        sqlite3_step(st); sqlite3_reset(st);
    }
    sqlite3_finalize(st);

    sqlite3_prepare_v2(db, "INSERT INTO Mutations VALUES (?, ?, ?, ?);", -1, &st, 0);
    sqlite3_bind_int64(st, 1, wt); sqlite3_bind_int64(st, 2, mut);
    sqlite3_bind_double(st, 3, qual_sum_cnt ? sum_qual / qual_sum_cnt : 0);
    sqlite3_bind_double(st, 4, gc_reads ? sum_gc / gc_reads : 0);
    sqlite3_step(st); sqlite3_finalize(st);

    sqlite3_prepare_v2(db, "INSERT INTO QCMetrics VALUES (?, ?, ?, ?, ?, ?, ?, ?);", -1, &st, 0);
    long qc[8] = {total_reads, passed_reads, fragments, f_unmapped, f_secondary, f_supplementary, f_qcfail, f_dup};
    for (int i = 0; i < 8; i++) sqlite3_bind_int64(st, i + 1, qc[i]);
    sqlite3_step(st); sqlite3_finalize(st);

    clock_gettime(CLOCK_MONOTONIC, &end);
    double runtime = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    sqlite3_prepare_v2(db, "INSERT INTO ProcessingMetrics VALUES (?, ?, ?, ?);", -1, &st, 0);
    sqlite3_bind_int64(st, 1, total_reads); sqlite3_bind_double(st, 2, runtime);
    sqlite3_bind_double(st, 3, total_reads ? runtime * 1000.0 / total_reads : 0);
    sqlite3_bind_int(st, 4, n_threads);
    sqlite3_step(st); sqlite3_finalize(st);

    sqlite3_exec(db, "COMMIT;", 0, 0, 0);
    sqlite3_close(db);
    bam_destroy1(aln); sam_hdr_destroy(header); sam_close(fp);
    return 0;
}
