import random, os

def create_million_reads_bam(filename, peak_size, mut_rate):
    sam_name = filename.replace(".bam", ".sam")
    bases = ['A', 'C', 'G', 'T']
    print(f"Generating 1,000,000 reads for {filename}...")

    with open(sam_name, "w") as f:
        f.write("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:250000000\n")
        # Εδώ είναι το 1 εκατομμύριο (1,000,000)
        for i in range(1000000):
            isize = int(random.gauss(peak_size, 12))
            pos = 10000 + (i % 100000) # Reusing positions for speed

            prefix = "".join(random.choices(bases, k=4))
            is_mutant = random.random() < mut_rate
            base_at_10 = 'T' if is_mutant else 'G'
            rest_seq = "".join(random.choices(bases, weights=[0.3, 0.3, 0.2, 0.2], k=25))

            full_seq = prefix + "AAAAAA" + base_at_10 + rest_seq
            qual = "I" * len(full_seq)

            f.write(f"read_{i}\t99\tchr1\t{pos}\t60\t{len(full_seq)}M\t=\t{pos+isize}\t{isize}\t{full_seq}\t{qual}\n")

    os.system(f"samtools view -b {sam_name} > {filename}")
    os.remove(sam_name)
    print(f"✅ Successfully created {filename} with 1,000,000 entries.")

create_million_reads_bam("healthy.bam", 167, 0.01)
create_million_reads_bam("cancer.bam", 145, 0.18)
