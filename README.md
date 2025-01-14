# FastBWA
Fast alignment tool based on bwa-mem

# 1. Compile the source code
make -j 8

# 2. Build the FMT-Index
./fastbwa index reference.fasta

# 3. Run sequence alignment with fastBWA
./fastbwa mem -t 64 -2 -M -R @RG\\tID:normal\\tSM:normal\\tPL:illumina\\tLB:normal\\tPG:fastbwa reference.fasta r1.fq.gz r2.fq.gz
