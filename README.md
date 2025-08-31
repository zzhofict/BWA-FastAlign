# FastAlign
Fast alignment tool based on bwa-mem

# 1. Compile the source code
make -j 8

# 2. Build the FMT-Index
./fastalign index reference.fasta

# 3. Run sequence alignment with FastAlign
./fastalign mem -t 64 -2 -M -R @RG\\tID:normal\\tSM:normal\\tPL:illumina\\tLB:normal\\tPG:fastbwa reference.fasta r1.fq.gz r2.fq.gz
