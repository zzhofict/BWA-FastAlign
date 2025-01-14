thread=128

n1=~/data/dataset/real/D1/n1.fq.gz
n2=~/data/dataset/real/D1/n2.fq.gz

reference=~/data/fmt_ref/human_g1k_v37_decoy.fasta

out=./out.sam

time ./fastbwa mem -t $thread -M -R @RG\\tID:normal\\tSM:normal\\tPL:illumina\\tLB:normal\\tPG:bwa \
    $reference \
    $n1 \
    $n2 \
    -o $out -2
