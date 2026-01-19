thread=64

n1=~/data/dataset/D2/n1.fq.gz
n2=~/data/dataset/D2/n2.fq.gz

reference=~/data/reference/fmt/human_g1k_v37_decoy.fasta

out=./out.sam

time ./fastalign mem -t $thread -M -R @RG\\tID:normal\\tSM:normal\\tPL:illumina\\tLB:normal\\tPG:bwa \
    $reference \
    $n1 \
    $n2 \
    -o $out -2
