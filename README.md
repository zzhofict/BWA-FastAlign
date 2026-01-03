# BWA-FastAlign: Faster and Cheaper Sequence Alignment on Commercial CPUs

**BWA-FastAlign** is a high-performance, cost-efficient software package for mapping low-divergent sequences against a large reference genome, such as the human genome.

It is designed as a drop-in replacement for the de facto standard **BWA-MEM**, offering **2.27× ∼ 3.28× throughput speedup** and **2.54× ∼ 5.65× cost reductions** on standard CPU servers, while guaranteeing **100% identical output** (SAM/BAM) to BWA-MEM.

## 🚀 Key Features

*   **High Throughput:** Achieves ~2.85× average speedup over BWA-MEM by optimizing both the seeding and extension phases.
*   **Cost Efficient:** Delivers 2.54× ∼ 5.65× cost reduction compared to state-of-the-art CPU and GPU baselines (including BWA-MEM2 and BWA-GPU).
*   **Identical Output:** Guarantees 100% output compatibility with BWA-MEM. You can swap it into your existing pipelines without changing downstream analysis results.
*   **Low Memory Footprint:** Uses a novel Multi-stage Seeding strategy (Hybrid Index) that improves search performance without the massive memory overhead seen in hash-based or learned-index aligners (e.g., ERT-BWA-MEM2).
*   **Optimized for Modern CPUs:** Features an Intra-query Parallel algorithm for the seed-extension phase, utilizing AVX2 instructions to eliminate computation bubbles caused by varying read lengths.

## 🔧 Technical Innovations

BWA-FastAlign revitalizes the traditional alignment pipeline with two core algorithmic contributions:

1.  **Multi-Stage Seeding (Hybrid Index)**
    *   Combines **Kmer-Index**, **FMT-Index** (Enhanced FM-Index with prefetching), and **Direct-Index**.
    *   Dynamically switches strategies based on seed length and match density.
    *   Achieves an **18.92× improvement in memory efficiency** (bases processed per GB per second).

2.  **Intra-Query Parallel Seed-Extension**
    *   Unlike BWA-MEM2 (which uses inter-query parallelism and suffers from load imbalance), BWA-FastAlign parallelizes the Smith-Waterman alignment *within* a single query.
    *   Includes **Dynamic Pruning** to skip zero-alignment scores.
    *   Implements a **Sliding Window** mechanism to reduce costly memory gather operations.
    *   Achieves **3.45× higher SIMD utilization**, performing consistently well on both WGS (Whole Genome Sequencing) and WES (Whole Exome Sequencing) data.

## 📥 Installation

### Prerequisites
*   Linux operating system (tested on Ubuntu 22.04).
*   GCC compiler (version 11.4 or higher recommended).
*   CPU supporting **AVX2** instructions (most modern Intel/AMD CPUs).
*   zlib development files.

### Compilation
```bash
git clone https://github.com/your-username/BWA-FastAlign.git
cd BWA-FastAlign
make
```

## 📖 Usage
BWA-FastAlign follows the same command-line interface as BWA-MEM.

1.  **Index the Reference.** Before alignment, you must index your reference genome (e.g., human_g1k_v37.fasta).
```bash
# This will generate the Hybrid Index files
./fastalign index ref.fa
```

2.  **Align Reads (Mem).** Map single-end or paired-end reads to the reference.
```bash
# Single-end alignment
./fastalign mem ref.fa reads.fq > aln.sam

# Paired-end alignment
./fastalign mem ref.fa read1.fq read2.fq > aln.sam

# Using multiple threads (Recommended: 32-128 threads for high throughput)
./fastalign mem -t 64 ref.fa read1.fq read2.fq > aln.sam
```

3.  **Options.** BWA-FastAlign supports the standard BWA-MEM options. Run ./fastalign mem to see the full list.

## 📜 Citation
If you find BWA-FastAlign is useful in your research, please cite our paper:
```bibtex
@inproceedings{fastalign2026,
  title={Faster and Cheaper: Pushing the Sequence Alignment Throughput with Commercial CPUs},
  author={Zhonghai Zhang, Yewen Li, Ke Meng, Chunming Zhang, Guangming Tan},
  booktitle={Proceedings of the 31st ACM SIGPLAN Symposium on Principles and Practice of Parallel Programming (PPoPP '26)},
  year={2026}
}
```
