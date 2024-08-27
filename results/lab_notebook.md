---
title: "Lab Notebook v0.1"
---

# Haplotype Threading Lab Notebook

The goal of this project is to determine if parts of the human pangenome are associated with certain traits. Kind of like a GWAS except instead of SNPs we will use paths of a pangenome.

We have to answer several questions first though:

1. Can we associate each individual with a specific path through a pangenome via alignment.
2. Can we break up those paths into chunks that make sense and can be associated with a trait via GWAS?

For this directory I'm copying the project structure reccomended by [this paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424) roughly.

Log syntax
```bash
function TIMESTAMP () { echo $(TZ=US/Pacific; date '+%Y-%m-%d_%H:%M:%S'); }
[command] 2>&1 | tee results/[out_date]/[name]_$(TIMESTAMP).log
```

All commands in this file will be from the root of the git repo haplotype-threading. Also, `./bin` is prepended to my PATH. 

## 2024-08-27 Compare minigraph and graphaligner

Working on tscc.

```bash
OUTPUT=results/2024-08-27
```

We have two graph alignment tools: `minigraph` and `graphaligner` and I want to test them on a trivial base case to see if one performs better than the other. So I'm going to test aligning GRCh38 to v1.1 of the hprc minigraph-cactus graph as it was used to generate the graph and should have a near perfect alignment.

We will compare the tools: runtime and alignment accuracy and the output files types

### Subset graph and grch38 to a loci
This is required to have any reasonable runtime. We'll choose two different regions so that we can have one control.

I highly complex region from [here](https://github.com/CAST-genomics/snakemake-gfa-complexity/blob/main/highly_complex_regions) and a simple region from another analysis of the same size.

- Complex Region: chr1:16478848-16878847
- Simple Region:  chr4:107201506-107601505

Linear reference genome at simple region is entirely N's.

```bash
# Subset graph
gfabase sub data/hprc-v1.1-mc-grch38.gfab --range GRCh38#chr1:16478848-16878847 --view -o $OUTPUT/chr1:16478848-16878847.gfa

gfabase sub data/hprc-v1.1-mc-grch38.gfab --range GRCh38#chr4:107201506-107601505 --view -o $OUTPUT/chr4:107201506-107601505.gfa

# Subset linear reference genome
samtools faidx data/hg38.fa chr1:16478848-16878847 > $OUTPUT/chr1:16478848-16878847.fa

samtools faidx data/hg38.fa chr4:107201506-107601505 > $OUTPUT/chr4:107201506-107601505.fa
```

We can try subsampling sequences as well using pbsim, if it is neccessary.

#### graphaligner:
Graphs must be in gfa or vg file formats. I'm skeptical that this will be efficient. 

Has several required parameters:
- `precise-clipping` is the minimum sequence similarity for an alignment before it's tossed, has suggested values
- `-x vg` is variation graph defaults, `-x dbg` for de Brujin graphs. In general the de Brujin graphs are more complex so may fit a large pangenome use case better.

```bash
conda create -n haplotype-threading
conda activate haplotype-threading
conda install -c bioconda graphaligner

GraphAligner -g $OUTPUT/chr1:16478848-16878847.gfa -f $OUTPUT/chr1:16478848-16878847.fa -t 2 -a $OUTPUT/hg38_v1.1-mc-grch38.gam --precise-clipping 0.95 -x dbg 2>&1 | tee $OUTPUT/logs/chr1tochr1-graphaligner_$(TIMESTAMP).log

GraphAligner -g $OUTPUT/chr4:107201506-107601505.gfa -f $OUTPUT/chr4:107201506-107601505.fa -t 2 -a $OUTPUT/hg38_v1.1-mc-grch38.gam --precise-clipping 0.95 -x dbg 2>&1 | tee $OUTPUT/logs/chr4-correct-alignment_$(TIMESTAMP).log

GraphAligner -g $OUTPUT/chr4:107201506-107601505.gfa -f $OUTPUT/chr1:16478848-16878847.fa -t 2 -a $OUTPUT/hg38_v1.1-mc-grch38.gam --precise-clipping 0.95 -x dbg 2>&1 | tee $OUTPUT/logs/chr4tochr1-alignment_$(TIMESTAMP).log
# 24,747 / 40,000 bases aligned --> no alignment

GraphAligner -g $OUTPUT/chr1:16478848-16878847.gfa -f $OUTPUT/chr4:107201506-107601505.fa -t 2 -a $OUTPUT/hg38_v1.1-mc-grch38.gam --precise-clipping 0.95 -x dbg 2>&1 | tee $OUTPUT/logs/chr1tochr4-alignment_$(TIMESTAMP).log
# 16,201 / 40,000 bases aligned --> no alignment

conda deactivate
```
All process instantly.

#### minigraph:

Default `-x asm` exists for assemblies. Should also test `lr` for long reads.

```bash
git submodule add https://github.com/lh3/minigraph src/minigraph
cd src/minigraph; make; 
cp minigraph ../../bin; cd ../..

minigraph -x asm $OUTPUT/chr1:16478848-16878847.gfa $OUTPUT/chr1:16478848-16878847.fa -t 4 > $OUTPUT/chr1tochr1-minigraph.gaf 2> >(tee $OUTPUT/logs/chr1tochr1-minigraph_$(TIMESTAMP).log >&2)
# Takes about 24 seconds, 3 different alignments, 1 the total length + 2 others 50% length. Found more than graphaligner due to lower threshold

minigraph -x asm $OUTPUT/chr4:107201506-107601505.gfa $OUTPUT/chr4:107201506-107601505.fa -t 4 > $OUTPUT/chr4tochr4-minigraph.gaf 2> >(tee $OUTPUT/logs/chr4tochr4-minigraph_$(TIMESTAMP).log >&2)
# Takes 5 seconds, 1 sequence aligned, almost all bases mapped

minigraph -x asm $OUTPUT/chr1:16478848-16878847.gfa $OUTPUT/chr4:107201506-107601505.fa -t 4 > $OUTPUT/chr1tochr4-minigraph.gaf 2> >(tee $OUTPUT/logs/chr1tochr4-minigraph_$(TIMESTAMP).log >&2)
# Instant, no alignments found

minigraph -x asm $OUTPUT/chr4:107201506-107601505.gfa $OUTPUT/chr1:16478848-16878847.fa -t 4 > $OUTPUT/chr4tochr1-minigraph.gaf 2> >(tee $OUTPUT/logs/chr4tochr1-minigraph_$(TIMESTAMP).log >&2)
# Instant, no alignments found
```

### Test Alignments with simulated sequences

Use pbsim3 for simulating long reads from a reference genome.

```bash
git submodule add https://github.com/yukiteruono/pbsim3 src/pbsim3
cd src/pbsim3; ./configure; make
cp src/pbsim ../../bin/; cd ../..

# Long reads are low coverage in All of Us
pbsim --prefix $OUTPUT/simulated_chr4 \
      --seed 1 \
      --strategy wgs \
      --genome $OUTPUT/chr4:107201506-107601505.fa \
      --depth 8 \
      --errhmm src/pbsim3/data/ERRHMM-RSII.model \
      --method errhmm src/pbsim3/data/ERRHMM-RSII.model \
      2>&1 | tee $OUTPUT/logs/simulated_chr4_$(TIMESTAMP).log

pbsim --prefix $OUTPUT/simulated_chr1 \
      --seed 1 \
      --strategy wgs \
      --genome $OUTPUT/chr1:16478848-16878847.fa \
      --depth 8 \
      --errhmm src/pbsim3/data/ERRHMM-RSII.model \
      --method errhmm src/pbsim3/data/ERRHMM-RSII.model \
      2>&1 | tee $OUTPUT/logs/simulated_chr1_$(TIMESTAMP).log

```

***#Idea*** 
This produces .fastq, .maf, and .ref files that can be used with a tool like quast once you generate a .bam to quantify the accuracy of the alignment. But no such tool exists for .gaf/.gam files yet. This is something we could build.

Alignment tools take in fasta files not fastqs so we have to convert them

```bash
sed -n '1~4s/^@/>/p;2~4p' $OUTPUT/simulated_chr4_0001.fastq > $OUTPUT/simulated_chr4.fa
sed -n '1~4s/^@/>/p;2~4p' $OUTPUT/simulated_chr1_0001.fastq > $OUTPUT/simulated_chr1.fa
```

#### graphaligner
Need vg to view graphaligner .gam files

```bash
wget https://github.com/vgteam/vg/releases/download/v1.59.0/vg
chmod +x vg; mv vg bin

vg view -a [path_to_.gam]
```

```bash
conda activate haplotype-threading

GraphAligner -g $OUTPUT/chr1:16478848-16878847.gfa -f $OUTPUT/simulated_chr1.fa -t 2 -a $OUTPUT/hg38_v1.1-mc-grch38.gam --precise-clipping 0.75 -x dbg -a $OUTPUT/simalignment_graphaligner_chr1tochr1.gam 2>&1 | tee $OUTPUT/logs/simalignment_graphaligner_chr1tochr1_$(TIMESTAMP).log
# 291/353 reads with an alignment
# 3395901/3200047 input bps aligned.
# 88 end to end alignments (851714bp)
# 2 seconds

GraphAligner -g $OUTPUT/chr4:107201506-107601505.gfa -f $OUTPUT/simulated_chr4.fa -t 2 -a $OUTPUT/hg38_v1.1-mc-grch38.gam --precise-clipping 0.75 -x dbg -a $OUTPUT/simalignment_graphaligner_chr4tochr4.gam 2>&1 | tee $OUTPUT/logs/simalignment_graphaligner_chr4tochr4_$(TIMESTAMP).log
# 317/342 reads with an alignment
# 4568702/3200052 input bps aligned
# 77 end to end alignments
# 2 seconds

GraphAligner -g $OUTPUT/chr4:107201506-107601505.gfa -f $OUTPUT/simulated_chr1.fa -t 2 -a $OUTPUT/hg38_v1.1-mc-grch38.gam --precise-clipping 0.75 -x dbg -a $OUTPUT/simalignment_graphaligner_chr1tochr4.gam 2>&1 | tee $OUTPUT/logs/simalignment_graphaligner_chr1tochr4_$(TIMESTAMP).log
# 209/353 reads with an alignment
# 98457/3200047 input bps aligned
# 0 end to end alignments
# instant

GraphAligner -g $OUTPUT/chr1:16478848-16878847.gfa -f $OUTPUT/simulated_chr4.fa -t 2 -a $OUTPUT/hg38_v1.1-mc-grch38.gam --precise-clipping 0.75 -x dbg -a $OUTPUT/simalignment_graphaligner_chr4tochr1.gam 2>&1 | tee $OUTPUT/logs/simalignment_graphaligner_chr4tochr1_$(TIMESTAMP).log
# 186/342 read with alignment
# 87838/3200052 bps aligned
# 0 end to end alignments
# instant

conda deactivate
```
Correctly makes 0 alignments in the wrong region. 

But only correctly maps a small portion of the reads. Seems like it might struggle with low quality assemblies or reads.

If some reads can pass through multiple paths in a bubble, we have repeated alignments when one read passes through one location.

#### minigraph

```bash
minigraph -x lr $OUTPUT/chr1:16478848-16878847.gfa $OUTPUT/simulated_chr1.fa -t 4 > $OUTPUT/simalignment_minigraph_chr1tochr1.gaf 2> >(tee $OUTPUT/logs/simalignment_minigraph_chr1tochr1_$(TIMESTAMP).log >&2)
# 1.5 seconds
# mapped 265/353 reads

minigraph -x lr $OUTPUT/chr4:107201506-107601505.gfa $OUTPUT/simulated_chr4.fa -t 4 > $OUTPUT/simalignment_minigraph_chr4tochr4.gaf 2> >(tee $OUTPUT/logs/simalignment_minigraph_chr4tochr4_$(TIMESTAMP).log >&2)
# 1 second
# mapped 302/342 reads

minigraph -x lr $OUTPUT/chr1:16478848-16878847.gfa $OUTPUT/simulated_chr4.fa -t 4 > $OUTPUT/simalignment_minigraph_chr4tochr1.gaf 2> >(tee $OUTPUT/logs/simalignment_minigraph_chr4tochr1_$(TIMESTAMP).log >&2)
# instant
# mapped 3/342 reads

minigraph -x lr $OUTPUT/chr4:107201506-107601505.gfa $OUTPUT/simulated_chr1.fa -t 4 > $OUTPUT/simalignment_minigraph_chr1tochr4.gaf 2> >(tee $OUTPUT/logs/simalignment_minigraph_chr1tochr4_$(TIMESTAMP).log >&2)
# instant
# 17/353 reads
```

Minigraph doesn't actually calculate any base pair level alignments. It just finds seeds and then maps reads to gfa segments. According to their paper this means they are slightly worse at determining which path a read is on, but faster and better at determining the general loci a read is placed at.

### Summary
graphaligner allows partial alignments of reads and allows reads to align in multiple locations. It also gives per base level alignments. But it was slightly less accurate the minigraph at actually determining the correct loci for a read.

minigraph aligns each read at one location and only gives the segments that it believes the path travels. We might be over-confident, then, that some high accuracy alignment is actually the best.

### Conclusion
I'm inclined to favor graphaligner in a vacuum because we can use the multiple alignments to decide to drop certain reads and maintain a really high accuracy.

In an ideal world the best method might be to map all reads to the pangenome by minigraph to determine the correct loci. Then on a per-loci level use graphaligner to determine the specific paths favored.

Next steps will be to try using graphaligner on the All of Us portal.

## 2024-08-27

