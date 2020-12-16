# Argiope-bruennichi-Data-Analysis
Documentation of Data Analysis in the Argiope bruennichi sequencing project

# Contents

- [Data](#data)
- [Assembly](#assembly)
- [Repeat Masking](#repeat-masking)
- [Annotation](#annotation)

# Data

- 21.8X PacBio long-reads

- 30X paired-end Illumina short-reads

- Hi-C data

- RNA-seq data


# Assembly

## Assembling using PacBio long-reads

```
# Step 1a - assemble with wtdbg2
wtdbg2 -x sq -t 30 -i pacbio_raw.fq -fo argiope_wtdbg2
# Step 1b - derive consensus
wtpoa - cns -t 30 -i argiope_wtdbg2.ctg.lay.gz -fo argiope_wtdbg2.ctg.fa
```

## Polishing using paired-end Illumina short-reads

```
# Step 2a - index reference
bwa index -a bwtsw argiope_wtdbg2.ctg.fa
# Step 2b - map the PE reads to the wtdbg2 FASTA file
bwa mem -t 10 argiope_wtdbg2.ctg.fa PE_data_R1 .fq.gz PE_data_R2 .fq.gz | samtools view -bS | samtools sort -o PE_data.sorted.bam
# Step 2c - index the mapped file
samtools index PE_data.sorted.bam
# Step 2d - Polish the wtdbg2 FASTA with the Illumina data
java -Xmx900G -jar pilon -1.23.jar -- genome argiope_wtdbg2.ctg.fa --frags PE_data.sorted.bam -- output argiope_wtdbg2_p1
# Repeat steps 2a – 2d twice (the resulting FASTA file is named argiope_wtdbg2_p3.fa)
```

## Scaffolding using Hi-C data


# Repeat Masking

```
# Repeat masking the contig assembly
BuildDatabase -name argiope -engine ncbi contig_assembly.fa
nohup RepeatModeler -engine ncbi -database argiope >& run.out &
RepeatMasker contig_assembly.fa -lib argiope-families.fa -s -xsmall
# Repeat masking the scaffolded assembly with the same library
RepeatMasker scaffolded_genome.fa -lib argiope-families.fa -s -xsmall
```

# Annotation

## Generating intron hints from RNA-Seq data


## Running AUGUSTUS


## Converting the GFF output file to GTF format


## Replacing genes with in-frame stop codon by newly predicted genes


## Extracting protein sequences of predicted genes


## Running the functional annotation with InterProScan

