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

A first draft assembly was generated by assembling the 21.8X PacBio long-reads using the wtdbg2/Readbean <b id="f1">[1]</b>  assembler. 

This was done by executing the following commands:

```
# Step 1a - assemble with wtdbg2
wtdbg2 -x sq -t 30 -i pacbio_raw.fq -fo argiope_wtdbg2
# Step 1b - derive consensus
wtpoa -cns -t 30 -i argiope_wtdbg2.ctg.lay.gz -fo argiope_wtdbg2.ctg.fa
```

## Polishing using paired-end Illumina short-reads

The first draft assembly was polished by applying three runs of Pilon <b id="f2">[2]</b> using 30X paired-end Illumina short-reads.

This was done by executing the following commands:
```
# Step 2a - index reference
bwa index -a bwtsw argiope_wtdbg2.ctg.fa
# Step 2b - map the PE reads to the wtdbg2 FASTA file
bwa mem -t 10 argiope_wtdbg2.ctg.fa PE_data_R1.fq.gz PE_data_R2.fq.gz | samtools view -bS | samtools sort -o PE_data.sorted.bam
# Step 2c - index the mapped file
samtools index PE_data.sorted.bam
# Step 2d - Polish the wtdbg2 FASTA with the Illumina data
java -Xmx900G -jar pilon-1.23.jar --genome argiope_wtdbg2.ctg.fa --frags PE_data.sorted.bam --output argiope_wtdbg2_p1
# Repeat steps 2a – 2d twice (the resulting FASTA file is named argiope_wtdbg2_p3.fa)
```
The resulting assembly is also called "contig assembly".

## Scaffolding using Hi-C data

The polished assembly was scaffolded with HiRise <b id="f3">[3]</b> using Hi-C data.

This step was performed by Dovetail Genomics.

The resulting assembly is the final assembly and is also called "scaffolded assembly".


# Repeat Masking

The assembly was repeat-masked using the *de novo* repeat finder RepeatModeler (v. open-1.0.11) <b id="f4">[4]</b> and the homology-based repeat finder RepeatMasker (v. open-4.0.8) <b id="f5">[5,6]</b>.

First, RepeatModeler was used to generate a repeat library for the contig assembly. <br>
Since RepeatModeler needs a database as input, this was created first using ```BuildDatabase``` with the options ```-name argiope``` to specify the database name and ```-engine ncbi``` to specify the used search engine:

```
BuildDatabase -name argiope -engine ncbi contig_assembly.fa
```

RepeatModeler was run with the *de novo* repeat finder RECON (v. 1.08) <b id="f7">[7]</b>, the search engine RMBlast (v. 2.2.27), *de novo* prepeat finder RepeatScout (v. 1.0.5) <b id="f8">[8]</b> and tandem repeat finder TRF (v. 4.0.9) <b id="f9">[9]</b> with the options ```-database argiope``` to specify the database name and ```-engine ncbi``` to specify the used search engine:

```
nohup RepeatModeler -engine ncbi -database argiope >& run.out &
```

Then, RepeatMasker was used, to find and mask the repeats in the contig assembly as well as the scaffolded assembly using the previously generated repeat library.
RepeatMasker was run with the RepeatMasker combined database Dfam_3.0 and the search engine RMBlast (v. 2.9.0) with the options ```-lib argiope-families.fa``` to include the previously generated repeat library, ```-s``` for using a slow search which is 0-5 % more sensitive but 2-3 times slower than default search and ```-xsmall``` to generate a soft-masked assembly instead of hard-masked:

```
# Repeat masking the contig assembly
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


# References

<b id="f1">[1]</b> Ruan, J. and Li, H. (2020). Fast and accurate long-read assembly with wtdbg2. *Nature Methods*, 17(2):155-158.

<b id="f2">[2]</b> Walker, B. J., Abeel, T., Shea, T., Priest, M., Abouelliel, A., Sakthikumar,
S., Cuomo, C. A., Zeng, Q., Wortman, J., Young, S. K., et al. (2014). Pilon: an
integrated tool for comprehensive microbial variant detection and genome assembly
improvement. *PLoS ONE*, 9(11):e112963.

<b id="f3">[3]</b> Putnam, N. H., O'Connell, B. L., Stites, J. C., Rice, B. J., Blanchette, M.,
Calef, R., Troll, C. J., Fields, A., Hartley, P. D., Sugnet, C. W., et al. (2016).
Chromosome-scale shotgun assembly using an in vitro method for long-range linkage.
*Genome Research*, 26(3):342-350.

<b id="f4">[4]</b> Smit, A. F. and Hubley, R. (2008-2015). RepeatModeler Open-1.0
[http://www.repeatmasker.org].

<b id="f5">[5]</b> Chen, N. (2004). Using RepeatMasker to identify repetitive elements in genomic
sequences. *Current Protocols in Bioinformatics*, 5(1):4-10.

<b id="f6">[6]</b> Smit, A. F., Hubley, R., and Green, P. (2013-2015). RepeatMasker Open-4.0
[http://www.repeatmasker.org].

<b id="f7">[7]</b> 

<b id="f8">[8]</b> 

<b id="f9">[9]</b> 

<b id="f10">[10]</b> 
