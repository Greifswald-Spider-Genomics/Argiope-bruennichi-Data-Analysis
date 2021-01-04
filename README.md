# Argiope-bruennichi-Data-Analysis
Documentation of Data Analysis in the Argiope bruennichi sequencing project

An overview of the workflow used for assembling and annotating the genome of *Argiope bruennichi* is shown [here](workflow.pdf).

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

A first draft assembly was generated by assembling the 21.8X PacBio long-reads using the wtdbg2/Readbean (v. 2.3) <b id="f1">[1]</b>  assembler. 

This was done by executing the following commands:

```
# Step 1a - assemble with wtdbg2
wtdbg2 -x sq -t 30 -i pacbio_raw.fq -fo argiope_wtdbg2
# Step 1b - derive consensus
wtpoa-cns -t 30 -i argiope_wtdbg2.ctg.lay.gz -fo argiope_wtdbg2.ctg.fa
```
By the first command, wtdbg2 assembles the raw reads and generates the contig layout, where the option ```-x sq``` specifies the sequencing technology as being PacBio Sequel, ```-i``` specifies the raw sequencing reads and ```-fo``` specifies the prefix of the output files. By the second command, the consenser ```wtpoa-cns``` produces the final consensus in FASTA format, where the option ```-i``` specifies the raw sequencing reads and ```-fo``` specifies the name of output file.


## Polishing using paired-end Illumina short-reads

The first draft assembly was polished by applying three runs of Pilon (v. 1.23) <b id="f2">[2]</b> using 30X paired-end Illumina short-reads.

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

The polished assembly was scaffolded with HiRise (v. 2.1.7) <b id="f3">[3]</b> using Hi-C data.

This step was performed by Dovetail Genomics.

The resulting assembly is the final assembly and is also called "scaffolded assembly".


# Repeat Masking

The assembly was repeat-masked using the *de novo* repeat finder RepeatModeler (v. open-1.0.11) <b id="f4">[4]</b> and the homology-based repeat finder RepeatMasker (v. open-4.0.9) <b id="f5">[5,6]</b>.

First, RepeatModeler was used to generate a repeat library for the contig assembly. <br>
Since RepeatModeler needs a database as input, this was created first using ```BuildDatabase``` with the options ```-name argiope``` to specify the database name and ```-engine ncbi``` to specify the used search engine:

```
BuildDatabase -name argiope -engine ncbi contig_assembly.fa
```

RepeatModeler was run with the *de novo* repeat finder RECON (v. 1.08) <b id="f7">[7]</b>, the search engine RMBlast (v. 2.2.27), *de novo* prepeat finder RepeatScout (v. 1.0.5) <b id="f8">[8]</b> and tandem repeat finder TRF (v. 4.0.9) <b id="f9">[9]</b> with the options ```-database argiope``` to specify the database name and ```-engine ncbi``` to specify the used search engine:

```
nohup RepeatModeler -engine ncbi -database argiope >& run.out &
```

Then, RepeatMasker was used, to find and mask the repeats in the scaffolded assembly using the previously generated repeat library.
RepeatMasker was run with the RepeatMasker combined database Dfam_3.0 and the search engine RMBlast (v. 2.9.0) with the options ```-lib argiope-families.fa``` to include the previously generated repeat library, ```-s``` for using a slow search which is 0-5 % more sensitive but 2-3 times slower than default search and ```-xsmall``` to generate a soft-masked assembly instead of hard-masked:

```
RepeatMasker scaffolded_genome.fa -lib argiope-families.fa -s -xsmall
```

# Annotation

Predicting protein-coding genes was done with AUGUSTUS (v. 3.3.2) <b id="f10">[10,11,12]</b> including extrinsic evidence in form of RNA-seq reads to genome alignments.

## Generating intron hints from RNA-Seq data

To generate so called intron hints from RNA-seq data, which serve as extrinsic evidence for gene prediction, the first step was aligning RNA-seq reads to the genome assembly. <br>
This was done with HISAT2 (v. 2.1.0) <b id="f13">[13]</b>.

First, the genome assembly was indexed using ```hisat2-build``` with the following command:

```
hisat2-build -p 12 genome.fa genome.index
```
Where ```genome.index``` is the base of the HISAT2 index files.

Then, all of the four RNA extractions (from eggs, spiderlings, an adult male and an adult female: ```eggs_1.fq, eggs_2.fq, spiderlings_1.fq, spiderlings_2.fq, 
adultMale_1.fq, adultMale_2.fq, adultFemale_1.fq, adultFemale_2.fq```) were aligned against the genome assembly by executing the following commands:

```
files =" eggs spiderlings adultMale adultFemale "
for f in $files;
  do
  echo $f
  hisat2 -p 12 -x genome.index -1 $f_1.fq -2 $f_2.fq -S $f.sam
done
```
Where the option ```-x``` specifies the index file name prefix, ```-1``` and ```-2``` specify the paired-end RNA-seq FASTQ files and ```-S``` specifies the name of the SAM output file. This resulted in four SAM files.

Next, these SAM files were converted to intron hints, by first converting the SAM files to BAM files, then sorting the BAM files using SAMtools (v. 1.7) <b id="f14">[14]</b> and afterwards extracting the intron information from the sorted BAM files using the AUGUSTUS auxiliary tool ```bam2hints``` found in ```Augustus/auxprogs/bam2hints```. <br>
This was done by executing the following commands:

```
for f in $files;
  do
  echo $f
  samtools view -bS -o $f. bam $f.sam
  samtools sort $f.bam -o $f.s.bam
  samtools sort -n $f.s.bam -o $f.ss.bam
  bam2hints --intronsonly --in=$f.ss.bam --out=$f.intron.hints
done
```

Since for most RNA-seq data, it is unclear from which strand an aligned read stems, in the next step the correct strand was tried to be guessed from genomic splice site information and those reads without appropriate splice site information were filtered out using the BRAKER <b id="f15">[15]</b> script ```filterIntronsFindStrand.pl``` found in ```BRAKER/scripts``` with the option ```--score``` (to set the score column to the 'mult' entry) by executing the following command:

```
for f in $files;
  do
  filterIntronsFindStrand.pl genome.RMsoft.fa $f.intron.hints --score > $f.intron.f.score.hints
done
```
 Afterwards, the four generated hints files were merged into one file by executing the following command:
```
for f in $files;
  do
  cat $f.intron.f.score.hints >> hintsfile.tmp.gff
done
```

Then, the resulting hints file was sorted by executing:
```
cat hintsfile.tmp.gff | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -n -k 3,3 | sort -s -k 1,1 > hints.tmp.sort.gff
```
Finally, the multiple hints were joined using the AUGUSTUS script ```join_mult_hints.pl``` that can be found in ```Augustus/scripts```.
This was done by executing the following command:
```
perl join_mult_hints.pl < hints.tmp.sort.gff > hintsfile.gff
```

## Running AUGUSTUS

Protein-coding genes were predicted for the soft-masked genome assembly using AUGUSTUS with the *Parasteatoda* parameter set. <br>
In order to speed up gene prediction, the genome and the hints file were split into smaller tasks to be executed in parallel. 
This was done according to Support Protocol 9 of <b id="f10">[10]</b> using the following commands:

```
# Split the genome into smaller files
mkdir split
splitMfasta.pl genome.RMsoft.fa --outputpath=split --minsize=1000000

# Determine the number of split genome files
ls split/genome.RMsoft.split.*.fa | wc -l 
# 37

# Split hints file
for ((i=1; i<=37; i++));
  do
  fgrep ">" split/genome.RMsoft.split.$i.fa | perl -pe 's/^>//' > part.$i.lst
done
mkdir split_hints
for ((i=1; i<=37; i++));
  do
  cat hintsfile.gff | getLinesMatching.pl part.$i.lst 1 > split_hints/hints.split.$i.gff
done
```
```splitMfasta.pl``` and ```getLinesMatching.pl``` are AUGUSTUS scripts that can be found in ```Augustus/scripts```.

Having split the genome assembly and hints file into smaller parts, AUGUSTUS was run with the options ```--UTR=on``` to additionally predict untranslated regions (UTRs), ```--print_utr=on``` to print the predicted UTRs in the GFF output file, ```--species=parasteatoda``` to specifiy the used species parameters, ```--alternatives-from-evidence=true``` to report alternative transcripts when they are supported by hints, ```--hintsfile=split_hints/hints.split.{}.gff``` to specify the hints file, ```--extrinsicCfgFile=rnaseq.cfg``` to specify the file containing the list of used sources for the hints and their boni and mali (here, the default file used by BRAKER was used), ```--allow_hinted_splicesites=gcag,atac``` to allow AUGUSTUS to also predict those (rare) introns that start with GC and end with AG and those that start with AT and end with AC, ```--softmaskin=on``` to specify that the assembly is soft-masked, ```--codingseq=on``` to output the coding DNA sequences in the GFF output file and additionally produce a file augustus.codingseq with the complete coding sequences and ```--exonnames=on``` to, in addition to the exon-line in the GFF output file, print the same line again but with the exon name ('initial', 'internal', 'terminal' or 'single') instead of 'exon':

```
# Run AUGUSTUS
seq 1 37 | parallel -j 8 --bar --no-notice " nice augustus \
  --UTR=on --print_utr=on --species=parasteatoda \
  --alternatives-from-evidence=true \
  --hintsfile=split_hints/hints.split.{}.gff \
  --extrinsicCfgFile=rnaseq.cfg \
  --allow_hinted_splicesites=gcag,atac \
  --softmasking=on --codingseq=on --exonnames=on \
  split/genome.RMsoft.split.{}.fa \
  > out/augustus.{}.out"
```

Then, the AUGUSTUS output files were merged using the AUGUSTUS script ```join_aug_pred.pl``` by executing the following commands:
```
# Merge the AUGUSTUS outputs
for ((i=1; i<=37; i++));
  do
  cat out/augustus.$i.out | join_aug_pred.pl >> augustus.tmp.gff
done

join_aug_pred.pl < augustus.tmp.gff > augustus.gff
```

## Converting the GFF output file to GTF format

In the next step, the GFF output file was converted to GTF format using the AUGUSTUS script ```gtf2gff.pl``` by executing the following command:
```
cat augustus.gff | perl -ne 'if(m/\tAUGUSTUS\t/) { print $_ ;}' | perl gtf2gff.pl --printExon --out=augustus.gtf
```

## Replacing genes with in-frame stop codon by newly predicted genes

After that, additional AUGUSTUS scripts: ```getAnnoFastaFromJoinGenes.py``` and ```fix_in_frame_stop_codon_genes.py``` were used to find and replace predicted protein-coding genes containing an in-frame stop codon with newly predicted genes. This was done by executing the following commands:

```
getAnnoFastaFromJoingenes.py -g genome.RMsoft.fa -f augustus.gtf -o augustus
# WARNING: The GTF file contained 428 gene(s) with internal Stop codons.

fix_in_frame_stop_codon_genes.py -g genome.RMsoft.fa --gtf augustus.gtf -o augustus.hints -b bad_genes.lst
-H hintsfile.gff -s parasteatoda -e rnaseq.cfg --UTR on --print_utr on --softmasking on --additional_aug_args
" --allow_hinted_splicesites=gcag,atac --codingseq=on"
```

## Extracting protein sequences of predicted genes

Having prepared the final gene set, the next step was extracting protein sequences of the predicted protein-coding genes. This was done using the AUGUSTUS script ```getAnnoFastaFromJoinGenes.py``` with the following command:
```
getAnnoFastaFromJoingenes.py -g genome.RMsoft.fa -f augustus.hints.gtf -o augustus.hints
```

## Running the functional annotation with InterProScan

Finally, functional annotation was performed using InterProScan (v. 5.39-77.0) <b id="f16">[16,17]</b>.

Before running InterProScan, the stars (\*) at the end of each entry in the protein FASTA file ```augustus.hints.aa``` had to be removed, because InterProScan does
not accept sequences with the \* character. This was done by executing the following command:
```
cat augustus.hints.aa | perl -pe 's/\*//;' > augustus.hints_no_stars.aa
```
Then, InterProScan was run on the modified FASTA file with the following command:

```
interproscan.sh -i augustus.hints_no_stars.aa -pa &> iprsc.log
```

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

<b id="f7">[7]</b> Bao, Z. and Eddy, S. R. (2002). Automated de novo identification of repeat
sequence families in sequenced genomes. *Genome Research*, 12(8):1269-1276.

<b id="f8">[8]</b> Price, A. L., Jones, N. C., and Pevzner, P. A. (2005). *De novo* identification of
repeat families in large genomes. *Bioinformatics*, 21(suppl 1):i351-i358.

<b id="f9">[9]</b> Benson, G. (1999). Tandem repeats finder: a program to analyze DNA sequences.
*Nucleic Acids Research*, 27(2):573-580.

<b id="f10">[10]</b> Hoff, K. J. and Stanke, M. (2019). Predicting Genes in Single Genomes with
AUGUSTUS. *Current Protocols in Bioinformatics*, 65(1):e57.

<b id="f11">[11]</b> Stanke, M., Steinkamp, R., Waack, S., and Morgenstern, B. (2004).
AUGUSTUS: a web server for gene finding in eukaryotes. *Nucleic Acids Research*,
32(suppl 2):W309-W312.

<b id="f12">[12]</b> Stanke, M., Tzvetkova, A., and Morgenstern, B. (2006). AUGUSTUS at EGASP:
using EST, protein and genomic alignments for improved gene prediction in the
human genome. *Genome Biology*, 7(1):S11.

<b id="f13">[13]</b> Kim, D., Paggi, J. M., Park, C., Bennett, C., and Salzberg, S. L. (2019). Graphbased
genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology*, 37(8):907-915.

<b id="f14">[14]</b> Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth,
G., Abecasis, G., and Durbin, R. (2009). The sequence alignment/map format and
SAMtools. *Bioinformatics*, 25(16):2078-2079.

<b id="f15">[15]</b> Hoff, K. J., Lomsadze, A., Borodovsky, M., and Stanke, M. (2019). Whole-
Genome Annotation with BRAKER. In *Gene Prediction*, pages 65{95. Springer.

<b id="f16">[16]</b> Jones P, Binns D, Chang HY, et al. InterProScan 5: 
Genome-scale protein function classification. *Bioinformatics* 2014;30:1236–40.

<b id="f17">[17]</b> Quevillon E, Silventoinen V, Pillai S, et al. InterProScan: Protein
domains identifier. *Nucleic Acids Res* 2005;33:W116–20.
