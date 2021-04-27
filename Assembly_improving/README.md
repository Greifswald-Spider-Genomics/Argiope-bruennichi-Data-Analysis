# Improving the assembly

The previously published assembly was improved by additional polishing steps described in the following.

An overview of the workflow used for improving the assembly of the genome of *Argiope bruennichi* is shown [here](workflow_assembly_improving.pdf).

## Data

- Genome assembly

- 21.8X PacBio long-reads

- 30X paired-end Illumina short-reads

- Omni-C data

## Polishing with HyPo using paired-end Illumina short-reads and PacBio long-reads

The previously published assembly was first improved by polishing with HyPo using paired-end Illumina short-reads and PacBio long-reads.

This was done by executing the following commands:

```
# Step 1 - Map short-reads to the assembly
minimap2 --secondary=no --MD -ax sr -t 6 genome.fa PE_data_R1.fq.gz PE_data_R2.fq.gz | samtools view -bS - > mapped-sr.bam
samtools sort -@ 6 -o mapped-sr.s.bam mapped-sr.bam 
samtools index mapped-sr.s.bam
rm mapped-sr.bam

# Step 2 - Map long-reads to the assembly
minimap2 --secondary=no --MD -ax map-pb -t 8 genome.fa pacbio_raw.fq | samtools view -bS - > mapped-lg.bam
samtools sort -@ 6 -o mapped-lg.s.bam mapped-lg.bam
samtools index mapped-lg.s.bam
rm mapped-lg.bam

# Step 3 - Create a text file containing the names of the short-reads files 
echo -e "PE_data_R1.fq.gz\nPE_data_R2.fq.gz" > il_names.txt

# Step 4 - Run HyPo
hypo -r @il_names.txt -d genome.fa -b mapped-sr.s.bam -c 30 -s 1.7g -B mapped-lg.s.bam -p 16 -t 8 -o genome_HyPo.fa -i
```

## Polishing with Pilon using Omni-C data

The resulting assembly was polished again with Pilon using Omni-C data.

This was done by executing the following commands:

```
# Step 1 - Index reference
bwa index -a bwtsw genome_HyPo.fa

# Step 2a - Map the Omni-C reads to the assembly FASTA file
bwa mem -t 10 genome_HyPo.fa omnic_1.fq.gz omnic_2.fq.gz | samtools view -bS | samtools sort -o omnic.s.bam
# Step 2b - index the mapped file
samtools index omnic.s.bam

# Step 2c - Polish the assembly with the OmniC data
java -Xmx40G -jar pilon-1.23.jar --genome genome_HyPo.fa --frags omnic.s.bam --output genome_new
```

## Results of the polishing

| Assembly | unpolished | polished with HyPo | polished with HyPo and Pilon 
| :---     | ---:       | ---:               | ---:                        |
|Number of scaffolds |2,231 |2,231 |2,231 
| Total length (Gb) |1.67 | 1.67 | 1.67 
|Number of Ns | 1,169,100 | 624,300 | 571,613 
|Number of gaps | 11,691 | 6,243 | 6,231 
|<td colspan=4>BUSCOs (%)  
<&emsp; Complete | 91.0 | 91.9 | 92.5 
|&emsp;&emsp; and single-copy | 87.6 | 88.6 | 89.3 
|&emsp;&emsp; and duplicated | 3.4 | 3.3 | 3.2 
|&emsp; Fragmented | 2.8 | 2.3 | 2.1 
|&emsp; Missing | 6.2 | 5.8 | 5.4 
