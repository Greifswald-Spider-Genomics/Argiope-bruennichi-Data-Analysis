# Improving the assembly

The previously published assembly was improved by additional polishing steps described in the following.

An overview of the workflow used for improving the assembly of the genome of *Argiope bruennichi* is shown [here](workflow_assembly_improving.pdf).

## Data

- 21.8X PacBio long-reads

- 30X paired-end Illumina short-reads

- Omni-C data

## Polishing with HyPo using paired-end Illumina short-reads and PacBio long-reads

The previously published assembly was first improved by polishing with HyPo using paired-end Illumina short-reads and PacBio long-reads.

This was done by executing the following command:

```
# Step 1 - Mapping the short-reads to scaffolds
minimap2 --secondary=no --MD -ax sr -t 6 genome.fa PE_data_R1.fq.gz PE_data_R2.fq.gz | samtools view -bS - > mapped-sr.bam
samtools sort -@ 6 -o mapped-sr.sorted.bam mapped-sr.bam 
samtools index mapped-sr.sorted.bam
rm mapped-sr.bam

# Step 2 - Mapping the long-reads to scaffolds
minimap2 --secondary=no --MD -ax map-pb -t 8 genome.fa pacbio_raw.fq | samtools view -bS - > mapped-lg.bam
samtools sort -@ 6 -o mapped-lg.sorted.bam mapped-lg.bam
samtools index mapped-lg.sorted.bam
rm mapped-lg.bam

# Step 3 - Create a text file containing the names of the short reads files 
echo -e "PE_data_R1.fq.gz\nPE_data_R2.fq.gz" > il_names.txt

# Step 4 - Run HyPo
hypo -r @il_names.txt -d genome.fa -b mapped-sr.sorted.bam -c 30 -s 1.7g -B mapped-lg.sorted.bam -p 16 -t 8 -o genome_HyPo.fa -i
```

## Polishing with Pilon using Omni-C data

The resulting assembly was polished again with Pilon using Omni-C data.
