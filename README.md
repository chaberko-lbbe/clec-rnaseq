# clec-rnaseq : *Cimex lectularius* differential expression analysis after insecticide exposure through RNA-sequencing

Contact : chloe.haberkorn@univ-lyon1.fr / chloehbk@gmail.com

## Short summary

We used two strain of *Cimex lectularius*: "London Lab", a susceptible strain to insecticides, and "London Field", with moderate resistant to pyrethroids.

For each strain, 24 unfed, virgin adult females of 7-days old were sampled.
Four replicates of 8 insects were performed for each condition: treated lab, untreated lab, treated field, untreated field; i.e. a total of 4 conditions. 
Treated insects were exposed to 0.1 ng of deltamethrin powder diluted in 1 μL of acetone, whereas control insects received 1 μL of acetone only. This dose of insecticide should correspond to a lethal dose for 50% of London Lab strain.
For each replicate, three individuals were randomly sampled **out of survivors**, and immediately flash-frozen in liquid nitrogen and stored at -80°, for a total of 48 individuals (3 individuals x 4 replicates x 4 conditions).

<img src="https://github.com/chaberko-lbbe/clec-rnaseq/blob/main/design-rnaseq.png" width="40%" height="40%">

## Table of Contents

- **[Installing tools](#Installing-tools)**

- **[RNA-seq data processing](#RNA-seq-data-processing)**
	- [Getting the data](#Getting-the-data)
	- [Checking quality](#Checking-quality)
	- [Mapping with STAR](#Mapping-with-STAR)

- **[Differential Expression analyzes](#Overall-SNPs-analyzes)**
	- [Formatting the data for R/DESeq](#Formatting-the-data-for-R/DESeq)
	- [Running and interpreting R/DESeq](#Running-and-interpreting-R/DESeq)
 

## Installing tools

Here are the tools and versions used (on a cluster): 
- FastQC 
- STAR v 2.7.3a
- R v3.5.2

They will be store in /your-path/Tools.
We also used R on a computer with packages DESeq2 v 1.34.0.

## RNA-seq data processing

The goal is first to map *Cimex lectularius* RNA-seq samples (48 individuals) on reference genome.

### Getting the data

Raw sequences (fastq.gz files) are available on SRA: xx
```
mkdir /your-path/clec_rnaseq/Raw_Clec
```

We used the recent reference genome and annotation, avalaible here: https://www.ncbi.nlm.nih.gov/assembly/GCF_000648675.2
```
mkdir /your-path/clec_rnaseq/Ref_Clec
```

### Checking quality

We used Fastqc to check the quality:
```
cd /your-path/clec_rnaseq/RawData/
/your-path/Tools/FastQC/fastqc --java /your-path/Tools/jre1.8.0_202/bin/java -f fastq *.fastq -o /your-path/clec_rnaseq/Fastqc_Results/
```

Quality seems good, so no need to trim it. Macrogen, the compagny which performed the sequencing, also removed adapters (adapter content on Fastqc report = 0\%). We can therefore directly switch to mapping.

### Mapping with STAR

We followed STAR workflow:
1 - Generating genome indexes files (ref genome Clec_2.1 in FASTA and GTF annotation fie : GCF_000648675.2)
2 - Mapping reads to the genome

First, we checked some parameters according to guidelines:
log2(856686592)/2 - 1) = 13.84 > we change "--genomeSAindexNbases 14" to "13" for mapping
min(18, log2[max(856686592/1574,100)]) = min(18,19) > we stick to the default value of "--genomeChrBinNbits 18"

1 - Indexing
```
cd /your-path/clec_rnaseq/
mkdir /your-path/clec_rnaseq/index_STAR_Clec

/your-path/Tools/STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /your-path/clec_rnaseq/index_STAR_Clec --genomeFastaFiles /your-path/clec_rnaseq/Ref_Clec/ncbi-genomes-2020-12-15/GCF_000648675.2_Clec_2.1_genomic.fna --sjdbGTFfile /your-path/clec_rnaseq/Ref_Clec/ncbi-genomes-2022-03-16/GCF_000648675.2_Clec_2.1_genomic.gtf --sjdbOverhang 100
```

2 - Let's map our data !

We used --quantMode GeneCounts to "count number reads per gene while mapping", which requires GTF or GFF file provided with –sjdbGTFfile option.
A read is counted if it overlaps (1nt or more) one and only one gene. Both ends of the pairedend read are checked for overlaps. 

```
mkdir /your-path/clec_rnaseq/mapping_STAR_Clec/
mkdir /your-path/clec_rnaseq/mapping_STAR_Clec/results_mapping/

/your-path/Tools/STAR-2.7.3a/bin/Linux_x86_64/STAR --runThreadN 8 --readFilesCommand gunzip -c --readFilesIn /your-path/clec_rnaseq/RawData/samplename_1.fastq.gz /beegfs/data/chaberkorn/Data_Macrogen/2201KNO-0014/01.RawData/samplename_2.fastq.gz --genomeDir /your-path/clec_rnaseq/index_STAR_Clec --genomeSAindexNbases 13 --genomeChrBinNbits 18 --quantMode GeneCounts --sjdbOverhang 100 --outFileNamePrefix /your-path/clec_rnaseq/mapping_STAR_Clec/results_mapping/
```

STAR gave us some important informations about the mapping, as in Log.final.out files, with the mean of uniquely mapped reads = 84,09% among all samples.

STAR outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which correspond to different strandedness options.
We used the 4th column in following analysis, since the sum of 4th column counts is superior to 3rd one.

Example on R:
```
Sample_2E1LL <- read.table(file="2E1LL_ReadsPerGene.out.tab", sep="\t")
Sample_2E1LL <- Sample_2E1LL[grepl('LOC.*', Sample_2E1LL$V1),]

sum(Sample_2E1LL[,4]) # = 10707278
sum(Sample_2E1LL[,3]) # = 223975.
```

## Differential Expression analyzes

Our goal was to performed a differential expression analyzes between lab and field strain, whether they survived an insecticide exposure or not.
We used as input "ReadsPerGene.out.tab" files 

We needed for DESeq2 to build a "cts" matrix with gene count by gene for each sample, together with a "coldata" table with infos on data :
- untreated (T) or treated (E as "exposed")
- lab (LL as "London Lab") or field (LF) strain
- biological replicate number from one to four.

### Formatting the data for R/DESeq

As we briefly saw before, we first need to format the output data from STAR.
Load all the "ReadsPerGene.out.tab" files, and select only rows truly corresponding to gene counts:
```
Sample_1T1LL <- read.table(file="1T1LL_ReadsPerGene.out.tab", sep="\t") # 14,419 rows
colnames(Sample_1T1LL) <- c("GeneID", "Counts_Unstranded", "Counts_1stReadStrand","Counts_2ndReadStrand")
Sample_1T1LL <- Sample_1T1LL[grepl('LOC.*', Sample_1T1LL$GeneID),] # 13,208 rows
```

Build a "cts" matrix with gene count by gene for each sample:
```
cts_1 <- cbind(Sample_1T1LL[,c(1,4)],Sample_2T1LL[,4],Sample_3T1LL[,4],Sample_4T2LL[,4],Sample_5T2LL[,4],Sample_6T2LL[,4],
               Sample_7T3LL[,4],Sample_8T3LL[,4],Sample_9T3LL[,4],Sample_10T4LL[,4],Sample_11T4LL[,4],Sample_12T4LL[,4],
               Sample_1E1LL[,4],Sample_2E1LL[,4],Sample_3E1LL[,4],Sample_4E2LL[,4],Sample_5E2LL[,4],Sample_6E2LL[,4],
               Sample_7E3LL[,4],Sample_8E3LL[,4],Sample_9E3LL[,4],Sample_10E4LL[,4],Sample_11E4LL[,4],Sample_12E4LL[,4],
               Sample_1T1LF[,4],Sample_2T1LF[,4],Sample_3T1LF[,4],Sample_4T2LF[,4],Sample_5T2LF[,4],Sample_6T2LF[,4],
               Sample_7T3LF[,4],Sample_8T3LF[,4],Sample_9T3LF[,4],Sample_10T4LF[,4],Sample_11T4LF[,4],Sample_12T4LF[,4],
               Sample_1E1LF[,4],Sample_2E1LF[,4],Sample_3E1LF[,4],Sample_4E2LF[,4],Sample_5E2LF[,4],Sample_6E2LF[,4],
               Sample_7E3LF[,4],Sample_8E3LF[,4],Sample_9E3LF[,4],Sample_10E4LF[,4],Sample_11E4LF[,4],Sample_12E4LF[,4])
               
colnames(cts_1) <- c("GeneID","untreatedlab1","untreatedlab2","untreatedlab3","untreatedlab4","untreatedlab5","untreatedlab6",
                     "untreatedlab7","untreatedlab8","untreatedlab9","untreatedlab10","untreatedlab11","untreatedlab12",
                     "treatedlab1","treatedlab2","treatedlab3","treatedlab4","treatedlab5","treatedlab6",
                     "treatedlab7","treatedlab8","treatedlab9","treatedlab10","treatedlab11","treatedlab12",
                     "untreatedfield1","untreatedfield2","untreatedfield3","untreatedfield4","untreatedfield5","untreatedfield6",
                     "untreatedfield7","untreatedfield8","untreatedfield9","untreatedfield10","untreatedfield11","untreatedfield12",
                     "treatedfield1","treatedfield2","treatedfield3","treatedfield4","treatedfield5","treatedfield6",
                     "treatedfield7","treatedfield8","treatedfield9","treatedfield10","treatedfield11","treatedfield12")

rownames(cts_1) <- cts_1[,1]
cts_1 <- cts_1[,-1]
cts_1 <- as.matrix(cts_1)
```

Build a "coldata" table with infos on data:

```
coldata_1 <- data.frame(treatment=c(rep("untreated",12), rep("treated",12),rep("untreated",12), rep("treated",12)),
                        strain=c(rep("lab",24), rep("field",24)))

rownames(coldata_1) <- colnames(cts_1) # all samples names

coldata_1$treatment <- factor(coldata_1$treatment)
coldata_1$strain <- factor(coldata_1$strain)
```

### Running and interpreting R/DESeq

Now try DESeq package !
```
library("DESeq2")

ncol(cts_1) # = 48
nrow(coldata_1) # = 48 > check same value as above

dds_1 <- DESeqDataSetFromMatrix(countData = cts_1,
                                colData = coldata_1,
                                design = ~strain+treatment+strain:treatment) # add interaction
dds_1$treatment <- relevel(dds_1$treatment, ref = "untreated") # define reference
dds_1$strain <- relevel(dds_1$strain, ref = "lab")

dds_1 <- DESeq(dds_1)
```

