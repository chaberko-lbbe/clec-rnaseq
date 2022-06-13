# clec-rnaseq : *Cimex lectularius* differential expression analysis after insecticide exposure through RNA-sequencing

Contact : chloe.haberkorn@univ-lyon1.fr / chloehbk@gmail.com

## Short summary

We used four strains of *Cimex lectularius*: "London Lab" and "German Lab", susceptible strains to pyrethroid insecticides, together with "London Field" and Sweden Field, with moderate to high resistant to pyrethroids.

For each strain, 24 unfed, virgin adult females of 7-days old were sampled.
Four replicates of 8 insects were performed for each condition: treated London Lab, untreated London Lab, treated German Lab, untreated German Lab, treated London Field, untreated London Field, treated Sweden Field, untreated Sweden Field ; i.e. a total of 8 conditions. 
Treated insects were exposed to 0.1 ng of deltamethrin powder diluted in 1 μL of acetone, whereas control insects received 1 μL of acetone only. This dose of insecticide should correspond to a lethal dose for 50% of susceptible strains.
For each replicate, three individuals were randomly sampled **out of survivors**, and immediately flash-frozen in liquid nitrogen and stored at -80°, for a total of 96 individuals (3 individuals x 4 replicates x 8 conditions).

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
- lab or field strain
- strain origin (LL for London Lab, GL for German Lab, LF for London Field, SF for Sweden Field)
- biological replicate number from one to four (according to [DESeq2 Group condition](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#group-specific-condition-effects-individuals-nested-within-groups), we have to distinguish the individuals nested within a group).

Values in the matrix should be un-normalized counts, as DESeq2 model internally corrects for library size.


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
cts <- cbind(Sample_1T1LL[,c(1,4)],Sample_2T1LL[,4],Sample_3T1LL[,4],Sample_4T2LL[,4],Sample_5T2LL[,4],Sample_6T2LL[,4],
             Sample_7T3LL[,4],Sample_8T3LL[,4],Sample_9T3LL[,4],Sample_10T4LL[,4],Sample_11T4LL[,4],Sample_12T4LL[,4],
             Sample_1E1LL[,4],Sample_2E1LL[,4],Sample_3E1LL[,4],Sample_4E2LL[,4],Sample_5E2LL[,4],Sample_6E2LL[,4],
             Sample_7E3LL[,4],Sample_8E3LL[,4],Sample_9E3LL[,4],Sample_10E4LL[,4],Sample_11E4LL[,4],Sample_12E4LL[,4],
             Sample_1T1LF[,4],Sample_2T1LF[,4],Sample_3T1LF[,4],Sample_4T2LF[,4],Sample_5T2LF[,4],Sample_6T2LF[,4],
             Sample_7T3LF[,4],Sample_8T3LF[,4],Sample_9T3LF[,4],Sample_10T4LF[,4],Sample_11T4LF[,4],Sample_12T4LF[,4],
             Sample_1E1LF[,4],Sample_2E1LF[,4],Sample_3E1LF[,4],Sample_4E2LF[,4],Sample_5E2LF[,4],Sample_6E2LF[,4],
             Sample_7E3LF[,4],Sample_8E3LF[,4],Sample_9E3LF[,4],Sample_10E4LF[,4],Sample_11E4LF[,4],Sample_12E4LF[,4],
             Sample_1T1GL[,4],Sample_2T1GL[,4],Sample_3T1GL[,4],Sample_4T2GL[,4],Sample_5T2GL[,4],Sample_6T2GL[,4],
             Sample_7T3GL[,4],Sample_8T3GL[,4],Sample_9T3GL[,4],Sample_10T4GL[,4],Sample_11T4GL[,4],Sample_12T4GL[,4],
             Sample_1E1GL[,4],Sample_2E1GL[,4],Sample_3E1GL[,4],Sample_4E2GL[,4],Sample_5E2GL[,4],Sample_6E2GL[,4],
             Sample_7E3GL[,4],Sample_8E3GL[,4],Sample_9E3GL[,4],Sample_10E4GL[,4],Sample_11E4GL[,4],Sample_12E4GL[,4],
             Sample_1T1SF[,4],Sample_2T1SF[,4],Sample_3T1SF[,4],Sample_4T2SF[,4],Sample_5T2SF[,4],Sample_6T2SF[,4],
             Sample_7T3SF[,4],Sample_8T3SF[,4],Sample_9T3SF[,4],Sample_10T4SF[,4],Sample_11T4SF[,4],Sample_12T4SF[,4],
             Sample_1E1SF[,4],Sample_2E1SF[,4],Sample_3E1SF[,4],Sample_4E2SF[,4],Sample_5E2SF[,4],Sample_6E2SF[,4],
             Sample_7E3SF[,4],Sample_8E3SF[,4],Sample_9E3SF[,4],Sample_10E4SF[,4],Sample_11E4SF[,4],Sample_12E4SF[,4])

colnames(cts) <- c("GeneID","untreatedLLlabone1","untreatedLLlabone2","untreatedLLlabone3",
                   "untreatedLLlabtwo4","untreatedLLlabtwo5","untreatedLLlabtwo6",
                   "untreatedLLlabthree7","untreatedLLlabthree8","untreatedLLlabthree9",
                   "untreatedLLlab10four","untreatedLLlabfour11","untreatedLLlabfour12",
                   "treatedLLlabone1","treatedLLlabone2","treatedLLlabone3",
                   "treatedLLlabtwo4","treatedLLlabtwo5","treatedLLlabtwo6",
                   "treatedLLlabthree7","treatedLLlabthree8","treatedLLlabthree9",
                   "treatedLLlabfour10","treatedLLlabfour11","treatedLLlabfour12",
                   "untreatedLFfieldone1","untreatedLFfieldone2","untreatedLFfieldone3",
                   "untreatedLFfieldtwo4","untreatedLFfieldtwo5","untreatedLFfieldtwo6",
                   "untreatedLFfieldthree7","untreatedLFfieldthree8","untreatedLFfieldthree9",
                   "untreatedLFfieldfour10","untreatedLFfieldfour11","untreatedLFfieldfour12",
                   "treatedLFfieldone1","treatedLFfieldone2","treatedLFfieldone3",
                   "treatedLFfieldtwo4","treatedLFfieldtwo5","treatedLFfieldtwo6",
                   "treatedLFfieldthree7","treatedLFfieldthree8","treatedLFfieldthree9",
                   "treatedLFfieldfour10","treatedLFfieldfour11","treatedLFfieldfour12",
                   "untreatedGLlabone1","untreatedGLlabone2","untreatedGLlabone3",
                   "untreatedGLlabtwo4","untreatedGLlabtwo5","untreatedGLlabtwo6",
                   "untreatedGLlabthree7","untreatedGLlabthree8","untreatedGLlabthree9",
                   "untreatedGLlabfour10","untreatedGLlabfour11","untreatedGLlabfour12",
                   "treatedGLlabone1","treatedGLlabone2","treatedGLlabone3",
                   "treatedGLlabtwo4","treatedGLlabtwo5","treatedGLlabtwo6",
                   "treatedGLlabthree7","treatedGLlabthree8","treatedGLlabthree9",
                   "treatedGLlabfour10","treatedGLlabfour11","treatedGLlabfour12",
                   "untreatedSFfieldone1","untreatedSFfieldone2","untreatedSFfieldone3",
                   "untreatedSFfieldtwo4","untreatedSFfieldtwo5","untreatedSFfieldtwo6",
                   "untreatedSFfieldthree7","untreatedSFfieldthree8","untreatedSFfieldthree9",
                   "untreatedSFfieldfour10","untreatedSFfieldfour11","untreatedSFfieldfour12",
                   "treatedSFfieldone1","treatedSFfieldone2","treatedSFfieldone3",
                   "treatedSFfieldtwo4","treatedSFfieldtwo5","treatedSFfieldtwo6",
                   "treatedSFfieldthree7","treatedSFfieldthree8","treatedSFfieldthree9",
                   "treatedSFfieldfour10","treatedSFfieldfour11","treatedSFfieldfour12")

rownames(cts) <- cts[,1]
cts <- cts[,-1]
cts <- as.matrix(cts)
```

Build a "coldata" table with infos on data:

```
coldata <- data.frame(treatment=c(rep(c(rep("untreated",12), rep("treated",12)),4)),
                      location=c(rep("LL",24), rep("LF",24), rep("GL",24), rep("SF",24)),
                      strain=c(rep(c(rep("lab",24), rep("field",24)),2)),
                      replicate=c(rep(c(rep("one",3), rep("two",3), rep("three",3), rep("four",3)),8)))

rownames(coldata) <- colnames(cts) # all samples names

coldata$treatment <- factor(coldata$treatment)
coldata$location <- factor(coldata$location)
coldata$strain <- factor(coldata$strain)
coldata$replicate <- factor(coldata$replicate)
```

### Running and interpreting R/DESeq

Now try DESeq package !
```
library("DESeq2")

ncol(cts_1) # = 48
nrow(coldata_1) # = 48 > check same value as above

ncol(cts)
nrow(coldata) # check same value

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~treatment+strain+treatment:strain)

dds$treatment <- relevel(dds$treatment, ref = "untreated")
dds$strain <- relevel(dds$strain, ref = "lab")
```

We performed a minimal pre-filtering to keep only rows that have at least 10 reads total before running DESeq:
```
dds <- dds[ rowSums(counts(dds)) > 10, ] # 12 438
dds <- DESeq(dds)
```

In case we needed to try an other design:
```
dds$group <- factor(paste0(dds$treatment, dds$strain))
design(dds) <- ~ group
dds <- DESeq(dds)
```







