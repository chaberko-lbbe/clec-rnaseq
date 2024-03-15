![Image 1](bedbugs.png)

You can learn more about the context of this repository in the article from which it stems: LINK TO COME

Contact : chloe.haberkorn@zoologi.su.se / chloehbk@gmail.com

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

Raw sequences (fastq.gz files) are available on SRA under BioProject PRJNA832557.
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

## Differential Expression analysis

Our goal was to performed a differential expression analyzes between London Lab (LL) and London Field (LF) strains, whether they survived an insecticide exposure or not.
We used as input "ReadsPerGene.out.tab" files 

We needed for DESeq2 to build a "cts" matrix with gene count by gene for each sample, together with a "coldata" table with infos on data :
- untreated (T) or treated (E as "exposed")
- strain origin, correlated with insecticide resistance status (LL for London Lab, LF for London Field)

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
               Sample_7E3LF[,4],Sample_8E3LF[,4],Sample_9E3LF[,4],Sample_10E4LF[,4],Sample_11E4LF[,4],Sample_12E4LF[,4])

colnames(cts) <- c("GeneID","untreatedlab1","untreatedlab2","untreatedlab3","untreatedlab4","untreatedlab5","untreatedlab6",
                     "untreatedlab7","untreatedlab8","untreatedlab9","untreatedlab10","untreatedlab11","untreatedlab12",
                     "treatedlab1","treatedlab2","treatedlab3","treatedlab4","treatedlab5","treatedlab6",
                     "treatedlab7","treatedlab8","treatedlab9","treatedlab10","treatedlab11","treatedlab12",
                     "untreatedfield1","untreatedfield2","untreatedfield3","untreatedfield4","untreatedfield5","untreatedfield6",
                     "untreatedfield7","untreatedfield8","untreatedfield9","untreatedfield10","untreatedfield11","untreatedfield12",
                     "treatedfield1","treatedfield2","treatedfield3","treatedfield4","treatedfield5","treatedfield6",
                     "treatedfield7","treatedfield8","treatedfield9","treatedfield10","treatedfield11","treatedfield12")

rownames(cts) <- cts[,1]
cts <- cts[,-1]
cts <- as.matrix(cts)
```

Build a "coldata" table with infos on data:

```
coldata <- data.frame(treatment=c(rep("untreated",12), rep("treated",12),
                                    rep("untreated",12), rep("treated",12)),
                        strain=c(rep("lab",24), rep("field",24)))

rownames(coldata) <- colnames(cts) # all samples names

coldata$treatment <- factor(coldata$treatment)
coldata$strain <- factor(coldata$strain)
```

### Running R/DESeq

Now try DESeq package !
```
library("DESeq2")

ncol(cts_1) # = 48
nrow(coldata_1) # = 48 > check same value as above

ncol(cts)
nrow(coldata) # check same value

additive.model <- as.formula(~ treatment + strain)
interaction.model <- as.formula(~ treatment + strain + treatment:strain)

library("DESeq2")
dds_additive <- DESeqDataSetFromMatrix(countData = cts,
                                       colData = coldata,
                                       design = additive.model)

dds_interaction <- DESeqDataSetFromMatrix(countData = cts,
                                          colData = coldata,
                                          design = interaction.model)
```

R will choose a reference level for factors based on alphabetical order, so we need to tell the DESeq2 functions which level we want to compare against (e.g. which level represents the control group)
```
dds_additive$treatment <- relevel(dds_additive$treatment, ref = "untreated")
dds_additive$strain <- relevel(dds_additive$strain, ref = "lab")

dds_interaction$treatment <- relevel(dds_interaction$treatment, ref = "untreated")
dds_interaction$strain <- relevel(dds_interaction$strain, ref = "lab")
```

We performed a minimal pre-filtering to keep only rows that have at least 10 reads total before running DESeq:
```
dds_interaction <- dds_interaction[ rowSums(counts(dds_interaction)) > 10, ] # Remove lines with no count or a single count : 12,706
dds_additive <- dds_additive[ rowSums(counts(dds_additive)) > 10, ] # Remove lines with no count or a single count : 12,706

dds_additive <- DESeq(dds_additive)
dds_interaction <- DESeq(dds_interaction)
```

For how many genes is interaction model a better fit?

```
dds_LRT <- DESeq(dds_interaction, test="LRT", reduced=additive.model)
results.Interaction_v_Additive <- results(dds_LRT)
table(results.Interaction_v_Additive$padj < 0.05) # only 1 !
results.Interaction_v_Additive <- as.data.frame(results.Interaction_v_Additive)
results.Interaction_v_Additive[results.Interaction_v_Additive$padj<0.05,] # same gene : LOC106661213
```

### Samples overview

```
vsd <- vst(dds_interaction, blind=FALSE) # estimate dispersion trend and apply a variance stabilizing transformation
vsd <- vst(dds_additive, blind=FALSE) # estimate dispersion trend and apply a variance stabilizing transformation

pcaData <- plotPCA(vsd, intgroup=c("treatment","strain"), returnData=TRUE) # by default : ntop=500
percentVar <- round(100 * attr(pcaData, "percentVar"))

library(scales)
#extract hex color codes for a plot with three elements in ggplot2 
hex <- hue_pal()(4)

pcaData.san <- plotPCA.san(vsd, intgroup=c("treatment","strain"), returnData=TRUE)
percentVar.san <- round(100 * attr(pcaData.san, "percentVar"))

library(ggplot2)
ggplot(pcaData, aes(PC1, PC2, color=strain, shape=treatment)) + # pcaData.san
  geom_point(size=3) +
  stat_ellipse(geom = "polygon",
               aes(fill = strain,linetype=treatment), 
               alpha = 0.25)+
  scale_color_manual(values = c("lab" = "#9E3124",
                                "field" = "#779E61"))+
  ggtitle("PCA C. lectularius depending on strains and insecticide treatment") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) + # PC2, percentVar.san
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + # PC3, percentVar.san
  theme_bw()+
  coord_fixed()
```
### Detection of DE genes

```
resultsNames(dds_interaction) # lists the coefficients
resultsNames(dds_additive) # lists the coefficients

res_strain <- results(dds_interaction, name="strain_field_vs_lab") 

sum(res_strain$padj < 0.05 & res_strain$log2FoldChange>0.5849625, na.rm=TRUE) # 20 
sum(res_strain$padj < 0.05 & res_strain$log2FoldChange<(-0.5849625), na.rm=TRUE) # 4 
# log2(1.5) = 0.5849625

summary(res_strain,  alpha = 0.05) # by default adjusted p-value < 0.05
summary(res_strain,  alpha = 1) # no adj
```

Then, we can plot the results:
```
library(EnhancedVolcano)
EnhancedVolcano(res_strain,
                lab = NA, # You can also alternatively use rownames(res_strain)
                x = 'log2FoldChange', y = 'padj',
                pCutoff = 0.05, FCcutoff = 0.5849625,
                xlim = c(-3.5,6), ylim = c(0,5.5),
                pointSize = 2, labSize = 1)

plotCounts(dds, gene="LOC106663998", intgroup=c("strain","treatment"), normalized = T) # Example for a single gene
```








