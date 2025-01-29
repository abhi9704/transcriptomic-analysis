# ...........................Methods and methodology-......................................

# 1. Obtain RNA-seq data from colon cancer tissues and normal colon tissues. 
# You can use publicly available datasets (eg. TCGA, NCBI SRA).
    
# Download the data using sra-toolkit (using terminal)
# Download the data = 
prefetch prefetch SRR11190598 prefetch SRR11190599  SRR11190635 SRR11190636

# Convert sra file to fastq =  
fastq-dump --split file SRR11190635.sra
# Two output fastq files(R1,R2) = SRR11190635_1.fastq and SRR11190635_2.fastq

# 2. Quality control (QC) of raw RNA-seq data.
# FastQC - check quality of raw data.
FastQC SRR11190635_1.fastq

# Fastp - trimming and filtering raw data.
fastp -i SRR11190635 -o  control_1.fastq  -I SRR11190635_2.fastq -O control_2.fastq 

# 3. Align the reads to the human reference genome using a tool  HISAT2.
 
# Download reference genome using ENSEMBLE. (eg. Human GRCh38)
# Indexing of reference genome  
hisat2-build  reference_genome.fa  genome_index
        
# Alignement of sample to reference genome
hisat2 -x grch38_genome/grch38/genome  -1 control_1.fastq -2 control_2.fastq -S control_mapped.sam

# Convert SAM to BAM 
samtools view -S -b control_mapped.sam > control_mapped.bam 

# Sort and index the BAM file
samtools sort control_mapped.bam -o control_sorted.bam 

samtools index control_sorted.bam 

#Check mapping statistics 
samtools flagstat control_sorted.bam

# 4. Quantification of gene expression levels using featureCounts.
featureCounts -T 8 -p -t exon -g gene_id -a Homo_sapiens.GRCh38.112.gtf -o gene_counts.txt control_sorted.bam cancer_sorted.bam

# differential gene expression using DESeq2
# open r studio

# load the required libraries
library(DESeq2)
library(ggplot2)

# load countdata}
counts <- read.csv("/home/abhi/sw48.csv", row.names = 1, header =  TRUE, sep = ",")
counts

counts <- counts[which(rowSums(counts) > 10),]

# create metadata
condition <- factor(c("cancer","cancer","treatment","treatment"))
colData <- data.frame(row.names = colnames(counts),condition)
colData

# validate
all(rownames(colData) == colnames(counts))

colData $condition <- factor(colData $condition)
table(colData$condition)

# run DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts , colData = colData, design = ~condition)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

res <- results(dds)
head(results(dds, tidy=TRUE))

summary(res)
head(res)


