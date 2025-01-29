**Transcriptomic-analysis**


Transcriptomic analysis of colon cancer focuses on identifying differentially expressed genes between cancerous and normal tissues. By applying differential gene expression (DGE) methods such as RNA sequencing, researchers can pinpoint genes that are upregulated or downregulated in the tumor, providing insights into the molecular mechanisms of cancer. This information helps in understanding tumor biology, uncovering potential biomarkers for diagnosis, and identifying therapeutic targets for treatment strategies. The DGE pipeline typically involves quality control, read alignment, normalization, statistical analysis, and visualization of the results to interpret gene expression patterns effectively.

i am providing my pipeline here - https://github.com/abhi9704/transcriptomic-analysis/blob/main/transcriptomic_analysis.sh

**Required tools**
fastqc : https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

fastp : https://github.com/OpenGene/fastp

hisat2 : http://daehwankimlab.github.io/hisat2/

samtools : http://www.htslib.org/doc/

featurecounts : https://rnnh.github.io/bioinfo-notebook/docs/featureCounts.html

r studio : https://rstudio-education.github.io/hopr/starting.html

DESeq2 : https://lashlock.github.io/compbio/R_presentation.html


**1) Data Collection**
For this study, transcriptomic data were obtained from publicly available repositories
such as the GEO omnibus (GEO), The Cancer Genome Atlas (TCGA) and the Sequence
Read Archive (NCBI SRA). These repositories provide high-quality sequence data
from various studies. Obtained the 2 biological replicates means matched tumor as
cancer sample and adjacent normal tissue sample as control sample to ensure accurate
comparison. Obtained Metadata to get patient information including patients
information, sequencing platforms, and experimental condition. The data downloaded
in FASTQ Format.

**2) Pre-Processing**
**2.1) Quality control of Data:**
FastQC is used to check the quality of raw reads. Tools like fastp are used to remove
reads with a Phred quality score below 30, trim adapter sequences, and filter out short
reads. This step generates cleaned FASTQ files and QC reports, which visualize data
quality metrics like GC content and base quality distribution, ensuring reliable data for
downstream analysis.

**2.2) Alignment**
Align the cleaned reads mapped to human reference genome (eg. grch38) using
HISAT2. This process identifies the genomic origins of reads, creating Sequence
Alignment Map (SAM) files. SAM file sorted and converted into BAM file using
samtools and indexed. Aligned reads enable accurate gene quantification and analysis.

**2.3) Quantification**
Align reads quantified using Featurecounts, a tool that counts aligned read overlapping
with genomic features such as genes or exons. This step generates counts matrix which
is used as a input for differential gene expression analysis.

**3) Differential Gene Expression**
Differential gene expression was done using DESeq2, an R package designed for RNA
count data. This analysis identified significantly identified unregulated and
downregulated genes in tumor samples compared to normal samples.

**4) Functional Enrichment Analysis**
DAVID was used for for Gene Ontology (GO) and Kyoto Encyclopedia genes
and genomes (KEGG) pathway enrichment analysis


