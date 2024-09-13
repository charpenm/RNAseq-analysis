#raw read quality check with FastQC 0.11.8
source package fastqc-0.11.8
fastqc Sample_file_name.fq.gz

#raw read trimming with TrimGalore 0.5.0
source package trim_galore-0.5.0
trim_galore --paired --retain_unpaired Sample_read_1.fq.gz Sample_read_2.fq.gz -o Path_to_outputFolder

#generate STAR genome index
source package STAR-2.5.a
STAR --runThreadsN 8 --runMode genomeGenerate --genomeDir path_to_outputFiles --genomeFastaFiles organism_genome.fa --sjdbGTFfile annotation.gtf --sjdbOverhang 149


#align reads with STAR-2.5.a
source package STAR-2.5.a
gunzip Sample_read_1_trimmed.fq.gz
gunzip Sample_read_2_trimmed.fq.gz
STAR --genomeDir Path_to_genomeIndex_folder --readFilesIn Sample_read_1_trimmed.fq Sample_read_2_trimmed.fq -outSAMtype BAM Unsorted --runThreadN 8 --outFileNamePrefix Path_to_outputFolder --twopassMode Basic

#sort alignment files and get read counts with HTseq 0.6.1
source package htseq-0.6.1
source package python-3.7.2
source samtools 1.10
samtools sort -o Sample_sorted_fileName.bam Sample_alignment_file.bam
samtools view -q 3 Sample_sorted_fileName.bam | htseq-count --mode=intersection-strict --stranded=no --type=exon --idattr=gene_id - annotation.gtf > Sample_counts_file.out

#get read counts matrix for Deseq2 by running the following perl script 

source perl-5.22.1
perl Mergecounts.pl file.log Path_to_ReadCounts_Folder MatrixFileName Path_to_output

#using R get the differentially expressed genes 

		#load libraries
		library(gplots)
		library(ggplot2)
		library(DESeq2)

		#load data
		rawCountTable <- read.table("MatrixCounts_file.txt", header=TRUE, sep="\t", row.names=1) #load matrix file with read counts
		sampleInfo <- read.table("Sample_annotation_file.txt", header=TRUE) # load the sample annotation file 
		attach(sampleInfo)
		head(genotype) #check if data is correct
		head(sampleInfo) #check if data is correct

		#generate Deseq object and normalize reads
		meta<-data.frame(row.names=colnames(rawCountTable), line=factor(genotype), treatment=factor(condition), sample_names=factor(individuals))
		dds<-DESeqDataSetFromMatrix(countData=rawCountTable,meta,formula(~treatment))
		idx <- which(rowSums(counts(dds)) > 0) #removes genes with 0 reads counts in all conditions
		dds <- dds[idx,]
		dds <- estimateSizeFactors(dds) #get library size estimate
		temp2<-counts(dds, normalized=TRUE) #normalize by read counts

		#generate DEGS
		dds$treatment <- relevel(dds$treatment, ref = "Name_of_reference_condition") #set the type of comparison (treatment or line) and the reference expression (ex. WT or mock )
		dds2 <- DESeq(dds, parallel=T) # calculates DEGs values
		res <- results(dds2, alpha = 0.05) # loads result of DEseq2
		summary(res)
		res05 <- subset(res, padj < 0.05) #return only DEG value with FDR< 0.05
		write.table(res05, file= "DEGs_filename.txt", sep="\t", quote=F, row.names = TRUE, col.names = TRUE) #saves DEGs into txt file



