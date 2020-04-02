args = commandArgs(trailingOnly = TRUE)

library(DESeq2)
library(TCGAbiolinks)
rootDirectory = args[1]

sampleFilesPath <- list.files(path=rootDirectory, pattern="htseq.counts.gz", full.names=TRUE, recursive=TRUE)

tableSplit <- strsplit(sampleFilesPath,"/")
tableSplit <- matrix(unlist(tableSplit),ncol=13,byrow=TRUE)
htseq.path <- "/Users/bcrone/Documents/CLASS/Winter20/BIOINF545/Project/BIOINF545Project/data/htseq/"
xml.path <- "/Users/bcrone/Documents/CLASS/Winter20/BIOINF545/Project/BIOINF545Project/data/xml/"
updated.path <- paste(htseq.path, gsub(".gz","",tableSplit[,13]),sep="")

# Master List
masterTable <- data.frame(sampleName=gsub(".htseq.counts.gz","",tableSplit[,13]),fileName=gsub(".gz","",tableSplit[,13]),
                          sampleID=tableSplit[,12],subjectID=tableSplit[,11])

# Pull Clinical Data XML
query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  file.type = "xml", 
                  barcode = masterTable$subjectID)

GDCdownload(query,directory="/Users/bcrone/Documents/CLASS/Winter20/BIOINF545/Project/BIOINF545Project/data/xml/")

clinical <- GDCprepare_clinic(query, clinical.info = "patient")

# Remap Stage ID
tumorStage <- data.frame(subjectID=clinical$bcr_patient_barcode, stage=clinical$stage_event_pathologic_stage)

tumorStage$stage <- gsub("Stage IA","stage_1",tumorStage$stage)
tumorStage$stage <- gsub("Stage IIA|Stage IIB|Stage IIC","stage_2",tumorStage$stage)
tumorStage$stage <- gsub("Stage IIIA|Stage IIIB|Stage IIIC","stage_3",tumorStage$stage)
tumorStage$stage <- gsub("Stage IVA|Stage IVB","stage_4",tumorStage$stage)
tumorStage$stage <- gsub("Stage IV","stage_4",tumorStage$stage)
tumorStage$stage <- gsub("Stage III","stage_3",tumorStage$stage)
tumorStage$stage <- gsub("Stage II","stage_2",tumorStage$stage)
tumorStage$stage <- gsub("Stage I","stage_1",tumorStage$stage)
tumorStage$stage <- as.factor(tumorStage$stage)

# Merge with master
masterTable <- merge(masterTable,tumorStage,by="subjectID")

masterTable <- masterTable[masterTable$stage!="",]

# Build sample table for DESeq2
sampleTable <- data.frame(sampleName=masterTable$sampleName,fileName=masterTable$fileName,
                          stage=masterTable$stage,sampleID=masterTable$sampleID)
# Duplicate removal (need to revisit - where are these?)
sampleTable <- sampleTable[!duplicated(sampleTable$sampleName),]

# Build DESeq dataset
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=htseq.path,design=~stage)

ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq))>=dim(ddsHTSeq)[2],]

ddsHTSeq <- DESeq(ddsHTSeq)

# PCA
rld <- vst(ddsHTSeq)

plotPCA(rld,intgroup=c("sampleType"))
