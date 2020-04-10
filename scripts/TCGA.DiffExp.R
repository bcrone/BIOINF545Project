library(DESeq2)
library(TCGAbiolinks)
library(ggplot2)
library(ggrepel)

sampleFilesPath <- list.files(pattern="htseq.counts.gz", full.names=TRUE, recursive=TRUE)

tableSplit <- strsplit(sampleFilesPath,"/")
tableSplit <- matrix(unlist(tableSplit),ncol=5,byrow=TRUE)
htseq.path <- "data/htseq/"
xml.path <- "data/xml/"
updated.path <- paste(htseq.path, gsub(".gz","",tableSplit[,5]),sep="")

# Master List
masterTable <- data.frame(sampleName=gsub(".htseq.counts.gz","",tableSplit[,5]),fileName=gsub(".gz","",tableSplit[,5]),
                          sampleID=tableSplit[,4],subjectID=tableSplit[,3],sampleType=gsub("^.*-","",tableSplit[,4]))

# Pull Clinical Data XML
queryCOAD <- GDCquery(project = "TCGA-COAD", 
                      data.category = "Clinical", 
                      file.type = "xml", 
                      barcode = masterTable$subjectID)

queryREAD <- GDCquery(project = "TCGA-READ",
                      data.category = "Clinical",
                      file.type = "xml",
                      barcode = masterTable$subjectID)

GDCdownload(queryCOAD,directory=xml.path)

GDCdownload(queryREAD,directory=xml.path)

# Gather Clinical Data
clinicalCOAD <- GDCprepare_clinic(queryCOAD,directory=xml.path,
                                  clinical.info = "patient")

clinicalREAD <- GDCprepare_clinic(queryREAD,directory=xml.path,
                                  clinical.info = "patient")

# Remap Stage ID
tumorStageCOAD <- data.frame(subjectID=clinicalCOAD$bcr_patient_barcode, stage=clinicalCOAD$stage_event_pathologic_stage)
tumorStageCOAD <- tumorStageCOAD[!duplicated(tumorStageCOAD),]

tumorStageREAD<- data.frame(subjectID=clinicalREAD$bcr_patient_barcode, stage=clinicalREAD$stage_event_pathologic_stage)
tumorStageREAD <- tumorStageREAD[!duplicated(tumorStageREAD),]

tumorStageCOAD$stage <- gsub("Stage IA","Stage 1",tumorStageCOAD$stage)
tumorStageCOAD$stage <- gsub("Stage IIA|Stage IIB|Stage IIC","Stage 2",tumorStageCOAD$stage)
tumorStageCOAD$stage <- gsub("Stage IIIA|Stage IIIB|Stage IIIC","Stage 3",tumorStageCOAD$stage)
tumorStageCOAD$stage <- gsub("Stage IVA|Stage IVB","Stage 4",tumorStageCOAD$stage)
tumorStageCOAD$stage <- gsub("Stage IV","Stage 4",tumorStageCOAD$stage)
tumorStageCOAD$stage <- gsub("Stage III","Stage 3",tumorStageCOAD$stage)
tumorStageCOAD$stage <- gsub("Stage II","Stage 2",tumorStageCOAD$stage)
tumorStageCOAD$stage <- gsub("Stage I","Stage 1",tumorStageCOAD$stage)
tumorStageCOAD$stage <- as.factor(tumorStageCOAD$stage)

tumorStageREAD$stage <- gsub("Stage IA","Stage 1",tumorStageREAD$stage)
tumorStageREAD$stage <- gsub("Stage IIA|Stage IIB|Stage IIC","Stage 2",tumorStageREAD$stage)
tumorStageREAD$stage <- gsub("Stage IIIA|Stage IIIB|Stage IIIC","Stage 3",tumorStageREAD$stage)
tumorStageREAD$stage <- gsub("Stage IVA|Stage IVB","Stage 4",tumorStageREAD$stage)
tumorStageREAD$stage <- gsub("Stage IV","Stage 4",tumorStageREAD$stage)
tumorStageREAD$stage <- gsub("Stage III","Stage 3",tumorStageREAD$stage)
tumorStageREAD$stage <- gsub("Stage II","Stage 2",tumorStageREAD$stage)
tumorStageREAD$stage <- gsub("Stage I","Stage 1",tumorStageREAD$stage)
tumorStageREAD$stage <- as.factor(tumorStageREAD$stage)

tumorStage <- rbind(tumorStageCOAD,tumorStageREAD)

# Merge with master
masterTable <- merge(masterTable,tumorStage,by="subjectID",all.x=TRUE)

# Remove NA/no stage information (465 -> 441 samples)
masterTable <- masterTable[masterTable$stage!="",]
masterTable <- masterTable[!is.na(masterTable$stage),]

# Build sample table for DESeq2
sampleTable <- data.frame(sampleName=masterTable$sampleName,fileName=masterTable$fileName,
                          stage=masterTable$stage,sampleType=masterTable$sampleType,sampleID=masterTable$sampleID)

# Isolate just 01A (375 samples)
sampleTable <- sampleTable[sampleTable$sampleType == "01A",]

# Build DESeq datasets
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=htseq.path,design=~stage)

# QC: remove genes with avg read < 1
ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq))>=dim(ddsHTSeq)[2],]

# Note: This is the object needed for differential expression
ddsHTSeq <- DESeq(ddsHTSeq)
