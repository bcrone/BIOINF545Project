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

# Drop control samples
sampleTableNC <- sampleTable[sampleTable$sampleType != "11A",]

# Build DESeq datasets
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=htseq.path,design=~stage)
ddsHTSeqNC <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTableNC,directory=htseq.path,design=~stage)

# QC: remove genes with avg read < 1
ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq))>=dim(ddsHTSeq)[2],]
ddsHTSeqNC <- ddsHTSeqNC[rowSums(counts(ddsHTSeqNC))>=dim(ddsHTSeqNC)[2],]

ddsHTSeq <- DESeq(ddsHTSeq)
ddsHTSeqNC <- DESeq(ddsHTSeqNC)
# PCA
rld <- vst(ddsHTSeq)
rldNC <- vst(ddsHTSeqNC)

pca <- prcomp(t(assay(rld)),scale=TRUE)
pcaNC <- prcomp(t(assay(rldNC)),scale=TRUE)

percentVar <- pca$sdev^2/sum(pca$sdev^2)
percentVarNC <- pcaNC$sdev^2/sum(pcaNC$sdev^2)

intgroup <- c("stage")

intgroup.df <- as.data.frame(colData(rld)[, intgroup, drop = FALSE])
intgroup.dfNC <- as.data.frame(colData(rldNC)[, intgroup, drop = FALSE])

group <- if (length(intgroup) > 1) {
  factor(apply(intgroup.df, 1, paste, collapse = " : "))
} else {
    colData(rld)[[intgroup]]
}

groupNC <- if (length(intgroup) > 1) {
  factor(apply(intgroup.dfNC, 1, paste, collapse = " : "))
} else {
  colData(rldNC)[[intgroup]]
}

# 1v2
d12 <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                intgroup.df, name = colData(rld)[,1])

dNC12 <- data.frame(PC1 = pcaNC$x[, 1], PC2 = pcaNC$x[, 2], group = groupNC, 
                intgroup.dfNC, name = colData(rldNC)[,1])

# 1v3
d13 <- data.frame(PC1 = pca$x[, 1], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colData(rld)[,1])

dNC13 <- data.frame(PC1 = pcaNC$x[, 1], PC3 = pcaNC$x[, 3], group = groupNC, 
                    intgroup.dfNC, name = colData(rldNC)[,1])

# 1v4
d14 <- data.frame(PC1 = pca$x[, 1], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, name = colData(rld)[,1])

dNC14 <- data.frame(PC1 = pcaNC$x[, 1], PC4 = pcaNC$x[, 4], group = groupNC, 
                    intgroup.dfNC, name = colData(rldNC)[,1])

# 2v3 
d23 <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colData(rld)[,1])

dNC23 <- data.frame(PC2 = pcaNC$x[, 2], PC3 = pcaNC$x[, 3], group = groupNC, 
                    intgroup.dfNC, name = colData(rldNC)[,1])

# 2v4
d24 <- data.frame(PC2 = pca$x[, 2], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, name = colData(rld)[,1])

dNC24 <- data.frame(PC2 = pcaNC$x[, 2], PC4 = pcaNC$x[, 4], group = groupNC, 
                    intgroup.dfNC, name = colData(rldNC)[,1])

# 3v4
d34 <- data.frame(PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, name = colData(rld)[,1])

dNC34 <- data.frame(PC3 = pcaNC$x[, 3], PC4 = pcaNC$x[, 4], group = groupNC, 
                    intgroup.dfNC, name = colData(rldNC)[,1])


# 1v2
ggplot(data = d12, aes_string(x = "PC1", y = "PC2", color = "group", label = "name")) + geom_point(size = 3) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + 
  ggtitle("PC1 vs PC2 - Stage (With Controls)")

ggplot(data = dNC12, aes_string(x = "PC1", y = "PC2", color = "group", label = "name")) + geom_point(size = 3) + 
  xlab(paste0("PC1: ", round(percentVarNC[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(percentVarNC[2] * 100), "% variance")) + coord_fixed() +
  ggtitle("PC1 vs PC2 - Stage (No Controls)")

# 1v3
ggplot(data = d13, aes_string(x = "PC1", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + coord_fixed() + 
  ggtitle("PC1 vs PC3 - Stage (With Controls)")

ggplot(data = dNC13, aes_string(x = "PC1", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + 
  xlab(paste0("PC1: ", round(percentVarNC[1] * 100), "% variance")) +
  ylab(paste0("PC3: ", round(percentVarNC[3] * 100), "% variance")) + coord_fixed() +
  ggtitle("PC1 vs PC3 - Stage (No Controls)")

# 1v4
ggplot(data = d14, aes_string(x = "PC1", y = "PC4", color = "group", label = "name")) + geom_point(size = 3) + 
  xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
  ylab(paste0("PC4: ", round(percentVar[4] * 100), "% variance")) + coord_fixed() + 
  ggtitle("PC1 vs PC4 - Stage (With Controls)")

ggplot(data = dNC14, aes_string(x = "PC1", y = "PC4", color = "group", label = "name")) + geom_point(size = 3) + 
  xlab(paste0("PC1: ", round(percentVarNC[1] * 100), "% variance")) +
  ylab(paste0("PC4: ", round(percentVarNC[4] * 100), "% variance")) + coord_fixed() +
  ggtitle("PC1 vs PC4 - Stage (No Controls)")

# 2v3
ggplot(data = d23, aes_string(x = "PC2", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + 
  xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + coord_fixed() + 
  ggtitle("PC2 vs PC3 - Stage (With Controls)")

ggplot(data = dNC23, aes_string(x = "PC2", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + 
  xlab(paste0("PC2: ", round(percentVarNC[2] * 100), "% variance")) +
  ylab(paste0("PC3: ", round(percentVarNC[3] * 100), "% variance")) + coord_fixed() +
  ggtitle("PC2 vs PC3 - Stage (No Controls)")

# 2v4
ggplot(data = d24, aes_string(x = "PC2", y = "PC4", color = "group", label = "name")) + geom_point(size = 3) + 
  xlab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) +
  ylab(paste0("PC4: ", round(percentVar[4] * 100), "% variance")) + coord_fixed() + 
  ggtitle("PC2 vs PC4 - Stage (With Controls)")

ggplot(data = dNC24, aes_string(x = "PC2", y = "PC4", color = "group", label = "name")) + geom_point(size = 3) + 
  xlab(paste0("PC2: ", round(percentVarNC[2] * 100), "% variance")) +
  ylab(paste0("PC4: ", round(percentVarNC[4] * 100), "% variance")) + coord_fixed() +
  ggtitle("PC2 vs PC4 - Stage (No Controls)")

# 3v4
ggplot(data = d34, aes_string(x = "PC3", y = "PC4", color = "group", label = "name")) + geom_point(size = 3) + 
  xlab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) +
  ylab(paste0("PC4: ", round(percentVar[4] * 100), "% variance")) + coord_fixed() +
  ggtitle("PC3 vs PC4 - Stage (With Controls)")

ggplot(data = dNC34, aes_string(x = "PC3", y = "PC4", color = "group", label = "name")) + geom_point(size = 3) + 
  xlab(paste0("PC3: ", round(percentVarNC[3] * 100), "% variance")) +
  ylab(paste0("PC4: ", round(percentVarNC[4] * 100), "% variance")) + coord_fixed() +
  ggtitle("PC3 vs PC4 - Stage (No Controls)")
