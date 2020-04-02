args = commandArgs(trailingOnly = TRUE)

library(DESeq2)
rootDirectory = args[1]

sampleFilesPath <- list.files(path=rootDirectory, pattern="htseq.counts.gz", full.names=TRUE, recursive=TRUE)

tableSplit <- strsplit(sampleFilesPath,"/")
tableSplit <- matrix(unlist(tableSplit),ncol=13,byrow=TRUE)
htseq.path <- "/Users/bcrone/Documents/CLASS/Winter20/BIOINF545/Project/BIOINF545Project/data/htseq/"
updated.path <- paste(htseq.path, gsub(".gz","",tableSplit[,13]),sep="")

# Master List
masterTable <- data.frame(sampleName=gsub(".htseq.counts.gz","",tableSplit[,13]),fileName=gsub(".gz","",tableSplit[,13]),
                          sampleID=tableSplit[,12],subjectID=tableSplit[,11],sampleType=gsub("^.*-","",tableSplit[,12]))

n_occur <- data.frame(table(masterTable$subjectID))
multiples <- masterTable[masterTable$subjectID %in% n_occur$Var1[n_occur$Freq>1],]

# All Samples (465)
sampleTable <- data.frame(sampleName=masterTable$sampleName,fileName=masterTable$fileName,
                         sampleType=masterTable$sampleType,sampleID=masterTable$sampleID)

# Only 01A/11A sample types (448)
twoCatTable <- sampleTable[sampleTable$sampleType=="01A" | sampleTable$sampleType=="11A",]

# Only subjects with > 1 sample (125)
multiplesTable <- data.frame(sampleName=multiples$sampleName,fileName=multiples$fileName,
                             sampleType=multiples$sampleType,sampleID=multiples$sampleID)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,directory=htseq.path,design=~sampleType)
ddsHTSeq.twoCat <- DESeqDataSetFromHTSeqCount(sampleTable=twoCatTable,directory=htseq.path,design=~sampleType)
ddsHTSeq.multiple <- DESeqDataSetFromHTSeqCount(sampleTable=multiplesTable,directory=htseq.path,design=~sampleType)

ddsHTSeq <- ddsHTSeq[rowSums(counts(ddsHTSeq))>=dim(ddsHTSeq)[2],]
ddsHTSeq.twoCat <- ddsHTSeq.twoCat[rowSums(counts(ddsHTSeq.twoCat))>=dim(ddsHTSeq.twoCat)[2],]
ddsHTSeq.multiple <- ddsHTSeq.multiple[rowSums(counts(ddsHTSeq.multiple))>=dim(ddsHTSeq.multiple)[2],]

ddsHTSeq <- DESeq(ddsHTSeq)
ddsHTSeq.twoCat <- DESeq(ddsHTSeq.twoCat)
ddsHTSeq.multiple <- DESeq(ddsHTSeq.multiple)

rld <- vst(ddsHTSeq)
rld.twoCat <- vst(ddsHTSeq.twoCat)
rld.multiple <- vst(ddsHTSeq.multiple)

plotPCA(rld,intgroup=c("sampleType"))
plotPCA(rld.twoCat,intgroup=c("sampleType"))
plotPCA(rld.multiple,intgroup=c("sampleType"))

res <- results(ddsHTSeq.twoCat,contrast=c("sampleType","01A","11A"),pAdjustMethod="fdr")

query <- GDCquery(project = "TCGA-COAD", 
                  data.category = "Clinical", 
                  file.type = "xml", 
                  barcode = masterTable$subjectID)

