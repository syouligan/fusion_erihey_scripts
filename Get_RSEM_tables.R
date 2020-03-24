#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Makes tables from RSEM output
# --------------------------------------------------------------------------

library(data.table)

#Setup directory structure
projectDir <-'/share/ScratchGeneral/scoyou/sarah_projects/fusion_erihey/'
biodataDir <- '/share/ScratchGeneral/biodata/contrib/scoyou/'
tool <- 'RSEM'
genome <- 'GRCh38_w_ercc'

#Practice dataset file paths
#sampleDir <- 'DGE/'
#resultsDir <- 'practice_results/'
#QCDir <- 'logs/practice/'

#Real dataset file paths
sampleDir <- 'RSEM/'
resultsDir <- 'project_results/'
QCDir <- 'logs/'

#Set input/ouput paths
samplePath <- paste(projectDir, resultsDir, sampleDir, sep="")
outPath <- paste(projectDir, resultsDir, tool, "_tables/", sep="")
dir.create(outPath)

ERCC_gtf <- paste(biodataDir, 'ercc92', '/ERCC92.gtf', sep="")

# Make tables containing RSEM expected gene counts, TPM, FPKM 
sampleFiles <- list.files(path=samplePath, pattern= "[0-9]\\.genes\\.results$", recursive=TRUE, full.names=TRUE)
sampleNames <- list.files(path=samplePath)
sampleFiles
sampleNames

RSEM_FPKM <- data.frame(fread(sampleFiles[1], header=TRUE, select=c(1)), row.names=1)
for (i in sampleFiles) {
  RSEM_FPKM <- cbind(RSEM_FPKM, fread(i, header=TRUE, select=c(7)))
}
colnames(RSEM_FPKM) <- sampleNames

RSEM_TPM <- data.frame(fread(sampleFiles[1], header=TRUE, select=c(1)), row.names=1)
for (i in sampleFiles) {
  RSEM_TPM <- cbind(RSEM_TPM, fread(i, header=TRUE, select=c(6)))
}
colnames(RSEM_TPM) <- sampleNames

RSEM_count <- data.frame(fread(sampleFiles[1], header=TRUE, select=c(1)), row.names=1)
for (i in sampleFiles) {
  RSEM_count <- cbind(RSEM_count, fread(i, header=TRUE, select=c(5)))
}
colnames(RSEM_count) <- sampleNames
RSEM_all_count <- RSEM_count

# Subset to ERCCs to different object and write to .csv files
ERCCFilter <- grep("^ERCC-", rownames(RSEM_FPKM), perl=TRUE, value=FALSE)

if (length(ERCCFilter) != 0) {
	RSEM_ERCC_FPKM <- RSEM_FPKM[ERCCFilter,]
	RSEM_ERCC_count <- RSEM_count[ERCCFilter,]
	RSEM_FPKM <- RSEM_FPKM[-ERCCFilter,]
	RSEM_TPM <- RSEM_TPM[-ERCCFilter,]
	RSEM_count <- RSEM_count[-ERCCFilter,]

	write.csv(RSEM_ERCC_FPKM, file=paste0(outPath, "RSEM_ERCC_FPKM.csv", sep=""))
	write.csv(RSEM_ERCC_count, file=paste0(outPath, "RSEM_ERCC_count.csv", sep=""))
	write.csv(RSEM_all_count, file=paste0(outPath, "RSEM_all_count.csv", sep=""))
	write.csv(RSEM_FPKM, file=paste0(outPath, "RSEM_FPKM.csv", sep=""))
	write.csv(RSEM_TPM, file=paste0(outPath, "RSEM_TPM.csv", sep=""))
	write.csv(RSEM_count, file=paste0(outPath, "RSEM_count.csv", sep=""))
} else {
	print("No ERCCs")
	write.csv(RSEM_FPKM, file=paste0(outPath, "RSEM_FPKM.csv", sep=""))
	write.csv(RSEM_TPM, file=paste0(outPath, "RSEM_TPM.csv", sep=""))
	write.csv(RSEM_count, file=paste0(outPath, "RSEM_count.csv", sep=""))
}
