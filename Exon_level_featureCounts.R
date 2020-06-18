#!/usr/bin/Rscript

# --------------------------------------------------------------------------
# Quantify exon counts with RSubread
# --------------------------------------------------------------------------

setwd("/share/ScratchGeneral/scoyou/sarah_projects/fusion_erihey/project_results/STAR_ENCODE")

dir.create("/share/ScratchGeneral/scoyou/sarah_projects/fusion_erihey/project_results/featureCounts")

ncores <- 2

library('Rsubread')
library('rtracklayer')
library('dplyr')

# List BAM files using all STAR reads (genome mapped)
gtffile <- file.path("/share/ClusterShare/biodata/contrib/scoyou/genomes/human/GRCh38_w_ercc/gencode.v32.chr_patch_hapl_scaff.annotation_ERCC92.gtf")
gtf <- rtracklayer::import(gtffile)
gtf_df <- as.data.frame(gtf)

# Map gene names to gene ids
gene_names <- gtf_df[, c("gene_id", "gene_name")]
gene_names <- distinct(gene_names)

samples <- list.files()
filenames <- file.path(paste0(samples, "/", samples, ".Aligned.sortedByCoord.out.dupmarked.bamProcessed.out.bam"))
filenames <- filenames[! grepl("merged", filenames)]
file.exists(filenames)

# Count all reads aligning to exons
fc_all <- featureCounts(files=filenames, annot.ext=gtffile, isGTFAnnotationFile=TRUE, GTF.featureType = "exon", GTF.attrType = "gene_id", GTF.attrType.extra = "exon_id", useMetaFeatures = FALSE, strandSpecific = 2, nthreads = ncores, isPairedEnd=TRUE, requireBothEndsMapped = TRUE, allowMultiOverlap = TRUE, minOverlap = 16, countMultiMappingReads = TRUE,  ignoreDup = FALSE )
write.csv(fc_all$counts, "/share/ScratchGeneral/scoyou/sarah_projects/fusion_erihey/project_results/featureCounts/Exon_level_counts_all_reads.csv")

idx <- match(fc_all$annotation$GeneID, gene_names$gene_id)
fc_all$annotation[ , 'GeneSymbol'] <- gene_names$gene_name [idx]
write.csv(fc_all$annotation, "/share/ScratchGeneral/scoyou/sarah_projects/fusion_erihey/project_results/featureCounts/Exon_level_annotation_all_reads.csv")

write.csv(fc_all$stat, "/share/ScratchGeneral/scoyou/sarah_projects/fusion_erihey/project_results/featureCounts/Exon_level_stats_all_reads.csv")

# Count dedup reads aligning to exons
fc_dedup <- featureCounts(files=filenames, annot.ext=gtffile, isGTFAnnotationFile=TRUE, GTF.featureType = "exon", GTF.attrType = "gene_id",GTF.attrType.extra = "exon_id", useMetaFeatures = FALSE, strandSpecific = 2, nthreads = ncores,  isPairedEnd=TRUE, requireBothEndsMapped = TRUE, allowMultiOverlap = TRUE, minOverlap = 16, countMultiMappingReads = TRUE, ignoreDup = TRUE)
write.csv(fc_dedup$counts, "/share/ScratchGeneral/scoyou/sarah_projects/fusion_erihey/project_results/featureCounts/Exon_level_counts_dedup_reads.csv")

idx <- match(fc_dedup$annotation$GeneID, gene_names$gene_id)
fc_dedup$annotation[ , 'GeneSymbol'] <- gene_names$gene_name [idx]
write.csv(fc_dedup$annotation, "/share/ScratchGeneral/scoyou/sarah_projects/fusion_erihey/project_results/featureCounts/Exon_level_annotation_dedup_reads.csv")

write.csv(fc_dedup$stat, "/share/ScratchGeneral/scoyou/sarah_projects/fusion_erihey/project_results/featureCounts/Exon_level_stats_dedup_reads.csv")

# Count dedup and uniquely mapping reads aligning to exons
fc_dedup_unique <- featureCounts(files=filenames, annot.ext=gtffile, isGTFAnnotationFile=TRUE, GTF.featureType = "exon", GTF.attrType = "gene_id",GTF.attrType.extra = "exon_id", useMetaFeatures = FALSE, strandSpecific = 2, nthreads = ncores,  isPairedEnd=TRUE, requireBothEndsMapped = TRUE, allowMultiOverlap = TRUE, minOverlap = 16, countMultiMappingReads = FALSE, ignoreDup = TRUE)
write.csv(fc_dedup_unique$counts, "/share/ScratchGeneral/scoyou/sarah_projects/fusion_erihey/project_results/featureCounts/Exon_level_counts_dedup_unique_reads.csv")

idx <- match(fc_dedup_unique$annotation$GeneID, gene_names$gene_id)
fc_dedup_unique$annotation[ , 'GeneSymbol'] <- gene_names$gene_name [idx]
write.csv(fc_dedup_unique$annotation, "/share/ScratchGeneral/scoyou/sarah_projects/fusion_erihey/project_results/featureCounts/Exon_level_annotation_dedup_unique_reads.csv")

write.csv(fc_dedup_unique$stat, "/share/ScratchGeneral/scoyou/sarah_projects/fusion_erihey/project_results/featureCounts/Exon_level_stats_dedup_unique_reads.csv")




