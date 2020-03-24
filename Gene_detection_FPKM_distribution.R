#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Identify novel binding partners in capture data
# --------------------------------------------------------------------------

# Working directory
setwd("/Users/mac/cloudstor/sarah_projects/fusions_erihey/project_results/gene_activity/")

library('ggplot2')
library('rtracklayer')
library('biomaRt')
library('org.Hs.eg.db')
library('dplyr')

# FPKM file of total read FPKMs
RSEM_FPKM <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/project_results/RSEM/RSEM_FPKM.csv", header = TRUE, row.names = 1)
RSEM_FPKM <- RSEM_FPKM[! grepl("_", rownames(RSEM_FPKM)), ]
rownames(RSEM_FPKM) <- sub("(.*?)\\..*", "\\1", rownames(RSEM_FPKM))

# Genes in panel
blood_panel <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/data_info/Blood_panel_genes.csv", header = TRUE, row.names = 1)
solid_panel <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/data_info/Solid_panel_genes.csv", header = TRUE, row.names = 1)

# Human gene information
human_orths <- read.csv("/Users/mac/cloudstor/sarah_projects/MDA231_bulk_chrcha/annotation_data/Ensembl_GRCh38_99_annotation.csv")

# count file of total read counts
RSEM_count <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/project_results/RSEM/RSEM_count.csv", header = TRUE, row.names = 1)
RSEM_count <- RSEM_count[! grepl("_", rownames(RSEM_count)), ]
rownames(RSEM_count) <- sub("(.*?)\\..*", "\\1", rownames(RSEM_count))

# count file of unique (HQ) read counts
RSEM_count_unique <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/project_results/RSEM/RSEM_count_unique.csv", header = TRUE, row.names = 1)
RSEM_count_unique <- RSEM_count_unique[! grepl("_", rownames(RSEM_count_unique)), ]
rownames(RSEM_count_unique) <- sub("(.*?)\\..*", "\\1", rownames(RSEM_count_unique))

# Find FPKM threshold of detection for each sample
# --------------------------------------------------------------------------

rm(activity_summary)
for (i in colnames(RSEM_FPKM)) {
  
  # Subset to individual sample
  Sample_FPKM <- RSEM_FPKM[, i, drop = FALSE]
  colnames(Sample_FPKM) <- "Sample"
  
  # Remove genes with zero FPKMs and log tranform FPKMs
  Sample_FPKM_non0 <- Sample_FPKM[Sample_FPKM > 0, , drop = FALSE] 
  Sample_FPKM_log <- log2(Sample_FPKM_non0)
  
  # Determine threshold of activity
  sample_name <- i
  se_num_vec <- Sample_FPKM_log[ , "Sample"]
  c <- density(se_num_vec, bw = "nrd")
  cGT0 <- subset(c$y, c$x > 8)
  MKDE <- subset(c$x, c$y == max(cGT0)) 
  U <- mean(subset(Sample_FPKM_log[ , "Sample"], Sample_FPKM_log[ , "Sample"] > MKDE))
  sd <- (U-MKDE)*(sqrt(pi/2))
  zFPKM <- data.frame((Sample_FPKM_log[, "Sample"] - MKDE)/sd, row.names = rownames(Sample_FPKM_log))
  colnames(zFPKM) <- "Sample"
  active <- zFPKM[zFPKM$Sample >= -2, , drop = FALSE]
  
  # Identify genes with minimum zFPKM threshold
  min_active_genes <- rownames(active[active$Sample == min(active$Sample),, drop = FALSE])
  
  # Retrieve FPKM value of genes with minimum zFPKM threshold
  log2_FPKM_threshold <- unique(Sample_FPKM_log[min_active_genes,])
  FPKM_threshold <- 2^log2_FPKM_threshold
  
  if (exists("activity_summary") == TRUE) {
    activity_summary <- rbind(activity_summary, c(MKDE, sd, FPKM_threshold))
  } else {
    activity_summary <- c(MKDE, sd, FPKM_threshold)
  }
  
  # Plot distribution of total FPKMs
  ggplot(Sample_FPKM_log) +
    geom_density(aes(x = Sample), color = "black", fill = "#1f78b4", size = 0.5) +
    geom_vline(aes(xintercept = log2_FPKM_threshold), color="red", linetype="dashed", size=0.5) +
    # geom_text(aes(x=FPKM_threshold, label=FPKM_threshold, y=0.15), colour="black", angle=0) +
    labs(title = i, x = "log2(FPKM)", y = "Density") +
    theme_classic() +
    ggsave(paste0(i, "_total_FPKM_distribution.pdf"))
}

colnames(activity_summary) <- c("Mean of distribution", "SD of distribution", "FPKM_threshold")
rownames(activity_summary) <- colnames(RSEM_FPKM)
activity_summary <- data.frame(activity_summary)
write.csv(activity_summary, "Activity_threshold_summary.csv")

# Make table summarising per sample activity. A gene was considered active if it was above the activity threshold.
# --------------------------------------------------------------------------

# Find active genes
active_matrix <- RSEM_FPKM >= activity_summary$FPKM_threshold
write.csv(active_matrix, "Activity_matrix.csv")
colnames(active_matrix) <- paste0("Detected_", colnames(active_matrix))
colnames(RSEM_FPKM) <- paste0("FPKM_", colnames(RSEM_FPKM))

# Calculate the percent of total reads that are uniquely aligning
RSEM_count_unique_FC <- RSEM_count_unique/RSEM_count * 100
colnames(RSEM_count_unique_FC) <- paste0("Percent_unique_", colnames(RSEM_count_unique_FC))

# Merge columns from three matrices by samples
Active_FPKM <- cbind(active_matrix, RSEM_FPKM, RSEM_count_unique_FC) # combine
Active_FPKM <- Active_FPKM[, c(matrix(1:ncol(Active_FPKM), nrow = 3, byrow = T))] # then reorder

Activity_table <- data.frame(Active_FPKM)
Activity_table$Detected_all <- rowSums(active_matrix) == ncol(active_matrix)
Activity_table$Detected_any <- rowSums(active_matrix) > 0
Activity_table$Detected_number <- rowSums(active_matrix > 0)
Activity_table$Detected_UMI_only <- rowSums(active_matrix > 0) == rowSums(active_matrix[,c(14, 32:36)] > 0)

# Add mean expression in detected samples
active_matrix_central <- active_matrix
active_matrix_central[active_matrix_central == FALSE] <- NA
RSEM_central <- as.matrix(RSEM_FPKM) * active_matrix_central
Activity_table$Detected_median_FPKM <- rowMedians(RSEM_central[ ,c(1:13, 15:31)], na.rm=TRUE)

# Add mean percent unique in detected samples
unique_central <- as.matrix(RSEM_count_unique_FC) * active_matrix_central
Activity_table$Detected_percent_unique <- rowMedians(unique_central[ ,c(1:13, 15:31)], na.rm=TRUE)
write.csv(Activity_table, "Total_gene_detection.csv", row.names = FALSE)

# Annotate with gene information
Activity_table_annot <- merge(human_orths, Activity_table, by.x = "EnsemblID", by.y = 0)
Activity_table_annot$Solid <- is.element(Activity_table_annot$EnsemblID, solid_panel$Gene.ID)
Activity_table_annot$Blood <- is.element(Activity_table_annot$EnsemblID, blood_panel$Gene.ID)
Activity_table_annot$Both <- Activity_table_annot$Solid & Activity_table_annot$Blood
Activity_table_annot <- Activity_table_annot %>%
  mutate(
    Type = case_when(
      Solid == TRUE & Both == FALSE ~ "Solid",
      Blood == TRUE & Both == FALSE ~ "Blood",
      Both == TRUE ~ "Both",
      TRUE ~ "Other"
    )
  )

write.csv(Activity_table_annot, "Total_gene_detection_annotated.csv", row.names = FALSE)

# Test with ERCCS
# --------------------------------------------------------------------------

# Add ERCC FPKM data to help calibrate
ERCC_FPKM <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/project_results/RSEM/RSEM_ERCC_FPKM.csv", header = TRUE, row.names = 1)

# ERCC inclusion
ERCC_panel <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/data_info/ERCC_included_excluded.csv", header = TRUE, row.names = 1)
ERCC_included <- ERCC_panel[ERCC_panel$Captured == "Captured", ]
ERCC_excluded <- ERCC_panel[ERCC_panel$Captured == "Excluded", ]

rm(activity_summary)
for (i in colnames(ERCC_FPKM)) {
  
  # Subset to individual sample
  Sample_FPKM <- ERCC_FPKM[, i, drop = FALSE]
  colnames(Sample_FPKM) <- "Sample"
  
  # Remove genes with zero FPKMs and log tranform FPKMs
  Sample_FPKM_non0 <- Sample_FPKM[Sample_FPKM > 0, , drop = FALSE] 
  Sample_FPKM_log <- log2(Sample_FPKM_non0)
  
  # Determine threshold of activity
  sample_name <- i
  se_num_vec <- Sample_FPKM_log[ , "Sample"]
  c <- density(se_num_vec, bw = "nrd")
  cGT0 <- subset(c$y, c$x > 8)
  MKDE <- subset(c$x, c$y == max(cGT0)) 
  U <- mean(subset(Sample_FPKM_log[ , "Sample"], Sample_FPKM_log[ , "Sample"] > MKDE))
  sd <- (U-MKDE)*(sqrt(pi/2))
  zFPKM <- data.frame((Sample_FPKM_log[, "Sample"] - MKDE)/sd, row.names = rownames(Sample_FPKM_log))
  colnames(zFPKM) <- "Sample"
  active <- zFPKM[zFPKM$Sample >= -2, , drop = FALSE]
  
  # Identify genes with minimum zFPKM threshold
  min_active_genes <- rownames(active[active$Sample == min(active$Sample),, drop = FALSE])
  
  # Retrieve FPKM value of genes with minimum zFPKM threshold
  log2_FPKM_threshold <- unique(Sample_FPKM_log[min_active_genes,])
  FPKM_threshold <- 2^log2_FPKM_threshold
  
  if (exists("activity_summary") == TRUE) {
    activity_summary <- rbind(activity_summary, c(MKDE, sd, FPKM_threshold))
  } else {
    activity_summary <- c(MKDE, sd, FPKM_threshold)
  }
  
  # Plot distribution of total FPKMs
  ggplot(Sample_FPKM_log) +
    geom_density(aes(x = Sample), color = "black", fill = "#1f78b4", size = 0.5) +
    geom_vline(aes(xintercept = log2_FPKM_threshold), color="red", linetype="dashed", size=0.5) +
    # geom_text(aes(x=FPKM_threshold, label=FPKM_threshold, y=0.15), colour="black", angle=0) +
    labs(title = i, x = "log2(FPKM)", y = "Density") +
    theme_classic() +
    ggsave(paste0(i, "_ERCC_FPKM_distribution.pdf"))
}

colnames(activity_summary) <- c("Mean of distribution", "SD of distribution", "FPKM_threshold")
rownames(activity_summary) <- colnames(RSEM_FPKM)
activity_summary <- data.frame(activity_summary)
write.csv(activity_summary, "Activity_ERCC_threshold_summary.csv")

active_matrix <- ERCC_FPKM >= activity_summary$FPKM_threshold
colnames(active_matrix) <- paste0("Detected_", colnames(active_matrix))
colnames(ERCC_FPKM) <- paste0("FPKM_", colnames(ERCC_FPKM))
Active_FPKM <- cbind(active_matrix, ERCC_FPKM)                   # combine
Active_FPKM <- Active_FPKM[, c(matrix(1:ncol(Active_FPKM), nrow = 2, byrow = T))]   # then reorder

Activity_table <- data.frame(Active_FPKM)
Activity_table$Detected_all <- rowSums(active_matrix) == ncol(active_matrix)
Activity_table$Detected_any <- rowSums(active_matrix) > 0
Activity_table$Detected_number <- rowSums(active_matrix > 0)

# Add mean expression in detected samples
active_matrix_central <- active_matrix
active_matrix_central[active_matrix_central == FALSE] <- NA
ERCC_central <- as.matrix(ERCC_FPKM) * active_matrix_central
Activity_table$Detected_mean_FPKM <- rowMeans(ERCC_central, na.rm=TRUE)

Activity_table$Excluded <- is.element(rownames(Activity_table), rownames(ERCC_excluded))
Activity_table$Captured <- is.element(rownames(Activity_table), rownames(ERCC_included))
Activity_table <- Activity_table %>%
  mutate(
    Type = case_when(
      Captured == TRUE ~ "Captured",
      Excluded == TRUE ~ "Excluded",
      TRUE ~ "Other"
    )
  )
rownames(Activity_table) <- rownames(ERCC_FPKM)

write.csv(Activity_table, "Total_ERCC_detection.csv", row.names = TRUE)



