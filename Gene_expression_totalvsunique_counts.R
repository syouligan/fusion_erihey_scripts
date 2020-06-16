#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Identify novel binding partners in capture data
# --------------------------------------------------------------------------

# Working directory
setwd("/Users/mac/cloudstor/sarah_projects/fusions_erihey/project_results/novel_partners/")

library('ggplot2')
library('rtracklayer')
library('biomaRt')
library('org.Hs.eg.db')
library('dplyr')


# count file of total read counts
RSEM_count <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/project_results/RSEM/RSEM_count.csv", header = TRUE, row.names = 1)
RSEM_count <- RSEM_count[! grepl("_", rownames(RSEM_count)), ]
rownames(RSEM_count) <- sub("(.*?)\\..*", "\\1", rownames(RSEM_count))

# count file of unique (HQ) read counts
RSEM_count_unique <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/project_results/RSEM/RSEM_count_unique.csv", header = TRUE, row.names = 1)
RSEM_count_unique <- RSEM_count_unique[! grepl("_", rownames(RSEM_count_unique)), ]
rownames(RSEM_count_unique) <- sub("(.*?)\\..*", "\\1", rownames(RSEM_count_unique))

# Genes in panel
blood_panel <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/data_info/Blood_panel_genes.csv", header = TRUE, row.names = 1)
solid_panel <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/data_info/Solid_panel_genes.csv", header = TRUE, row.names = 1)

# Plot count distribution for each sample
# --------------------------------------------------------------------------

for (i in colnames(RSEM_count)) {
  # Create folder for sample
  dir.create(i)
  
  # Subset to individual sample
  Sample_count <- RSEM_count[, i, drop = FALSE]
  colnames(Sample_count) <- "Sample"
  
  # Remove genes with zero counts and log tranform counts
  Sample_count_non0 <- Sample_count[Sample_count > 0, , drop = FALSE] 
  Sample_count_log <- log2(Sample_count_non0)
  
  # Plot distribution of total counts
  ggplot(Sample_count_log) +
    geom_density(aes(x = Sample), color = "black", fill = "#1f78b4", size = 0.5) +
    # geom_vline(aes(xintercept=count_threshold), color="red", linetype="dashed", size=0.5) +
    # geom_text(aes(x=count_threshold, label=count_threshold, y=0.15), colour="black", angle=0) +
    labs(title = i, x = "log2(count)", y = "Density") +
    theme_classic() +
    ggsave(paste0(i, "/", i, "_total_count_distribution.pdf"))
  
  # Subset to individual sample
  Sample_count <- RSEM_count_unique[, i, drop = FALSE]
  colnames(Sample_count) <- "Sample"
  
  # Remove genes with zero counts and log tranform counts
  Sample_count_non0 <- Sample_count[Sample_count > 0, , drop = FALSE] 
  Sample_count_log <- log2(Sample_count_non0)

  # Plot distribution of unique counts
  ggplot(Sample_count_log) +
    geom_density(aes(x = Sample), color = "black", fill = "#1f78b4", size = 0.5) +
    # geom_vline(aes(xintercept=count_threshold), color="red", linetype="dashed", size=0.5) +
    # geom_text(aes(x=count_threshold, label=count_threshold, y=0.15), colour="black", angle=0) +
    labs(title = i, x = "log2(count)", y = "Density") +
    theme_classic() +
    ggsave(paste0(i, "/", i, "_unique_count_distribution.pdf"))
}

# Plot LFC betwwen total and unique mapping reads
# --------------------------------------------------------------------------

RSEM_count_unique_FC <- RSEM_count_unique/RSEM_count

for (i in colnames(RSEM_count_unique_FC)) {
  # Subset to individual sample
  Sample_count <- RSEM_count_unique_FC[, i, drop = FALSE]
  colnames(Sample_count) <- "Sample"
  
  # Remove genes with zero counts and log tranform counts
  Sample_count_non0 <- Sample_count[Sample_count > 0 & is.na(Sample_count) == FALSE, , drop = FALSE] 
  Sample_count_log <- log2(Sample_count_non0)
  
  # Plot distribution of log fold change between unique counts and total counts
  ggplot(Sample_count_log) +
    geom_density(aes(x = Sample), color = "black", fill = "#1f78b4", size = 0.5) +
    # geom_vline(aes(xintercept=count_threshold), color="red", linetype="dashed", size=0.5) +
    # geom_text(aes(x=count_threshold, label=count_threshold, y=0.15), colour="black", angle=0) +
    labs(title = i, x = "log2(fold-change)", y = "Density") +
    theme_classic() +
    ggsave(paste0(i, "/", i, "_LFC_distribution.pdf"))
}

# Plot count distribution vs LFC between total and unique mapping reads
# --------------------------------------------------------------------------

# Add ERCC data to help calibrate
ERCC_count <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/project_results/RSEM/RSEM_ERCC_count.csv", header = TRUE, row.names = 1)
RSEM_count <- rbind(RSEM_count, ERCC_count)
ERCC_count_unique <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/project_results/RSEM/RSEM_ERCC_count_unique.csv", header = TRUE, row.names = 1)
RSEM_count_unique <- rbind(RSEM_count_unique, ERCC_count_unique)
RSEM_count_unique_FC <- RSEM_count_unique/RSEM_count * 100

# ERCC inclusion
ERCC_panel <- read.csv("/Users/mac/cloudstor/sarah_projects/fusions_erihey/data_info/ERCC_included_excluded.csv", header = TRUE, row.names = 1)
ERCC_included <- ERCC_panel[ERCC_panel$Captured == "Captured", ]
ERCC_excluded <- ERCC_panel[ERCC_panel$Captured == "Excluded", ]

for (i in colnames(RSEM_count)) {
  # Subset to individual sample
  Sample_count <- RSEM_count[, i, drop = FALSE]
  colnames(Sample_count) <- "Total_count"
  Sample_count$Unique_count <- RSEM_count_unique[, i]
  
  # Remove genes with zero counts and log tranform counts
  Sample_count_non0 <- Sample_count[Sample_count$Total_count > 0 & Sample_count$Unique_count > 0, , drop = FALSE] 
  Sample_count_log <- log2(Sample_count_non0)
  
  # Add % of unique
  idx <- match(rownames(Sample_count_log), rownames(RSEM_count_unique_FC))
  Sample_count_log$Percent_unique <- RSEM_count_unique_FC[, i] [idx]
  Sample_count_log <- Sample_count_log[!(is.na(Sample_count_log$Percent_unique) | Sample_count_log$Percent_unique == 0 | Sample_count_log$Percent_unique == Inf), ]
  
  # Label with blood or solid panel
  Sample_count_log$Solid <- is.element(rownames(Sample_count_log), solid_panel$Gene.ID)
  Sample_count_log$Blood <- is.element(rownames(Sample_count_log), blood_panel$Gene.ID)
  Sample_count_log$Both <- Sample_count_log$Solid & Sample_count_log$Blood
  Sample_count_log$Excluded <- is.element(rownames(Sample_count_log), rownames(ERCC_excluded))
  Sample_count_log$Captured <- is.element(rownames(Sample_count_log), rownames(ERCC_included))
  
  # Label with gene type
  Sample_count_log <- Sample_count_log %>%
    mutate(
      Type = case_when(
        Solid == TRUE & Both == FALSE ~ "Solid",
        Blood == TRUE & Both == FALSE ~ "Blood",
        Both == TRUE ~ "Both",
        Captured == TRUE ~ "Captured",
        Excluded == TRUE ~ "Excluded",
        TRUE ~ "Other"
      )
    )
  rownames(Sample_count_log) <- rownames(Sample_count_non0)
  Sample_count_log$Type <- factor(Sample_count_log$Type, levels = c("Solid", "Both", "Blood", "Captured", "Excluded", "Other"))
  
  # Plot distribution of total counts relative to the percent unique
  ggplot(Sample_count_log) +
    geom_point(aes(x = Total_count, y = Percent_unique, color = Type), size = 1) +
    # geom_vline(aes(xintercept=count_threshold), color="red", linetype="dashed", size=0.5) +
    labs(title = i, x = "Total_count", y = "% Unique") +
    theme_classic() +
    ggsave(paste0(i, "/", i, "_Total_count_percent_unique_distribution.pdf"))

  # Plot distribution of total counts relative to the percent unique
  ggplot(Sample_count_log) +
    geom_point(aes(x = Unique_count, y = Total_count, color = Type), size = 1) +
    # geom_vline(aes(xintercept = 1.5), color = "red", linetype = "dashed", size=0.5) +
    labs(title = i, x = "Unique_count", y = "Total_count") +
    theme_classic() +
    ggsave(paste0(i, "/", i, "_Unique_Total_count_distribution.pdf"))
}

