library(recount3)
library(edgeR)
library(recount)
library(ggplot2)

# Load datasets
rse_brain <- readRDS("rse_brain.RDS")
rse_liver <- readRDS("rse_liver.RDS")
rse_pancreas <- readRDS("rse_pancreas.RDS")

# Transform raw counts
assays(rse_brain)$counts <- transform_counts(rse_brain)
assays(rse_liver)$counts <- transform_counts(rse_liver)
assays(rse_pancreas)$counts <- transform_counts(rse_pancreas)

# Compute TPM for complete datasets
assays(rse_brain)$TPM    <- recount::getTPM(rse_brain)
assays(rse_liver)$TPM    <- recount::getTPM(rse_liver)
assays(rse_pancreas)$TPM <- recount::getTPM(rse_pancreas)

# Function to build dataframe with selected vs non-selected genes
make_plot_df <- function(gene_names) {
  idx_sel <- which(rowData(rse_brain)$gene_name %in% gene_names)
  idx_rest <- setdiff(seq_len(nrow(rse_brain)), idx_sel)

  sel_mat <- list(
    Brain = assays(rse_brain)$TPM[idx_sel, ],
    Liver = assays(rse_liver)$TPM[idx_sel, ],
    Pancreas = assays(rse_pancreas)$TPM[idx_sel, ]
  )

  rest_mat <- list(
    Brain = assays(rse_brain)$TPM[idx_rest, ],
    Liver = assays(rse_liver)$TPM[idx_rest, ],
    Pancreas = assays(rse_pancreas)$TPM[idx_rest, ]
  )

  sel_means <- lapply(sel_mat, colMeans)
  rest_means <- lapply(rest_mat, colMeans)

  df_sel <- data.frame(
    Tissue = rep(names(sel_means), lengths(sel_means)),
    AvgTPM = unlist(sel_means, use.names = FALSE),
    Set = "Selected"
  )
  df_rest <- data.frame(
    Tissue = rep(names(rest_means), lengths(rest_means)),
    AvgTPM = unlist(rest_means, use.names = FALSE),
    Set = "NonSelected"
  )
  rbind(df_sel, df_rest)
}

# Example lists of DE genes
# result_genes_BoLP, result_genes_LoBP, result_genes_PoBL should
# contain the genes overexpressed in Brain, Liver, and Pancreas respectively.

# Plot for Brain overexpressed genes
plot_df_BoLP <- make_plot_df(result_genes_BoLP)
plot_df_BoLP$Tissue <- factor(plot_df_BoLP$Tissue, levels = c("Brain", "Liver", "Pancreas"))

ggplot(plot_df_BoLP, aes(x = Tissue, y = AvgTPM, fill = Set)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.8), alpha = 0.5) +
  geom_jitter(aes(color = Set), position = position_dodge(width = 0.8),
              width = 0.1, size = 2, alpha = 0.5) +
  ggtitle("Brain overexpressed genes") +
  ylab("Mean TPM") +
  theme_minimal()

# Plot for Liver overexpressed genes
plot_df_LoBP <- make_plot_df(result_genes_LoBP)
plot_df_LoBP$Tissue <- factor(plot_df_LoBP$Tissue, levels = c("Brain", "Liver", "Pancreas"))

ggplot(plot_df_LoBP, aes(x = Tissue, y = AvgTPM, fill = Set)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.8), alpha = 0.5) +
  geom_jitter(aes(color = Set), position = position_dodge(width = 0.8),
              width = 0.1, size = 2, alpha = 0.5) +
  ggtitle("Liver overexpressed genes") +
  ylab("Mean TPM") +
  theme_minimal()

# Plot for Pancreas overexpressed genes
plot_df_PoBL <- make_plot_df(result_genes_PoBL)
plot_df_PoBL$Tissue <- factor(plot_df_PoBL$Tissue, levels = c("Brain", "Liver", "Pancreas"))

ggplot(plot_df_PoBL, aes(x = Tissue, y = AvgTPM, fill = Set)) +
  geom_violin(trim = FALSE, position = position_dodge(width = 0.8), alpha = 0.5) +
  geom_jitter(aes(color = Set), position = position_dodge(width = 0.8),
              width = 0.1, size = 2, alpha = 0.5) +
  ggtitle("Pancreas overexpressed genes") +
  ylab("Mean TPM") +
  theme_minimal()

