# Henry Lock Final Project 2 - Multi-Tissue Analysis

# These all load in the libraries that are used later in the code for cleaning,
# modeling, and plotting
library(ggplot2)
library(readxl)
library(dplyr)
library(tidyr)
library(broom)
library(mclust)


# This reads all the metadata, both phenotype and attributes from GTEx portal 
gtex_attr <- read_xlsx("/Users/henrylock/Library/Mobile Documents/com~apple~CloudDocs/R Stuff/Final Project/GTEx_Analysis_v10_Annotations_SampleAttributesDD.xlsx")

url <- "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
destfile <- "/Users/henrylock/Library/Mobile Documents/com~apple~CloudDocs/R Stuff/Final Project/GTEx_SampleAttributesDS.txt"
download.file(url, destfile, method = "curl")  
gtex_attr_metadata <- read.delim(destfile, sep = "\t", header = TRUE)

gtex_phen <- read_xlsx("/Users/henrylock/Library/Mobile Documents/com~apple~CloudDocs/R Stuff/Final Project/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDD.xlsx")

url <- "https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
destfile <- "/Users/henrylock/Library/Mobile Documents/com~apple~CloudDocs/R Stuff/Final Project/GTEx_SubjectPhenotypesDS.txt"
download.file(url, destfile, method = "curl")  
gtex_phen_metadata <- read.delim(destfile, sep = "\t", header = TRUE)


#Functions

# Creates a function to read each gct file from GTEx
read_gct <- function(file) {
  lines <- readLines(file)
  version <- lines[1]
  dims <- strsplit(lines[2], "\t")[[1]]
  df <- read.delim(file, skip = 2, header = TRUE, stringsAsFactors = FALSE)
  expr <- as.matrix(df[, -(1:2)])
  rownames(expr) <- df$Name
  list(version = version, data = df, expression = expr)
}

# Function to get the top 10 most variable genes from each tissue type, but also
# filters for ones above threshold.
get_top_variable_genes <- function(expr_matrix, sex_labels, top_n = 10, threshold = 20, min_samples = 20) {
  expr_filtered <- expr_matrix[rowSums(expr_matrix > threshold) >= min_samples, ]
  
  gene_vars <- apply(expr_filtered, 1, var)
  
  top_genes <- sort(gene_vars, decreasing = TRUE)[1:top_n]
  
  return(top_genes)
}

# Function to create regression plot for each gene based on age x expression level,
# with a line for each sex per gene.
plot_gene_regressions <- function(expr_long_df, tissue_name = "") {
  ggplot(expr_long_df, aes(x = AGE_numeric, y = TPM, color = as.factor(SEX))) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "lm", se = FALSE) +
    facet_wrap(~ Gene, scales = "free_y") +
    labs(
      title = paste("Gene Expression vs Age by Sex -", tissue_name),
      x = "Age (Midpoint of Bin)",
      y = "TPM (Transcripts Per Million)",
      color = "Sex"
    ) +
    scale_color_manual(
      values = c("1" = "red", "2" = "blue"),
      labels = c("1" = "Male", "2" = "Female")
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Makes the bar graph for the top ten variable genes in each tissue type.
plot_top_gene_variance <- function(top_genes, tissue_name = "") {
  df <- data.frame(Gene = names(top_genes), Variance = top_genes)
  
  ggplot(df, aes(x = reorder(Gene, -Variance), y = Variance)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(
      title = paste("Top Variable Genes -", tissue_name),
      x = "Gene (Ensembl ID)",
      y = "Variance of TPM"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Function to make and fit the linear model to test if age x sex are correlated
# for any of the top variable genes.
analyze_gene_expression <- function(expr_matrix, top_genes, sample_metadata) {
  genes_to_plot <- names(top_genes)
  expr_long <- as.data.frame(t(expr_matrix[genes_to_plot, , drop = FALSE]), check.names = FALSE)
  expr_long$SAMPID <- rownames(expr_long)
  expr_merged <- left_join(expr_long, sample_metadata, by = "SAMPID")
  expr_long_df <- expr_merged |> pivot_longer(cols = all_of(genes_to_plot), names_to = "Gene", values_to = "TPM") |>
    filter(!is.na(AGE), !is.na(SEX), !is.na(TPM))
  expr_long_df$AGE <- factor(expr_long_df$AGE, levels = c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79"), ordered = TRUE)
  expr_long_df$AGE_numeric <- recode(expr_long_df$AGE, "20-29" = 25, "30-39" = 35, "40-49" = 45, "50-59" = 55, "60-69" = 65, "70-79" = 75)
  model_results <- expr_long_df |> group_by(Gene) |>
    do({ model <- lm(TPM ~ AGE_numeric * SEX, data = .); tidy(model) }) |>
    ungroup()
  interaction_effects <- model_results |>
    filter(term == "AGE_numeric:SEX") |>
    mutate(signif = case_when(is.na(p.value) ~ "Missing", p.value < 0.05 ~ "Significant", TRUE ~ "Not significant"))
  list(expression_data = expr_long_df, interaction_results = interaction_effects)
}

# Creates a volcano plot that shows if the difference between the regression slopes
# are significant or not.
plot_interaction_volcano <- function(interaction_effects, tissue_name = "") {
  ggplot(interaction_effects, aes(x = estimate, y = -log10(p.value))) +
    geom_point(aes(color = signif)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = c("Significant" = "red", "Not significant" = "gray", "Missing" = "black"), name = "Significance") +
    labs(
      title = paste("Volcano Plot: Age × Sex Interaction -", tissue_name),
      x = "Interaction Effect Size (Slope Difference)",
      y = "-log10(P-value)"
    ) +
    theme_minimal()
}

# Does the logistic regression classification model to predict the sex of each
# tissue sample without the metadata label
evaluate_logistic_model <- function(expr_matrix, top_genes, metadata) {
  selected_genes <- names(top_genes)
  expr_df <- as.data.frame(t(expr_matrix[selected_genes, , drop = FALSE]))
  expr_df$SAMPID <- rownames(expr_df)
  logit_data <- left_join(expr_df, metadata, by = "SAMPID") |> filter(!is.na(SEX))
  logit_data$SEX <- factor(ifelse(logit_data$SEX == 1, "Male", "Female"))
  formula <- as.formula(paste("SEX ~", paste(selected_genes, collapse = " + ")))
  model <- glm(formula, data = logit_data, family = binomial)
  logit_data$predicted_prob <- predict(model, type = "response")
  logit_data$predicted_class <- ifelse(logit_data$predicted_prob > 0.5, "Male", "Female")
  confusion <- table(Predicted = logit_data$predicted_class, Actual = logit_data$SEX)
  list(model = model, data = logit_data, confusion = confusion)
}

# Makes the confusion matrix for the logistic regression showing which ones
# it predicted the correct sex for.
plot_confusion_matrix <- function(conf_mat, tissue_name = "") {
  conf_df <- as.data.frame(conf_mat)
  ggplot(conf_df, aes(x = Actual, y = Predicted, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 6, fontface = "bold") +
    scale_fill_gradient(low = "lightgray", high = "steelblue") +
    labs(
      title = paste("Confusion Matrix -", tissue_name),
      x = "Actual Sex",
      y = "Predicted Sex",
      fill = "Count"
    ) +
    theme_minimal(base_size = 14)
}


# Main

# Makes a list assigning the correct GTEx file to the name of the tissue type.
tissue_files <- list(
  "Brain - Frontal cortex ba9" = "~/Library/Mobile Documents/com~apple~CloudDocs/R Stuff/Final Project/gene_tpm_v10_brain_frontal_cortex_ba9.gct.gz",
  "Heart - Left Ventricle" = "~/Library/Mobile Documents/com~apple~CloudDocs/R Stuff/Final Project/gene_tpm_v10_heart_left_ventricle.gct.gz",
  "Lung" = "~/Library/Mobile Documents/com~apple~CloudDocs/R Stuff/Final Project/gene_tpm_v10_lung.gct.gz",
  "Skeletal Muscle" = "~/Library/Mobile Documents/com~apple~CloudDocs/R Stuff/Final Project/gene_tpm_v10_muscle_skeletal.gct.gz",
  "Pituitary Gland" = "~/Library/Mobile Documents/com~apple~CloudDocs/R Stuff/Final Project/gene_tpm_v10_pituitary.gct.gz"
)

# Actually reads in each file so the data is associated with the name
tissue_data_list <- lapply(tissue_files, read_gct)
names(tissue_data_list) <- names(tissue_files)

# Makes each data file go through all of the previoiusly defined functions to
# be analyzed in all three ways and make the plots to visualize them.
results_list <- lapply(names(tissue_data_list), function(tissue_name) {
  expr_all <- tissue_data_list[[tissue_name]]$expression
  
  sample_subject_full <- data.frame(SAMPID = colnames(expr_all)) |>
    mutate(SUBJID = sapply(strsplit(as.character(SAMPID), "\\."), function(x) paste(x[1:2], collapse = "-"))) |>
    left_join(gtex_phen_metadata |> select(SUBJID, SEX, AGE), by = "SUBJID") |>
    left_join(gtex_attr_metadata |> select(SAMPID, SMTS, SMTSD), by = "SAMPID")
  
  filtered_subject_samples <- sample_subject_full |> filter(!is.na(SEX))
  
  top_genes <- get_top_variable_genes(expr_all, filtered_subject_samples)
  p1 <- plot_top_gene_variance(top_genes, tissue_name)
  analysis <- analyze_gene_expression(expr_all, top_genes, sample_subject_full)
  p4 <- plot_gene_regressions(analysis$expression_data, tissue_name)
  p2 <- plot_interaction_volcano(analysis$interaction_results, tissue_name)
  logit <- evaluate_logistic_model(expr_all, top_genes, sample_subject_full)
  p3 <- plot_confusion_matrix(logit$confusion, tissue_name)
  
  list(
    tissue = tissue_name,
    top_genes = top_genes,
    interaction = analysis$interaction_results,
    logistic_model = logit$model,
    plots = list(top_variance = p1, volcano = p2, confusion = p3, regression = p4)
  )
})


# Prints all of the plots created
for (res in results_list) {
  cat("\n====", res$tissue, "====\n")
  
  print(res$plots$top_variance)
  print(res$plots$volcano)
  print(res$plots$confusion)
  print(res$plots$regression)
}

# Combines all the top variable genes from each tissue type to be classified
# all together. Predicts sex for everything to have a combined confusion matrix.
all_top_genes <- unique(unlist(lapply(results_list, function(x) names(x$top_genes))))

combined_expr_matrix <- do.call(cbind, lapply(results_list, function(x) {
  x_data <- tissue_data_list[[x$tissue]]$expression
  x_data[all_top_genes, , drop = FALSE]
}))

combined_metadata <- do.call(rbind, lapply(names(tissue_data_list), function(tissue_name) {
  expr_all <- tissue_data_list[[tissue_name]]$expression
  data.frame(SAMPID = colnames(expr_all)) |>
    mutate(SUBJID = sapply(strsplit(as.character(SAMPID), "\\."), function(x) paste(x[1:2], collapse = "-"))) |>
    left_join(gtex_phen_metadata |> select(SUBJID, SEX, AGE), by = "SUBJID") |>
    mutate(Tissue = tissue_name)
}))

expr_df <- as.data.frame(t(combined_expr_matrix))
expr_df$SAMPID <- rownames(expr_df)

logit_data_all <- left_join(expr_df, combined_metadata, by = "SAMPID") |>
  filter(!is.na(SEX))

logit_data_all$SEX <- factor(ifelse(logit_data_all$SEX == 1, "Male", "Female"))

formula <- as.formula(paste("SEX ~", paste(all_top_genes, collapse = " + ")))
logit_model_all <- glm(formula, data = logit_data_all, family = binomial)

logit_data_all$predicted_prob <- predict(logit_model_all, type = "response")
logit_data_all$predicted_class <- ifelse(logit_data_all$predicted_prob > 0.5, "Male", "Female")

conf_all <- table(Predicted = logit_data_all$predicted_class, Actual = logit_data_all$SEX)
print(conf_all)

# Actually plots the matrix formed before
plot_confusion_matrix(conf_all, tissue_name = "All Tissues Combined")


# Re-established all of the tissue files as one data matrix and cleaned to run
# clustering
tissue_data_list <- lapply(tissue_files, read_gct)
names(tissue_data_list) <- names(tissue_files)

expr_matrices <- lapply(tissue_data_list, function(x) x$expression)

common_genes <- Reduce(intersect, lapply(expr_matrices, rownames))

expr_all_tissues <- do.call(cbind, lapply(expr_matrices, function(mat) mat[common_genes, ]))

all_sample_ids <- colnames(expr_all_tissues)

sample_subject_full <- data.frame(SAMPID = all_sample_ids) |>
  mutate(SUBJID = sapply(strsplit(as.character(SAMPID), "\\."), function(x) paste(x[1:2], collapse = "-"))) |>
  left_join(gtex_phen_metadata |> select(SUBJID, SEX, AGE), by = "SUBJID") |>
  left_join(gtex_attr_metadata |> select(SAMPID, SMTS, SMTSD), by = "SAMPID")

# Mega-function that cleans, filters, transposes, and scales all of the samples
# to run a pca and k-means of the top 1000 variable genes.
run_clean_kmeans_pca <- function(expr_matrix, metadata, k_values = c(2, 3, 5), top_n = 1000) {
  # 1. Filter for most variable genes
  gene_vars <- apply(expr_matrix, 1, var, na.rm = TRUE)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(top_n, length(gene_vars))]
  expr_filtered <- expr_matrix[top_genes, , drop = FALSE]
  
  # 2. Transpose: samples as rows
  expr_t <- t(expr_filtered)
  
  # 3. Remove columns (genes) with zero variance or NA-only values
  good_cols <- apply(expr_t, 2, function(x) sd(x, na.rm = TRUE) > 0 & !all(is.na(x)))
  expr_t <- expr_t[, good_cols, drop = FALSE]
  
  # 4. Remove rows (samples) with NA/NaN/Inf
  is_finite_row <- apply(expr_t, 1, function(row) all(is.finite(row)))
  expr_t_clean <- expr_t[is_finite_row, , drop = FALSE]
  
  # 5. Scale and PCA
  expr_scaled <- scale(expr_t_clean)
  pca_result <- prcomp(expr_scaled, center = TRUE, scale. = FALSE)
  pca_df <- as.data.frame(pca_result$x[, 1:2])
  pca_df$Sample <- rownames(expr_scaled)
  
  # 6. Merge metadata
  metadata_renamed <- metadata
  if (!"Sample" %in% colnames(metadata)) {
    metadata_renamed <- metadata_renamed |> dplyr::rename(Sample = SAMPID)
  }
  pca_df <- left_join(pca_df, metadata_renamed, by = "Sample")
  
  # 7. Cluster and plot
  plot_list <- list()
  cluster_assignments_list <- list()
  
  for (k in k_values) {
    kmeans_result <- kmeans(expr_scaled, centers = k)
    cluster_df <- data.frame(Sample = rownames(expr_scaled),
                             Cluster = as.factor(kmeans_result$cluster))
    pca_df$Cluster <- cluster_df$Cluster
    
    # Save cluster assignments and plot
    cluster_assignments_list[[paste0("k", k)]] <- kmeans_result$cluster
    
    plot_list[[paste0("k", k)]] <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster)) +
      geom_point(size = 3, alpha = 0.7) +
      theme_minimal() +
      labs(
        title = paste("K-means Clustering (k =", k, ") on Expression PCA"),
        x = "PC1", y = "PC2", color = "Cluster"
      )
  }
  
  return(list(plots = plot_list, clusters = cluster_assignments_list))
}
  

# Function to compare the formed clusters form pca/k-means to the metadata
# labels.
compare_clusters_to_all_labels <- function(cluster_assignments, metadata, label_columns = c("SEX", "AGE", "Tissue")) {
  if (!"Sample" %in% colnames(metadata)) {
    metadata <- metadata |> dplyr::rename(Sample = SAMPID)
  }
  
  df <- data.frame(Sample = names(cluster_assignments), Cluster = cluster_assignments)
  df <- left_join(df, metadata, by = "Sample")
  
  results <- list()
  
  for (label_col in label_columns) {
    df_label <- df[!is.na(df[[label_col]]), ]
    ari <- mclust::adjustedRandIndex(df_label$Cluster, df_label[[label_col]])
    cat("\n====", label_col, "====\n")
    cat("Adjusted Rand Index (ARI):", round(ari, 3), "\n")
    print(table(Cluster = df_label$Cluster, Label = df_label[[label_col]]))
    results[[label_col]] <- list(ARI = ari, Table = table(Cluster = df_label$Cluster, Label = df_label[[label_col]]))
  }
  
  invisible(results)
}

# Actually runs the k-means
kmeans_result <- run_clean_kmeans_pca(expr_all_tissues, combined_metadata)

# Compares the clusters to metadata using the function
compare_clusters_to_all_labels(kmeans_result$clusters$k2, combined_metadata)
compare_clusters_to_all_labels(kmeans_result$clusters$k3, combined_metadata)
compare_clusters_to_all_labels(kmeans_result$clusters$k5, combined_metadata)

# Prints the k-means plots for 2,3, and 5 clusters
print(kmeans_result$plots$k2)
print(kmeans_result$plots$k3)
print(kmeans_result$plots$k5)

# Use the clusters/labels comparison to make confusion matrices to show how
# effective/ineffective the clustering was.
visualize_cluster_vs_label <- function(cluster_assignments, metadata, label_col) {
  if (!"Sample" %in% colnames(metadata)) {
    metadata <- metadata |> dplyr::rename(Sample = SAMPID)
  }
  
  df <- data.frame(Sample = names(cluster_assignments), Cluster = as.factor(cluster_assignments)) |>
    left_join(metadata, by = "Sample")
  
  df_label <- df[!is.na(df[[label_col]]), ]
  count_table <- as.data.frame(table(Cluster = df_label$Cluster, Label = df_label[[label_col]]))
  
  ggplot(count_table, aes(x = Label, y = Cluster, fill = Freq)) +
    geom_tile(color = "white") +
    geom_text(aes(label = Freq), size = 5, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = paste("Cluster vs", label_col),
         x = label_col, y = "K-means Cluster", fill = "Count") +
    theme_minimal(base_size = 14)
}

# Prints out all of the plots to see the comparison
visualize_cluster_vs_label(kmeans_result$clusters$k3, combined_metadata, "Tissue")
visualize_cluster_vs_label(kmeans_result$clusters$k3, combined_metadata, "SEX")
visualize_cluster_vs_label(kmeans_result$clusters$k3, combined_metadata, "AGE")

visualize_cluster_vs_label(kmeans_result$clusters$k2, combined_metadata, "Tissue")
visualize_cluster_vs_label(kmeans_result$clusters$k2, combined_metadata, "SEX")
visualize_cluster_vs_label(kmeans_result$clusters$k2, combined_metadata, "AGE")

visualize_cluster_vs_label(kmeans_result$clusters$k5, combined_metadata, "Tissue")
visualize_cluster_vs_label(kmeans_result$clusters$k5, combined_metadata, "SEX")
visualize_cluster_vs_label(kmeans_result$clusters$k5, combined_metadata, "AGE")

