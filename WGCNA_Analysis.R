#############################################################
###          Enhanced WGCNA Co-expression Analysis       ###
###              Complete Pipeline - Clean Version       ###
#############################################################

#-----------------------------
# Part 1: Initial Setup
#-----------------------------
original_wd <- "D:/WGCNA"
setwd(original_wd)

output_dir <- file.path(original_wd, format(Sys.time(), "WGCNA_Results_%Y%m%d_%H%M%S"))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(">> Created output directory:", output_dir, "\n")
}

subdirs <- c(
  "01_Sample_Clustering",
  "02_Soft_Threshold", 
  "03_Network_Construction",
  "04_Module_Visualization",
  "05_Module_Trait_Relationships",
  "06_Gene_Significance",
  "07_GS_MM_Plots",
  "08_Module_Membership_Distribution", 
  "09_Module_Correlation",
  "10_Cytoscape_Export",
  "11_Edge_Summary"
)

cat(">> Creating analysis subdirectories:\n")
for (dir in subdirs) {
  dir_path <- file.path(output_dir, dir)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    cat("   - Created:", dir_path, "\n")
  }
}

cat(">> Installing/loading required R packages...\n")
required_packages <- c("WGCNA", "ggplot2", "reshape2", "ggrepel", "dendextend", "gplots")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap")
}
if (!requireNamespace("circlize", quietly = TRUE)) {
  install.packages("circlize")
}

library(WGCNA)
library(ComplexHeatmap)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(dendextend)
library(circlize)

enableWGCNAThreads()

color_breaks <- c(-1, 0, 1)
color_palette <- c("blue", "white", "red")
n_colors <- 50

safe_setwd <- function(path) {
  if (dir.exists(path)) {
    setwd(path)
    cat("Working directory set to:", path, "\n")
  } else {
    stop("Directory does not exist: ", path)
  }
}

#-----------------------------
# Part 2: Data Import and Preprocessing
#-----------------------------
cat("\n===== Reading Data =====\n")

isCountData <- FALSE
topGeneFilter <- FALSE
filterMethod <- "MAD"
topGeneCount <- 10000

gene_file <- file.path(original_wd, "GENEFPKM.txt")
if (!file.exists(gene_file)) stop("Expression data file not found: ", gene_file)
cat(">> Reading expression data file:", gene_file, "\n")
geneData <- read.delim(gene_file, header = TRUE, row.names = 1)
original_gene_count <- nrow(geneData)
original_sample_count <- ncol(geneData)
cat(">> Original data dimensions: Genes=", original_gene_count, "Samples=", original_sample_count, "\n")

cat("\n>> Executing data preprocessing pipeline\n")

if (isCountData) {
  keep_genes <- apply(geneData, 1, function(gene) {
    sum(gene >= 10) > 0.1 * original_sample_count
  })
  geneData <- geneData[keep_genes, ]
  cat(">> Gene filtering result: Retained", sum(keep_genes), "genes (filtered", 
      original_gene_count - sum(keep_genes), ")\n")
} else {
  keep_genes <- apply(geneData, 1, function(gene) {
    sum(gene >= 0) > 0.1 * original_sample_count
  })
  geneData <- geneData[keep_genes, ]
  cat(">> Gene filtering result: Retained", sum(keep_genes), "genes (filtered", 
      original_gene_count - sum(keep_genes), ")\n")
}

if (isCountData) {
  cat(">> Performing DESeq2 vst transformation\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
  library(DESeq2)
  
  if(!all(geneData%%1 == 0)) warning("Count data contains non-integer values")
  
  colData <- data.frame(row.names = colnames(geneData))
  dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(geneData),
    colData = colData,
    design = ~ 1
  )
  vsd <- vst(dds, blind = TRUE)
  vstData <- signif(assay(vsd), digits = 6)
  geneData <- as.data.frame(vstData)
} else {
  cat(">> No log transformation, only handling negative values\n")
  geneData[geneData < 0] <- 0
  geneData <- signif(geneData, digits = 6)
}

if (topGeneFilter) {
  cat("\n>> Performing advanced gene filtering: Based on", filterMethod, "retaining top", topGeneCount, "highly variable genes\n")
  
  if (filterMethod == "MAD") {
    gene_variability <- apply(geneData, 1, mad)
    cat(">> Using MAD as variability measure\n")
  } else if (filterMethod == "Var") {
    gene_variability <- apply(geneData, 1, var)
    cat(">> Using Variance as variability measure\n")
  } else {
    stop("Invalid filtering method: ", filterMethod, ". Please choose 'MAD' or 'Var'")
  }
  
  if (nrow(geneData) > topGeneCount) {
    top_indices <- order(gene_variability, decreasing = TRUE)[1:topGeneCount]
    geneData <- geneData[top_indices, ]
    cat(">> Advanced gene filtering result: Retained", topGeneCount, "highly variable genes (filtered", 
        nrow(geneData) - topGeneCount, ")\n")
  } else {
    cat(">> Gene count (", nrow(geneData), ") is less than target (", topGeneCount, "), skipping this step\n")
  }
}

datExpr <- as.data.frame(t(geneData))
cat(">> Processed data dimensions: Samples=", nrow(datExpr), "Genes=", ncol(datExpr), "\n")

write.table(datExpr, file.path(output_dir, "Processed_Expression_Data.txt"), 
            sep = "\t", quote = FALSE)
cat(">> Preprocessed expression data saved: Processed_Expression_Data.txt\n")

cat("\n===== Trait Data Import =====\n")
trait_file <- file.path(original_wd, "biaoxing.txt")
if (file.exists(trait_file)) {
  datTraits <- read.delim(trait_file, header = TRUE, row.names = 1)
  cat(">> Successfully read trait data: ", trait_file, "\n")
  cat(">> Trait data dimensions: Samples=", nrow(datTraits), "Traits=", ncol(datTraits), "\n")
  
  if(all(rownames(datExpr) %in% rownames(datTraits))) {
    datTraits <- datTraits[match(rownames(datExpr), rownames(datTraits)), , drop = FALSE]
    cat(">> Aligned trait data with expression data sample order\n")
  } else {
    warning("!! Sample names mismatch between expression and trait data")
  }
} else {
  datTraits <- data.frame(
    row.names = rownames(datExpr),
    TissueType = rep(NA, nrow(datExpr)),
    Group = rep(NA, nrow(datExpr))
  )
  cat("!! Trait file not found, creating empty trait data frame for basic analysis\n")
}

#-----------------------------
# Part 3: Sample Quality Check
#-----------------------------
target_dir <- file.path(output_dir, subdirs[1])
safe_setwd(target_dir)

cat("\n>> Performing sample and gene quality check\n")

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    cat(">> Removing", sum(!gsg$goodGenes), "low-quality genes\n")
    datExpr <- datExpr[, gsg$goodGenes]
  }
  if (sum(!gsg$goodSamples) > 0) {
    cat(">> Removing", sum(!gsg$goodSamples), "low-quality samples\n")
    datExpr <- datExpr[gsg$goodSamples, ]
    
    if (exists("datTraits") && nrow(datTraits) > 0) {
      datTraits <- datTraits[gsg$goodSamples, , drop = FALSE]
    }
  }
  cat(">> Final data dimensions: Samples=", nrow(datExpr), "Genes=", ncol(datExpr), "\n")
}

save(datExpr, file = file.path(output_dir, "Processed_Data.RData"))
cat(">> Final analysis dataset saved: Processed_Data.RData\n")

#-----------------------------
# Part 4: Outlier Sample Detection
#-----------------------------
cat("\n===== Outlier Sample Detection =====\n")

sampleTree <- hclust(dist(datExpr), method = "average")
pdf("SampleClustering.pdf", width = 12, height = 6)
par(cex = 0.6)
plot(sampleTree, 
     main = "Sample Clustering to Detect Outliers", 
     sub = "", 
     xlab = "")
dev.off()
cat(">> Sample clustering plot saved: SampleClustering.pdf\n")

#-----------------------------
# Part 5: Soft Threshold Selection
#-----------------------------
cat("\n===== Soft Threshold Selection =====\n")

target_dir <- file.path(output_dir, subdirs[2])
safe_setwd(target_dir)

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

write.table(sft$fitIndices, "SoftThreshold_FitIndices.txt", sep = "\t", quote = FALSE)
cat(">> Soft threshold indices file saved: SoftThreshold_FitIndices.txt\n")

useAdjustedR2 <- TRUE

if (useAdjustedR2) {
  r2_values <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
  ylab_text <- "Adjusted R²"
  cat(">> Using adjusted R² for soft threshold selection\n")
} else {
  r2_values <- sft$fitIndices[, 2]
  ylab_text <- "Scale Free Topology Fit"
  cat(">> Using raw R² for soft threshold selection\n")
}

pdf("SoftThreshold.pdf", width = 9, height = 5)
par(mfrow = c(1, 2))

plot(sft$fitIndices[, 1], r2_values,
     xlab = "Soft Threshold (power)", 
     ylab = ylab_text,
     main = "Scale Independence", 
     type = "n")
text(sft$fitIndices[, 1], r2_values, labels = powers, col = "blue")
abline(h = 0.8, col = "red", lty = 2)
abline(h = 0.9, col = "darkgreen", lty = 2)

plot(sft$fitIndices[, 1], sft$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     main = "Mean Connectivity", 
     type = "n")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "blue")
dev.off()
cat(">> Soft threshold plot saved: SoftThreshold.pdf\n")

candidate_powers <- sft$fitIndices$Power[r2_values > 0.8]

if (length(candidate_powers) > 0) {
  best_power <- min(candidate_powers)
  cat(">> Recommended soft threshold power:", best_power, " (", ylab_text, "=", 
      round(r2_values[sft$fitIndices$Power == best_power], 3), 
      ", Mean connectivity =", sft$fitIndices$mean.k.[sft$fitIndices$Power == best_power], ")\n")
} else {
  best_idx <- which.max(r2_values)
  best_power <- sft$fitIndices$Power[best_idx]
  best_r2 <- r2_values[best_idx]
  
  mean_conn <- sft$fitIndices$mean.k.[best_idx]
  
  if (mean_conn > 100) {
    viable_idx <- which(sft$fitIndices$mean.k. < 100)
    
    if (length(viable_idx) > 0) {
      best_viable_idx <- viable_idx[which.max(r2_values[viable_idx])]
      alt_power <- sft$fitIndices$Power[best_viable_idx]
      alt_r2 <- r2_values[best_viable_idx]
      
      cat("!! Warning: Initial power", best_power, "has high mean connectivity (", round(mean_conn, 1), ")\n")
      cat(">> Using alternative power:", alt_power, " (", ylab_text, "=", round(alt_r2, 3), 
          ", Mean connectivity =", sft$fitIndices$mean.k.[best_viable_idx], ")\n")
      best_power <- alt_power
    } else {
      cat("!! Warning: No threshold with R²>0.8 found, using highest", ylab_text, 
          "(", round(best_r2, 3), ") with power: ", best_power, 
          ", Mean connectivity =", round(mean_conn, 1), "\n")
    }
  } else {
    cat("!! Warning: No threshold with R²>0.8 found, using highest", ylab_text, 
        "(", round(best_r2, 3), ") with power: ", best_power, 
        ", Mean connectivity =", round(mean_conn, 1), "\n")
  }
}

cat(">> Final selected soft threshold power:", best_power, "\n")

manual_power <- 11
if (!is.na(manual_power) && manual_power %in% powers) {
  best_power <- manual_power
  cat(">> User manually overrode soft threshold power:", best_power, "\n")
}

#-----------------------------
# Part 6: Network Construction
#-----------------------------
cat("\n===== Network Construction =====\n")
target_dir <- file.path(output_dir, subdirs[3])
safe_setwd(target_dir)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

cat(">> Building network (this may take a while)...\n")
net <- blockwiseModules(
  datExpr,
  power = best_power,
  TOMType = "unsigned",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "WGCNATOM",
  verbose = 3,
  maxBlockSize = 16000
)

module_counts <- as.data.frame(table(net$colors))
colnames(module_counts) <- c("ModuleColor", "GeneCount")
cat(">> Total modules identified:", length(unique(net$colors)), "\n")
print(module_counts)

write.csv(module_counts, "Module_Gene_Counts.csv")
cat(">> Module gene counts saved: Module_Gene_Counts.csv\n")

save(net, file = "Network_Object.RData")
cat(">> Network object saved: Network_Object.RData\n")

#-----------------------------
# Part 7: Module Visualization
#-----------------------------
cat("\n===== Module Visualization =====\n")
target_dir <- file.path(output_dir, subdirs[4])
safe_setwd(target_dir)

moduleColors <- labels2colors(net$colors)

MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
rownames(MEs) <- rownames(datExpr)

cat("\n>> Calculating and saving KME values...\n")
kME <- signedKME(datExpr, MEs, outputColumnName = "kME_")

kme_data <- data.frame(
  GeneID = colnames(datExpr),
  ModuleColor = moduleColors,
  kME = kME
)

kme_file_path <- "Gene_Module_Membership_KME_Values.txt"
write.table(kme_data, kme_file_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat(">> KME values file saved:", kme_file_path, "\n")

cat("\n>> Plotting module KME distributions...\n")

kme_cols <- grep("^kME\\.kME_", colnames(kme_data), value = TRUE)
for (col in kme_cols) {
  kme_data[[col]] <- as.numeric(as.character(kme_data[[col]]))
}

unique_colors <- unique(moduleColors)
unique_colors <- unique_colors[unique_colors != "grey"]

for (color in unique_colors) {
  kme_column <- paste0("kME.kME_", color)
  
  if (!(kme_column %in% colnames(kme_data))) {
    cat(">> Warning: Module", color, "missing KME column (expected:", kme_column, "), skipping\n")
    next
  }
  
  module_data <- kme_data[kme_data$ModuleColor == color, ]
  module_kme <- module_data[[kme_column]]
  
  module_kme <- as.numeric(module_kme)
  module_kme <- module_kme[!is.na(module_kme) & is.finite(module_kme) & module_kme >= -1 & module_kme <= 1]
  
  if (length(module_kme) == 0) {
    cat(">> Skipping module", color, "(no valid genes)\n")
    next
  }
  
  mean_kme <- mean(module_kme, na.rm = TRUE)
  median_kme <- median(module_kme, na.rm = TRUE)
  n_genes <- length(module_kme)
  
  cat(">> Module", color, "data range:", round(range(module_kme), 3), 
      "valid values:", n_genes, "\n")
  
  plot_df <- data.frame(kME = module_kme)
  
  p <- ggplot(plot_df, aes(x = kME)) +
    geom_histogram(
      binwidth = 0.05,
      boundary = -1,
      fill = color,
      color = "black",
      alpha = 0.7
    ) +
    aes(y = after_stat(density)) +
    geom_density(
      color = "red",
      linewidth = 1,
      alpha = 0.5
    ) +
    geom_vline(
      xintercept = mean_kme,
      color = "blue",
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_vline(
      xintercept = median_kme,
      color = "green",
      linetype = "dashed",
      linewidth = 1
    ) +
    labs(
      title = paste("KME Distribution in", color, "Module"),
      subtitle = paste(
        "n =", n_genes, "genes | ",
        "Mean =", round(mean_kme, 3), "| ",
        "Median =", round(median_kme, 3)
      ),
      x = paste("kME (Module Membership in", color, "module)"),
      y = "Density"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    ) +
    coord_cartesian(xlim = c(-1, 1)) +
    annotate(
      "text",
      x = -0.95,
      y = Inf,
      label = paste("Mean:", round(mean_kme, 3)),
      color = "blue",
      size = 4,
      hjust = 0,
      vjust = 1.5
    ) +
    annotate(
      "text",
      x =  -0.95,
      y = Inf,
      label = paste("Median:", round(median_kme, 3)),
      color = "green",
      size = 4,
      hjust = 0,
      vjust = 3.5
    )
  
  pdf_file <- paste0("KME_Distribution_", color, ".pdf")
  ggsave(
    pdf_file,
    plot = p,
    width = 8,
    height = 6,
    device = "pdf"
  )
  
  png_file <- paste0("KME_Distribution_", color, ".png")
  ggsave(
    png_file,
    plot = p,
    width = 8,
    height = 6,
    dpi = 300,
    device = "png"
  )
  
  cat(">> Saved:", color, "module KME distribution (genes =", n_genes, ")\n")
}

pdf("ModuleDendrogram.pdf", width = 12, height = 8)
plotDendroAndColors(
  net$dendrograms[[1]], 
  moduleColors, 
  dendroLabels = FALSE, 
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene Clustering and Module Colors",
  marAll = c(1, 5, 3, 1)
)
dev.off()
cat(">> Module dendrogram saved: ModuleDendrogram.pdf\n")

#-----------------------------
# Module Correlation Analysis
#-----------------------------
cat("\n===== Module Correlation Analysis =====\n")

target_dir <- file.path(output_dir, subdirs[9])
safe_setwd(target_dir)

rownames(MEs) <- rownames(datExpr)

cat(">> Calculating inter-module correlations...\n")
ME_cor <- cor(MEs, use = "pairwise.complete.obs")

ME_pvalue <- matrix(NA, ncol(MEs), ncol(MEs))
for (i in 1:ncol(MEs)) {
  for (j in 1:ncol(MEs)) {
    ME_pvalue[i, j] <- cor.test(MEs[,i], MEs[,j])$p.value
  }
}

module_names <- colnames(MEs)
module_colors <- gsub("ME", "", module_names)

ha_col <- HeatmapAnnotation(
  Module = module_names,
  col = list(Module = structure(module_colors, names = module_names)),
  show_legend = FALSE,
  annotation_name_gp = gpar(col = "transparent")
)

ha_row <- rowAnnotation(
  Module = module_names,
  col = list(Module = structure(module_colors, names = module_names)),
  show_legend = FALSE,
  annotation_name_gp = gpar(col = "transparent")
)

n_mods <- ncol(MEs)
base_font_size <- 10 - n_mods * 0.15
if (base_font_size < 4) base_font_size <- 4

pdf("ModuleCorrelationHeatmap.pdf", width = max(10, n_mods*0.8), height = max(8, n_mods*0.7))
Heatmap(
  matrix = ME_cor,
  name = "Correlation",
  col = colorRamp2(seq(-1, 1, length = n_colors), colorRampPalette(color_palette)(n_colors)),
  top_annotation = ha_col,
  left_annotation = ha_row,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", ME_cor[i, j]), x, y + height*0.15, 
              gp = gpar(fontsize = base_font_size, col = "black"))
    p_val_text <- ifelse(ME_pvalue[i, j] < 0.001, "(<0.001)", 
                         sprintf("(%.3f)", ME_pvalue[i, j]))
    grid.text(p_val_text, x, y - height*0.15, 
              gp = gpar(fontsize = base_font_size*0.8, col = "black"))
  },
  show_row_names = FALSE,
  show_column_names = FALSE,
  row_dend_reorder = TRUE,
  column_dend_reorder = TRUE,
  rect_gp = gpar(col = "white", lwd = 1),
  row_names_gp = gpar(fontsize = max(8, 12 - n_mods*0.1)),
  column_names_gp = gpar(fontsize = max(8, 12 - n_mods*0.1))
)
dev.off()
cat(">> Module correlation heatmap saved: ModuleCorrelationHeatmap.pdf\n")

write.table(ME_cor, "Module_Correlation_Values.txt", sep = "\t", quote = FALSE)
write.table(ME_pvalue, "Module_Correlation_Pvalues.txt", sep = "\t", quote = FALSE)
cat(">> Module correlation numerical files saved\n")

#-----------------------------
# Part 8: Module-Trait Relationships
#-----------------------------
cat("\n===== Module-Trait Relationship Analysis =====\n")

target_dir <- file.path(output_dir, subdirs[5])
safe_setwd(target_dir)

tissueTypes <- c("Bark", "Flower", "Fruit", "Leaf", "Root", "Wood")

traitData <- data.frame(row.names = rownames(datExpr))

if ("TissueType" %in% colnames(datTraits)) {
  datTraits$TissueType <- factor(datTraits$TissueType, levels = tissueTypes)
  tissueData <- model.matrix(~0 + datTraits$TissueType)
  colnames(tissueData) <- levels(datTraits$TissueType)
  traitData <- cbind(traitData, tissueData)
}

numericTraits <- datTraits[, sapply(datTraits, is.numeric), drop = FALSE]
if (ncol(numericTraits) > 0) {
  traitData <- cbind(traitData, numericTraits)
}

if (ncol(traitData) > 0 && nrow(traitData) > 0) {
  cat(">> Analyzing the following traits:\n")
  cat("   - Categorical traits:", if ("TissueType" %in% colnames(datTraits)) "TissueType" else "None", "\n")
  cat("   - Numerical traits:", if (ncol(numericTraits) > 0) paste(colnames(numericTraits), collapse = ", ") else "None", "\n")
  
  cat(">> Calculating module-trait correlations...\n")
  moduleTraitCor <- cor(MEs, traitData, use = "pairwise.complete.obs")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  n_traits <- ncol(traitData)
  n_modules <- ncol(MEs)
  
  pdf_width <- max(8, n_traits * 0.8)
  pdf_height <- max(8, n_modules * 0.7)
  base_font_size <- 10 - max(n_traits, n_modules) * 0.1
  if (base_font_size < 5) base_font_size <- 5
  
  pdf("ModuleTraitCorrelation.pdf", width = pdf_width, height = pdf_height)
  
  textMatrix <- paste0(
    signif(moduleTraitCor, 2), "\n(", 
    signif(moduleTraitPvalue, 1), ")"
  )
  dim(textMatrix) <- dim(moduleTraitCor)
  
  par(mar = c(
    max(6, n_traits * 0.5),
    max(8, n_modules * 0.5),
    3,
    1
  ))
  
  labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = colnames(traitData),
    yLabels = colnames(MEs),
    ySymbols = colnames(MEs),
    colorLabels = FALSE,
    colors = colorRampPalette(color_palette)(n_colors),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = base_font_size * 0.1,
    zlim = c(-1, 1),
    main = "Module-Trait Relationships"
  )
  dev.off()
  cat(">> Module-trait heatmap saved: ModuleTraitCorrelation.pdf\n")
  
  write.table(moduleTraitCor, "ModuleTrait_Correlation.txt", sep = "\t", quote = FALSE)
  write.table(moduleTraitPvalue, "ModuleTrait_Pvalues.txt", sep = "\t", quote = FALSE)
  cat(">> Module-trait correlation values saved\n")
} else {
  warning(">> No valid trait data available, skipping module-trait analysis")
}

#============================================================
# Gene Significance (GS) and Module Membership (MM/kME) Analysis
#============================================================
#-----------------------------
# Part 9: Gene Significance Analysis
#-----------------------------
cat("\n===== Gene Significance Analysis =====\n")

target_dir <- file.path(output_dir, subdirs[6])
safe_setwd(target_dir)

if (exists("traitData") && ncol(traitData) > 0 && nrow(traitData) > 0) {
  geneTraitSignificance <- as.data.frame(cor(datExpr, traitData, use = "p"))
  GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr)))
  
  colnames(geneTraitSignificance) <- paste("GS.", colnames(traitData), sep = "")
  colnames(GSPvalue) <- paste("p.GS.", colnames(traitData), sep = "")
  
  geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
  MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))
  
  colnames(geneModuleMembership) <- paste("MM.", substring(colnames(MEs), 3), sep = "")
  colnames(MMPvalue) <- paste("p.MM.", substring(colnames(MEs), 3), sep = "")
  
  annotation_file <- file.path(original_wd, "gene_annotation.txt")
  if (file.exists(annotation_file)) {
    cat(">> Detected gene annotation file:", annotation_file, "\n")
    geneAnnotation <- read.delim(annotation_file, header = TRUE, stringsAsFactors = FALSE)
  } else {
    geneAnnotation <- data.frame(GeneID = colnames(datExpr), GeneSymbol = colnames(datExpr))
    cat(">> Gene annotation file not found, using GeneID as gene symbol\n")
  }
  
  geneStats <- data.frame(
    GeneID = colnames(datExpr),
    Module = moduleColors
  )
  
  if ("GeneSymbol" %in% colnames(geneAnnotation)) {
    geneStats$GeneSymbol <- geneAnnotation$GeneSymbol[match(geneStats$GeneID, geneAnnotation$GeneID)]
  } else {
    geneStats$GeneSymbol <- geneStats$GeneID
  }
  
  geneStats <- data.frame(
    geneStats,
    geneTraitSignificance,
    GSPvalue,
    geneModuleMembership,
    MMPvalue
  )
  
  write.table(geneStats, "GeneSignificance_ModuleMembership.txt", sep = "\t", quote = FALSE, row.names = FALSE)
  cat(">> Gene analysis results saved: GeneSignificance_ModuleMembership.txt\n")
} else {
  warning(">> No valid trait data, skipping gene significance analysis")
}

#-----------------------------
# Part 10: GS vs MM Visualization
#-----------------------------
cat("\n===== GS vs MM Relationship Visualization =====\n")

target_dir <- file.path(output_dir, subdirs[7])
safe_setwd(target_dir)

if (exists("traitData") && ncol(traitData) > 0 && nrow(traitData) > 0) {
  for (trait in colnames(traitData)) {
    traitColumn <- paste("GS.", trait, sep = "")
    
    sigModules <- rownames(moduleTraitCor)[moduleTraitPvalue[, trait] < 0.05]
    
    if (length(sigModules) > 0) {
      for (module in sigModules) {
        moduleColor <- substring(module, 3)
        mmColumn <- paste("MM.", moduleColor, sep = "")
        
        moduleGenes <- (moduleColors == moduleColor)
        nGenes <- sum(moduleGenes)
        
        if (nGenes < 10) {
          cat(">> Skipping module", moduleColor, "(genes <10)\n")
          next
        }
        
        plot_data <- data.frame(
          MM = geneModuleMembership[moduleGenes, mmColumn]^2,
          GS = abs(geneTraitSignificance[moduleGenes, traitColumn]),
          Gene = geneStats$GeneSymbol[moduleGenes]
        )
        
        cor_test <- cor.test(plot_data$MM, plot_data$GS)
        cor_value <- round(cor_test$estimate, 3)
        p_value <- ifelse(cor_test$p.value < 0.001, "< 0.001", 
                          paste0("= ", round(cor_test$p.value, 3)))
        
        p_base <- ggplot(plot_data, aes(x = MM, y = GS)) +
          geom_point(shape = 21, fill = moduleColor, color = "black", size = 2.5, stroke = 0.5, alpha = 0.7) +
          geom_smooth(
            method = "lm", 
            formula = y ~ x,
            se = TRUE, 
            color = "red", 
            fill = "#40404060", 
            alpha = 0.4, 
            linewidth = 0.8,
            fullrange = TRUE,
            level = 0.95
          ) +
          labs(
            x = paste("Module Membership (kME²) in", moduleColor, "module"),
            y = paste("Gene Significance (|GS|) for", trait),
            title = paste(trait, "vs", moduleColor, "module"),
            subtitle = paste("n =", nGenes, "genes | r =", cor_value, "| p", p_value)
          ) +
          theme_classic(base_size = 12) +
          theme(
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 12, hjust = 0.5),
            axis.title = element_text(size = 12, face = "bold"),
            axis.text = element_text(size = 10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line.x = element_line(color = "black", linewidth = 0.5),
            axis.line.y = element_line(color = "black", linewidth = 0.5)
          ) +
          coord_cartesian(xlim = c(min(plot_data$MM)*0.95, max(plot_data$MM)*1.05), 
                          ylim = c(min(plot_data$GS)*0.95, max(plot_data$GS)*1.05))
        
        ggsave(paste0("GSvsMM_", trait, "_", moduleColor, ".pdf"), 
               p_base, width = 8, height = 6)
        
        p_labeled <- p_base + 
          ggrepel::geom_text_repel(
            aes(label = Gene),
            data = subset(plot_data, GS > quantile(GS, 0.9) | MM > quantile(MM, 0.9)),
            size = 3,
            max.overlaps = 100,
            segment.color = "gray30",
            box.padding = 0.5,
            force = 0.5
          )
        
        ggsave(paste0("GSvsMM_", trait, "_", moduleColor, "_labeled.png"), 
               p_labeled, width = 8, height = 6, dpi = 300)
        
        cat(">> Saved plot:", trait, "vs", moduleColor, "\n")
      }
    }
  }
} else {
  warning(">> No valid trait data, skipping GS-MM analysis")
}

#-----------------------------
# Part 11: Cytoscape Network Export
#-----------------------------
cat("\n===== Cytoscape Network Export =====\n")

target_dir <- file.path(output_dir, subdirs[10])
if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
setwd(target_dir)
cat(">> Export directory set to:", target_dir, "\n")

annotation_file <- file.path(original_wd, "gene_annotation.txt")
if (file.exists(annotation_file)) {
  geneAnnotation <- read.delim(annotation_file, header = TRUE, stringsAsFactors = FALSE)
  cat(">> Detected gene annotation file:", annotation_file, "\n")
} else {
  geneAnnotation <- NULL
  warning("Gene annotation file not found, using GeneID as display name")
}

geneModule <- data.frame(
  GeneID = colnames(datExpr), 
  Module = moduleColors,
  stringsAsFactors = FALSE
)

if (!is.null(geneAnnotation)) {
  geneModule <- merge(geneModule, geneAnnotation, by = "GeneID", all.x = TRUE)
  geneModule$GeneSymbol <- ifelse(
    is.na(geneModule$GeneSymbol), 
    geneModule$GeneID, 
    geneModule$GeneSymbol
  )
} else {
  geneModule$GeneSymbol <- geneModule$GeneID
}

write.table(geneModule, "CytoscapeNodes_AllModules.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(">> Global node file saved: CytoscapeNodes_AllModules.txt\n")

tom_file <- file.path(output_dir, subdirs[3], "WGCNATOM-block.1.RData")
if (file.exists(tom_file)) {
  load(tom_file)
  TOM <- as.matrix(TOM)
  cat(">> Loaded TOM matrix:", tom_file, "\n")
} else {
  cat(">> Warning: TOM file not found, recalculating...\n")
  TOM <- TOMsimilarityFromExpr(
    datExpr,
    power = best_power,       
    networkType = "unsigned",
    verbose = 0
  )
}
dimnames(TOM) <- list(colnames(datExpr), colnames(datExpr))

modules <- setdiff(unique(moduleColors), "grey")  

for (module in modules) {
  moduleGenes <- (moduleColors == module)
  subGenes <- colnames(datExpr)[moduleGenes]
  
  if (length(subGenes) < 3) {
    cat(">> Skipping module ", module, " (genes=", length(subGenes), "<3)\n", sep = "")
    next
  }
  
  geneModule_sub <- geneModule[moduleGenes, ]
  node_filename <- paste0("CytoscapeNodes_", module, ".txt")
  write.table(geneModule_sub, node_filename, sep = "\t", quote = FALSE, row.names = FALSE)
  
  subTOM <- TOM[moduleGenes, moduleGenes]
  dimnames(subTOM) <- list(subGenes, subGenes)
  
  edges <- reshape2::melt(subTOM)
  colnames(edges) <- c("fromGene", "toGene", "weight")
  edges <- edges[edges$fromGene != edges$toGene, ]
  edges <- edges[order(-edges$weight), ]
  
  threshold <- quantile(edges$weight, probs = 0.95, na.rm = TRUE)
  
  if (threshold > 0.99) {
    threshold <- quantile(edges$weight, probs = 0.90, na.rm = TRUE)
    cat(">> Module", module, "using dynamic threshold: ", round(threshold, 3), "\n")
  }
  
  valid_edges <- edges[edges$weight > threshold, ]
  if (nrow(valid_edges) < 20) {
    valid_edges <- edges[1:min(20, nrow(edges)), ]
    cat(">> Module", module, "retaining top 20 high-weight edges\n")
  }
  
  valid_edges$direction <- "undirected"
  
  if (!is.null(geneAnnotation)) {
    valid_edges <- merge(valid_edges, geneAnnotation, 
                         by.x = "fromGene", by.y = "GeneID", all.x = TRUE)
    valid_edges <- merge(valid_edges, geneAnnotation, 
                         by.x = "toGene", by.y = "GeneID", all.x = TRUE,
                         suffixes = c("_from", "_to"))
    valid_edges$fromAltName <- valid_edges$GeneSymbol_from
    valid_edges$toAltName <- valid_edges$GeneSymbol_to
  } else {
    valid_edges$fromAltName <- valid_edges$fromGene
    valid_edges$toAltName <- valid_edges$toGene
  }
  
  output_edges <- valid_edges[, c("fromGene", "toGene", "weight", 
                                  "direction", "fromAltName", "toAltName")]
  colnames(output_edges) <- c("fromNode", "toNode", "weight", 
                              "direction", "fromAltName", "toAltName")
  
  edge_filename <- paste0("CytoscapeEdges_", module, ".txt")
  write.table(output_edges, edge_filename, sep = "\t", 
              quote = FALSE, row.names = FALSE)
  
  cat(">> Module ", module, "export completed:\n",
      "  - Node file: ", node_filename, " (genes: ", nrow(geneModule_sub), ")\n",
      "  - Edge list: ", edge_filename, " (edges: ", nrow(output_edges), ")\n\n", sep = "")
}

#-----------------------------
# Part 12: Edge List Summary
#-----------------------------
cat("\n===== Edge List Summary =====\n")

target_dir <- file.path(output_dir, subdirs[11])
if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)
setwd(target_dir)
cat(">> Summary directory set to:", target_dir, "\n")

module_files <- list.files(
  path = file.path(output_dir, subdirs[10]), 
  pattern = "CytoscapeEdges_.*\\.txt$", 
  full.names = TRUE
)

all_edges <- data.frame(
  fromNode = character(),
  toNode = character(),
  weight = numeric(),
  direction = character(),
  fromAltName = character(),
  toAltName = character(),
  Module = character(),
  stringsAsFactors = FALSE
)

if (length(module_files) > 0) {
  cat(">> Found", length(module_files), "module edge list files\n")
  
  for (file in module_files) {
    module <- gsub("CytoscapeEdges_", "", gsub("\\.txt$", "", basename(file)))
    cat(">> Processing module:", module, "\n")
    
    edges <- tryCatch({
      read.delim(file, stringsAsFactors = FALSE)
    }, error = function(e) {
      warning("File reading failed: ", file, " - ", e$message)
      NULL
    })
    
    if (!is.null(edges) && nrow(edges) > 0) {
      edges$Module <- module
      all_edges <- rbind(all_edges, edges)
    }
  }
  
  if (nrow(all_edges) > 0) {
    write.table(all_edges, "AllModules_EdgesSummary.txt", 
                sep = "\t", quote = FALSE, row.names = FALSE)
    cat(">> Edge list summary completed:\n",
        "  - File: AllModules_EdgesSummary.txt\n",
        "  - Total modules: ", length(unique(all_edges$Module)), "\n",
        "  - Total edges: ", nrow(all_edges), "\n")
  } else {
    warning("No valid edge lists found, summary file not generated")
  }
} else {
  warning("No module edge list files found")
}

#-----------------------------
# Final Output Summary
#-----------------------------
setwd(output_dir)
cat("\n===== WGCNA Analysis Completed =====")
cat("\n>> Analysis completion time:", date(), "\n")

cat("\n>> Saving analysis parameter file...\n")
params_list <- list(
  "Output Directory" = output_dir,
  "Soft Threshold Power" = best_power,
  "TOM Type" = "unsigned",
  "Minimum Module Size" = 30,
  "Module Merge Threshold" = 0.25,
  "Gene Filtering Method" = if(topGeneFilter) paste(filterMethod, "top", topGeneCount) else "None",
  "Genes Retained" = if(topGeneFilter) topGeneCount else "None",
  "Is Count Data" = isCountData,
  "Use Adjusted R²" = useAdjustedR2,
  "Trait Data Types" = if(ncol(traitData) > 0) paste(colnames(traitData), collapse = ", ") else "No trait data",
  "Sample Count" = nrow(datExpr),
  "Gene Count" = ncol(datExpr),
  "Expression File" = basename(gene_file),
  "Trait File" = if(file.exists(trait_file)) basename(trait_file) else "Not used",
  "Run Start Time" = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

params_data <- data.frame(
  Parameter = names(params_list),
  Value = as.character(params_list),
  stringsAsFactors = FALSE
)
params_file_path <- "Analysis_Parameters_Summary.txt"
write.table(params_data, params_file_path, sep = ":\t", quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
cat(">> Parameter file saved:", params_file_path, "\n")

cat("\n>> All results saved in:", output_dir)
cat("\n>> Core Results Summary:")
cat("\n  1. Sample Clustering: ", file.path(subdirs[1], "SampleClustering.pdf"))
cat("\n  2. Soft Threshold Analysis: ", file.path(subdirs[2], "SoftThreshold.pdf"))
cat("\n  3. Module Statistics: ", file.path(subdirs[3], "Module_Gene_Counts.csv"))
cat("\n  4. Module Dendrogram: ", file.path(subdirs[4], "ModuleDendrogram.pdf"))
cat("\n  5. Module KME Values: ", file.path(subdirs[4], "Gene_Module_Membership_KME_Values.txt")) 
cat("\n  6. Module-Trait Relationships: ", file.path(subdirs[5], "ModuleTraitCorrelation.pdf"))
cat("\n  7. Gene Analysis Results: ", file.path(subdirs[6], "GeneSignificance_ModuleMembership.txt"))
cat("\n  8. GS-MM Relationship Plots: ", subdirs[7])
cat("\n  9. Module Membership Distribution: ", file.path(subdirs[8], "ModuleMembership_Distribution.pdf"))
cat("\n  10. Inter-module Correlation Heatmap: ", file.path(subdirs[9], "ModuleCorrelationHeatmap.pdf"))
cat("\n  11. Cytoscape Data: ", subdirs[10])
cat("\n  12. Edge List Summary: ", file.path(subdirs[11], "AllModules_EdgesSummary.txt"))
cat("\n  13. Parameter File: ", params_file_path)
cat("\n\n>> Analysis completed, thank you for using the enhanced WGCNA pipeline!\n")