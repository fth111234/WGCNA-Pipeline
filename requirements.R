# Required R packages for WGCNA analysis
# Install using: source("requirements.R")

# CRAN packages
cran_packages <- c(
  "WGCNA",
  "ggplot2", 
  "reshape2",
  "ggrepel",
  "dendextend",
  "gplots",
  "circlize"
)

# Bioconductor packages
bioc_packages <- c(
  "ComplexHeatmap"
)

# Install CRAN packages
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

cat("All required packages installed successfully!\n")