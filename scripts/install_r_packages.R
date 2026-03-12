# Install required R packages not available via conda

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Function to install if not present
install_if_missing <- function(pkg, bioc = FALSE) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("Installing", pkg, "..."))
        if (bioc) {
            BiocManager::install(pkg, ask = FALSE, update = FALSE)
        } else {
            install.packages(pkg)
        }
    } else {
        message(paste(pkg, "already installed"))
    }
}

# Bioconductor packages
install_if_missing("ANCOMBC", bioc = TRUE)
install_if_missing("microbiomeMarker", bioc = TRUE)
install_if_missing("ALDEx2", bioc = TRUE)
install_if_missing("phyloseq", bioc = TRUE)
install_if_missing("microbiome", bioc = TRUE)
install_if_missing("ComplexHeatmap", bioc = TRUE)
install_if_missing("DESeq2", bioc = TRUE)

# GitHub packages
if (!requireNamespace("qiime2R", quietly = TRUE)) {
    message("Installing qiime2R from GitHub...")
    if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
    devtools::install_github("jbisanz/qiime2R")
}

if (!requireNamespace("microViz", quietly = TRUE)) {
    message("Installing microViz from GitHub...")
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github("david-barnett/microViz")
}

message("\nR package installation complete!")
