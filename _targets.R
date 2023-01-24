library(targets)

source("code/functions.R")

tar_option_set(
  packages = c(
    "Signac",
    "Seurat",
    "GenomeInfoDb",
    "EnsDb.Hsapiens.v75",
    "patchwork",
    "tidyverse",
    "here"
  )
)

list(
  # target raw data, metadata, fragments
  tar_target(raw_data, here("data", "atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"), format = "file"),
  tar_target(meta, here("data", "atac_v1_pbmc_10k_singlecell.csv"), format = "file"),
  tar_target(frags, here("data", "atac_v1_pbmc_10k_fragments.tsv.gz"), format = "file"),
  
  # Load data
  tar_target(counts, get_data(raw_data)),
  tar_target(metadata, get_metadata(meta)),
  
  # Create Chromatin Assay, Seurat object and preprocess object
  tar_target(chrom_assay, create_chrom_assay(counts, frags)),
  tar_target(
    pbmc, 
    create_seurat(chrom_assay, metadata) %>% 
      get_annotation() %>% 
      compute_qc() %>% 
      modify_metadata()
  ),
  
  # Plot some qc
  tar_target(qc_plots, plot_qc(pbmc)),
  tar_target(qc_violins, plot_qc_violins(pbmc))
)

