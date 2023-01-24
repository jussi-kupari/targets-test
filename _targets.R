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
  # download data files
  tar_target(peak_matrix_h5, get_peak_matrix_h5(), format = "file"),
  tar_target(metadata_csv, get_metadata_csv(), format = "file"),
  tar_target(frags_file, get_frags(), format = "file"),
  tar_target(frags_indx, get_frags_indx(), format = "file"),
  tar_target(proc_rna_data, get_proc_rna_data() , format = "file"),
  
  # Load peak matrix and metadata
  tar_target(counts, data_load(peak_matrix_h5)),
  tar_target(metadata, metadata_load(metadata_csv)),
  
  # Create Chromatin Assay, Seurat object and preprocess object
  tar_target(chrom_assay, create_chrom_assay(counts, frags_file)),
  tar_target(pbmc, 
             create_seurat(chrom_assay, metadata) %>% 
               get_annotation() %>% 
               compute_qc() %>% 
               modify_metadata()
  ),
  
  # Plot some qc
  tar_target(qc_plots, plot_qc(pbmc)),
  tar_target(qc_violins, plot_qc_violins(pbmc)),
  
  # Create final Seurat object
  tar_target(rna_data, load_rna_data(proc_rna_data)),
  tar_target(pbmc_final,
             qc_filter(pbmc) %>% 
               normalize_and_reduce_dims() %>% 
               cluster() %>% 
               add_gene_activity() %>% 
               label_transfer_from_rna_data(rna_data)
  ),
  
  # Plot UMAPS of RNA and ATAC
  tar_target(rna_atac_umaps, plot_rna_atac_umaps(pbmc_final, rna_data))
)