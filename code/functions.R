get_peak_matrix_h5 <- function() {
  system(
    "cd ./data; \
    wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
    )
  fs::dir_ls("data") %>% 
    str_subset(".h5")
}

get_metadata_csv <- function() {
  system(
    "cd ./data; \
    wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_singlecell.csv"
  )
  fs::dir_ls("data") %>% 
    str_subset(".csv")
}

get_frags <- function() {
  system(
    "cd ./data; \
    wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz"
  )
  fs::dir_ls("data") %>% 
    str_subset("fragments.tsv.gz$")
}

get_frags_indx <- function() {
  system(
    "cd ./data; wget wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi"
  )
  fs::dir_ls("data") %>% 
    str_subset("tsv.gz.tbi")
}

get_proc_rna_data <- function() {
  system(
    "cd ./data; https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds"
  )
  fs::dir_ls("data") %>% 
    str_subset("10k_v3.rds")
}

data_load <- function(peak_matrix_h5) {
  counts <- Read10X_h5(peak_matrix_h5)
}

metadata_load <- function(metadata_csv) {
  metadata <- read.csv(metadata_csv, header = TRUE, row.names = 1)
}

create_chrom_assay <- function(counts, frags) {
  chrom_assay <- 
    CreateChromatinAssay(
      counts = counts,
      sep = c(":", "-"),
      genome = "hg19",
      fragments = frags,
      min.cells = 10,
      min.features = 200
    )
}

create_seurat <- function(chrom_assay, metadata) {
  pbmc <- 
    CreateSeuratObject(
      counts = chrom_assay, 
      assay = "peaks", 
      meta.data = metadata
    )
}

get_annotation <- function(pbmc) {
  annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
  GenomeInfoDb::seqlevelsStyle(annotations) <- "UCSC"
  Annotation(pbmc) <- annotations
  pbmc
}

compute_qc <- function(pbmc) {
  pbmc <-NucleosomeSignal(pbmc)
  pbmc <- TSSEnrichment(pbmc, fast = FALSE)
}

modify_metadata <- function(pbmc) {
  pbmc@meta.data <-
    mutate(
      pbmc@meta.data,
      pct_reads_in_peaks = peak_region_fragments / passed_filters * 100,
      blacklist_ratio = blacklist_region_fragments / peak_region_fragments,
      high_tss = if_else(TSS.enrichment > 2, "High", "Low"),
      nucleosome_group = if_else(nucleosome_signal > 4, "NS > 4", "NS < 4")
    )
  pbmc
}

plot_qc <- function(pbmc) {
  qc_plots <-
    list(
      TSSPlot(pbmc, group.by = "high_tss"),
      FragmentHistogram(pbmc, group.by = "nucleosome_group")
    )
}

plot_qc_violins <- function(pbmc) {
  VlnPlot(
    object = pbmc,
    features = c(
      'pct_reads_in_peaks', 
      'peak_region_fragments',
      'TSS.enrichment', 
      'blacklist_ratio', 
      'nucleosome_signal'
    ),
    pt.size = 0.1,
    ncol = 5
  )
}