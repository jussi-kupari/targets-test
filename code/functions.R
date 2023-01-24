
get_data <- function(raw_data) {
  counts <- Read10X_h5(raw_data)
}

get_metadata <- function(metadata) {
  meta_data <- 
    read.csv(metadata, header = TRUE, row.names = 1)
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
  pbmc
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
  pbmc
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













