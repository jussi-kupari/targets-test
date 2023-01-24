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
    "cd ./data; \
    wget https://cf.10xgenomics.com/samples/cell-atac/1.0.1/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz.tbi"
  )
  fs::dir_ls("data") %>% 
    str_subset("tsv.gz.tbi")
}


get_proc_rna_data <- function() {
  system(
    "cd ./data; wget https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds"
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
  pbmc <- NucleosomeSignal(pbmc)
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


qc_filter <- function(pbmc) {
  pbmc <- 
    subset(pbmc,
      peak_region_fragments > 3000 &
        peak_region_fragments < 20000 &
        pct_reads_in_peaks > 15 &
        blacklist_ratio < 0.05 &
        nucleosome_signal < 4 &
        TSS.enrichment > 2
    )
} 


normalize_and_reduce_dims <- 
  function(pbmc) {
    pbmc <-
      RunTFIDF(pbmc) %>% 
      FindTopFeatures(min.cutoff = 'q0') %>% 
      RunSVD()
  } 


cluster <- function(pbmc) {
  pbmc <-
    RunUMAP(pbmc, reduction = 'lsi', dims = 2:30) %>% 
    FindNeighbors(reduction = 'lsi', dims = 2:30) %>% 
    FindClusters(verbose = FALSE, algorithm = 3)
}


add_gene_activity <- 
  function(pbmc) {
    gene.activities <- GeneActivity(pbmc)
    pbmc[['RNA']] <- CreateAssayObject(gene.activities)
    pbmc <- 
      NormalizeData(
        pbmc,
        assay = 'RNA',
        normalization.method = 'LogNormalize',
        scale.factor = median(pbmc$nCount_RNA)
      )
    pbmc
  }


load_rna_data <- 
  function(proc_rna_data) {
    rna_data <-
      readRDS(proc_rna_data) %>% 
      UpdateSeuratObject()
}
  

label_transfer_from_rna_data <- 
  function(pbmc, rna_data) {
    
    # Label transfer
    DefaultAssay(pbmc) <- 'RNA'
    
    transfer.anchors <- 
      FindTransferAnchors(
        reference = rna_data,
        query = pbmc,
        reduction = 'cca'
      )
    
    predicted.labels <- 
      TransferData(
        anchorset = transfer.anchors,
        refdata = rna_data$celltype,
        weight.reduction = pbmc[['lsi']],
        dims = 2:30
      )

    # Add predicted labels to scATAC pbmc metadata
    pbmc <- AddMetaData(pbmc, metadata = predicted.labels)
    
    # Drop cluster 14 (likely low quality cells) and rename clusters
    pbmc <- 
      subset(
        pbmc, 
        idents = 14, 
        invert = TRUE
      ) %>%  
      RenameIdents(
        '0' = 'CD14 Mono',
        '1' = 'CD4 Memory',
        '2' = 'CD8 Effector',
        '3' = 'CD4 Naive',
        '4' = 'CD14 Mono',
        '5' = 'DN T',
        '6' = 'CD8 Naive',
        '7' = 'NK CD56Dim',
        '8' = 'pre-B',
        '9' = 'CD16 Mono',
        '10' = 'pro-B',
        '11' = 'DC',
        '12' = 'NK CD56bright',
        '13' = 'pDC'
      )
    pbmc
  }


plot_rna_atac_umaps <- 
  function(pbmc_final, rna_data) {
  (DimPlot(
    object = rna_data,
    group.by = 'celltype',
    label = TRUE,
    repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')) +
    
    (DimPlot(
      object = pbmc_final,
      group.by = 'predicted.id',
      label = TRUE,
      repel = TRUE) + NoLegend() + ggtitle('scATAC-seq'))
}