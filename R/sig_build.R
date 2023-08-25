#modify sig_vec_list data
if(F){
  meta$sig_vec_list = list()
}

if(F){
  # 1) T cell inflammaed
  sig_vec_list$T_cell_inflamed = c(
    "IRF1", "CD8A", "CCL2", "CCL3", "CCL4", "CXCL9", "CXCL10", "ICOS",
    "GZMK", "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB")
  meta$sig_vec_list$T_cell_inflamed = "Melanoma-intrinsic β-catenin signalling prevents anti-tumour immunity"

  # 2) IFNγ-related gene signature
  sig_vec_list$IFN_Gamma_related_gene = c(
    "CD8A", "CCL5", "CD27", "CD274", "PDCD1LG2", "CD276", "CMKLR1",
    "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PSMB10", "STAT1", "TIGIT" )
  meta$sig_vec_list$IFN_Gamma_related_gene = "IFNgamma-related mRNA profile predicts clinical response to PD-1 blockade"

  # 3) T effector signature
  sig_vec_list$T_effector = c("GZMA", "GZMB", "PRF1", "EOMES", "IFNG", "TNF", "CXCL9", "CXCL10", "CD8A", "CD4", "FOXP3", "ICOS", "CTLA4")
  meta$sig_vec_list$T_effector = "Predictive correlates of response to the anti-PD-L1 antibody MPDL3280A in cancer patients"

  # 4) Immune_cytolytic_activity
  sig_vec_list$Immune_cytolytic_activity = c("GZMA", "PRF")
  meta$sig_vec_list$Immune_cytolytic_activity = "Molecular and genetic properties of tumors associated with local immune cytolytic activity"

  # 5)Martinez_Gordon_M1/2
  sig_vec_list$Martinez_Gordon_M1 = c("CD64", "IDO1", "SOCS1", "CXCL10")
  sig_vec_list$Martinez_Gordon_M2 = c("MRC1", "TGM2", "CD23", "CCL22")
  meta$sig_vec_list$Martinez_Gordon_M1 = "The M1 and M2 paradigm of macrophage activation: time for reassessment"
  meta$sig_vec_list$Martinez_Gordon_M2 = "The M1 and M2 paradigm of macrophage activation: time for reassessment"

  # 6)
  sig_vec_list$Murray_M1 = c("IL23A", "IDO1","PTGS2","COX2","IL12B","NOS2","SOCS3")
  sig_vec_list$Murray_M2 = c("KLF4", "CCL24", "CCL12", "CXCL13","CHIA","IRF4","SOCS2","RETNLB","CHI3L1","CHI3L2","CHI3L3")

  meta$sig_vec_list$Murray_M1 = "The M1 and M2 paradigm of macrophage activation: time for reassessment"
  meta$sig_vec_list$Murray_M2 = "The M1 and M2 paradigm of macrophage activation: time for reassessment"

}


if(F){


  db_kegg_metabolism = readRDS("~/data/project/db/genes/kegg_metabolism.rds")
  use_data(db_kegg_metabolism)

  db_reactome_metabolism = readRDS("~/data/project/db/genes/reactome_metabolism.rds")
  use_data(db_reactome_metabolism)
}



# tcga data
if(F){
  tcga_tide_tidepy = data.table::fread("~/data/tcga_tide_tidepy.csv") %>%
    as.data.frame() %>%
    set_rownames(.$sample)
  use_data(tcga_tide_tidepy)

  tcga_tide = readRDS("~/data/tmp_tide.rds")
  use_data(tcga_tide,overwrite = T)

}
