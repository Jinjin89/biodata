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
  # 1) get MHC
  mhc_string = c("B2M MHC-classical, class-I TAP1 MHC-classical, class-I TAP2 MHC-classical, class-I TAPBP MHC-classical, class-I HLA-A MHC-classical, class-I HLA-B MHC-classical, class-I HLA-C MHC-classical, class-I HLA-DPA1 MHC-classical, class-II HLA-DPB1 MHC-classical, class-II HLA-DQA1 MHC-classical, class-II HLA-DQA2 MHC-classical, class-II HLA-DQB1 MHC-classical, class-II HLA-DRA MHC-classical, class-II HLA-DRB1 MHC-classical, class-II HLA-E MHC-non-class, class-I HLA-F MHC-non-class, class-I HLA-G MHC-non-class, class-I HLA-DMA MHC-non-class, class-II HLA-DMB MHC-non-class, class-II HLA-DOA MHC-non-class, class-II HLA-DOB MHC-non-class, class-II")
  mhc_string =  stringr::str_split(mhc_string,"MHC-") %>%
    unlist() %>%
    purrr::map_chr(\(x){
      stringr::str_remove_all(x,"non-class, class-I(I){0,1}|classical, class-I(I){0,1}| ")
    })
  mhc_string = mhc_string[mhc_string!= ""]
  sig_vec_list$MHC =mhc_string
  meta$sig_vec_list$MHC = "Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade"

  #2) simmulator genes
  si_string = "MICA Immunostimulator MICB Immunostimulator CD27 Immunostimulator CD274 Immunostimulator CD28 Immunostimulator CD40 Immunostimulator CD40LG Immunostimulator CD70 Immunostimulator CD80 Immunostimulator CD86 Immunostimulator ICOS Immunostimulator ICOSLG Immunostimulator IL6 Immunostimulator IL6R Immunostimulator PDCD1LG2 Immunostimulator TMEM173 Immunostimulator TNFRSF13B Immunostimulator TNFRSF13C Immunostimulator TNFRSF14 Immunostimulator TNFRSF17 Immunostimulator TNFRSF18 Immunostimulator TNFRSF4 Immunostimulator TNFRSF9 Immunostimulator TNFSF13 Immunostimulator TNFSF13B Immunostimulator TNFSF18 Immunostimulator TNFSF4 Immunostimulator TNFSF9 Immunostimulator TNFSF15 Immunostimulator TNFRSF25 Immunostimulator HHLA2 Immunostimulator TMIGD2 Immunostimulator C10orf54 Immunostimulator BTNL2 Immunostimulator CD276 Immunostimulator CD48 Immunostimulator TNFSF14 Immunostimulator TNFRSF8 Immunostimulator CD70 Immunostimulator PVR Immunostimulator LTA Immunostimulator MICA Immunostimulator MICB Immunostimulator IL2RA Immunostimulator ENTPD1 Immunostimulator NT5E Immunostimulator CXCR4 Immunostimulator CXCL12 Immunostimulator KLRK1 Immunostimulator NKG2A Immunostimulator RAET1E Immunostimulator ULBP1 Immunostimulator"
  si_string = stringr::str_split(si_string,"Immunostimulator") %>%
    unlist() %>%
    purrr::map(\(x)stringr::str_remove_all(x," ")) %>%
    unlist() %>%
    {.[. != ""]}

  sig_vec_list$Immunostimulator =si_string
  meta$sig_vec_list$Immunostimulator = "Pan-cancer Immunogenomic Analyses Reveal Genotype-Immunophenotype Relationships and Predictors of Response to Checkpoint Blockade"

  # 3) chemokines
  sig_vec_list$Chemokine =
    "CCL1	CCL11	CCL13	CCL14	CCL15	CCL16	CCL17	CCL18	CCL19	CCL2	CCL20	CCL21	CCL22	CCL23	CCL24	CCL25	CCL26	CCL27	CCL28	CCL3	CCL4	CCL5	CCL7	CCL8	CX3CL1	CXCL1	CXCL10	CXCL11	CXCL12	CXCL13	CXCL14	CXCL16	CXCL17	CXCL2	CXCL3	CXCL5	CXCL6	CXCL8	CXCL9" %>%
    str_split("\t") %>%
    unlist()

  sig_vec_list$Chemokine_receptor <-
    "CCR1	CCR10	CCR2	CCR3	CCR4	CCR5	CCR6	CCR7	CCR8	CCR9	CXCR1	CXCR2	CXCR3	CXCR4	CXCR5	CXCR6	CX3CR1" %>%
    str_split("\t") %>%
    unlist()
  meta$sig_vec_list$Chemokine = "Siglec15 shapes a non-inflamed tumor microenvironment and predicts the molecular subtype in bladder cancer"
  meta$sig_vec_list$Chemokine_receptor = "Siglec15 shapes a non-inflamed tumor microenvironment and predicts the molecular subtype in bladder cancer"


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

# tcga immune cellai

if(F){
  tcga_immuneCellAI = readRDS("~/data/project/db/tcga/others/immuneCellAI.rds")
  use_data(tcga_immuneCellAI)
}

# tcga tip data

if(F){
  all_files <- list.files("~/data/tmp/download/",full.names = T)
  purrr::map(all_files,\(each_file){
    fread(each_file)%>%
      as.data.frame() %>%
      set_rownames(.[[1]]) %>%
      select(-1) %>%
      t %>%
      as.data.frame() %>%
      mutate(sample = str_sub(rownames(.),1,15)) %>%
      filter(!duplicated(sample)) %>%
      set_rownames(.$sample)
  }) %>%
    list_c() ->tcga_tip
}

# gene data

if(F){
  if(F){
    # 1) get all the data into data.frame
    df <- fread("~/data/tmp/hgnc_complete_set_2023-10-01.txt") %>%
      mutate(alias = paste(symbol,alias_symbol,prev_symbol,sep = "|")) %>%
      mutate(alias = str_remove(alias,"^\\|")) %>%
      mutate(alias = str_remove(alias,"\\|$"))


    # 2) get Genes column
    gene_set <- tibble::tibble(symbol = df$symbol,
                               hgnc_id = df$hgnc_id,
                               name = df$name,
                               entrez_id = df$entrez_id,
                               ensembl = df$ensembl_gene_id,
                               refseq = df$refseq_accession,
                               alias = purrr::map(str_split(df$alias,"\\|"),\(x) unique(x) %>% setdiff("")),
                               gene_family = df$gene_family,

                               # meta data
                               locus_type = df$locus_type,
                               location = df$location
    )
    gene_set %<>%
      as.data.frame() %>%
      set_rownames(.$symbol)
    use_data(gene_set,overwrite = T)

    genes_used <- gene_set %>%
      dplyr::select(symbol,alias) %>%
      as.data.frame() %>%
      set_rownames(.$symbol)

    tictoc::tic()
    alias2symbol <- purrr::map(genes_used$symbol,\(each_symbol){
      current_alias <- genes_used[each_symbol,"alias"] %>% unlist()
      data.frame(alias = current_alias) %>%
        mutate(symbol = each_symbol)
    }) %>%
      purrr::list_rbind()
    gene_alias2symbol = alias2symbol
    tictoc::toc()

    # save data
    all_symbols = gene_alias2symbol$symbol %>% unique
    symbol_test <-
      gene_alias2symbol %>%
      group_by(alias) %>%
      summarise(symbol_new = paste0(symbol,collapse = "/")) %>%
      ungroup() %>%
      mutate(symbol_new = ifelse(alias %in% all_symbols,alias,symbol_new)) %>%
      as.data.frame() %>%
      set_rownames(.$alias) %>%
      set_colnames(c("alias","symbol"))
    all_symbols %>% setdiff(symbol_test %>% filter(alias %in% symbol) %>% pull(symbol))

    gene_alias2symbol <- symbol_test
    use_data(gene_alias2symbol,overwrite = T)

  }

}

# pancancer Immune metagenes
