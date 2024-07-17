
# biodata

<!-- badges: start -->
<!-- badges: end -->

The goal of biodata is to ...

## Installation

You can install the development version of biodata from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Jinjin89/biodata")
```

# data

## tcga-data

data are collected from

1. xena
2. cbioportal
3. tcia

## gene database data

data are collected from
1. go
2. kegg
3. hallmark

## genes info

[hgnc gene data](https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/archive/monthly/tsv/hgnc_complete_set_2021-03-01.txt)

gene info

  * name (contains alias,symbol,previous symbol)
  * symbol
  * ensg
  * entrez
  * locus group: "protein-coding gene" "non-coding RNA"      "pseudogene"          "other" 
  * location
  



## signatures 

* sig_df_list: storing data in data.frame/matrix-like format

* sig_vec_list: stroing data in vector-like format


# meta

* meta_sig_vec_list
