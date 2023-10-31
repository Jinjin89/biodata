fun_check_list_of_vector <- function(input_list){
  if(is.data.frame(input_list)){
    stop("Input should be a list of vectors, found data.frame-like obj")
  }else if(is.list(input_list)){
    if(all(purrr::map_lgl(input_list,purrr::is_atomic))){
      return(T)
    }else{
      stop("Not all the elements in the list are atomic")
    }
  }else{
    stop("The input is not a list")
  }
}
fun_list2df_each <- function(input_list,colname){
  df = data.frame(row.names = unique(input_list[[1]]))
  df[[colname]] = if(names(input_list) == "") "Unknown" else names(input_list)
  df
}

#' change a list into dataframe format
#'
#' @param input_list the list, each element is is_atomic format
#' @param multiple_items_solver how to deal with items appear multiple time, accept drop and keep
#' @param colname what the colname to show, default is pathway
#'
#' @return df
#' @export
#'
#' @examples
#' fun_list2df(sig_vec_list)
fun_list2df <- function(input_list,multiple_items_solver="keep",
                        colname = "pathway"){
  if(fun_check_list_of_vector(input_list)){
    if(length(input_list) == 1){
      return(fun_list2df_each(input_list,colname = colname))
    }else{
      common_items = table(unname(unlist(input_list)))
      common_items = names(common_items)[common_items > 1]
      # 1) if found common genes
      if(length(common_items) == 0 ){
        message("no common genes found")
        purrr::map_df(seq_along(input_list),function(i){
          fun_list2df_each(input_list[i],colname = colname)
          }) -> df
        return(df)
      }else{
        purrr::map_df(seq_along(input_list),function(i){
          list_new = input_list[i]
          list_new[[1]] = setdiff(list_new[[1]],common_items)
          fun_list2df_each(list_new,colname = colname)
        }) -> df_unique

        # solve the common items
        if(multiple_items_solver == "drop"){
          message(paste0("Drop items appeared in multile categories: ",length(common_items)))
          return(df_unique)
        }else{
          message(paste0("keep items appeared in multile categories: ",length(common_items)))
          # 1) final df
          df = data.frame(row.names = common_items)
          df[[colname]] = NA

          # 2) merge data
          for(each_item in common_items){
            for(each_cate in names(input_list)){
              if(each_item %in% input_list[[each_cate]]){
                item_now = c(df[each_item,colname],each_cate)
                item_now = stats::na.omit(item_now)
                df[each_item,colname] = paste0(item_now,collapse = "&")
              }
            }
          }
          return(rbind(df_unique,df))
        }
      }
    }
  }
}

#' change the data.frame into list
#' @param input_df data.frame
#' @param input_name_regex regex pattern to search in data.frame column data
#' @param name_col which column to search(list names in results)
#' @param value_col the values(list values in results)
#'
#' @return list
#' @export
#'
#' @examples
#'  tmp = fun_df2list(input_df=db_hallmark,input_name_regex = "NF")
fun_df2list <- function(input_df,input_name_regex=NULL,name_col="term",
                        value_col = "gene"){
  if(is.null(input_name_regex)){
    message("translate all data into df")
    db_filtered = input_df
  }else{
    db_filtered =
      input_df %>%
      dplyr::filter(stringr::str_detect(.[[name_col]],input_name_regex))
  }

  name_found = unique(db_filtered[[name_col]])
  if(length(name_found) == 0){
    warning("No items found!")
    return(NULL)
  }else{
    message(paste0("Found: ",length(name_found)))
    found_list = split(db_filtered[[value_col]],db_filtered[[name_col]])
    purrr::map(found_list,unique)
  }
}


#' convert gene data into symbol
#'
#' @param input_genes input_genes
#' @param from from type
#' @param alias_not_found_keep_or_drop whether keep notfound data
#'
#' @return data.frame
#' @export
#'
#' @examples
#' fun_gene2symbol(c("ACTB","CT",'afdfs',"MT1"))
fun_gene2symbol <- function(input_genes,from = "alias",
                            # alias params
                            alias_not_found_keep_or_drop = c("drop","keep")

                            # entrez params

                            # engs params
                            ){
  message("This script use `hgnc_complete_set_2023-10-01.txt`\n",
          "which stores 43736 genes ",
          "accept: alias,entrezid,ensg")
  stopifnot("genes should be unique" = (length(input_genes) == length(unique(input_genes))))
  from = match.arg(from,c("alias","entrezid","ensg"))

  # 1) check alias
  if(from == "alias"){
    message("conert genes from alias into symbol")
    # 0) params prepare
    alias_not_found_keep_or_drop = match.arg(alias_not_found_keep_or_drop,c("drop","keep"))
    gene_info_found = input_genes %>% dplyr::intersect(gene_alias2symbol$alias)
    gene_info_lost = input_genes %>% dplyr::setdiff(gene_alias2symbol$alias)

    res = data.frame(
      row.names = input_genes,
      alias = input_genes) %>%
      mutate(symbol = NA_character_)
    res[gene_info_found,"symbol"] = gene_alias2symbol[gene_info_found,"symbol"]
    if(alias_not_found_keep_or_drop == "drop"){
      res[gene_info_lost,"symbol"] = NA_character_
    }else{
      res[gene_info_lost,"symbol"] = gene_info_lost
    }
    message(
      "Keep or drop for alias with symbol: ",alias_not_found_keep_or_drop,
      "\nThe input genes list: ",length(input_genes),
      "\nHow many alias not found in database: ",length(gene_info_lost),
      "\nThe symbols found: ", length(unique(res$symbol)),
      "\nof which, symbol with multiple map: ", sum(table(res$symbol)>1),
      "\nand, symbol with single map: ", sum(table(res$symbol)==1)

    )
    res %>%
      dplyr::mutate(alias_is_symbol = .$alias %in% gene_set$symbol) %>%
      dplyr::arrange(dplyr::desc(alias_is_symbol)) %>%
      return()

  }else if(from == "entrezid"){
    stop("TODO")

  }else if(from == "ensg"){
    stop("TODO")

  }else{
    stop("wrong input")
  }
  # 2)

}
