#' @author A-BN / A-AB
#' @date 2018/07/10
#' @description gd_merge() call gd_read() on gd files and merge them in an unique data frame with a gd_origin column
#' @param gd_list a list of path to gd file
#' @return a data frame
#' @export

library(dplyr)
library(stringr)
gd_merge <-
  function(gd_list) {
    gd_merge_df <- list()
    for (i in 1:length(gd_list)) {
      curr_name <- 
        str_replace(string = gd_list[[i]], pattern = ".*/(.*)\\.gd", replacement = "\\1")
      gd_merge_df[[i]] <- gd_load(gd_list[[i]])
      gd_merge_df[[i]]$origin <- curr_name
    }
    gd_merge_df <- do.call(what = rbind, args = gd_merge_df)
    return(gd_merge_df)
  }


gd_list <- c('example_gd/output.gd', 'example_gd/output_b.gd', 'example_gd/output_c.gd')
bla <- gd_merge(gd_list)
