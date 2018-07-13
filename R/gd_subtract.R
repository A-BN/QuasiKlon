#' @author A-BN / A-AB
#' @date 2018/07/10
#' @description gd_subtract() remove mutation present from reference in other gd file an output it as data.frame
#' @param gd_list a list of path to gd file
#' @return a data frame
#' @export

library(dplyr)
library(tidyr)

gd_subtract <-
	function(gd_merged, ref_name, filtered_out = TRUE) {
	  ref_name <- stringr::str_remove(string = ref_name, pattern = "\\.gd")
	  invert_thre <- length(unique(gd_merged$origin)) / 2
	  gd_merged <-
	  gd_merged %>%
	    # select(-evidence_id, -parent_ids) %>%
	    mutate(unique_id = paste(type, seq_id, position, sep = "|")) %>%
	    mutate(in_ref = unique_id %in% unique_id[origin == ref_name]) %>%
	    group_by(unique_id) %>%
	    mutate(mutated_ref = if_else(condition = !(in_ref) & n() > invert_thre, 
	                                 true = TRUE, 
	                                 false = FALSE)) %>%
	    select(-unique_id)
	  
	  if(! filtered_out) return(gd_merged)
	  
    gd_out <- 
      gd_merged %>%
      filter(! in_ref) %>%
      mutate(origin = if_else(condition = mutated_ref, true = ref_name, false = origin)) %>%
      distinct() %>% # An 'inverted' mutation is present n times w/ n being the number of gd it was present in.
      select(-in_ref, -mutated_ref)
   return(gd_out)   
	}

