#' gd_subtract() remove mutation present from reference in other gd file an output it as data frame
#'
#' @param gd_merged a dataframe from gd_merged()
#' @param ref_name name (with .gd) of the reference gd file
#' @param filtered_out remove row in ref. Default = TRUE
#'
#' @return a dataframe
#' @export

gd_subtract <-
	function(gd_merged, ref_name, filtered_out = TRUE) {
		 ref_name <- stringr::str_remove(string = ref_name, pattern = "\\.gd")
		 invert_thre <- length(unique(gd_merged$origin)) / 2
		 gd_merged <-
		 gd_merged %>%
		 	dplyr::select(-evidence_id, -parent_ids) %>%
		 	dplyr::mutate(unique_id = paste(type, seq_id, position, sep = "|")) %>%
		 	dplyr::mutate(in_ref = unique_id %in% unique_id[origin == ref_name]) %>%
		 	dplyr::group_by(unique_id) %>%
		 	dplyr::mutate(mutated_ref = dplyr::if_else(condition = !(in_ref) & n() > invert_thre, 
			                                 true = TRUE, 
			                                 false = FALSE)) %>%
		 	dplyr::select(-unique_id)
	 
	 if(! filtered_out) return(gd_merged)
	 
	gd_out <- 
		gd_merged %>%
			dplyr::filter(! in_ref) %>%
			dplyr::mutate(origin = dplyr::if_else(condition = mutated_ref, true = ref_name, false = origin)) %>%
			dplyr::distinct() %>% # An 'inverted' mutation is present n times w/ n being the number of gd it was present in.
			dplyr::select(-in_ref, -mutated_ref)
	
	return(gd_out)   
	}
