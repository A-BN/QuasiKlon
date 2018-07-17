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

		 # If mutation is in more than half gd but not present in ref, we consider the mutation to be in the ref
		 #
		 invert_thre <- length(unique(gd_merged$origin)) / 2
		 gd_merged <-
			 gd_merged %>%
			 	dplyr::select(-evidence_id, -parent_ids) %>%
			 	dplyr::mutate(unique_id = dplyr::if_else(condition = (!is.na(type) & !is.na(seq_id) & !is.na(position)),
			 											 true = paste(type, seq_id, position, sep = "|"),
			 											 false = as.character(NA))) %>%
		 		# Next line return true if the mutation is present in ref
			 	dplyr::mutate(in_ref = unique_id %in% unique_id[origin == ref_name]) %>%
			 	dplyr::group_by(unique_id) %>%
			 	dplyr::mutate(mutated_ref = dplyr::if_else(condition = !(in_ref) & n() > invert_thre,
				                                 true = TRUE,
				                                 false = FALSE))

	if(! filtered_out) return(gd_merged)



	gd_out <-
		gd_merged %>%
			dplyr::filter(! in_ref) %>%
			dplyr::mutate(origin = dplyr::if_else(condition = mutated_ref, true = ref_name, false = origin)) %>%
			dplyr::distinct() %>% # An 'inverted' mutation is present n times w/ n being the number of gd it was present in.
			dplyr::select(-in_ref, -mutated_ref)

	# absent_ori_df <- as_tibble(t(c(rep(NA,5),unique(gd_merged$origin)[!(unique(gd_merged$origin) %in% gd_out$origin)], NA)))
	# names(absent_ori_df) <- names(gd_out)
	# base::rbind(gd_out, absent_ori_df) %>% View

	for(origin_in_merged in unique(gd_merged$origin)) {
		if(! origin_in_merged %in% unique(gd_out$origin)) {
			gd_out[nrow(gd_out) + 1, ] <- NA
			gd_out$origin[nrow(gd_out)] <- origin_in_merged
		}
	}

	return(gd_out)
	}
