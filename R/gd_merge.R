#' @author A-BN / A-AB
#' @date 2018/07/10
#' @description gd_merge() call gd_read() on gd files and merge them in an unique data frame with a gd_origin column
#' @param gd_list a list of path to gd file
#' @return a data frame
#' @export

library(dplyr)

gd_merge <-
	function(gd_list) {
		# gd_df <- data.frame(
		# 		type=character(),
		# 		evidence_id=integer(),
		# 		parent_ids=character(), # character because the comma is for spliting
		# 		seq_id=character(),
		# 		position=integer(),
		# 		new_seq=character(),
		# 		size=character(),
		# 		repeat_name=character(),
		# 		strand=character(),
		# 		new_copy_number=character(),
		# 		region=character()
		# 	)
		
		tmp = list()
		
		for (gd_file in gd_list) {
			tmp_df <- gd_load(gd_file)
			tmp_df <-
				tmp_df %>%
					dplyr::mutate(gd_origin = 
						strsplit(basename(gd_file), '.')[1] # Remove extension (if there is one and path)
					) 
			# tmp <- list(c(unlist(tmp, recursive = F), tmp_df))
			
			# rbind(gd_df, tmp_df)
		}
		return(tmp)
	}

bla <- gd_merge(list('R/output.gd', 'R/output_b.gd'))
