#' @author A-BN / A-AB
#' @date 2018/07/10
#' @description gd_subtract() remove mutation present from reference in other gd file an output it as data.frame
#' @param gd_list a list of path to gd file
#' @return a data frame
#' @export

library(dplyr)

gd_subtract <-
	function(gd_list, reference = 0) {
		# Get merged list TODO = as argument?
		gd_df <- gd_merge(gd_list)
		gd_df <-
			gd_df %>% group_by(gd_origin) # Why?
		
		# parcours chaque gd
			# si tu rencontres une mutation (= position + type + new_seq + size ?)
			# qui est aussi presente chez la reference, supprime la
		
		# exporte ce nouveau dataframe (voire output-le)
		
		if(reference == 0) {
			# Don't use it
			reference = gd_df[gd_origin==unique(gd_origin)[1]]
		}
		
		for (gd_file in gd_list) {
			tmp_df <-
				gd_load(gd_file)
			tmp_df <-
				tmp_df %>%
					dplyr::mutate(gd_origin = 
						strsplit(basename(gd_file), '.')[1] # Remove extension (if there is one and path)
					) 
			
			rbind(gd_df, tmp_df)
		}
	}

gd_list = list('R/output.gd', 'R/output_b.gd')
gd_merge(gd_list)
