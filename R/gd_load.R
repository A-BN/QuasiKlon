#' @author A-BN / A-AB
#' @date 2018/07/10
#' @description gd_load() import a gd (genome diff) file in a data frame
#' @param gd_file path to the gd file
#' @return a data frame
#' @export

library(dplyr)

gd_load <-
	function(gd_file) {
		output_ <- read_delim(gd_file, 
			delim = "\t",
			comment='#',
			col_names=c('type',
				'evidence-id',
				'parent-ids',
				'seq_id',
				'position',
				'a',  # temp column name for type-based calling
				'b')) 
		
		output <- 
			output_ %>%
				filter(type != 'RA', type != 'MC', type != 'JC', type != 'UN') %>%
			
				mutate(new_seq = 
					case_when(
					  	type == 'SNP' ~ a,
					   	type == 'SUB' ~ b,
					   	type == 'DEL' ~ '0'
					 )
				) %>%
				mutate(size = 
					case_when(
						type == 'SNP' ~ '0',
						type == 'SUB' ~ a,
						type == 'DEL' ~ a
					)
				) %>%
				select(-a, -b)
	}
