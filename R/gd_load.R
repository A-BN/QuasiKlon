#' @author A-BN / A-AB
#' @date 2018/07/10
#' @description gd_load() import a gd (genome diff) file in a data frame
#' @param gd_file path to the gd file
#' @return a data frame
#' @export

library(dplyr)

gd_load <-
	function(gd_file) {
		gd_tmp <- read_delim(gd_file, # Raw data frame with few columns called wrongly
			delim = "\t",
			comment='#',
			col_names=c('type',
				'evidence-id',
				'parent-ids',
				'seq_id',
				'position',
				'a',  # temp column name for rename column name
				'b',
				'c')) 
		
		gd_df <- # Arranged data frame
			gd_tmp %>%
				filter(type != 'RA', type != 'MC', type != 'JC', type != 'UN') %>%
				mutate(new_seq = 
					case_when(
					  	type == 'SNP' ~ a,
					   	type == 'SUB' ~ b,
					   	type == 'DEL' ~ 'NA',
					  	type == 'INS' ~ a
					 )
				) %>%
				mutate(size = 
					case_when(
						type == 'SNP' ~ 'NA',
						type == 'SUB' ~ a,
						type == 'DEL' ~ a,
						type == 'MOB' ~ c,
						type == 'AMP' ~ a,
						type == 'CON' ~ a,
						type == 'INV' ~ a
					)
				) %>%
				mutate(repeat_name = if_else(type == 'MOB', a, 'NA')) %>%
				mutate(strand = if_else(type == 'MOB', b, 'NA')) %>%
				mutate(size = if_else(type == 'MOB', b, 'NA')) %>%
				mutate(new_copy_number = if_else(type == 'AMP', b, 'NA')) %>%
				mutate(region = if_else(type == 'region', b, 'NA')) %>%
				select(-a, -b, -c)
		
		return(gd_df)
	}
