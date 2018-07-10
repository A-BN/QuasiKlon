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
				'c')
		) 
		
		gd_df <- data.frame(type=character(),
			evidence_id=integer(),
			parent_ids=character(), # character because the comma is for spliting
			seq_id=character(),
			position=integer(),
			new_seq=character(),
			size=character(),
			repeat_name=character(),
			strand=integer(),
			region=character(),
			new_copy_number=integer()
		)
		
		gd_df <- # Arranged data frame
			gd_tmp %>%
				dplyr::filter(type != 'RA', type != 'MC', type != 'JC', type != 'UN') %>%
			
				dplyr::mutate(new_seq = case_when(
					type=='SNP' ~ a,
					type=='SUB' ~ b,
					type=='DEL' ~ '0',
					type=='INS' ~ a)) %>%

				dplyr::mutate(size = case_when(
					type == 'SNP' ~ '0',
					type == 'SUB' ~ a,
					type == 'DEL' ~ a, # Do not work
					type == 'MOB' ~ c,
					type == 'AMP' ~ a,
					type == 'CON' ~ a,
					type == 'INV' ~ a)) %>%
			
				dplyr::mutate(repeat_name = if_else(type == 'MOB', a, '0')) %>%
			
				dplyr::mutate(strand = if_else(type == 'MOB', b, '0')) %>%
			
				dplyr::mutate(new_copy_number = if_else(type == 'AMP', b, '0')) %>%
			
				dplyr::mutate(region = if_else(type == 'MOB', b, '0')) %>%

				dplyr::select(-a, -b, -c)
		
		return(gd_df)
	}

gd_file="R/output.gd"

test = gd_load(gd_file = gd_file)
