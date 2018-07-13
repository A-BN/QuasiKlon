#' gd_load() import a gd (genome diff) file in a data frame
#'
#' @param gd_file path to the gd file
#'
#' @return a dataframe
#' @export

gd_load <-
	function(gd_file) {
		gd_colnames <- 
			c('type', 'evidence_id', 'parent_ids', 'seq_id', 'position', 'NewOrSize', 'SizeOrNew')
		gd_coltypes <-
			cols(type = col_character(),
					evidence_id = col_integer(),
			         parent_ids = col_character(),
					seq_id = col_character(),
					position = col_integer(),
					NewOrSize = col_character(),
					SizeOrNew = col_character() 
				)
		unwanted_types <-
			c("MOB","RA", "MC", "JC", "UN")

		gd_in <- readr::read_delim(gd_file, # Raw data frame with few columns called wrongly    
						 delim = "\t",
                         comment ='#',
                         col_names = gd_colnames,
                         col_types = gd_coltypes)
    
		if(nrow(gd_in) == 0){
			gd_in <- data.frame(type = NA_character_, 
                          evidence_id = NA_character_, 
                          parent_ids = NA_character_, 
                          seq_id = NA_character_, 
                          position = NA_integer_, 
                          NewOrSize = NA_character_, 
                          SizeOrNew = NA_character_, 
                          stringsAsFactors = FALSE)
    }
	gd_df <- # Arranged data frame
		gd_in %>%
	        dplyr::filter(! type %in% unwanted_types) %>%
	        dplyr::mutate(new_seq = as.character(dplyr::case_when(
	          type == 'SNP' ~ NewOrSize,
	          type == 'SUB' ~ SizeOrNew,
	          type == 'DEL' ~ 'none',
	          type == 'INS' ~ NewOrSize,
	          TRUE ~ NA_character_))) %>%
	        dplyr::mutate(size = as.integer(dplyr::case_when(
	          type == 'SNP' ~ '0',
	          type == 'SUB' ~ NewOrSize,
	          type == 'DEL' ~ NewOrSize, # Do not work
	          type == 'AMP' ~ NewOrSize,
	          type == 'CON' ~ NewOrSize,
	          type == 'INV' ~ NewOrSize,
	          TRUE ~ NA_character_))) %>%
	        dplyr::select(-NewOrSize, -SizeOrNew)
	  if(nrow(gd_df) == 0) gd_df[1, ] <- NA
	
	return(gd_df)
  }
