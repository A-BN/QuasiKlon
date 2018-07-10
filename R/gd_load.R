#' @author A-BN / A-AB
#' @date 2018/07/10
#' @description gd_load() import a gd (genome diff) file in a data frame
#' @param gd_file path to the gd file
#' @return a data frame
#' @export

library(dplyr)
library(readr)

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

gd_load <-
  function(gd_file) {
    gd_in <- read_delim(gd_file, # Raw data frame with few columns called wrongly
                         delim = "\t",
                         comment='#',
                         col_names= gd_colnames,
                         col_types = gd_coltypes) 
   gd_df <- # Arranged data frame
      gd_in %>%
        dplyr::filter(! type %in% unwanted_types) %>%
      
        dplyr::mutate(new_seq = case_when(
          type=='SNP' ~ NewOrSize,
          type=='SUB' ~ SizeOrNew,
          type=='DEL' ~ 'N',
          type=='INS' ~ NewOrSize)) %>%
        dplyr::mutate(size = as.integer(case_when(
          type == 'SNP' ~ '0',
          type == 'SUB' ~ NewOrSize,
          type == 'DEL' ~ NewOrSize, # Do not work
          type == 'AMP' ~ NewOrSize,
          type == 'CON' ~ NewOrSize,
          type == 'INV' ~ NewOrSize))) %>%
        dplyr::select(-NewOrSize, -SizeOrNew)
    return(gd_df)
  }

gd_file <- "R/output.gd"

test <- gd_load(gd_file = gd_file)
