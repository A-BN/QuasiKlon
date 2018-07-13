#' gff_parser() returns annotation from a list of mutation from df
#'
#' @param gff_file path to gff3 reference file from breseq/data
#' @param subtracted come from gd_substract()
#'
#' @return a data frame containing for each mutation a row with anotation
#' @export

gd_anotate <-
	function(gff_file, subtracted) {
		# Cut gff3 file from breseq in two: the second one is the fasta
		gff_total <- read_file(gff_file) %>%
			stringr::str_split(string = ., pattern = "##FASTA")
		
		write(gff_total[[1]][2], file="/tmp/fasta.delete")
		
		# Getting metadata from the first part of gff file
		my_ref <- rtracklayer::import.gff3(gff_file)
		
		# Open fasta
		my_fasta <- Biostrings::readDNAStringSet(filepath = "/tmp/fasta.delete")
		
		subtracted <-
			subtracted %>%
				dplyr::filter(unique_id != "NA|NA|NA")
		
		
		# Left join metadata from gff on subtracted data frame
		my_gd <- makeGRangesFromDataFrame(df = subtracted, keep.extra.columns = TRUE, 
			start.field = "position", end.field = "position", seqnames.field = "seq_id")
		
		nearest_feat <- IRanges::nearest(x = my_gd, subject = my_ref, ignore.strand = TRUE)
		dist_feat <- IRanges::distanceToNearest(x = my_gd, subject = my_ref, ignore.strand = TRUE)
		my_gd_df <- as.data.frame(my_gd)
		my_gd_df$dist_to_feat <- 0
		my_gd_df$dist_to_feat[queryHits(dist_feat)] <- elementMetadata(dist_feat)$distance
		my_gd_df <- cbind(my_gd_df, ID = my_ref@elementMetadata[nearest_feat, c("ID")])
		
		## Preparing the result dataframe
		my_gd_df <- dplyr::left_join(x = as.data.frame(my_gd_df), 
							  y = as.data.frame(my_ref), 
							  by = 'ID', 
							  suffix = c("_mut", "_ref"))
		
		# Reset size for SNP, INS et SUB		
		my_gd_df <-
			my_gd_df %>%
				dplyr::mutate(size = dplyr::case_when(type_mut == "SNP" ~ as.double(1),
										type_mut == "INS" | type_mut == "SUB" ~ as.double(nchar(new_seq)),
										TRUE ~ as.double(size))) %>%
				dplyr::rowwise() %>%
				dplyr::mutate(ref_seq = as.character(XVector::subseq(my_fasta[(names(my_fasta) == seqnames_mut)], 
													 start = start_mut, 
													 width = size))) %>%
				dplyr::select(ref_seq, new_seq, size, type_mut, everything(), -unique_id, -seqnames_ref, -width_mut, -evidence_id, -parent_ids, -score, -phase, -source, -locus_tag, -note)
		
		# Get codon pos
		my_gd_df <-
			my_gd_df %>%
				dplyr::mutate(codon_start = dplyr::if_else(condition = (type_mut == "SNP"),
											true = start_mut - start_ref,
											false = NA_integer_)) %>% 
				dplyr::mutate(codon_pos = dplyr::case_when(codon_start %% 3 == 2 ~ 1, # In a 1-based world, pos %% 3 + 1
											 codon_start %% 3 == 0 ~ 2, # don't gave us 1 2 3 but 3 1 2,
											 codon_start %% 3 == 1 ~ 3, # this mutate rearrange that (without + 1)
											 TRUE ~ as.double(NA))) %>% # If not a SNP
									  							
				dplyr::mutate(codon_start = dplyr::case_when(codon_pos == 1 ~ codon_start,
											   codon_pos == 2 ~ as.integer(codon_start - 1),
											   codon_pos == 3 ~ as.integer(codon_start - 2),
											   TRUE ~ NA_integer_)) %>% # If not a SNP
				dplyr::mutate(codon_ref = dplyr::if_else(condition = (type_mut == "SNP"),
										   true = as.character(XVector::subseq(my_fasta[(names(my_fasta) == seqnames_mut)], 
										   						   start = start_mut,
										   						   width = 3)),
										   false = NA_character_)) %>%
				dplyr::mutate(codon_mut = dplyr::if_else(condition = (type_mut == "SNP"),
										   true = stringr::str_replace(string = codon_ref,
											   pattern =  paste0("(.{", codon_pos - 1, "}).{1}(.{", 3 - codon_pos, "})"),
											   replacement =  paste0("\\1", new_seq, "\\2")),
										   false = NA_character_))
		
		my_gd_df <-
			my_gd_df %>%
				dplyr::mutate(aa_ref = dplyr::ifelse(test = (type_mut == "SNP"), # if_else don't work : it looks like it evaluated the DNAStringSet with NA
										yes = as.character(
											Biostrings::translate(x = Biostrings::DNAStringSet(codon_ref))),
										no = NA_character_)) %>%
				dplyr::mutate(aa_mut = dplyr::ifelse(test = (type_mut == "SNP"),
										yes = as.character(
											Biostrings::translate(x = Biostrings::DNAStringSet(codon_mut))),
										no = NA_character_)) %>%
				dplyr::mutate(is_syn = dplyr::if_else(condition = (type_mut == "SNP"),
										true = (aa_ref == aa_mut),
					   					false = NA)) %>%
				dplyr::select(origin, gene, product, type_mut, is_syn, aa_ref, aa_mut, codon_ref, codon_mut, codon_pos, ref_seq, new_seq, size, everything()) %>%
			
			return(my_gd_df)
	}
