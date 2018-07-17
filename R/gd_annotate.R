#' gff_annotate() returns annotation from a list of mutation from df
#'
#' @param gff_file path to gff3 reference file from breseq/data
#' @param subtracted comes from gd_substract()
#'
#' @return a data frame containing for each mutation a row with anotation
#' @export

gd_annotate <-
	function(gff_file, subtracted) {
		# Cut gff3 file from breseq in two: the second one is the fasta
		gff_total <- read_file(gff_file) %>%
			stringr::str_split(string = ., pattern = "##FASTA")

		write(gff_total[[1]][2], file="/tmp/fasta.delete")

		# Getting metadata from the first part of gff file
		my_ref <- rtracklayer::import.gff3(gff_file)

		# Open fasta
		my_fasta <- Biostrings::readDNAStringSet(filepath = "/tmp/fasta.delete")

		# They will be reinjected in the final df
		na_backup <- subtracted$origin[is.na(subtracted$type)]

		subtracted <-
			subtracted %>%
			ungroup() %>%
			dplyr::filter(!is.na(unique_id))


		# Left join metadata from gff on subtracted data frame
		my_gd <- GenomicRanges::makeGRangesFromDataFrame(df = subtracted, keep.extra.columns = TRUE,
														 start.field = "position", end.field = "position",
														 seqnames.field = "seq_id")

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

		# # Reset size for SNP, INS et SUB
		# my_gd_df <-
		# 	my_gd_df %>%
		# 	dplyr::mutate(size = dplyr::case_when(type_mut == "SNP" ~ as.double(1),
		# 										  type_mut == "INS" | type_mut == "SUB" ~ as.double(nchar(new_seq)),
		# 										  TRUE ~ as.double(size))) %>%
		# 	dplyr::rowwise() %>%
		# 	dplyr::mutate(ref_seq = as.character(XVector::subseq(my_fasta[(names(my_fasta) == seqnames_mut)],
		# 														 start = start_mut,
		# 														 width = size))) %>%
		# 	dplyr::select(ref_seq, new_seq, size, type_mut, everything())

		# Get codon position (inside the contig) and the mutation position (inside the codon)
		# With the codon position we can get the reference sequence, with the mutation position, we exchange
		# the old nucleotide with the new one at the good position. For - strand, we need to inverse the direction
		my_gd_df <-
			my_gd_df %>%
			# Not elegant at all, the sequence is 1-based (position 1), %%3 on each position gives
			# 1 2 0 1 2 0 and (+ 1) 2 3 1 2 (one position decalage), this block is needed to hot fix
			# this. Futhermore, if the strand is negative, you need to inverse codon 1 and 3 (if_else part)
			# In gff3, negative strands have there start_ref and end_ref inversed (start_ref < end_ref)
			# We're looking for an elegant and math-based solution, please help us save this part of the code
			dplyr::mutate(codon_start = dplyr::if_else(condition = (type_mut == "SNP"),
													   true = dplyr::case_when(strand_ref == "+" ~ start_mut - start_ref,
													   				 	strand_ref == "-" ~ end_ref - start_mut,
													   				 	TRUE ~ NA_integer_),
													   false = NA_integer_)) %>%
			dplyr::mutate(codon_pos = dplyr::case_when(codon_start %% 3 == 2 ~ dplyr::if_else(condition = (strand_ref == "+"),
																					   true = 1,
																					   false = 3),
													   codon_start %% 3 == 0 ~ 2,
													   codon_start %% 3 == 1 ~ dplyr::if_else(condition = (strand_ref == "+"),
													   								true = 3,
													   								false = 1),
													   TRUE ~ as.double(NA))) %>% # If not a SNP
			# Update the begining of the codon knowing the position of the mutation inside the codon
			# If you have mutation on pos 2, you need to take 3 nucleotides beginning one position before the mutation
			dplyr::mutate(codon_start = dplyr::case_when(codon_pos == 1 ~ dplyr::if_else(condition = strand_ref == "+",
																				  true = start_mut,
																				  false = as.integer(start_mut - 2)),
														 codon_pos == 2 ~ as.integer(start_mut - 1),
														 codon_pos == 3 ~ dplyr::if_else(condition = strand_ref == "+",
														 						 true = as.integer(start_mut - 2),
														 						 false = start_mut),
														 TRUE ~ NA_integer_)) %>% # If not a SNP
			dplyr::mutate(codon_ref = dplyr::if_else(condition = (type_mut == "SNP"),
													 true = dplyr::if_else(condition = (strand_ref == "+"),
													 			   true = as.character(XVector::subseq(my_fasta[(names(my_fasta) == seqnames_mut)],
													 									start = codon_start,
													 									width = 3)),
													 			   false = IRanges::reverse(as.character(XVector::subseq(my_fasta[(names(my_fasta) == seqnames_mut)],
													 			   									 start = codon_start,
													 			   									 width = 3)))),
													 false = NA_character_)) %>%
			dplyr::mutate(codon_mut = dplyr::if_else(condition = (type_mut == "SNP"),
													 true = stringr::str_replace(string = codon_ref,
													 							pattern =  paste0("(.{", codon_pos - 1, "}).{1}(.{", 3 - codon_pos, "})"),
													 							replacement =  paste0("\\1", new_seq, "\\2")),
													 false = NA_character_))

		my_gd_df <-
			my_gd_df %>%
			dplyr::mutate(aa_ref = ifelse(test = (type_mut == "SNP"), # if_else don't work : it looks like it evaluated the DNAStringSet with NA
										  yes = as.character(
										  	Biostrings::translate(x = Biostrings::DNAStringSet(codon_ref))),
										  no = NA_character_)) %>%
			dplyr::mutate(aa_mut = ifelse(test = (type_mut == "SNP"),
										  yes = as.character(
										  	Biostrings::translate(x = Biostrings::DNAStringSet(codon_mut))),
										  no = NA_character_)) %>%
			dplyr::mutate(is_syn = dplyr::if_else(condition = (type_mut == "SNP"),
												  true = (aa_ref == aa_mut),
												  false = NA)) %>%
			dplyr::select(origin, type_mut, is_syn, aa_ref, aa_mut, codon_ref, codon_mut, codon_pos, ref_seq, new_seq, codon_start, ID, everything())

		for(to_insert in na_backup) {
			my_gd_df[nrow(my_gd_df) + 1, ] <- NA
			my_gd_df$origin[nrow(my_gd_df)] <- to_insert
		}

		return(my_gd_df)
	}
