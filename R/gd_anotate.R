#' @author A-BN / A-AB
#' @date 2018/07/11
#' @description gff_parser() returns annotation from a list of mutation from df
#' @param gff_file is gff3 reference from breseq/data and subtracted comes from the gd_substract() function
#' @return a data frame containing for each mutation a row with anotation
#' @export

library(readr)
library(ape)
library(stringr)
library(tidyr)
library(rtracklayer)
library(fastaUtils)
library(seqinr)
library(Biostrings)
library(XVector) # subseq needs it


gd_anotate <-
	function(gff_file, subtracted) {
		# Cut gff3 file from breseq in two: the second one is the fasta
		gff_total <- read_file(gff_file) %>%
			str_split(string = ., pattern = "##FASTA")
		
		write(gff_total[[1]][2], file="/tmp/fasta.delete")
		
		# Getting metadata from the first part of gff file
		my_ref <- rtracklayer::import.gff3(gff_file)
		
		# Open fasta
		#my_fasta <- seqinr::read.fasta(file = "/tmp/fasta.delete", as.string=TRUE)
		my_fasta <- Biostrings::readDNAStringSet(filepath = "/tmp/fasta.delete")
		
		subtracted <-
			subtracted %>%
				filter(unique_id != "NA|NA|NA")
		
		
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
		my_gd_df <- left_join(x = as.data.frame(my_gd_df), 
							  y = as.data.frame(my_ref), 
							  by = 'ID', 
							  suffix = c("_mut", "_ref"))
		
		# Reset size for SNP, INS et SUB		
		my_gd_df <-
			my_gd_df %>%
				mutate(size = case_when(type_mut == "SNP" ~ as.double(1),
										type_mut == "INS" | type_mut == "SUB" ~ as.double(nchar(new_seq)),
										TRUE ~ as.double(size))) %>%
				rowwise() %>%
				mutate(ref_seq = as.character(subseq(my_fasta[(names(my_fasta) == seqnames_mut)], 
													 start = start_mut, 
													 width = size))) %>%
				select(ref_seq, new_seq, size, type_mut, everything(), -unique_id, -seqnames_ref, -width_mut, -evidence_id, -parent_ids, -score, -phase, -source, -locus_tag, -note)
		
		# Get codon pos
		my_gd_df <-
			my_gd_df %>%
				mutate(codon_start = if_else(condition = (type_mut == "SNP"),
											true = start_mut - start_ref,
											false = NA_integer_)) %>% 
				mutate(codon_pos = case_when(codon_start %% 3 == 2 ~ 1, # In a 1-based world, pos %% 3 + 1
											 codon_start %% 3 == 0 ~ 2, # don't gave us 1 2 3 but 3 1 2,
											 codon_start %% 3 == 1 ~ 3, # this mutate rearrange that (without + 1)
											 TRUE ~ as.double(NA))) %>% # If not a SNP
									  							
				mutate(codon_start = case_when(codon_pos == 1 ~ codon_start,
											   codon_pos == 2 ~ as.integer(codon_start - 1),
											   codon_pos == 3 ~ as.integer(codon_start - 2),
											   TRUE ~ NA_integer_)) %>% # If not a SNP
				mutate(codon_ref = if_else(condition = (type_mut == "SNP"),
										   true = as.character(subseq(my_fasta[(names(my_fasta) == seqnames_mut)], 
										   						   start = start_mut,
										   						   width = 3)),
										   false = NA_character_)) %>%
				mutate(codon_mut = if_else(condition = (type_mut == "SNP"),
										   true = str_replace(string = codon_ref,
											   pattern =  paste0("(.{", codon_pos - 1, "}).{1}(.{", 3 - codon_pos, "})"),
											   replacement =  paste0("\\1", new_seq, "\\2")),
										   false = NA_character_))
		
		#my_gd_df <-
			my_gd_df %>%
				mutate(aa_ref = ifelse(test = (type_mut == "SNP"), # if_else don't work : it looks like it evaluated the DNAStringSet with NA
										yes = as.character(
											Biostrings::translate(x = Biostrings::DNAStringSet(codon_ref))),
										no = NA_character_)) %>%
				mutate(aa_mut = ifelse(test = (type_mut == "SNP"),
										yes = as.character(
											Biostrings::translate(x = Biostrings::DNAStringSet(codon_mut))),
										no = NA_character_)) %>%
				mutate(is_syn = if_else(condition = (type_mut == "SNP"),
										true = (aa_ref == aa_mut),
					   					false = NA)) %>%
				select(origin, gene, product, is_syn, aa_ref, aa_mut, codon_ref, codon_mut, codon_pos, ref_seq, new_seq, size, type_mut, everything())	%>%
				identity()
			
			return(my_gd_df)
		
		# Todo get the type of mutation (syn or not?)
		
	}

gd_list <- list.files(path = "~/Desktop", pattern = ".*.gd", full.names = TRUE)
subtracted <- gd_subtract(gd_merged = gd_merge(gd_list), ref_name = "output_d.gd", filtered_out = TRUE)
gff_file <- '~/Desktop/reference.gff3'

tada <- gd_anotate(gff_file = gff_file, subtracted = subtracted)
tada %>% View()
