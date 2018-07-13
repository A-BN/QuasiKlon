#' Perform neighbour joining (perform a distance matrix and sample matrix, then use APE)
#'
#' @param annotated_mutations result from gd_annotate()
#'
#' @return a phylo object
#' @export

tree_generator <-
	function(annotated_mutations) {
		## Matrix where row = sample and col = mut with TRUE if mutation appears in a tupple row/mut
		m_muts <- matrix(nrow = length(unique(annotated_mutations$origin)),
						 ncol = length(unique(annotated_mutations$ID)))
		rownames(m_muts) <- unique(annotated_mutations$origin)[1:length(unique(annotated_mutations$origin))]
		colnames(m_muts) <- unique(annotated_mutations$ID)[1:length(unique(annotated_mutations$ID))]
		
		for(sample in unique(annotated_mutations$origin)) 
			for(mut in unique(annotated_mutations$ID))
				m_muts[sample,mut] <- sum(annotated_mutations$ID[annotated_mutations$origin==sample]==mut) > 0
		
		## Matrix between each sample with 1 point if there is one difference of mutation
		widthus <- length(unique(annotated_mutations$origin))
		m_subklon <- matrix(nrow = widthus, ncol = widthus)
		rownames(m_subklon) <- unique(annotated_mutations$origin)[1:widthus]
		colnames(m_subklon) <- unique(annotated_mutations$origin)[1:widthus]
		
		for(rowus in 1:widthus)
			for(colus in 1:widthus)
				m_subklon[rowus,colus] <- sum(m_muts[rowus,] != m_muts[colus,])
		
		phylo_data <- ape::nj(m_subklon)
		plot(phylo_data)
		
		return(phylo_data)
	}
