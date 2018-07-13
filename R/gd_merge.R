#' gd_merge() call gd_read() on gd files and merge them in an unique data frame with an origin column
#'
#' @param gd_list a list of path to gd file
#'
#' @return a data frame containing each mutation from each gd files indexed by a new origin col
#' @export

gd_merge <-
	function(gd_list) {
		gd_merge_df <- list()
		for (i in 1:length(gd_list)) {
			curr_name <- 
			stringr::str_replace(string = gd_list[[i]], pattern = ".*/(.*)\\.gd", replacement = "\\1")
			gd_merge_df[[i]] <- gd_load(gd_list[[i]])
			gd_merge_df[[i]]$origin <- curr_name
		}
		gd_merge_df <- do.call(what = rbind, args = gd_merge_df)
		return(gd_merge_df)
	}
