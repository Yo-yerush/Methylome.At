edit_TE_file <- function(TE_df) {
  
  TE_df <- TE_df %>%
    mutate(seqnames = NA) %>%  # Add a new column with NA values
    dplyr::select(seqnames,Transposon_min_Start,Transposon_max_End,orientation_is_5prime, everything()) %>%
    dplyr::rename(gene_id = Transposon_Name)
  
  for (i in 1:5) {
    TE_df$seqnames[grep(paste0("AT",i,"TE"),TE_df$gene_id)] = paste0("Chr",i)
  }
  TE_df$orientation_is_5prime = gsub("true","+",TE_df$orientation_is_5prime)
  TE_df$orientation_is_5prime = gsub("false","-",TE_df$orientation_is_5prime)
  
  names(TE_df)[1:4] = c("seqnames","start","end","strand")
  TE_gr = makeGRangesFromDataFrame(TE_df, keep.extra.columns = T)
  
  return(TE_gr)
}