trimm_and_rename <- function(gr_obj) {
  
  #### change 'seqnames' to 'TAIR10' if needed
  if (any(grepl("NC_0", seqlevels(gr_obj)))) {
    # change 'TAIR10.1'
    gr_obj = renameSeqlevels(gr_obj, gsub("NC_003070.9","Chr1",seqlevels(gr_obj)))
    gr_obj = renameSeqlevels(gr_obj, gsub("NC_003071.7","Chr2",seqlevels(gr_obj)))
    gr_obj = renameSeqlevels(gr_obj, gsub("NC_003074.8","Chr3",seqlevels(gr_obj)))
    gr_obj = renameSeqlevels(gr_obj, gsub("NC_003075.7","Chr4",seqlevels(gr_obj)))
    gr_obj = renameSeqlevels(gr_obj, gsub("NC_003076.8","Chr5",seqlevels(gr_obj)))
    gr_obj = renameSeqlevels(gr_obj, gsub("NC_037304.1","ChrM",seqlevels(gr_obj))) # MT
    gr_obj = renameSeqlevels(gr_obj, gsub("NC_000932.1","ChrC",seqlevels(gr_obj))) # Pltd
    names(mcols(gr_obj)) = gsub("gene_biotype","gene_model_type", names(mcols(gr_obj)))
    names(mcols(gr_obj)) = gsub("locus_tag","gene_id", names(mcols(gr_obj)))
    
  } else if (all(seqlevels(gr_obj)[1:5] == as.character(1:5)) & any(grepl("Pt", seqlevels(gr_obj)))) {
    # change 'Ensembl'
    gr_obj = renameSeqlevels(gr_obj, gsub("Mt","M",seqlevels(gr_obj))) # MT
    gr_obj = renameSeqlevels(gr_obj, gsub("Pt","C",seqlevels(gr_obj))) # Pltd
    gr_obj = renameSeqlevels(gr_obj, paste0("Chr",seqlevels(gr_obj))) # Pltd
    names(mcols(gr_obj)) = gsub("biotype","gene_model_type", names(mcols(gr_obj)))
  }
  
  ### trimm 'ChrM' and 'ChrC' from the genome
  remove_seqnames = c("ChrM","ChrC")
  gr_obj <- gr_obj[!seqnames(gr_obj) %in% remove_seqnames]
  seqlevels(gr_obj) = as.character(unique(seqnames(gr_obj))) #setdiff(seqlevels(gr_obj), remove_seqnames)
  
  return(gr_obj)
}
