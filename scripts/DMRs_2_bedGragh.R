DMRs_2_bedGragh <- function(var1,var2,context,ann.res) {
  
  if (ann.res == "all") {
    DMRsReplicates = read.csv(paste0("DMRs_",context,"_",var2,"_vs_",var1,".csv"))
    file_name = paste(context,"all",var2,"vs",var1, sep = "_")
    
  } else {
    DMRsReplicates = read.csv(paste0("genome_annotation/",context,"/",ann.res))
    ann.name <- sub(paste0("_",context,".*$"), "", ann.res)
    file_name = paste(context,ann.name,var2,"vs",var1, sep = "_")
  }
  
  # keep genome location and methylation proportion
  DMRsReplicates_bed = DMRsReplicates[,c("seqnames","start","end")]
  DMRsReplicates_bed$proportion = round(log2(DMRsReplicates$proportion2/DMRsReplicates$proportion1),2)
  
  # save as bedGragh file
  names(DMRsReplicates_bed) = c("chrom", "chromStart", "chromEnd", "count")
  bedgragh_file = PeakSegDisk::writeBedGraph(DMRsReplicates_bed, paste0("DMRs_bedGragh/",file_name,".bedGragh"))
}