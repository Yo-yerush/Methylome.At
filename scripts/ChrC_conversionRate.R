conversionRate <- function(methylationData, var) {
  x = data.frame(a = NA, b = NA)
  
  methylationData = methylationData[grepl("ChrC|NC_000932.1|Pt|chloroplast", methylationData@seqnames)]
  methylation_df = methylationData@elementMetadata
  readsM_pos = grep("readsM",names(methylation_df))
  ii=1
  for (i in readsM_pos) {
    convRate.0 = (methylation_df[,i]) / methylation_df[,i+1] # readsM / readsN
    convRate = (1-mean(convRate.0[!is.na(convRate.0)]))*100
    message(time_msg(), var,":\t",round(convRate,2),"%")
    
    x[ii,] = c(var,round(convRate,2))
    ii=ii+1
  }
  names(x) = c("sample","conversion_rate")
  return(x)
}