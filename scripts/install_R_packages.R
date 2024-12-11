pkg_name = c("dplyr","tidyr","ggplot2","lattice","PeakSegDisk","geomtextpath","parallel","BiocManager","circlize")
pkg_biocond = c("DMRcaller","rtracklayer","topGO","KEGGREST","Rgraphviz","org.At.tair.db","GenomicFeatures","plyranges")

install.packages(pkg_name, repos = "https://mirror.howtolearnalanguage.info/")
BiocManager::install(pkg_biocond)

pkg = .packages(all.available = TRUE)
message("\n")
for (i in c(pkg_name,pkg_biocond)) {
  tryCatch({
    if (pkg[grep(paste0("^",i,"$"), pkg)] == i) {
      message(paste0("*\tinstall ",i,": yes"))
    }
  }, error = function(cond) {message(paste0("*\tinstall ",i,": no"))})
}

