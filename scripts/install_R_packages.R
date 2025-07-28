# Update 'textshaping' version to '0.4.1' if required
if (packageVersion("textshaping") < "0.4.1") {
  tryCatch({
    suppressMessages(devtools::install_github("r-lib/textshaping", ref = "v0.4.1", quiet = T))
    message(paste0("Update 'textshaping' to v", packageVersion("textshaping")))
  }, error = function(e) {
    message("Error installing 'textshaping'")
  })
}

pkg_name = c("dplyr", "tidyr", "ggplot2", "lattice", "PeakSegDisk", "geomtextpath", "parallel", "BiocManager", "RColorBrewer", "circlize", "cowplot", "knitr")
pkg_biocond = c("DMRcaller","rtracklayer","topGO","KEGGREST","Rgraphviz","org.At.tair.db","GenomicFeatures","plyranges")

# install packages if required
installed_packages = rownames(installed.packages())
i=1
for (pkg in pkg_name) {
  if (!(pkg %in% installed_packages)) {
    try({
      suppressMessages(install.packages(pkg, repos = "http://cran.r-project.org", quiet = T))
    }, silent = T)
  }
  perc_val <- (i / length(c(pkg_name, pkg_biocond))) * 100
  cat(paste0("\rInstall R packages...\t", round(perc_val, 1), "% "))
  i=i+1
}
    
for (pkg in pkg_biocond) {
  if (!(pkg %in% installed_packages)) {
    try({
      suppressMessages(BiocManager::install(pkg, quiet = T, update = F, ask = F))
    }, silent = T)
  }
  perc_val <- (i / length(c(pkg_name, pkg_biocond))) * 100
  cat(paste0("\rInstall R packages...\t", round(perc_val, 1), "% "))
  i=i+1
}


cat("\n\n")

# Check each package if installed
c.pkg = .packages(all.available = TRUE)
message("\n")
for (i in c("textshaping", pkg_name, pkg_biocond)) {
  tryCatch(
    {
      if (c.pkg[grep(paste0("^", i, "$"), c.pkg)] == i) {
        cat(paste0("*\tinstalled ", i, ": yes\n"))
        message(paste0("*\tinstalled ", i, ": yes"))
      }
    },
    error = function(cond) {
      cat(paste0("*\tinstalled ", i, ": no\n"))
      message(paste0("*\tinstalled ", i, ": no"))
    }
  )
}