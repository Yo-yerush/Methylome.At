chromosome_plot <- function(var1_trimmed,var2_trimmed,var3_trimmed=NULL,var1_name,var2_name=NULL,var3_name=NULL,
                            chr.n,cntx,difference=F,y.new.scale=F,y_cntx=0.5) {
  chr.vec = var1_trimmed[var1_trimmed@seqnames == chr.n]
  regions <- GRanges(seqnames = Rle(chr.n), ranges = IRanges(chr.vec@ranges@start[1],chr.vec@ranges@start[length(chr.vec@ranges@start)]))
  
  profile_var1 <- computeMethylationProfile(var1_trimmed,
                                            regions,
                                            windowSize = 150000,
                                            context = cntx)
  
  profile_var2 <- computeMethylationProfile(var2_trimmed,
                                            regions,
                                            windowSize = 150000,
                                            context = cntx)
  if (difference) {
    methylationProfiles = profile_var1[,0]
    methylationProfiles$Proportion = profile_var2$Proportion - profile_var1$Proportion
    methylationProfiles <- GRangesList("var1" = methylationProfiles)
    legend_text = paste0(var2_name," vs ",var1_name)
    
    if (chr.n == 1) {ylab = "methylation difference"} else {ylab=''}
    col = "#3d3333"
    pch = rep(26,length(methylationProfiles))
    lty = rep(1,length(methylationProfiles))
    col_legend = rep("#ffffff",length(methylationProfiles))
    lty_legend = rep(0,length(methylationProfiles))
    
  } else {
    if (!is.null(var3_trimmed)) {
      profile_var3 <- computeMethylationProfile(var3_trimmed,
                                                regions,
                                                windowSize = 150000,
                                                context = cntx)
      methylationProfiles <- GRangesList("var1" = profile_var1, "var2" = profile_var2, "var3" = profile_var3)
      legend_text = c(var1_name,var2_name,var3_name)
    } else {
      methylationProfiles <- GRangesList("var1" = profile_var1, "var2" = profile_var2)
      legend_text = c(var1_name,var2_name)
    }
    
    if (chr.n == 1) {ylab = "methylation"} else {ylab=''}
    col = c("gray50","#bf6828","#0072B2")
    pch = rep(26,length(methylationProfiles))
    lty = rep(1,length(methylationProfiles))
    col_legend = col[1:length(methylationProfiles)]
    lty_legend = lty[1:length(methylationProfiles)]
  }
  
  if (y.new.scale) {
    if (difference) {
      ymax = max(methylationProfiles[[1]]$Proportion)*1.1
      ymin = min(methylationProfiles[[1]]$Proportion)*1.1
      if (ymax < y_cntx) {ymax = y_cntx}
      if (ymin > -y_cntx) {ymin = -y_cntx}
    } else {
      ymax = y_cntx
      ymin = 0
    } 
  } else {
    if (difference) {
      ymax = 1
      ymin = -1
    } else {
      ymax = 1
      ymin = 0
    }
  }
  
  pos = (start(methylationProfiles[[1]]) + end(methylationProfiles[[1]]))/2
  plot(pos, methylationProfiles[[1]]$Proportion, type = "o", 
       ylim = c(ymin, ymax), xlab = "genomic coordinate", ylab = ylab, 
       col = col[1], pch = pch[1], lty = lty[1], 
       yaxt = "n", xaxt = "n", 
       
       main = ""
  )
  if (length(methylationProfiles) > 1) {
    for (i in 2:length(methylationProfiles)) {
      lines(pos, methylationProfiles[[i]]$Proportion, 
            type = "o", col = col[i], lty = lty[i], pch = pch[i])
    }
  }
  
  mtext(paste0("Chr ",chr.n), side = 1, line = 0, cex = 1)#adj = 0,
}

###################################################################

ChrPlots_CX <- function(comparison_name, meth_var1, meth_var2, var1, var2, scripts_dir, n.cores, y_cntx_diff=0.1) {
  #### Low resolution profiles plot #### 

  #### change seqnames for plot
  meth_var1_trimmed = renameSeqlevels(meth_var1, gsub("Chr","",seqlevels(meth_var1)))
  meth_var2_trimmed = renameSeqlevels(meth_var2, gsub("Chr","",seqlevels(meth_var2)))

  
  par_function <- function(cntx) {
    #### two conditions plot
    if (cntx == "CG") {
      y_cntx = 1
    } else if (cntx == "CHG") {
      y_cntx = 0.5
    } else if (cntx == "CHH") {
      y_cntx = 0.2
    }
    svg(paste0("ChrPlot_",cntx,"_",comparison_name,".svg"), width = 8, height = 2, family = "serif")
    par(mar = c(1,4,2,0))
    par(fig=c(0,2,0,10)/10)
    plot(runif(10),runif(10),xlim=c(0, 0.01), 
         ylim=c(0,y_cntx),axes=FALSE,type="n",ylab=paste(cntx," methylation"),xlab="")
    axis(2, c(0, y_cntx/2 ,y_cntx), lty = 1, labels = c(0, y_cntx/2 ,y_cntx))
    par(new=T)
    
    i=1
    u = (10-i)/6
    for (chr_number in 1:5) {
      par(mar = c(1,0,2,0))
      par(fig=c(i,i+u,0,10)/10)
      chromosome_plot(meth_var1_trimmed,meth_var2_trimmed,NULL,var1,var2,NULL,chr_number,cntx,difference=F,y.new.scale=T,y_cntx=y_cntx)
      par(new=T)
      i=i+u
    }
    par(mar = c(1,0,2,0))
    par(fig=c(8,10,0,10)/10)
    legend("topright", legend = c(var1,var2), 
           lty = 0, col = c("gray50","#bf6828"), 
           pch = 16, bty = "n")
    dev.off()

    #### difference (delta) plot
    svg(paste0("ChrPlot_difference_",cntx,"_",comparison_name,".svg"), width = 8, height = 2, family = "serif")
    par(mar = c(1,4,2,0))
    par(fig=c(0,2,0,10)/10)
    plot(runif(10),runif(10),xlim=c(0, 0.01), 
         ylim=c(-y_cntx_diff,y_cntx_diff),axes=FALSE,type="n",ylab=paste(cntx," methylation (Δ)"),xlab="")
    axis(2, c(-y_cntx_diff, 0, y_cntx_diff), lty = 1, labels = c(-y_cntx_diff, 0, y_cntx_diff))
    par(new=T)
    
    i=1
    u = (10-i)/6
    for (chr_number in 1:5) {
      par(mar = c(1,0,2,0))
      par(fig=c(i,i+u,0,10)/10)
      chromosome_plot(meth_var1_trimmed,meth_var2_trimmed,NULL,var1,var2,NULL,chr_number,cntx,
                      difference=T,y.new.scale=T,y_cntx=y_cntx_diff)
      par(new=T)
      i=i+u
    }
    par(mar = c(1,0,2,0))
    par(fig=c(8,10,0,10)/10)
    legend("topright", legend = c(paste0(var2," vs ",var1)), 
           lty = 0, col = "#3d3333", 
           pch = 16, bty = "n")
    dev.off()
  }
  
  if (n.cores > 2) {n.cores.plot = 3} else {n.cores.plot = n.cores}
  mclapply(c("CG","CHG","CHH"), par_function, mc.cores = n.cores.plot)
}