gainORloss = function(DMRsReplicates, context) {
  
  gainORloss_vec = DMRsReplicates$regionType
  
  gain_DMRs = length(grep("gain",gainORloss_vec))
  loss_DMRs = length(grep("loss",gainORloss_vec))
  total = gain_DMRs+loss_DMRs
  pres_gain = round((gain_DMRs/total)*100, 1)
  pres_loss = round((loss_DMRs/total)*100, 1)
  
  pie_data = data.frame(group = c("", paste0("Loss (",pres_loss,"%)"), paste0("Gain (",pres_gain,"%)")),
                        value = c(0,loss_DMRs, gain_DMRs))
  
  svg(file = paste0("pie_",context,"_gainORloss.svg"), width = 4.81, height = 3.93, family = "serif")
  pie_plot = pie(pie_data$value, pie_data$group, main = paste0(total," DMRs in ",context," context"), radius = 1,clockwise = T)
  #plot(pie_plot)
  pie_plot
  dev.off()
}

########################################################

ratio.distribution = function(DMRsReplicates, var1, var2, context, comparison_name) {
  
  dmrs_ratio = DMRsReplicates$proportion2/DMRsReplicates$proportion1
  
  svg(paste0("ratio.distribution_",context,"_gainORloss.svg"), width = 3.58, height = 3.3, family = "serif")
  
  h <- hist(log(dmrs_ratio), breaks=1000, plot=FALSE)
  hist(log(dmrs_ratio),
       main = paste0(gsub("_"," ",comparison_name)," - ",context),
       xlab = "Ratio of methylation proportions (log scale)",
       ylab = "Frequency",
       #xlim = c(-10,10),
       #col = c(rep("blue",length(h$counts)/2),rep("red",length(h$counts)/2)),
       #border = c(rep("blue",length(h$counts)/2),rep("red",length(h$counts)/2)),
       border = ifelse((h$mids > 0), "#b36b74", "#9396c2"),
       breaks = 1000)
  # if((length(h$breaks) %% 2) != 0) {bar_colors = append(bar_colors,"red",after = T)} 
  
  dev.off()
}