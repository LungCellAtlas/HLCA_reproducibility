# Look at cibersort results on actual bulk data
rm(list=ls())
gc()

library(ggplot2)
library(tidyr)
library(ggpubr)
library(plyr)
library(dplyr)
library(introdataviz)

setwd("~/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/bulk_analysis/")
tissue = "airway"
withNORM = TRUE # T for airway, F for parenchyma and nose

ctrl = FALSE # v2, Control, or FALSE (for nose, parenchyma and airway, resp.)
other = TRUE # v4, GOLD 3/4, or TRUE

##### proportions ground truth: 
# for now leaving this in to match colours with earlier plots
gtruth <- read.table(paste0("../pseudobulk_analysis/cell_counts_pseudobulk_samples/cell_counts_pseudobulk_", tissue, ".csv"), sep=",", header = T)

count = 0
for (s in unique(gtruth$sample)) {
  subset_gt <- gtruth[gtruth$sample == s,]
  subset_gt$sample_proportions <- subset_gt$sample.1 / sum(subset_gt$sample.1)
  if (count == 0) {
    gtruth_withProportions = subset_gt
  } else {
    gtruth_withProportions = rbind(gtruth_withProportions, subset_gt)
  }
  count = count + 1
} # result: gtruth_withProportions

#### proportions deconvolution:
if (withNORM == TRUE) {
  deconv <- read.table(paste0("CIBERSORTx_output_", tissue, "_withNORM/CIBERSORTx_Results.txt"), sep = "\t", header = T, row.names = 1)
} else {
  deconv <- read.table(paste0("CIBERSORTx_output_", tissue, "/CIBERSORTx_Results.txt"), sep = "\t", header = T, row.names = 1)
}
num_columns = length(colnames(deconv)) -3 # columns to include in pivot longer
deconv$sample = rownames(deconv)
deconv_plottable <- pivot_longer(deconv, 1:num_columns, names_to = "custom_labels", values_to = "sample_proportions")

# prepare for plotting
deconv_plottable$custom_labels <- gsub("\\.", " ", deconv_plottable$custom_labels)
gtruth_withProportions$which_plot <- "Ground truth"
deconv_plottable$which_plot <- "Deconvolution"
names(gtruth_withProportions)[names(gtruth_withProportions) == "custom_label"] <- "custom_labels"
df_to_plot <- rbind(gtruth_withProportions[, c("sample", "custom_labels", "sample_proportions", "which_plot")], 
                    deconv_plottable[, c("sample", "custom_labels", "sample_proportions", "which_plot")])
df_to_plot <- df_to_plot[df_to_plot$which_plot == "Deconvolution",]

if (tissue == "airway") {
  # add status info
  df_status <- read.table(paste0(tissue,'_status_withNORM.txt'))
  status_info = vector()
  
  for (row in seq(df_to_plot[,1])) {
    new_status = df_status[df_to_plot[row,]$sample,]$asthma_status
    status_info = c(status_info , c(new_status))
  }
  df_to_plot$status = status_info
}

if (tissue == "nose") {
  # add status info
  df_status <- read.table(paste0(tissue,'_status.txt'))
  rownames(df_status) = df_status$Sample # only if nose!!
  status_info = vector()
  for (row in seq(df_to_plot[,1])) {
    new_status = df_status[df_to_plot[row,]$sample,]$Visit
    status_info = c(status_info , c(new_status))
  }
  df_to_plot$status = status_info
}

if (tissue == "parenchyma") {
  # add status info
  df_status <- read.table(paste0(tissue,'_status.txt'))
  status_info = vector()
  for (row in seq(df_to_plot[,1])) {
    new_status = df_status[df_to_plot[row,]$sample,]$copdgolddeff # fencopddef or copdgolddeff
    status_info = c(status_info , c(new_status))
  }
  df_to_plot$status = ifelse(status_info %in% c(3,4), "GOLD 3/4", ifelse(status_info == 9, "Control", 
                                                                         ifelse(status_info == 2, "GOLD 2", 1)))
  print(table(df_to_plot$status, status_info))
  df_to_plot = df_to_plot[df_to_plot$status != 1,]
}

if (tissue == "parenchyma") {
  df_to_plot$status = factor(df_to_plot$status, levels =c("Control", "GOLD 2", "GOLD 3/4", 1))
}

# plotting colours
library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# plot stacked bars composition
plot_proportions <- function(df, plot_title, col_vector) {
  ggplot(df, aes(x = factor(sample) ,y = sample_proportions, fill = factor(custom_labels))) + 
    geom_bar(stat = "identity") +
    theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
          legend.title = element_text(size=14), #change legend title font size
          legend.text = element_text(size=8), #change legend text font size
          axis.text.x = element_text(size=8, angle=90)) +
    labs(title=plot_title, x ="Sample", y = "Cell type proportions", fill = "annotation") +
    facet_grid(cols = vars(status), scales="free", space = "free") +
    scale_fill_manual(values = col_vector)
}

plot_proportions(df_to_plot, paste("Deconvolution", tissue), col_vector)#_airway)

### chuck out >60% zeros
accepted_celltypes = c()
chucked_celltypes = data.frame(cell_type = character(), median_difference = numeric(), p = numeric(), adj_p = numeric())
for (celltp in unique(df_to_plot$custom_labels)) {
  check_subset = df_to_plot[df_to_plot$custom_labels == celltp,]
  total_zeros = sum(check_subset$sample_proportions == 0)

  if (total_zeros < 0.6*length(check_subset$sample_proportions)) {
    accepted_celltypes = c(accepted_celltypes, celltp)
  } else {
    print(paste("removing", celltp, total_zeros, "zero values out of", length(check_subset$sample_proportions), "values"))
    
    med_other = median(df_to_plot[df_to_plot$custom_labels == celltp & df_to_plot$status == other,]$sample_proportions)
    med_ctrl = median(df_to_plot[df_to_plot$custom_labels == celltp & df_to_plot$status == ctrl,]$sample_proportions)
    
    exclusion = "excluded: >60% zero values"
    chucked_celltypes = rbind(chucked_celltypes, data.frame(cell_type = celltp, median_difference = med_other - med_ctrl, p = exclusion, adj_p = NA))
  }
}

# remove unreliable cell types
if (tissue == 'airway') {
  df_to_plot = df_to_plot[!(df_to_plot$custom_labels %in% c("EC capillary", "Interstitial macrophages")),]
  exclusion = "excluded: unreliable results in benchmark"
  chucked_celltypes = rbind(chucked_celltypes, data.frame(cell_type = "EC capillary", median_difference = exclusion, p = NA, adj_p = NA))
  chucked_celltypes = rbind(chucked_celltypes, data.frame(cell_type = "Interstitial macrophages", median_difference = exclusion, p = NA, adj_p = NA))
}

#### check differences
wilcoxon_results = data.frame(cell_type = character(), median_difference = numeric(), p = numeric())
for (celltp in unique(df_to_plot[df_to_plot$custom_labels %in% accepted_celltypes,]$custom_labels)) {

  wc <- wilcox.test(df_to_plot[df_to_plot$custom_labels == celltp & df_to_plot$status == other,]$sample_proportions,
                    df_to_plot[df_to_plot$custom_labels == celltp & df_to_plot$status == ctrl,]$sample_proportions)
  
  med_other = median(df_to_plot[df_to_plot$custom_labels == celltp & df_to_plot$status == other,]$sample_proportions)
  med_ctrl = median(df_to_plot[df_to_plot$custom_labels == celltp & df_to_plot$status == ctrl,]$sample_proportions)
  
  wilcoxon_results = rbind(wilcoxon_results, data.frame(cell_type = celltp, median_difference = med_other - med_ctrl, p = wc$p.value))
}

# save results
wilcoxon_results$adj_p = p.adjust(wilcoxon_results$p, method= "fdr")
wilcoxon_results = rbind(wilcoxon_results, chucked_celltypes)
write.table(wilcoxon_results, file = paste0("results_", tissue, ".txt"))

### plot results: split box plot with jitter

# prettify labels
if (tissue == "nose") {
  df_to_plot$status <- ifelse(df_to_plot$status == "v2", "Baseline", "After ICS")
}
if (tissue == "airway") {
  df_to_plot$status <- ifelse(df_to_plot$status == TRUE, "Asthma", "Control")
}
if (tissue == "parenchyma") {
  df_to_plot = df_to_plot[df_to_plot$status %in% c("Control", "GOLD 3/4"),]
}
df_to_plot$custom_labels <- gsub('   ', ' & ', df_to_plot$custom_labels)

# plot
ggplot(df_to_plot, aes(x = custom_labels, y = sample_proportions, fill = status)) + 
  geom_boxplot(outlier.shape=NA) + theme_classic() +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.title = element_text(size=14), 
        legend.text = element_text(size=8),
        axis.text.x = element_text(size=8, angle=90)) +
  labs(x ="Cell type", y = "Predicted proportion", fill = "Asthma status") + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2), size = 0.1, alpha = 1)  +
  scale_x_discrete(labels = function(x) 
    stringr::str_wrap(x, width = 15)) 