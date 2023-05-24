# process results CIBERSORT on pseudobulk
# Compare cibersort pseudobulk results to ground truth.

rm(list=ls())
gc()

library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(RColorBrewer)
library(ggrepel)

setwd("~/Documents/20220627_HLCAdeconvolutionReference/restart_CIBERSORT_analyses/pseudobulk_analysis/")

tissues = c("airway", "nose", "parenchyma")

for (tissue in tissues) {
  ##### get proportions ground truth:
  gtruth <- read.table(paste0("cell_counts_pseudobulk_samples/cell_counts_pseudobulk_", tissue, ".csv"), sep=",", header = T)
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
  
  # get proportions deconvolution:
  deconv <- read.table(paste0("CIBERSORTx_output_", tissue, "/CIBERSORTx_Results.txt"), sep = "\t", header = T, row.names = 1)
  num_columns = length(colnames(deconv)) -3 # columns to include in pivot longer
  deconv$sample = rownames(deconv)
  deconv_plottable <- pivot_longer(deconv, 1:num_columns, names_to = "custom_labels", values_to = "sample_proportions")
  
  # prepare to plot
  deconv_plottable$custom_labels <- gsub("\\.\\.\\.", " & ", deconv_plottable$custom_labels) # match labels to ground truth formatting
  deconv_plottable$custom_labels <- gsub("\\.", " ", deconv_plottable$custom_labels) # match labels to ground truth formatting
  deconv_plottable$custom_labels <- gsub("Hillock like", "Hillock-like", deconv_plottable$custom_labels) # match labels to ground truth formatting
  gtruth_withProportions$which_plot <- "Ground truth"
  deconv_plottable$which_plot <- "Deconvolution"
  names(gtruth_withProportions)[names(gtruth_withProportions) == "custom_label"] <- "custom_labels"
  df_to_plot <- rbind(gtruth_withProportions[, c("sample", "custom_labels", "sample_proportions", "which_plot")], 
                      deconv_plottable[, c("sample", "custom_labels", "sample_proportions", "which_plot")])
  
  # colours to plot
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
  
  # plot ground truth and deconvolution proportions (stacked bar chart)
  plot_proportions <- function(df, plot_title) {
    ggplot(df, aes(x = factor(sample) ,y = sample_proportions, fill = factor(custom_labels))) + 
      geom_bar(stat = "identity") +
      theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
            legend.title = element_text(size=14), #change legend title font size
            legend.text = element_text(size=8), #change legend text font size
            axis.text.x = element_text(size=8, angle=90)) +
      labs(title=plot_title, x ="Sample", y = "Cell type proportions", fill = "annotation") +
      facet_wrap(~which_plot, nrow =2) + scale_fill_manual(values = col_vector)
  }
  print(plot_proportions(df_to_plot, paste("Pseudobulk composition", tissue)))
  
  # get proportional error
  gtruth_withProportions$sample_label <- paste0(gtruth_withProportions$sample, "_", gtruth_withProportions$custom_labels)
  deconv_plottable$sample_label <- paste0(deconv_plottable$sample, "_", deconv_plottable$custom_labels)
  prop_error_df <- full_join(gtruth_withProportions, deconv_plottable, by = "sample_label")
  prop_error_df <- prop_error_df[!(is.na(prop_error_df$P.value)),]
  prop_error_df <- prop_error_df[, c("sample.x", "custom_labels.x", "sample_proportions.x", "sample_proportions.y")]
  prop_error_df$prop_error <- (prop_error_df$sample_proportions.x - prop_error_df$sample_proportions.y) / prop_error_df$sample_proportions.x
  prop_error_df$tissue <- tissue
  
  # combine proportional error data over the multiple tissues
  if (tissue == tissue[1]) {
    prop_error_all_tissues = prop_error_df
    prop_error_all_tissues_backup = prop_error_df
  } else {
    prop_error_all_tissues = rbind(prop_error_all_tissues, prop_error_df)
    prop_error_all_tissues_backup = rbind(prop_error_all_tissues, prop_error_df)
  }
  
  # get correlations of deconvoluted cell type proportions with ground truth
  df_corr_results = data.frame('cell type' = character(), 'estimate' = numeric(), 'p-value' = numeric(),
                               'avg % ground truth' = numeric(), 'stdev % ground truth' = numeric())
  
  for (celltp in unique(gtruth_withProportions$custom_labels)) {
    truth = gtruth_withProportions[gtruth_withProportions$custom_labels == celltp,c("sample", "sample_proportions")]
    pred = deconv_plottable[deconv_plottable$custom_labels == celltp,c("sample", "sample_proportions")]
    
    corr_df <- full_join(truth, pred, by = "sample")
    corr_df[["sample_proportions.y"]][is.na(corr_df[["sample_proportions.y"]])] <- 0
    corr_test <- cor.test(corr_df$sample_proportions.x, corr_df$sample_proportions.y)
    
    res <- data.frame(celltp, corr_test$estimate, corr_test$p.value, mean(corr_df$sample_proportions.x), sd(corr_df$sample_proportions.x))
    df_corr_results <- rbind(df_corr_results, res)
  }
  df_corr_results$tissue = tissue
  
  # combine correlation data over the multiple tissues
  if (tissue == tissues[1]) {
    all_dfcorr = df_corr_results
  } else {
    all_dfcorr = rbind(all_dfcorr, df_corr_results)
  }
  
}

# plot prop error per cell type per tissue
ggplot(prop_error_all_tissues, aes(x = prop_error, y = sample_proportions.x, color = custom_labels.x)) + geom_point() +
  labs(title=paste0("Proportional error"), x ="Error", y = "Ground truth proportion", color = "annotation") +
  scale_color_manual(values = col_vector) +
  facet_wrap(~tissue)+xlim(-10, NA)

# plot correlation data, per tissue
ggplot(all_dfcorr, aes(x = corr_test.estimate, y=mean.corr_df.sample_proportions.x., label = celltp, color = ifelse(corr_test.estimate < 0.3, F, T))) + 
  geom_point(size=3, alpha = 0.7) +
  labs(x="Correlation estimate", y = "Avgerage ground truth sample proportion", color = "Inclusion analyses") + theme_classic() +
  facet_wrap(~tissue) + geom_text_repel(mapping = aes(text.size = 0.01))

# inspect & save results
all_dfcorr[all_dfcorr$corr_test.estimate < 0.5,][,c(1,2,6)]
write.table(all_dfcorr, file = 'results_correlations.txt')
