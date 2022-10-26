setwd("/Users/lisa.sikkema/Documents/LungCellAtlas/scib_in_r_updated/")

# OPTION 1: plot all methods:
source("./plotSingleTaskRNA.R")
## test data:
# plotSingleTaskRNA(csv_metrics_path = "./data/metrics_RNA_allTasks.csv")
# real data
plotSingleTaskRNA(csv_metrics_path = "./metrics_scgen_added.csv")

# OPTION 2: plot only best-performing pre-processing for every output:
source("./plotSingleTaskRNA_best_preprocessing_only.R")
plotSingleTaskRNA_best_preprocessing_only(csv_metrics_path = "./metrics_scgen_added.csv", summary_scores_file_path = "./lung_atlas_fixed_summary_scores.csv")
