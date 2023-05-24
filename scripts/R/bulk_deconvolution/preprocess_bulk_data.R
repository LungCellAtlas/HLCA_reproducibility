rm(list=ls())
gc()

library('biomaRt')
library(stringr)

setwd("~/Documents/20220627_HLCAdeconvolutionReference/bulk_data_as_received/")

################ GET OLIVIA DATA ################
olivia_counts <- read.delim("Olivia/DGE/2020-02-10_Count Data for DGE Analysis.txt")
olivia_info <- read.delim("Olivia/DGE/2020-02-10_Clinical Data for DGE Analysis.txt")

# rownames(olivia_info) == colnames(olivia_counts) # check: all TRUE

# 55 samples have study == olivia, visit == 2 (baseline after ICS wash out)
olivia_info_subset_2 <- olivia_info[which(olivia_info$Study == "olivia" & olivia_info$Visit == "v2"),]
olivia_info_subset_4 <- olivia_info[which(olivia_info$Study == "olivia" & olivia_info$Visit == "v4"),]
olivia_info_subset <- rbind(olivia_info_subset_2, olivia_info_subset_4)
olivia_counts_subset_2 <- olivia_counts[,which(olivia_info$Study == "olivia" & olivia_info$Visit == "v2")]
olivia_counts_subset_4 <- olivia_counts[,which(olivia_info$Study == "olivia" & olivia_info$Visit == "v4")]
olivia_counts_subset <- cbind(olivia_counts_subset_2, olivia_counts_subset_4)

# rownames(olivia_info_subset) == colnames(olivia_counts_subset) # check: all TRUE

# get CPM, then write to file
olivia_CPM <- t(t(olivia_counts_subset) /  colSums(olivia_counts_subset))
write.table(olivia_CPM, file = "../bulk_analysis/nose_bulk.txt", quote = F, sep = "\t")

# write status info to file
olivia_status = data.frame(Sample = rownames(olivia_info_subset), Visit = olivia_info_subset$Visit)
write.table(olivia_status, file = "../bulk_analysis/nose_status.txt", quote = F, sep = "\t")


################ GET INDURAIN DATA ################ using INDURAIN *and* NORM!
indurain_counts <- read.delim("Indurain/NORM & INDURAIN - mRNA Gene Expression (RNA Seq)/1504_VanenBerge.expression.genelevel.v75.htseq.txt", row.names = 1)
indurain_info <- read.delim("Indurain/NORM & INDURAIN - Patients/master_table-RNA-Seq.txt")

# reformat counts data colnames to match info sample names
colnames_indurain = colnames(indurain_counts)
colnames_indurain = substr(colnames_indurain, 2, nchar(colnames_indurain)-10)
colnames_indurain = gsub("_", "-", colnames_indurain)
colnames(indurain_counts) = colnames_indurain

# select count data that is in the samples in info
indurain_counts = indurain_counts[, colnames_indurain %in% indurain_info$rnaseq.id]

# do CPM, then write to file
indurain_CPM <- t(t(indurain_counts) /  colSums(indurain_counts))
write.table(indurain_CPM, file = "../bulk_analysis/airway_bulk_withNORM.txt", quote = F, sep = "\t")

# write status info to file
rownames(indurain_info) = indurain_info$rnaseq.id
indurain_info$asthma_status <- substr(indurain_info$link.id, 1, 4) != 'NORM'
indurain_status = indurain_info[colnames(indurain_CPM),][,c('rnaseq.id', 'currentsmoking', 'asthma_status')]
# indurain_status$rnaseq.id == colnames(indurain_CPM) # check!
write.table(indurain_status, file = "../bulk_analysis/airway_status_withNORM.txt", quote = F, sep = "\t")


################ GET lung tissue db DATA ################

ltdb_counts <- read.delim("lung_tissue_database/merckgeno22oct2014.txt")
ltdb_info <- read.delim("lung_tissue_database/fengeno22oct2014.txt")
ltdb_labels <- read.delim("lung_tissue_database/correct_labels_from_Maarten.txt")

# narrow down to rows with gene info
ltdb_genes <- read.delim("lung_tissue_database/orig.txt")
ltdb_genes_include <- read.delim("lung_tissue_database/list_of_probes_Alen.txt", sep=" ")
ltdb_counts = ltdb_counts[ltdb_genes$Probe %in% ltdb_genes_include$Probe,] 
ltdb_genes = ltdb_genes[ltdb_genes$Probe %in% ltdb_genes_include$Probe,]

# subset to samples in ltdb_labels
ltdb_counts <- ltdb_counts[, colnames(ltdb_counts) %in% ltdb_labels$Sample_ID]

# write to file
ltdb_counts_2 <- cbind(ltdb_genes$HUGO, ltdb_counts)
write.table(ltdb_counts_2, file = "../bulk_analysis/parenchyma_bulk.txt", quote = F, sep = "\t", row.names = F )

# get relevant status & write to file
rownames(ltdb_labels) <- ltdb_labels$Sample_ID
ltdb_labels <- ltdb_labels[colnames(ltdb_counts_2[2:length(colnames(ltdb_counts_2))]),]
ltdb_status <- ltdb_labels[, c("Sample_ID", "fencopddef", "copdgolddeff")]
write.table(ltdb_status, file = "../bulk_analysis/parenchyma_status.txt")
