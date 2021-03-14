# All data sets should have format: pts-rows, fts-cols
# Save as csv, one per dataset subtype and one per dataset
path='/Users/senosam/Documents/Massion_lab/data_integration/models'
###############################################################################
# CyTOF
###############################################################################
## 1. Frequencies table
# Main cell types frequencies
# Normalized immune subtype frequencies 
# Normalized cluster frequencies 
## 2. Median bulk protein expression


source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/51_supervised_analysis_viz.R")
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/20_ClustAnnot_functions.R")

# Load data
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/cellclusters.RData")
ref <- ref[-42,] # sample 13376 has been collected twice and will be merged, only one ref is needed
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/percent_pt.RData")
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/protein_correlations.RData")

# function to edit percent tables
rc_prcnt <- function(dt, col_rm){
  dt[,col_rm]<- NULL
  tt <- apply(dt, 1, sum)
  dt <- dt/tt
}

# Get frequencies

# Immune
prcnt_imm <- rc_prcnt(prcnt_subtypes, col_rm = c('Epithelial', 'Endothelial', 'Fib_Mesenchymal'))
# Epi clusters
col_rm <- colnames(prcnt_clusters)[-c(grep("ECC", colnames(prcnt_clusters)))]
prcnt_epi <- rc_prcnt(prcnt_clusters, col_rm = col_rm)
# Endo clusters
col_rm <- colnames(prcnt_clusters)[-c(grep("endo", colnames(prcnt_clusters)))]
prcnt_endo <- rc_prcnt(prcnt_clusters, col_rm = col_rm)
# Fmes clusters
col_rm <- colnames(prcnt_clusters)[-c(grep("fmes", colnames(prcnt_clusters)))]
prcnt_fmes <- rc_prcnt(prcnt_clusters, col_rm = col_rm)
# Tc clusters
col_rm <- colnames(prcnt_clusters)[-c(grep("Tc", colnames(prcnt_clusters)))]
prcnt_tc <- rc_prcnt(prcnt_clusters, col_rm = col_rm)
# Th clusters
col_rm <- colnames(prcnt_clusters)[-c(grep("Th", colnames(prcnt_clusters)))]
prcnt_th <- rc_prcnt(prcnt_clusters, col_rm = col_rm)
# DNT clusters
col_rm <- colnames(prcnt_clusters)[-c(grep("DNT", colnames(prcnt_clusters)))]
prcnt_dnt <- rc_prcnt(prcnt_clusters, col_rm = col_rm)
# Mye clusters
col_rm <- colnames(prcnt_clusters)[-c(grep("Mye", colnames(prcnt_clusters)))]
prcnt_mye <- rc_prcnt(prcnt_clusters, col_rm = col_rm)
# Other clusters
col_rm <- colnames(prcnt_clusters)[-c(grep("Other", colnames(prcnt_clusters)))]
prcnt_oimm <- rc_prcnt(prcnt_clusters, col_rm = col_rm)

freq_cytof <- cbind(prcnt_celltypes, prcnt_imm, prcnt_epi, prcnt_endo, prcnt_fmes, 
    prcnt_tc, prcnt_th, prcnt_dnt, prcnt_mye, prcnt_oimm)
write.csv(freq_cytof, file.path(path, 'data', 'cytof_freq.csv'))


# Median bulk protein expression
ft_cols <- c(15,17:31, 33:35, 37:51)
df_cluster <- denoisingCTF::t_asinh(cbind(annot_df[,ft_cols], 
                                          'pt_ID'= annot_df[,'pt_ID']))

cl_median <- aggregate(df_cluster[,-ncol(df_cluster)], list(df_cluster$pt_ID), median)
cl_median <- cbind(cl_median[match(CDE$pt_ID, cl_median$Group.1),], 'n_op2'=CDE$n_op2, 'SILA'=CDE$SILA)
rownames(cl_median)<- cl_median$Group.1
cl_median <- cl_median[,-1]
colnames(cl_median) <- gsub(".*_",  "",colnames(cl_median)) 
write.csv(cl_median, file.path(path, 'data', 'cytof_medianprot.csv'))


###############################################################################
# RNA Seq
###############################################################################
## 1. Top expressed genes (25K) vst normed
## 2. REACTOME eigen genes
## 3. TF activity



###############################################################################
# WES 
###############################################################################
## 1. Binary matrix

###############################################################################
# Radiomics
###############################################################################
## 1. Healthmyne



###############################################################################
# Clinical data
###############################################################################

## read all datasets and determine which patients have complete data
## Select patients that have complete data (n=57)
## split 75% training 25% test (43, 14) (by pt_id)

## SMOTE for training (?)



## DATASETS (csv: 

# cytof_freq
# cytof_exp
# cytof_all
# rna_top
# rna_reac
# rna_viper
# rna_all
# wes_bin
# rad_hm
# all (cytof_all, rna_reac(?), wes_bin, rad_hm)

## ALGORITHMS (regression)
# PLS
# elastic net
# random forest
# xgboost


## Script plan
# Create a list csv file names/path
# Per file:
    # create directory
    # ALG 1
        # train
        # test
        # evaluate metrics
        # feature importance, latent factors
        # Save RData
    # ALG 2 (...)

# One folder for all datasets (csv files)
# One folder per data set to save RData per algorithm
# One R markdown per data set (one master)


