# All data sets should have format: pts-rows, fts-cols
# Save as csv, one per dataset subtype and one per dataset

###############################################################################
# CyTOF
###############################################################################
## 1. Frequencies table
# Main cell types frequencies
# Normalized immune subtype frequencies 
# Normalized cluster frequencies 
## 2. Median bulk protein expression

path='/Users/senosam/Documents/Massion_lab/data_integration/models'
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
## 1. Top expressed genes (20K) vst normed
## 2. REACTOME eigen genes
## 3. TF activity
path='/Users/senosam/Documents/Massion_lab/data_integration/models'

source('/Users/senosam/Documents/Repositories/Research/data_analysis_rnaseq/R/30_DEGanalysis.R')
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/51_supervised_analysis_viz.R")
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/20_ClustAnnot_functions.R")
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData")
environment_set()

ls_preprocessed <- preprocess_rna(path_rnaseq = '/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData', 
                                  correct_batch = F, lowvargenesrm = T, prot_coding_only = F, xychr_rm = T)

CDE <- read.csv(file = '/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2021JAN13_SA_MF.csv')
pData_rnaseq <- CDE[match(ls_preprocessed$pData_rnaseq$pt_ID, CDE$pt_ID),]
ls_preprocessed$pData_rnaseq <- pData_rnaseq

## 1. Top expressed genes (20K) vst normed
vsd_mat <- t(ls_preprocessed$vsd_mat)
rownames(vsd_mat) <- as.character(ls_preprocessed$p_all$pt_ID)

genes_info <- ls_preprocessed$rna_all[,1:7]

write.csv(vsd_mat, file.path(path, 'data', 'rna_topgenes.csv'))
write.csv(genes_info, file.path(path, 'data', 'rna_genesinfo.csv'))


## 2. REACTOME eigen genes
library(biomaRt)
# Reactome pathways

reapaths <- read.table('/Users/senosam/Documents/Repositories/Research/DH_project01/data/Reactome2020/NCBI2Reactome_All_Levels.txt')
genetopaths <- read.delim('/Users/senosam/Documents/Repositories/Research/DH_project01/data/Reactome2020/Ensembl2Reactome_PE_Pathway.txt', header=F, stringsAsFactors = F)

dt <- ls_preprocessed$vsd_mat
rownames(dt) <- sapply(strsplit(rownames(dt), "\\."), "[[", 1)

genetopaths=genetopaths[which(genetopaths$V8=='Homo sapiens'),]
genetopaths$V3=sapply(genetopaths$V3,function(x){return(unlist(strsplit(x, split=' \\['))[1])})

allgenes=rownames(dt)
ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
allgenes.with.id=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id_version", "entrezgene_id"),values=allgenes, mart= ensembl)

##summarise expression per pathway
paths=unique(genetopaths$V6)

pathexpl=list()

for (i in 1:length(paths)){
  #print(i)
  ens=unique(genetopaths$V1[which(genetopaths$V6==paths[i])])
  genname=unique(allgenes.with.id[which(allgenes.with.id$ensembl_gene_id %in% ens),1])
  genname_com=intersect(genname, rownames(dt))
  if (length(genname_com)>1){
  pathexpl[[paths[i]]]=colMeans(dt[genname_com,])
  }
  if (length(genname_com)<1){
    pathexpl[[paths[i]]]=NA
  } 
  if (length(genname_com)==1){
    pathexpl[[paths[i]]]=dt[genname_com,]
  } 
  
}

pathexpdf=Reduce(rbind, pathexpl)
rownames(pathexpdf)=names(pathexpl)

pathexpdf=pathexpdf[complete.cases(pathexpdf),] 
colnames(pathexpdf) <- ls_preprocessed$p_all$pt_ID
path_nm <- rownames(pathexpdf)
pathexpdf <- t(pathexpdf)
colnames(pathexpdf) <- rep(paste0("reactome_", 1:ncol(pathexpdf)))

path_nm <- data.frame(path_ID = colnames(pathexpdf), path_names = path_nm)

write.csv(pathexpdf, file.path(path, 'data', 'rna_reactome.csv'))
write.csv(path_nm, file.path(path, 'data', 'rna_reactomeinfo.csv'))

## 3. TF activity

library(viper)
dorothea2viper_regulons <- function(df) {
  regulon_list <- split(df, df$tf)
  viper_regulons <- lapply(regulon_list, function(regulon) {
    tfmode <- stats::setNames(regulon$mor, regulon$target)
    list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
  })
  return(viper_regulons)
}


result <- WGCNA::collapseRows(ls_preprocessed$vsd_mat,
                          rowGroup=ls_preprocessed$rna_all$Feature_gene_name,
                          rowID=rownames(ls_preprocessed$vsd_mat),
                          method="MaxMean")

data <- data.frame(result$datETcollapsed)

regulons = dorothea::dorothea_hs %>%
  filter(confidence %in% c("A", "B"))

regu <- dorothea2viper_regulons(regulons)
vpres_25 <- t(viper(data, regu, verbose = FALSE, minsize = 25))
rownames(vpres_25) <- ls_preprocessed$p_all$pt_ID

vpres_4 <- t(viper(data, regu, verbose = FALSE, minsize = 4))
rownames(vpres_4) <- ls_preprocessed$p_all$pt_ID

write.csv(vpres_25, file.path(path, 'data', 'rna_viper25.csv'))
write.csv(vpres_4, file.path(path, 'data', 'rna_viper4.csv'))


###############################################################################
# WES 
###############################################################################
## 1. Binary matrix

library(maftools)
library(dplyr)
library(tidyr)
maf_dir <- "/Users/senosam/Documents/Massion_lab/WES_summary/TwistWES_Tumor_pipeline.freq0.01.filtered.tsv_020921.maf"
mut = read.maf(maf = read.delim(maf_dir))
x = getClinicalData(mut)
x$Tumor_Sample_Barcode = as.character(x$Tumor_Sample_Barcode)
pt_ID <- sapply(strsplit(x$Tumor_Sample_Barcode, "pt"), "[[", 2)
pt_ID <- sapply(strsplit(pt_ID, "_"), "[[", 1)
clin_dt <- read.csv('/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2020FEB25_SA.csv')
clin_dt <- clin_dt[match(pt_ID, clin_dt$pt_ID),]
clin_dt <- cbind('Tumor_Sample_Barcode'=x$Tumor_Sample_Barcode, clin_dt)
mut = read.maf(maf = read.delim(maf_dir), clinicalData = clin_dt)


mut_df <- mut@data

h_sym <- unique(mut_df$Hugo_Symbol)
pt_ID <- unique(mut_df$Tumor_Sample_Barcode)

mut_dt <- matrix(0, nrow=length(h_sym), ncol=length(pt_ID))
rownames(mut_dt) <- h_sym
colnames(mut_dt) <- pt_ID


for (i in 1:length(h_sym)){
    k <- which(mut_df$Hugo_Symbol == h_sym[i])
    pt <- mut_df$Tumor_Sample_Barcode[k]
    mut_dt[i,pt] <- 1
}

colnames(mut_dt) <- sapply(strsplit(as.character(pt_ID), "pt"), "[[", 2)
colnames(mut_dt) <- sapply(strsplit(colnames(mut_dt), "_"), "[[", 1)
mut_dt <- t(mut_dt)

# Remove genes that are mutated in less than X% of samples
prcnt <- 0.1 # cutoff: genes must be mutated in prcnt*100 % of the samples
mut_dt_sum <- data.frame(t(mut_dt)) %>% 
  mutate(sum = rowSums(., na.rm = TRUE), 
         genes = rownames(.)) %>%
  arrange(., desc(sum)) %>%
  dplyr::select(genes, sum) %>%
  filter(., sum > nrow(mut_dt)*prcnt)

mut_dt <- mut_dt[,mut_dt_sum$genes]

write.csv(mut_dt, file.path(path, 'data', 'wes_binary.csv'))


###############################################################################
# Radiomics
###############################################################################
## 1. Healthmyne
m_HM <- read.csv('/Users/senosam/Documents/Massion_lab/radiomics_summary/TMA36_HM.csv')
rownames(m_HM) <- m_HM[,1]
m_HM <- m_HM[-which(is.na(m_HM[,2])),]
m_HM <- m_HM[,15:ncol(m_HM)]
nona <- colnames(m_HM)[apply(m_HM, 2, anyNA)]
m_HM <- m_HM[,-which(colnames(m_HM) %in% nona)]

write.csv(m_HM, file.path(path, 'data', 'rad_healthmyne.csv'))


###############################################################################
# Clinical data
###############################################################################
## read all datasets and determine which patients have complete data
## Select patients that have complete data (n=57)
## split 75% training 25% test (43, 14) (by pt_id)

path='/Users/senosam/Documents/Massion_lab/data_integration/models'

cytof_freq <- read.csv(file.path(path, 'data', 'cytof_freq.csv'), row.names = 1)
cytof_medianprot <- read.csv(file.path(path, 'data', 'cytof_medianprot.csv'), row.names = 1)
rna_topgenes <- read.csv(file.path(path, 'data', 'rna_topgenes.csv'), row.names = 1)
rna_reactome <- read.csv(file.path(path, 'data', 'rna_reactome.csv'), row.names = 1)
rna_viper25 <- read.csv(file.path(path, 'data', 'rna_viper25.csv'), row.names = 1)
rna_viper4 <- read.csv(file.path(path, 'data', 'rna_viper4.csv'), row.names = 1)
wes_binary <- read.csv(file.path(path, 'data', 'wes_binary.csv'), row.names = 1)
rad_healthmyne <- read.csv(file.path(path, 'data', 'rad_healthmyne.csv'), row.names = 1)

# Selecting list of patients with complete data
k = which(rownames(cytof_freq) %in% rownames(rna_topgenes))
pt_ID <- rownames(cytof_freq)[k]
k = which(rownames(wes_binary) %in% pt_ID)
pt_ID <- rownames(wes_binary)[k]

CDE <- read.csv(file = '/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2021JAN13_SA_MF.csv')
CDE <- CDE[match(pt_ID, CDE$pt_ID),]
pts_y <- data.frame(pt_ID=CDE$pt_ID, SILA=CDE$SILA)

training_y <- sample_n(pts_y, size = round(nrow(pts_y)*0.75))
test_y <- pts_y[-which(pts_y$pt_ID %in% training_y$pt_ID),]

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


