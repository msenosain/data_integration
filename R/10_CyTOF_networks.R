####################-- LOAD LIBRARIES AND FUNCITONS--####################
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/20_ClustAnnot_functions.R")
library(RCy3)
library(igraph)


############################-- LOAD DATA--############################
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/cellsubtypes.RData")
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/analysis/subsets/subset_exprmat.RData")
pdata <- read.csv("/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2020FEB25_SA.csv")


############################-- ORGANIZING DATA--############################
# All cell/sub types B have Fib_Mesenchymal
ref <- ref[-42,] # sample 13376 has been collected twice and will be merged, only one ref is needed

# Changing un-identified T cells into "Other_immune"
k <- which(annot_df$subtype_B =='T_cells')
k2 <- which(annot_df$subtype_B =='NK_cells')

annot_df['subtype_B'] <- as.character(annot_df$subtype_B)
annot_df[c(k,k2), 'subtype_B'] <- 'Other_immune'
annot_df$subtype_B <- as.factor(annot_df$subtype_B)

annot_df['subtype_B4'] <- as.character(annot_df$subtype_B)
annot_df[which(annot_df$cell_type_B=='Epithelial'), 'subtype_B4'] <- 'Epithelial'
annot_df$subtype_B4 <- as.factor(annot_df$subtype_B4)
annot_df$subtype_B4 <- factor(annot_df$subtype_B4, levels = c("Epithelial", 
    "Endothelial", "Fib_Mesenchymal", "Other_immune", "Myeloid", "Tc_cells", 
    "Th_cells"))

prcnt_by_pt1 <- ClassAbundanceByPt(data=annot_df, ptID_col = 'pt_ID', 
    class_col = 'subtype_B')
prcnt_by_pt2 <- ClassAbundanceByPt(data=annot_df, ptID_col = 'pt_ID', 
    class_col = 'subtype_B4') # this is epi as one


############################-- BUILD NETWORK --############################
# Compute correlation matrix
corr_pt <- Hmisc::rcorr(t(as.matrix(prcnt_by_pt1)), type = 'spearman')
corr_mat <- corr_pt$r

# Select only significant correlations
corr_pt$P[which(is.na(corr_pt$P)==TRUE)] <- 0
k <- which(corr_pt$P<0.05)
corr_mat[-k] <- 0
#corr_mat[corr_mat<0.5] <- 0 # Keep only high correlations
 
# Make an Igraph object from this matrix:
network <- graph_from_adjacency_matrix(corr_mat, weighted=T, mode="undirected", 
    diag=F)
E(network)$weight

# Open Cytoscape and confirm connexion
cytoscapePing()
createNetworkFromIgraph(network,"cytof_ntw")


###############-- APPEND NODE ATTRIBUTES TO NTW IN CYTOSCAPE --###############
## Clincal data
pdata <- data.frame(pdata, row.names = pdata$pt_ID, stringsAsFactors = FALSE)
loadTableData(pdata)

## Cell type percentages
prcnt_by_pt1 <- data.frame(prcnt_by_pt1, row.names = row.names(prcnt_by_pt1), 
    stringsAsFactors = FALSE)
loadTableData(prcnt_by_pt1)

## Median protein expression per cell type per patient
sbst_exp <- data.frame(sbst_exp, row.names = sbst_exp$pt_ID, 
    stringsAsFactors = FALSE)
loadTableData(sbst_exp)
#sbst_exp_2 # epithelial by subset








