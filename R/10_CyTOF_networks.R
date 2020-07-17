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

# Changing un-identified T cells and NK into "Other_immune"
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
    class_col = 'subtype_B4') # this is epi as one
prcnt_by_pt2 <- ClassAbundanceByPt(data=annot_df, ptID_col = 'pt_ID', 
    class_col = 'subtype_B')

# Changing un-identified T cells and NK into "Other_immune" in protein exp data
sbst_exp <- sbst_exp[-grep('NK|Tcells', colnames(sbst_exp))]
sbst_exp_2 <- sbst_exp_2[-grep('NK|Tcells', colnames(sbst_exp_2))]

##########################-- BUILD NODE ATTRIBUTES --##########################
## Clinical data
pdata <- data.frame(pdata, row.names = pdata$pt_ID, stringsAsFactors = FALSE)

## Cell type percentages
prcnt <- data.frame(prcnt_by_pt1, row.names = row.names(prcnt_by_pt1), 
    stringsAsFactors = FALSE)
prcnt2 <- data.frame(prcnt_by_pt2, row.names = row.names(prcnt_by_pt2), 
    stringsAsFactors = FALSE)

## Median protein expression per cell type per patient
sbst <- data.frame(sbst_exp, row.names = sbst_exp$pt_ID, 
    stringsAsFactors = FALSE)
sbst2 <- data.frame(sbst_exp_2, row.names = sbst_exp_2$pt_ID, 
    stringsAsFactors = FALSE)


###############################################################################
# NETWORK 1:  BASED ON CELL TYPE PERCENTAGES (EPI AS ONE)
###############################################################################

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
createNetworkFromIgraph(network,"cytof_ntw1")


###############-- APPEND NODE ATTRIBUTES TO NTW1 IN CYTOSCAPE --###############
## Clinical data
loadTableData(pdata)

## Cell type percentages
loadTableData(prcnt)
loadTableData(prcnt2)

## Median protein expression per cell type per patient
loadTableData(sbst)
loadTableData(sbst2)


###############################################################################
# NETWORK 1.2:  BASED ON CELL TYPE PERCENTAGES (EPI BY SUBSETS) use this always
###############################################################################

############################-- BUILD NETWORK --############################
# Compute correlation matrix
corr_pt <- Hmisc::rcorr(t(as.matrix(prcnt_by_pt2)), type = 'spearman')
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
#createNetworkFromIgraph(network,"cytof_ntw1.2")
createNetworkFromIgraph(network,"cytof_ntw")


###############-- APPEND NODE ATTRIBUTES TO NTW1 IN CYTOSCAPE --###############
## Clinical data
loadTableData(pdata)

## Cell type percentages
loadTableData(prcnt)
loadTableData(prcnt2)

## Median protein expression per cell type per patient
loadTableData(sbst)
loadTableData(sbst2)


###############################################################################
# NETWORK 2: BASED ON PROTEIN EXPRESSION BY CELL TYPE
###############################################################################

############################-- BUILD NETWORK --############################
sbst_A <- data.frame(sbst_exp[,3:ncol(sbst_exp)], row.names=sbst_exp$pt_ID)
# Compute correlation matrix
corr_pt <- Hmisc::rcorr(t(as.matrix(sbst_A)), type = 'spearman')
corr_mat <- corr_pt$r

# Select only significant correlations
corr_pt$P[which(is.na(corr_pt$P)==TRUE)] <- 0
k <- which(corr_pt$P<0.05)
corr_mat[-k] <- 0
 
# Make an Igraph object from this matrix:
network <- graph_from_adjacency_matrix(corr_mat, weighted=T, mode="undirected", 
    diag=F)
E(network)$weight

# Open Cytoscape and confirm connexion
cytoscapePing()
createNetworkFromIgraph(network,"cytof_ntw2")


###############-- APPEND NODE ATTRIBUTES TO NTW1 IN CYTOSCAPE --###############
## Clinical data
loadTableData(pdata)

## Cell type percentages
loadTableData(prcnt)
loadTableData(prcnt2)

## Median protein expression per cell type per patient
loadTableData(sbst)
loadTableData(sbst2)



###############################################################################
# NETWORK 2.2: BASED ON PROTEIN EXPRESSION BY CELL TYPE
###############################################################################

############################-- BUILD NETWORK --############################
sbst_A <- data.frame(sbst_exp_2[,3:ncol(sbst_exp_2)], row.names=sbst_exp_2$pt_ID)
# Compute correlation matrix
corr_pt <- Hmisc::rcorr(t(as.matrix(sbst_A)), type = 'spearman')
corr_mat <- corr_pt$r

# Select only significant correlations
corr_pt$P[which(is.na(corr_pt$P)==TRUE)] <- 0
k <- which(corr_pt$P<0.05)
corr_mat[-k] <- 0
 
# Make an Igraph object from this matrix:
network <- graph_from_adjacency_matrix(corr_mat, weighted=T, mode="undirected", 
    diag=F)
E(network)$weight

# Open Cytoscape and confirm connexion
cytoscapePing()
createNetworkFromIgraph(network,"cytof_ntw2.2")


###############-- APPEND NODE ATTRIBUTES TO NTW1 IN CYTOSCAPE --###############
## Clinical data
loadTableData(pdata)

## Cell type percentages
loadTableData(prcnt)
loadTableData(prcnt2)

## Median protein expression per cell type per patient
loadTableData(sbst)
loadTableData(sbst2)








