####################-- LOAD LIBRARIES AND FUNCITONS--####################
source('/Users/senosam/Documents/Repositories/Research/data_analysis_rnaseq/R/30_DEGanalysis.R')
environment_set()
library(genefilter)

############################-- LOAD DATA--############################
ls_preprocessed <- preprocess_rna(path_rnaseq = '/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData', 
    correct_batch = T, correct_gender = T)


############################-- ORGANIZING DATA--############################
vsd_mat <- ls_preprocessed$vsd_mat
colnames(vsd_mat) <- ref$pt_ID
ref <- ls_preprocessed$p_all
pData <- ls_preprocessed$pData_rnaseq
rna_raw <- ls_preprocessed$rna_all

# TOP 50%
vsd_mat_50 <- varFilter(vsd_mat, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
#rownames(vsd_mat_50) <- make.names(rna_raw$Feature_gene_name[which(rna_raw$Feature %in% rownames(vsd_mat_50))], unique = TRUE)

# TOP 25%
vsd_mat_25 <- varFilter(vsd_mat, var.func=IQR, var.cutoff=0.75, filterByQuantile=TRUE)
#rownames(vsd_mat_25) <- make.names(rna_raw$Feature_gene_name[which(rna_raw$Feature %in% rownames(vsd_mat_25))], unique = TRUE)

# TOP 10%
    
#rownames(vsd_mat_25) <- make.names(rna_raw$Feature_gene_name[which(rna_raw$Feature %in% rownames(vsd_mat_25))], unique = TRUE)

# TOP 0.5%
vsd_mat_2 <- varFilter(vsd_mat, var.func=IQR, var.cutoff=0.995, filterByQuantile=TRUE)
#rownames(vsd_mat_25) <- make.names(rna_raw$Feature_gene_name[which(rna_raw$Feature %in% rownames(vsd_mat_25))], unique = TRUE)

##########################-- BUILD NODE ATTRIBUTES --##########################

## Clinical data
pData <- data.frame(pData, row.names = pData$pt_ID, stringsAsFactors = FALSE)

## Cell type percentages (deconvolution) XCELL
xcell_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_only/XCELL/xCell_rnaseq_tpm_xCell_1212060320.txt", row.names=1)),
    row.names = ref$pt_ID, stringsAsFactors = FALSE)
mcp_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_only/MCP/mcp_dcv.txt", row.names=1)),
    row.names = ref$pt_ID, stringsAsFactors = FALSE)
qts_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_only/QTS/qts_dcv.txt", row.names=1)),
    row.names = ref$pt_ID, stringsAsFactors = FALSE)
cbs_dcv <- data.frame(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_only/CBS/CIBERSORT.rna_only_tpm.txt", row.names=1), 
    row.names = ref$pt_ID, stringsAsFactors = FALSE)
#table(rownames(cbs_dcv)==ref$Vantage_ID)

## Gene expression per patient, top x%
rownames(vsd_mat_2) <- make.names(rna_raw$Feature_gene_name[which(rna_raw$Feature %in% rownames(vsd_mat_2))], unique = TRUE)
vsd_mat_t <- t(vsd_mat_2)
vsd_mat_t <- data.frame(vsd_mat_t, row.names = colnames(vsd_mat), stringsAsFactors = FALSE)


###############################################################################
# NETWORK 1: USING ALL GENES
###############################################################################

############################-- BUILD NETWORK --############################
# Compute correlation matrix
corr_pt <- Hmisc::rcorr(vsd_mat, type = 'spearman')
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
createNetworkFromIgraph(network,"rna_ntw1")

###############-- APPEND NODE ATTRIBUTES TO NTW1 IN CYTOSCAPE --###############
## Clinical data
loadTableData(pData)

## Cell type percentages
loadTableData(xcell_dcv)
loadTableData(mcp_dcv)
loadTableData(qts_dcv)
loadTableData(cbs_dcv)

## Median protein expression per cell type per patient
loadTableData(vsd_mat_t)


###############################################################################
# NETWORK 2: TOP 50% VARIANT GENES
###############################################################################

############################-- BUILD NETWORK --############################
# Compute correlation matrix
corr_pt <- Hmisc::rcorr(vsd_mat_50, type = 'spearman')
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
createNetworkFromIgraph(network,"rna_ntw50")

###############-- APPEND NODE ATTRIBUTES TO NTW1 IN CYTOSCAPE --###############
## Clinical data
loadTableData(pData)

## Cell type percentages
loadTableData(xcell_dcv)
loadTableData(mcp_dcv)
loadTableData(qts_dcv)
loadTableData(cbs_dcv)

## Median protein expression per cell type per patient
loadTableData(vsd_mat_t)


###############################################################################
# NETWORK 3: TOP 25% VARIANT GENES
###############################################################################

############################-- BUILD NETWORK --############################
# Compute correlation matrix
corr_pt <- Hmisc::rcorr(vsd_mat_25, type = 'spearman')
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
createNetworkFromIgraph(network,"rna_ntw25")

###############-- APPEND NODE ATTRIBUTES TO NTW1 IN CYTOSCAPE --###############
## Clinical data
loadTableData(pData)

## Cell type percentages
loadTableData(xcell_dcv)
loadTableData(mcp_dcv)
loadTableData(qts_dcv)
loadTableData(cbs_dcv)

## Median protein expression per cell type per patient
loadTableData(vsd_mat_t)


###############################################################################
# NETWORK 3: TOP 10% VARIANT GENES
###############################################################################

############################-- BUILD NETWORK --############################
# Compute correlation matrix
corr_pt <- Hmisc::rcorr(vsd_mat_10, type = 'spearman')
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
createNetworkFromIgraph(network,"rna_ntw10")

###############-- APPEND NODE ATTRIBUTES TO NTW1 IN CYTOSCAPE --###############
## Clinical data
loadTableData(pData)

## Cell type percentages
loadTableData(xcell_dcv)
loadTableData(mcp_dcv)
loadTableData(qts_dcv)
loadTableData(cbs_dcv)

## Median protein expression per cell type per patient
loadTableData(vsd_mat_t)
