############################################################################
# RNA
############################################################################
source('/Users/senosam/Documents/Repositories/Research/data_analysis_rnaseq/R/30_DEGanalysis.R')
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData")
environment_set()

# Data preprocessing
ls_preprocessed <- preprocess_rna(path_rnaseq = '/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData', correct_batch = T, correct_gender = T)
k <- which(p_all$Vantage_ID %in% colnames(ls_preprocessed$vsd_mat))
p_all <- p_all[k,]
pData_rnaseq <- pData_rnaseq[k,]

# Top varaint genes
vsd_mat <- ls_preprocessed$vsd_mat
variances <- apply(vsd_mat, 1, var)
n_genes <- length(which(variances > 0.5))
#n_genes <- 200
top_genes <- data.frame(vsd_mat) %>%
   mutate(gene=rownames(.),
          symbol=ls_preprocessed$rna_all$Feature_gene_name,
          variances = variances) %>%
   arrange(desc(variances)) %>%
   dplyr::select(gene, symbol) %>%
   head(n_genes)
vsd_matTOP <- vsd_mat[top_genes$gene,]
vsd_matTOP_ENSEMBL <- vsd_matTOP
rownames(vsd_matTOP) <- top_genes$symbol

vsd_matTOP_ENSEMBL <- cbind(gene = rownames(vsd_matTOP_ENSEMBL), data.frame(vsd_matTOP_ENSEMBL))
rownames(vsd_matTOP_ENSEMBL) <- c(1:nrow(vsd_matTOP_ENSEMBL))


#---------------------------------------------------------------------------
# Matrix with genes extracted by Clust
#---------------------------------------------------------------------------
clust_genes <- c(as.matrix(data.frame(read_tsv('/Users/senosam/Documents/Massion_lab/RNASeq_summary/Results_23_Sep_20_2/Clusters_Objects.tsv', skip = 1))))
clust_genes <- na.omit(clust_genes)
vsd_matTOP_clust <- vsd_matTOP[which(vsd_matTOP_ENSEMBL$gene %in% clust_genes),]
vsd_matTOP_clust <- data.frame(t(vsd_matTOP_clust))
rownames(vsd_matTOP_clust) <- p_all$pt_ID

vsd_matTOP_clust_E <- vsd_matTOP_ENSEMBL[which(vsd_matTOP_ENSEMBL$gene %in% clust_genes),]
rownames(vsd_matTOP_clust_E) <- vsd_matTOP_clust_E$gene
vsd_matTOP_clust_E$gene <- NULL
vsd_matTOP_clust_E <- data.frame(t(vsd_matTOP_clust_E))
rownames(vsd_matTOP_clust_E) <- p_all$pt_ID

#---------------------------------------------------------------------------
# Matrix with top variant genes
#---------------------------------------------------------------------------
vsd_matTOP <- data.frame(t(vsd_matTOP))
rownames(vsd_matTOP) <- p_all$pt_ID

rownames(vsd_matTOP_ENSEMBL) <- vsd_matTOP_ENSEMBL$gene
vsd_matTOP_ENSEMBL$gene <- NULL
vsd_matTOP_ENSEMBL <- data.frame(t(vsd_matTOP_ENSEMBL))
rownames(vsd_matTOP_ENSEMBL) <- p_all$pt_ID

#---------------------------------------------------------------------------
# Clust transcriptional programs
#---------------------------------------------------------------------------
clust_eigen <- data.frame(readr::read_tsv('/Users/senosam/Documents/Massion_lab/RNASeq_summary/Results_23_Sep_20_2/Eigengenes.tsv'))
rownames(clust_eigen) <- clust_eigen$X1
clust_eigen$X1 <- NULL
clust_eigen <- data.frame(t(clust_eigen))
rownames(clust_eigen) <- p_all$pt_ID

#---------------------------------------------------------------------------
# Deconvolution results
#---------------------------------------------------------------------------
xcell_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_only/XCELL/xCell_rnaseq_fpkm_xCell_1132060320.txt", row.names=1)))
rownames(xcell_dcv) <- p_all$pt_ID


save(vsd_matTOP, vsd_matTOP_ENSEMBL, vsd_matTOP_clust, vsd_matTOP_clust_E, clust_eigen, xcell_dcv,
    file='/Users/senosam/Documents/Massion_lab/data_integration/RNA_data.Rdata')

