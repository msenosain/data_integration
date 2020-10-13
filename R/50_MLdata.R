############################################################################
# CyTOF
############################################################################

source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/20_ClustAnnot_functions.R")
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/30_DA_functions.R")
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/40_DE_functions.R")

# Load data
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/cellsubtypes.RData")

# Data pre-processing
## Changing un-identified T cells into "Other_immune"
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

## All cell/sub types B have Fib_Mesenchymal
ref <- ref[-42,] # sample 13376 has been collected twice and will be merged, only one ref is needed
grp <- 'B'


#---------------------------------------------------------------------------
# CELL TYPE PROPORTIONS
#---------------------------------------------------------------------------

## Major cell types (immune as one)
prcnt_by_pt <- ClassAbundanceByPt(data=annot_df, ptID_col = 'pt_ID', 
    class_col = paste0('cell_type_', grp))
imm_oa <- prcnt_by_pt$Immune

## All subtypes corrected by epithelial content
sbst <- subset(annot_df, cell_type_A != 'Epithelial')
sbst[,paste0('subtype_', grp)] <- factor(sbst[,paste0('subtype_', grp)])
prcnt_by_pt <- ClassAbundanceByPt(data=sbst, ptID_col = 'pt_ID', 
    class_col = paste0('subtype_', grp))
sbst <- subset(annot_df, cell_type_A == 'Epithelial')
sbst$subtype_A <- factor(sbst$subtype_A)
prcnt_by_ptE <- ClassAbundanceByPt(data=sbst, ptID_col = 'pt_ID', 
    class_col = 'subtype_A')
prcnt_by_pt <- cbind(prcnt_by_pt/2, prcnt_by_ptE/2, 'ALL_Immune' = imm_oa)


#---------------------------------------------------------------------------
# PROTEIN EXPRESSION PER CELL TYPE
#---------------------------------------------------------------------------

med_all <- median_by_pt(annot_df, ref, subset_celltype=F,
    ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)

match_all <- function(med_subtype, med_all){
    med_subtype <- med_subtype[match(med_all$pt_ID,med_subtype$pt_ID),]
    med_subtype$pt_ID <- med_all$pt_ID 
    med_subtype$CANARY <- med_all$CANARY

    med_subtype
}

med_endo <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B4',
    celltype_name='Endothelial', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_endo)[3:ncol(med_endo)] <- paste0('Endo_',colnames(med_endo)[3:ncol(med_endo)])
med_endo <- match_all(med_endo, med_all)


med_fibmes <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B4',
    celltype_name='Fib_Mesenchymal', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_fibmes)[3:ncol(med_fibmes)] <- paste0('FMes_',colnames(med_fibmes)[3:ncol(med_fibmes)])
med_fibmes <- match_all(med_fibmes, med_all)


med_CD8T <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B4',
    celltype_name='Tc_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_CD8T)[3:ncol(med_CD8T)] <- paste0('Tc_',colnames(med_CD8T)[3:ncol(med_CD8T)])
med_CD8T <- match_all(med_CD8T, med_all)


med_CD4T <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B4',
    celltype_name='Th_cells', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_CD4T)[3:ncol(med_CD4T)] <- paste0('Th_',colnames(med_CD4T)[3:ncol(med_CD4T)])
med_CD4T <- match_all(med_CD4T, med_all)


med_mye <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B4',
    celltype_name='Myeloid', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_mye)[3:ncol(med_mye)] <- paste0('Mye_',colnames(med_mye)[3:ncol(med_mye)])
med_mye <- match_all(med_mye, med_all)


# Epithelial as one
med_epi <- median_by_pt(annot_df, ref, subset_celltype=T, celltype_col='subtype_B4',
    celltype_name='Epithelial', ask_features=F, ft_idxs = c(15, 17:31, 33:35, 37:48, 50, 51), 
    ptID_col = 'pt_ID', compare_groups = F)
colnames(med_epi)[3:ncol(med_epi)] <- paste0('Epi_',colnames(med_epi)[3:ncol(med_epi)])
med_epi <- match_all(med_epi, med_all)


# Merging all
sbst_exp <- do.call('cbind', 
    list(med_epi, med_endo, med_fibmes, med_CD8T, med_CD4T, med_mye))

k <- c(grep('CANARY', colnames(sbst_exp)), grep('pt_ID', colnames(sbst_exp)))
sbst_exp <- sbst_exp[,-sort(k)[3:length(k)]]
sbst_exp[is.na(sbst_exp)] <- 0


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
# Clust transcriptional programs
#---------------------------------------------------------------------------
clust_eigen <- data.frame(readr::read_tsv('/Users/senosam/Documents/Massion_lab/RNASeq_summary/Results_23_Sep_20_2/Eigengenes.tsv'))
rownames(clust_eigen) <- clust_eigen$X1
clust_eigen$X1 <- NULL
clust_eigen <- data.frame(t(clust_eigen))
rownames(clust_eigen) <- p_all$pt_ID

#---------------------------------------------------------------------------
# Matrix with genes extracted by Clust
#---------------------------------------------------------------------------
clust_genes <- c(as.matrix(data.frame(read_tsv('/Users/senosam/Documents/Massion_lab/RNASeq_summary/Results_23_Sep_20_2/Clusters_Objects.tsv', skip = 1))))
clust_genes <- na.omit(clust_genes)
vsd_matTOP_clust <- vsd_matTOP[which(vsd_matTOP_ENSEMBL$gene %in% clust_genes),]
vsd_matTOP_clust <- data.frame(t(vsd_matTOP_clust))
rownames(vsd_matTOP_clust) <- p_all$pt_ID

#---------------------------------------------------------------------------
# Deconvolution results
#---------------------------------------------------------------------------
xcell_dcv <- data.frame(t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_only/XCELL/xCell_rnaseq_fpkm_xCell_1132060320.txt", row.names=1)))
rownames(xcell_dcv) <- p_all$pt_ID


############################################################################
# WES
############################################################################
# binary matrix for somatic mutations, remove genes that are mutated in less than 5% of samples
#...



############################################################################
# Building df for models
############################################################################
y <- read.csv('/Users/senosam/Documents/Massion_lab/radiomics_summary/TMA36_CANARY_khushbu.csv')

#---------------------------------------------------------------------------
# M1: CYTOF
#---------------------------------------------------------------------------
prcnt_by_pt <- na.omit(prcnt_by_pt[match(y$pt_ID, rownames(prcnt_by_pt)),])
sbst_exp <- na.omit(sbst_exp[match(y$pt_ID, sbst_exp$pt_ID),])
y_cytof <- y[which(y$pt_ID %in% sbst_exp$pt_ID),]
m_cytof <- cbind('SILA_S'=y_cytof$SILA_S, prcnt_by_pt, sbst_exp[,3:ncol(sbst_exp)])

write.csv(m_cytof, '/Users/senosam/Documents/Massion_lab/data_integration/ML/m_cytof.csv', row.names = T)

#---------------------------------------------------------------------------
# M2: RNA SEQ
#---------------------------------------------------------------------------
clust_eigen <- na.omit(clust_eigen[match(y$pt_ID, rownames(clust_eigen)),])
vsd_matTOP_clust <- na.omit(vsd_matTOP_clust[match(y$pt_ID, rownames(vsd_matTOP_clust)),])
xcell_dcv <- na.omit(xcell_dcv[match(y$pt_ID, rownames(xcell_dcv)),])
y_rna <- y[which(y$pt_ID %in% rownames(xcell_dcv)),]
m_rna <- cbind('SILA_S'=y_rna$SILA_S, clust_eigen, vsd_matTOP_clust, xcell_dcv)

write.csv(m_rna, '/Users/senosam/Documents/Massion_lab/data_integration/ML/m_rna.csv', row.names = T)

#---------------------------------------------------------------------------
# M3: MUTATION
#---------------------------------------------------------------------------
#...

#---------------------------------------------------------------------------
# M4: CYTOF & RNA
#---------------------------------------------------------------------------
k <- which(rownames(m_cytof) %in% rownames(m_rna))
m_cytof_rna <- cbind(m_cytof[which(rownames(m_cytof) %in% rownames(m_rna)),],
    m_rna[which(rownames(m_rna) %in% rownames(m_cytof)),2:ncol(m_rna)])

write.csv(m_cytof_rna, '/Users/senosam/Documents/Massion_lab/data_integration/ML/m_cytof_rna.csv', row.names = T)



#---------------------------------------------------------------------------
# M1: CYTOF & RNA & MUTATION
#---------------------------------------------------------------------------
#...
