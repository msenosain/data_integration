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


prcnt_by_pt

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

rownames(sbst_exp) <- sbst_exp$pt_ID
sbst_exp$pt_ID <- NULL
sbst_exp$CANARY <- NULL


save(prcnt_by_pt, sbst_exp,
    file='/Users/senosam/Documents/Massion_lab/data_integration/CyTOF_data.Rdata')

