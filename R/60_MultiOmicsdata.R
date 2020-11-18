############################################################################
# CyTOF
############################################################################
load('/Users/senosam/Documents/Massion_lab/data_integration/CyTOF_data.Rdata')
e <- c(15,22,25,33,36,38,41,52,56,61,64,68,80,83,94,103,106,111,124,127,131)
sbst_exp <- sbst_exp[-e,]

############################################################################
# RNA
############################################################################
load('/Users/senosam/Documents/Massion_lab/data_integration/RNA_data.Rdata')
rownames(vsd_matTOP_clust_E) <- sapply(strsplit(rownames(vsd_matTOP_clust_E), "\\."), "[[", 1)
rownames(vsd_matTOP_ENSEMBL) <- sapply(strsplit(rownames(vsd_matTOP_ENSEMBL), "\\."), "[[", 1)

############################################################################
# WES
############################################################################
load('/Users/senosam/Documents/Massion_lab/data_integration/WES_data.Rdata')

############################################################################
# CDE
############################################################################
CDE_TMA36 <- read.csv(file = '/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2020FEB25_SA.csv')

############################################################################
# Building df for multiplex
############################################################################
pt_ID <- as.character(sort(as.numeric(unique(c(rownames(sbst_exp), rownames(vsd_matTOP), rownames(mut_dt))))))

format_mo <- function(df, pt_ID){
    df <- df[pt_ID, ]
    rownames(df) <- pt_ID
    df <- t(df)
    df
}


CyTOF_prcnt <- format_mo(prcnt_by_pt, pt_ID)
CyTOF_exp <- format_mo(sbst_exp, pt_ID)
RNA_top12K <- format_mo(vsd_matTOP, pt_ID)
RNA_top12K_E <- format_mo(vsd_matTOP_ENSEMBL, pt_ID)
RNA_topclust <- format_mo(vsd_matTOP_clust, pt_ID)
RNA_topclust_E <- format_mo(vsd_matTOP_clust_E, pt_ID)
RNA_clusteigen <- format_mo(clust_eigen, pt_ID)
RNA_xcell <- format_mo(xcell_dcv, pt_ID)
mut_dt <- format_mo(mut_dt, pt_ID)
CDE_TMA36 <- format_mo(CDE_TMA36, pt_ID) #only including those pts with cytof, rna or wes data


save(CyTOF_prcnt, CyTOF_exp, RNA_top12K, RNA_top12K_E, RNA_topclust, 
    RNA_topclust_E, RNA_clusteigen, RNA_xcell, mut_dt, CDE_TMA36,
    file='/Users/senosam/Documents/Massion_lab/data_integration/MO_data.Rdata')


