###############################################################################
# CLINICAL DATA
###############################################################################
CDE <- read.csv(file = '/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2021JAN13_SA_MF.csv')

# Remove PHI info
k <- grep("Date|Age", colnames(CDE))
CDE <- CDE[-k,]

write.csv(CDE, "/Users/senosam/Documents/Massion_lab/data_integration/Multiplex/CDE_TMA36_DI.csv", row.names = FALSE)

###############################################################################
# CyTOF DATA
###############################################################################
load("/Users/senosam/Documents/Massion_lab/CyTOF_summary/both/percent_pt.RData")

corr_f <- function(data, rcorr_type = 'spearman', p.adjust_method = 'BH'){

    res <- Hmisc::rcorr(as.matrix(data), type = rcorr_type) #for corr plot
    # corrplot
    corrected_pvals <- p.adjust(res$P, method = p.adjust_method)
    corrected_pvals <- matrix(corrected_pvals, nrow = ncol(res$P), 
        ncol = ncol(res$P))
    colnames(corrected_pvals)<- colnames(res$P)
    rownames(corrected_pvals)<- rownames(res$P)

    res$P <- corrected_pvals

    res    
}


prcnt_subtypes_corr <- corr_f(t(prcnt_subtypes))
prcnt_clusters_corr <- corr_f(t(prcnt_clusters))

dir <- '/Users/senosam/Documents/Massion_lab/data_integration/Multiplex/CyTOF'
save(prcnt_subtypes, prcnt_subtypes_corr, prcnt_clusters, prcnt_clusters_corr,
    file = file.path(dir,'percent_pt_corr.RData'))