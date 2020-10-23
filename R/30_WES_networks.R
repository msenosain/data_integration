library(maftools)
library(dplyr)
maf_dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Data/WES/TwistWES_Tumor_pipeline.filter.allSamples.maf"
mut = read.maf(maf = maf_dir)
x = getClinicalData(mut)
x$Tumor_Sample_Barcode = as.character(x$Tumor_Sample_Barcode)
pt_ID <- sapply(strsplit(x$Tumor_Sample_Barcode, "pt"), "[[", 2)
pt_ID <- sapply(strsplit(pt_ID, "_"), "[[", 1)
clin_dt <- read.csv('/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2020FEB25_SA.csv')
clin_dt <- clin_dt[match(pt_ID, clin_dt$pt_ID),]
clin_dt <- cbind('Tumor_Sample_Barcode'=x$Tumor_Sample_Barcode, clin_dt)
mut = read.maf(maf = maf_dir, clinicalData = clin_dt)


mut_dt <- mut@data

h_sym <- unique(mut_dt$Hugo_Symbol)
pt_ID <- unique(mut_dt$Tumor_Sample_Barcode)

mut_df <- matrix(0, nrow=length(h_sym), ncol=length(pt_ID))
rownames(mut_df) <- h_sym
colnames(mut_df) <- pt_ID


for (i in 1:length(h_sym)){
    k <- which(mut_dt$Hugo_Symbol == h_sym[i])
    pt <- mut_dt$Tumor_Sample_Barcode[k]
    mut_df[i,pt] <- 1
}

colnames(mut_df) <- sapply(strsplit(as.character(pt_ID), "pt"), "[[", 2)
colnames(mut_df) <- sapply(strsplit(colnames(mut_df), "_"), "[[", 1)