library(maftools)
library(dplyr)
maf_dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Data/WES/TwistWES_Tumor_pipeline.filter.allSamples.maf"
mut = read.maf(maf = maf_dir)
x = getClinicalData(laml)
x$Tumor_Sample_Barcode = as.character(x$Tumor_Sample_Barcode)
pt_ID <- sapply(strsplit(x$Tumor_Sample_Barcode, "pt"), "[[", 2)
pt_ID <- sapply(strsplit(pt_ID, "_"), "[[", 1)
clin_dt <- read.csv('/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2020FEB25_SA.csv')
clin_dt <- clin_dt[match(pt_ID, clin_dt$pt_ID),]
clin_dt <- cbind('Tumor_Sample_Barcode'=x$Tumor_Sample_Barcode, clin_dt)
mut = read.maf(maf = maf_dir, clinicalData = clin_dt)