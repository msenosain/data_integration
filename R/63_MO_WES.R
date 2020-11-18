############################################################################
# WES
############################################################################

#---------------------------------------------------------------------------
# Generate binary matrix for somatic mutations
#---------------------------------------------------------------------------

library(maftools)
library(dplyr)
maf_dir <- "~/Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Data/WES/TwistWES_Tumor_pipeline.filter.allSamples.maf"
mut = read.maf(maf = read.delim(maf_dir))
x = getClinicalData(mut)
x$Tumor_Sample_Barcode = as.character(x$Tumor_Sample_Barcode)
pt_ID <- sapply(strsplit(x$Tumor_Sample_Barcode, "pt"), "[[", 2)
pt_ID <- sapply(strsplit(pt_ID, "_"), "[[", 1)
clin_dt <- read.csv('/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2020FEB25_SA.csv')
clin_dt <- clin_dt[match(pt_ID, clin_dt$pt_ID),]
clin_dt <- cbind('Tumor_Sample_Barcode'=x$Tumor_Sample_Barcode, clin_dt)
mut = read.maf(maf = read.delim(maf_dir), clinicalData = clin_dt)


mut_df <- mut@data

h_sym <- unique(mut_df$Hugo_Symbol)
pt_ID <- unique(mut_df$Tumor_Sample_Barcode)

mut_dt <- matrix(0, nrow=length(h_sym), ncol=length(pt_ID))
rownames(mut_dt) <- h_sym
colnames(mut_dt) <- pt_ID


for (i in 1:length(h_sym)){
    k <- which(mut_df$Hugo_Symbol == h_sym[i])
    pt <- mut_df$Tumor_Sample_Barcode[k]
    mut_dt[i,pt] <- 1
}

colnames(mut_dt) <- sapply(strsplit(as.character(pt_ID), "pt"), "[[", 2)
colnames(mut_dt) <- sapply(strsplit(colnames(mut_dt), "_"), "[[", 1)

mut_dt <- t(mut_dt)

#write.csv(mut_dt, file='/Users/senosam/Documents/Massion_lab/WES_summary/summary/binaryWES.csv')
#mut_dt <- read.csv(file='/Users/senosam/Documents/Massion_lab/WES_summary/summary/binaryWES.csv', row.names = 1)

#---------------------------------------------------------------------------
# Remove genes that are mutated in less than X% of samples
#---------------------------------------------------------------------------
library(dplyr)
library(tidyr)
prcnt <- 0.1 # cutoff: genes must be mutated in prcnt*100 % of the samples
mut_dt_sum <- as_tibble(mut_dt) %>% 
  mutate(sum = rowSums(., na.rm = TRUE), 
         genes = rownames(mut_dt)) %>%
  arrange(., desc(sum)) %>%
  select(genes, sum) %>%
  filter(., sum > nrow(mut_dt)*prcnt)

mut_dt <- mut_dt[mut_dt_sum$genes,]

save(mut_dt,
    file='/Users/senosam/Documents/Massion_lab/data_integration/WES_data.Rdata')