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

###############################################################################
# RNA Seq DATA
###############################################################################

#-------------------------------- TOP GENES --------------------------------#
source('/Users/senosam/Documents/Repositories/Research/data_analysis_rnaseq/R/30_DEGanalysis.R')
source("/Users/senosam/Documents/Repositories/Research/data_analysis_cytof/R/31_supervised_analysis_viz.R")
environment_set()
load("/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData")
ls_preprocessed <- preprocess_rna(path_rnaseq = '/Users/senosam/Documents/Massion_lab/RNASeq_summary/rnaseq.RData', correct_batch = T, correct_gender = T)
k <- which(p_all$Vantage_ID %in% colnames(ls_preprocessed$vsd_mat))
p_all <- p_all[k,]
pData_rnaseq <- pData_rnaseq[k,]
CDE <- read.csv(file = '/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2021JAN13_SA_MF.csv')
pData_rnaseq <- CDE[match(pData_rnaseq$pt_ID, CDE$pt_ID),]

# TOP VARIANT GENES
vsd_mat <- ls_preprocessed$vsd_mat
colnames(vsd_mat) <- p_all$pt_ID
variances <- apply(vsd_mat, 1, var)
n_genes <- length(which(variances > 0.5))
top_genes <- data.frame(vsd_mat) %>%
   mutate(gene=rownames(.),
          symbol=ls_preprocessed$rna_all$Feature_gene_name,
          variances = variances) %>%
   arrange(desc(variances)) %>%
   dplyr::select(gene, symbol) %>%
   head(n_genes)
vsd_matTOP<- vsd_mat[top_genes$gene,]
rownames(vsd_matTOP) <- top_genes$symbol

vsd_matTOP_corr <- corr_f(vsd_matTOP)

dir <- '/Users/senosam/Documents/Massion_lab/data_integration/Multiplex/RNA_Seq'
save(vsd_matTOP, vsd_matTOP_corr, file = file.path(dir,'topgenes_pt_corr.RData'))


#-------------------------------- REACTOME --------------------------------#
library(biomaRt)
# Reactome pathways

reapaths <- read.table('/Users/senosam/Documents/Repositories/Research/DH_project01/data/Reactome2020/NCBI2Reactome_All_Levels.txt')
genetopaths <- read.delim('/Users/senosam/Documents/Repositories/Research/DH_project01/data/Reactome2020/Ensembl2Reactome_PE_Pathway.txt', header=F, stringsAsFactors = F)

dt <- ls_preprocessed$vsd_mat
rownames(dt) <- sapply(strsplit(rownames(dt), "\\."), "[[", 1)

genetopaths=genetopaths[which(genetopaths$V8=='Homo sapiens'),]
genetopaths$V3=sapply(genetopaths$V3,function(x){return(unlist(strsplit(x, split=' \\['))[1])})

allgenes=rownames(dt)
ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
allgenes.with.id=getBM(attributes=c("ensembl_gene_id","hgnc_symbol","ensembl_gene_id_version", "entrezgene_id"),values=allgenes, mart= ensembl)

##summarise expression per pathway
paths=unique(genetopaths$V6)

pathexpl=list()

for (i in 1:length(paths)){
  #print(i)
  ens=unique(genetopaths$V1[which(genetopaths$V6==paths[i])])
  genname=unique(allgenes.with.id[which(allgenes.with.id$ensembl_gene_id %in% ens),1])
  genname_com=intersect(genname, rownames(dt))
  if (length(genname_com)>1){
  pathexpl[[paths[i]]]=colMeans(dt[genname_com,])
  }
  if (length(genname_com)<1){
    pathexpl[[paths[i]]]=NA
  } 
  if (length(genname_com)==1){
    pathexpl[[paths[i]]]=dt[genname_com,]
  } 
  
}

pathexpdf=Reduce(rbind, pathexpl)
rownames(pathexpdf)=names(pathexpl)

pathexpdf=pathexpdf[complete.cases(pathexpdf),]

#ls_preprocessed$p_all$Vantage_ID == colnames(pathexpdf)
colnames(pathexpdf) <- ls_preprocessed$p_all$pt_ID

pathexpdf_corr <- corr_f(pathexpdf)

dir <- '/Users/senosam/Documents/Massion_lab/data_integration/Multiplex/RNA_Seq'
save(pathexpdf, pathexpdf_corr, file = file.path(dir,'reactome_pt_corr.RData'))


#------------------------------ DECONVOLUTION -------------------------------#

mcp_dcv <- t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_only/MCP/mcp_fpkm_dcv.txt", row.names=1))
qts_dcv <- t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_only/QTS/qts_fpkm_dcv.txt", row.names=1))
cbs_dcv <- as.matrix(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_only/CBS/CIBERSORT.rna_only_fpkm.txt", row.names=1)[,1:22])
xcell_dcv <- t(read.delim("~/Documents/Massion_lab/RNASeq_summary/deconvolution/output/rna_only/XCELL/xCell_rnaseq_fpkm_xCell_1132060320.txt", row.names=1))

rownames(mcp_dcv) <- ls_preprocessed$p_all$pt_ID
rownames(qts_dcv) <- ls_preprocessed$p_all$pt_ID
rownames(cbs_dcv) <- ls_preprocessed$p_all$pt_ID
rownames(xcell_dcv) <- ls_preprocessed$p_all$pt_ID


mcp_dcv_corr <- corr_f(t(mcp_dcv))
qts_dcv_corr <- corr_f(t(qts_dcv))
cbs_dcv_corr <- corr_f(t(cbs_dcv))
xcell_dcv_corr <- corr_f(t(xcell_dcv))


dir <- '/Users/senosam/Documents/Massion_lab/data_integration/Multiplex/RNA_Seq'
save(mcp_dcv, mcp_dcv_corr, 
    qts_dcv, qts_dcv_corr,
    cbs_dcv, cbs_dcv_corr,
    xcell_dcv, xcell_dcv_corr, file = file.path(dir,'dcv_pt_corr.RData'))



############################################################################
# RADIOMICS - HEALTHMYNE
############################################################################
m_HM <- read.csv('/Users/senosam/Documents/Massion_lab/radiomics_summary/TMA36_HM.csv')
rownames(m_HM) <- m_HM[,1]
m_HM <- m_HM[-which(is.na(m_HM[,2])),]
m_HM <- m_HM[,15:ncol(m_HM)]
nona <- colnames(m_HM)[apply(m_HM, 2, anyNA)]
m_HM <- m_HM[,-which(colnames(m_HM) %in% nona)]

m_HM_corr <- corr_f(t(m_HM))

dir <- '/Users/senosam/Documents/Massion_lab/data_integration/Multiplex/Radiomics'
save(m_HM, m_HM_corr, file = file.path(dir,'radiomics_pt_corr.RData'))

