# Loading libraries
library(dplyr)
library(tidyverse)
library(tidyr)
library(MOFA)
library(MOFAdata)
library(reticulate)
library(MultiAssayExperiment)
use_python("/Users/senosam/opt/anaconda3/bin/python", required=TRUE) 


MOFA_model <- function(mofa_data, mae = F, CDE){
    set.seed(1234)

    if (mae){
        mae_x <- MultiAssayExperiment(
            experiments = mofa_data, 
            colData = CDE
            )
        MOFAobject <- createMOFAobject(mae_x)
    } else {
        MOFAobject <- createMOFAobject(mofa_data)
    }
    

    # Setting MOFA options
    DataOptions <- getDefaultDataOptions()
    ModelOptions <- getDefaultModelOptions(MOFAobject)
    TrainOptions <- getDefaultTrainOptions()
    TrainOptions$DropFactorThreshold <- 0.02
    TrainOptions$seed <- 2017

    # Prepare MOFA object
    MOFAobject <- prepareMOFA(
      MOFAobject, 
      DataOptions = DataOptions,
      ModelOptions = ModelOptions,
      TrainOptions = TrainOptions
    )

    # Run MOFA
    MOFAobject <- runMOFA(MOFAobject)

    return(MOFAobject)
}



# Loading data
load('/Users/senosam/Documents/Massion_lab/data_integration/MO_data.Rdata')
CDE_TMA36 <- read.csv(file = '/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2020FEB25_SA.csv')
rownames(CDE_TMA36) <- CDE_TMA36$pt_ID
CDE_TMA36 <- CDE_TMA36[colnames(CyTOF_exp),]


# Preparing MOFA data
mofa_data1 <- list('CyTOF_exp'   = CyTOF_exp, 
                  'CyTOF_prcnt' = CyTOF_prcnt, 
                  'RNA_topgenes'= RNA_top12K_E, 
                  'RNA_xcell'   = RNA_xcell,
                  'RNA_viper25' = RNA_viper25,
                  'mutation'    = mut_dt)

mofa_data2 <- list('CyTOF_exp'   = CyTOF_exp, 
                  'CyTOF_prcnt' = CyTOF_prcnt, 
                  'RNA_topgenes'= RNA_top12K_E, 
                  'RNA_epidish' = RNA_epidish,
                  'RNA_viper25' = RNA_viper25,
                  'mutation'    = mut_dt)

mofa_data3 <- list('CyTOF_exp'   = CyTOF_exp, 
                  'CyTOF_prcnt' = CyTOF_prcnt, 
                  'RNA_topgenes'= RNA_top12K_E, 
                  'RNA_xcell'   = RNA_xcell,
                  'RNA_viper4' = RNA_viper4,
                  'mutation'    = mut_dt)

mofa_data4 <- list('CyTOF_exp'   = CyTOF_exp, 
                  'CyTOF_prcnt' = CyTOF_prcnt, 
                  'RNA_topgenes'= RNA_top12K_E, 
                  'RNA_epidish' = RNA_epidish,
                  'RNA_viper4' = RNA_viper4,
                  'mutation'    = mut_dt)

mofa_simp <- list('CyTOF_exp'   = CyTOF_exp, 
                  'RNA_topgenes'= RNA_top12K_E, 
                  'mutation'    = mut_dt)

mofa_simp2 <- list('CyTOF_prcnt' = CyTOF_prcnt,
                  'RNA_topgenes'= RNA_top12K_E, 
                  'mutation'    = mut_dt)

mofa_simp3 <- list('CyTOF_exp'   = CyTOF_exp,
                  'CyTOF_prcnt' = CyTOF_prcnt, 
                  'RNA_topgenes'= RNA_top12K_E, 
                  'mutation'    = mut_dt)


MOFA_m_x25 <- MOFA_model(mofa_data1)
MOFA_m_e25 <- MOFA_model(mofa_data2)
MOFA_m_x4 <- MOFA_model(mofa_data3)
MOFA_m_e4 <- MOFA_model(mofa_data4)



save(MOFA_m_x25, MOFA_m_e25, MOFA_m_x4, MOFA_m_e4, 
    file='/Users/senosam/Documents/Massion_lab/data_integration/MOFA_model.RData')

MOFA_simple <- MOFA_model(mofa_simp, mae = F, CDE_TMA36)
MOFA_simple2 <- MOFA_model(mofa_simp2, mae = F, CDE_TMA36)
MOFA_simple3 <- MOFA_model(mofa_simp3, mae = F, CDE_TMA36)

save(MOFA_simple, MOFA_simple2, MOFA_simple3,
    file='/Users/senosam/Documents/Massion_lab/data_integration/MOFA_model_simple.RData')


# Model including radiomics
radiomics <- t(scale(read.csv('/Users/senosam/Documents/Massion_lab/data_integration/ML/m_HM.csv', row.names = 1)))
radiomics <- na.omit(radiomics)

add_pts <- function(dt, rad_df) {
  dt <- t(dt)
  rad_df <- t(rad_df)
  k<- match(rownames(rad_df), rownames(dt))
  dt <- dt[k,]
  rownames(dt) <- rownames(rad_df)
  dt <- t(dt)
  dt
}

x <- add_pts(CyTOF_prcnt, radiomics)

CyTOF_prcnt <- add_pts(CyTOF_prcnt, radiomics)
RNA_top12K_E <- add_pts(RNA_top12K_E, radiomics)
mut_dt <- add_pts(mut_dt, radiomics)


mofa_simp4 <- list('CyTOF_prcnt' = CyTOF_prcnt,
                  'RNA_topgenes'= RNA_top12K_E, 
                  'mutation'    = mut_dt,
                  'radiomics'   = radiomics)

MOFA_simple4 <- MOFA_model(mofa_simp4, mae = F)
save(MOFA_simple4,
    file='/Users/senosam/Documents/Massion_lab/data_integration/MOFA_model_simpleRAD.RData')

