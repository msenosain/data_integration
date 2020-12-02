# Loading libraries
library(dplyr)
library(tidyverse)
library(tidyr)
library(MOFA)
library(MOFAdata)
library(reticulate)
library(reticulate)
use_python("/Users/senosam/opt/anaconda3/bin/python", required=TRUE) 

# Loading data
load('/Users/senosam/Documents/Massion_lab/data_integration/MO_data.Rdata')

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


MOFA_model <- function(mofa_data){
    set.seed(1234)
    MOFAobject <- createMOFAobject(mofa_data)

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

MOFA_m_x25 <- MOFA_model(mofa_data1)
MOFA_m_e25 <- MOFA_model(mofa_data2)
MOFA_m_x4 <- MOFA_model(mofa_data3)
MOFA_m_e4 <- MOFA_model(mofa_data4)
MOFA_simple <- MOFA_model(mofa_simp)




save(MOFA_m_x25, MOFA_m_e25, MOFA_m_x4, MOFA_m_e4, 
    file='/Users/senosam/Documents/Massion_lab/data_integration/MOFA_model.RData')

