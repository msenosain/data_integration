# Data compilation HEALTHMYNE
healthmyne <- read.csv("Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Data/Radiomics/TMA36_RAD_RAW - Cleaned - linked PT ID.csv")
healthmyne <- healthmyne[-which(healthmyne$MRN == ''),]
aneri <- read_excel("Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Data/Radiomics/MergeFile.Mafe's Final List HM Raw Data AB.xlsx")
aneri$...1 <- NULL
colnames(aneri)[1] <- 'pt_ID'
CDE <- read.csv(file = '/Users/senosam/Documents/Massion_lab/CyTOF_summary/CDE_TMA36_2020FEB25_SA.csv')

k <- match(healthmyne$StudyInstanceUID, aneri$StudyInstanceUID)
healthmyne <- cbind('pt_ID'= aneri$pt_ID[k], healthmyne)
healthmyne$PT.ID <- NULL

# Matching with CDE
k <- match(CDE$pt_ID, healthmyne$pt_ID)
healthmyne <- healthmyne[k,]
healthmyne$pt_ID <- CDE$pt_ID

write.csv(healthmyne, file = '/Users/senosam/Documents/Massion_lab/radiomics_summary/TMA36_HM.csv', row.names = F)


# CANARY
khushbu <- read_excel("Library/Mobile Documents/com~apple~CloudDocs/Documents/Vanderbilt/Massion_Lab/Projects/CyTOF_ADC_project/Data/Radiomics/ADC_TMA 36_ khushbu 10-02-2020.xlsx")
khushbu <- khushbu[match(CDE$pt_ID, khushbu$PatientID),]
CANARY <- khushbu[,c(1,7:10)]
colnames(CANARY) <- c('pt_ID', 'SILA_P', 'SILA_S', 'CANARY_P', 'CANARY_S')

write.csv(CANARY, file = '/Users/senosam/Documents/Massion_lab/radiomics_summary/TMA36_CANARY_khushbu.csv', row.names = F)
