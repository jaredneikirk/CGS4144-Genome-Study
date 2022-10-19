rm(list = ls()) 

BiocManager::install("DESeq2")
library(DESeq2)

exprsData <- read.csv(file = 'e_data.csv')

# Drop first column
exprsData <- subset(exprsData, select=-c(SYMBOL))
# Drop rows with genes less than ten
# exprsData <- ExprsData %>%
  # dplyr::filter(rowSums(.) >= 10)
# Round the data
exprsData <- round(exprsData)
# Convert NA to 0 if needed
exprsData[is.na(exprsData)] <- 0
# Convert negative values to 0 if needed
exprsData[exprsData<0] <- 0

condition <- read.csv(file = 'f_data.csv')
# Drop first column
condition <- subset(feature_data, select=-c(PatientID))
dds <- DESeqDataSetFromMatrix(exprsData, condition, ~Symptoms)

vst_data <- vst(dds)
plotPCA(vst_data, intgroup="Symptoms")
