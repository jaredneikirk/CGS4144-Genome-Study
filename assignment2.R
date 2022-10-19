rm(list = ls()) 
options(stringsAsFactors = FALSE)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")
BiocManager::install("GEOquery")
library(GEOquery)
library(org.Hs.eg.db)

# retrieve meta data
# gset <- getGEO('GSE198449', GSEMatrix = FALSE, destdir = ".")
# class(gset)
# metaData <- Meta(gset)[c("title", "type", "platform_id", "summary",
            "supplementary_file")]

# Download supplement files
supp_info <- getGEOSuppFiles(GEO = "GSE198449")
gunzip('GSE198449/GSE198449_featureCounts.txt.gz')

# Store feature counts data into a table
data1 <- read.table('GSE198449/GSE198449_featureCounts.txt', sep="", header = TRUE)
# Keep only the Ensembl id without decimal
data1$Geneid <- gsub("\\..*","", data1$Geneid)
# Rename Geneid to ENSEMBL
colnames(data1)[1] <- "ENSEMBL"

# Map ENSEMBL ID to corresponding gene name
x <- data1$ENSEMBL
get <- c("SYMBOL")
result <- select(org.Hs.eg.db, keys = x, columns = get, keytype = "ENSEMBL")
result1 <- merge(result, data1, by = "ENSEMBL")
# Drop N/A value
result1 <- na.omit(result1)
# Get final data
data <- subset(result1, select=-c(ENSEMBL))
write.csv(data,"Exprs.csv", row.names = TRUE)

# log scale data
data[, 2:1859] <- log2(data[, 2:1859])
# data2 <- data[!is.infinite(rowSums(data[, 2:1859])),]
# Replace Inf in data by NA
data <- do.call(data.frame,lapply(data, function(x) replace(x, is.infinite(x), NA)))


# Plot density plot
ranges = apply(data[,2:1859], 1, max, na.rm=T) - apply(data[,2:1859], 1, min, na.rm=T)
plot(density(ranges))
