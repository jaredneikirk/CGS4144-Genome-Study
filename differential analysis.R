rm(list = ls()) 

if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("apeglm", update = FALSE)
}
if (!("devtools" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("devtools", update = FALSE)
}
if (!("ComplexHeatmap" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("ComplexHeatmap", update = FALSE)
}

library(DESeq2)
library(ggplot2)
library(magrittr)

mydf <- read.csv(file = 'e_data.csv')
exprsData <- mydf[-1]
row.names(exprsData) <- mydf$SYMBOL

# Drop rows with genes less than ten
exprsData <- exprsData %>%
 dplyr::filter(rowSums(.) >= 10)
# Round the data
exprsData <- round(exprsData)
# Convert NA to 0 if needed
exprsData[is.na(exprsData)] <- 0
# Convert negative values to 0 if needed
exprsData[exprsData<0] <- 0

mydf2 <- read.csv(file = 'f_data.csv')
feature_data <- mydf2[-1]
row.names(feature_data) <- mydf2$PatientID

ddset <- DESeqDataSetFromMatrix(exprsData, feature_data, ~Symptoms)

deseq_object <- DESeq(ddset)
deseq_results <- results(deseq_object)
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 2, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)
# this is of class DESeqResults -- we want a data frame
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))

readr::write_tsv(
  deseq_df,
  file.path("diff_expr_results.tsv"
  )
)

volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)

ggsave(
  plot = volcano_plot,
  file.path("volcano_plot.png")
)

diffsig <- deseq_df[(deseq_df$pvalue < 0.5 & abs(deseq_df$log2FoldChange) > 0),]
write.csv(diffsig, "significant_difference.csv")
# Standardize
vsd <- vst(ddset, blind = FALSE)
normalizeExp <- assay(vsd)

library(ComplexHeatmap)
# HeatMap
df <- read.csv("significant_difference.csv", header = T)
df02 <- as.character(df$Gene)
diff_expr <- normalizeExp[df02,]
ac <- subset(mydf2, select=-c(PatientID))
rownames(ac) <- colnames(diff_expr)

pdf("heatmap.pdf", height = 8, width = 12)
p <- pheatmap(diff_expr,
              annotation_col = ac,
              color = colorRampPalette(c("green","black","red"))(50),
              cluster_cols = F,
              show_rownames = F,
              show_colnames = F,
              scale = "row", ## none, row, column
              fontsize = 12,
              fontsize_row = 12,
              fontsize_col = 6,
              border = FALSE)
print(p)
dev.off()




