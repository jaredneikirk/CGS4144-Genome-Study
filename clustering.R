rm(list = ls()) 

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/factoextra")

library(factoextra)
library(cluster)

# Read original expression data
df <- read.csv(file = 'e_data.csv')
expr_Data <- df[-1]
row.names(expr_Data) <- df$SYMBOL

# Read variable expression data
variable_data <- read.csv(file = 'significant_difference.csv')
# Sort by pvalue
variable_data <- variable_data[order(variable_data$pvalue, decreasing = FALSE), ]
# 5000 variable genes
var5k <- variable_data[0:5000,]
expr_data_5k <- expr_Data[var5k$Gene,]
write.csv(expr_data_5k, "5000 most variable genes.csv")
# 10K variable genes
var10k <- variable_data[0:10000,]
expr_data_10k <- expr_Data[var10k$Gene,]
# 1000 variable genes
var1k <- variable_data[0:1000,]
expr_data_1k <- expr_Data[var1k$Gene,]
# 100 variable genes
var100 <- variable_data[0:100,]
expr_data_100 <- expr_Data[var100$Gene,]
# 10 variable genes
var10 <- variable_data[0:10,]
expr_data_10 <- expr_Data[var10$Gene,]

# Scale and center data
scaledata <- scale(t(expr_data_5k))
scaledata1 <- scale(t(expr_data_10k))
scaledata2 <- scale(t(expr_data_1k))
scaledata3 <- scale(t(expr_data_100))
scaledata4 <- scale(t(expr_data_10))

# decide k value
fviz_nbclust(scaledata, kmeans, method = "wss", k.max=10)

# Perform k means when k=2
set.seed(1)
km <- kmeans(scaledata, centers = 2, nstart = 25)
#plot results of final k-means model
fviz_cluster(km, data = scaledata)

cat_5k <- km[[1]]
row <- rownames(scaledata)
cat_5k <- do.call(rbind, Map(data.frame, A=row, B=cat_5k))
rownames(cat_5k) <- 1:nrow(cat_5k)
colnames(cat_5k) <- c('Sample', 'Cluster')

# alluvial plot
install.packages("ggalluvial")
library(ggplot2)
library(ggalluvial)

ggplot(data = cat_5k,
       aes(axis1 = Sample, axis2 = Cluster, y = 1)) + 
  scale_x_discrete(limits = c("Sample", 'Cluster')) + 
  geom_alluvium(aes(fill = Cluster), show.legend = FALSE) + 
  geom_stratum() + aes(label = after_stat(stratum)) + 
  theme_void() + ggtitle("5000 genes clusters") + 
  theme(plot.title = element_text(hjust = 0.5, size = 20))


# Heatmap
install.packages('pheatmap')
library(pheatmap)
mydf2 <- read.csv(file = 'f_data.csv')
anno <- mydf2[-1]
row.names(anno) <- mydf2$PatientID
anno$Cluster <- cat_5k$Cluster


pheatmap(scaledata,
         cluster_cols = T, cluster_rows = T, scale = "none",
         treeheight_row = 40, 
         treeheight_col = 50,
         annotation_legend = T,
         display_numbers = F,
         show_rownames = F,
         show_colnames = F,
         annotation_row = anno,
         border_color = "black",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(40))













