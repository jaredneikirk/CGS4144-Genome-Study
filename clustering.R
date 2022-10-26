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
km1 <- kmeans(scaledata1, centers = 2, nstart = 25)
km2 <- kmeans(scaledata2, centers = 2, nstart = 25)
km3 <- kmeans(scaledata3, centers = 2, nstart = 25)
km4 <- kmeans(scaledata4, centers = 2, nstart = 25)
#plot results of final k-means model
fviz_cluster(km, data = scaledata)

cat_5k <- km[[1]]
row <- rownames(scaledata)
cat_5k <- do.call(rbind, Map(data.frame, A=row, B=cat_5k))
rownames(cat_5k) <- 1:nrow(cat_5k)
colnames(cat_5k) <- c('Sample', 'Cluster')

cat_10k <- km1[[1]]
row <- rownames(scaledata1)
cat_10k <- do.call(rbind, Map(data.frame, A=row, B=cat_10k))
rownames(cat_10k) <- 1:nrow(cat_10k)
colnames(cat_10k) <- c('Sample', 'Cluster')

cat_1k <- km2[[1]]
row <- rownames(scaledata2)
cat_1k <- do.call(rbind, Map(data.frame, A=row, B=cat_1k))
rownames(cat_1k) <- 1:nrow(cat_1k)
colnames(cat_1k) <- c('Sample', 'Cluster')

cat_100 <- km3[[1]]
row <- rownames(scaledata3)
cat_100 <- do.call(rbind, Map(data.frame, A=row, B=cat_100))
rownames(cat_100) <- 1:nrow(cat_100)
colnames(cat_100) <- c('Sample', 'Cluster')

cat_10 <- km4[[1]]
row <- rownames(scaledata4)
cat_10 <- do.call(rbind, Map(data.frame, A=row, B=cat_10))
rownames(cat_10) <- 1:nrow(cat_5k)
colnames(cat_10) <- c('Sample', 'Cluster')

mydf2 <- read.csv(file = 'f_data.csv')
anno <- mydf2[-1]
anno$Cluster <- cat_5k$Cluster
anno1 <- mydf2[-1]
anno1$Cluster <- cat_10k$Cluster
anno2 <- mydf2[-1]
anno2$Cluster <- cat_1k$Cluster
anno3 <- mydf2[-1]
anno3$Cluster <- cat_100$Cluster
anno4 <- mydf2[-1]
anno4$Cluster <- cat_10$Cluster

# alluvial plot
install.packages("ggalluvial")
library(ggplot2)
library(ggalluvial)

ggplot(data = anno4,
       aes(axis1 = Symptoms, axis2 = Cluster)) + 
  geom_alluvium(aes(fill = Cluster)) + 
  geom_stratum() + 
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) + 
  theme_void()


# Heatmap
install.packages('pheatmap')
library(pheatmap)
row.names(anno) <- mydf2$PatientID

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


# Statistics
chisq.test(anno$Symptoms, anno$Cluster, correct = FALSE)










