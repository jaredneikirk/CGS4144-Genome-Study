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

# Scale and center data
scaledata <- t(scale(t(expr_data_5k)))

# decide k value
fviz_nbclust(scaledata, kmeans, method = "wss", k.max=10)

# Perform k means when k=3
set.seed(1)
km <- kmeans(scaledata, centers = 3, nstart = 25)
#plot results of final k-means model
fviz_cluster(km, data = scaledata)
