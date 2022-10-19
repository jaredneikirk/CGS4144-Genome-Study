rm(list = ls()) 
options(stringsAsFactors = FALSE)

BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db) #人类注释数据库
library(clusterProfiler)#进行GO富集和KEGG富集
library(dplyr) #进行数据转换
library(ggplot2)#绘图

f <- read.csv('e_data.csv')
x <- f[,1]
hg <- bitr(x,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL","SYMBOL"),OrgDb="org.Hs.eg.db")

ggo <- groupGO(gene = hg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC",level = 3,readable = TRUE)
ego_ALL <- enrichGO(gene = hg$ENTREZID, 
                    OrgDb = org.Hs.eg.db, 
                    ont = "ALL", 
                    pAdjustMethod = "BH",
                    pvalueCutoff = 1,
                    qvalueCutoff = 1,
                    readable = TRUE)
ego_MF <- enrichGO(gene = hg$ENTREZID,OrgDb = org.Hs.eg.db,ont = "MF", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)
ego_MF1 <- setReadable(ego_MF, OrgDb = org.Hs.eg.db)
write.csv(summary(ego_ALL),"ALL-enrich.csv",row.names =FALSE)

barplot(ego_MF, showCategory=20,title="EnrichmentGO_MF")
dotplot(ego_MF,title="EnrichmentGO_MF_dot")

ego_ALL.sig <- ego_ALL[ego_ALL$pvalue <= 0.05]

