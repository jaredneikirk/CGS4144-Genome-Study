library(topGO)

# Annotate gene data and select subsample of signficiant genes
GOdata <- new("topGOdata",
              ontology='BP',
              allGenes=genes,
              geneSel=function(v){v<0.05}, # select statistically significant values
              nodeSize=10,
              annotationFun=annFUN.org,
              mapping="org.Hs.eg.db",
              ID="SYMBOL")

#  use Kolmogorov-Smirnov (K-S) test to make use of rank info:
# (suggested here: https://www.biostars.org/p/350710/)

results.ks <- runTest(GOdata, algorithm="classic", statistic="ks")
goEnrichment <- GenTable(GOdata, KS=results.ks, orderBy="KS", topNodes=20)
goEnrichment$KS <- as.numeric(goEnrichment$KS)
goEnrichment <- goEnrichment[goEnrichment$KS<0.05,]
goEnrichment <- goEnrichment[,c("GO.ID","Term","KS")]
goEnrichment$Term <- gsub(" [a-z]*\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- gsub("\\.\\.\\.$", "", goEnrichment$Term)
goEnrichment$Term <- paste(goEnrichment$GO.ID, goEnrichment$Term, sep=", ")
goEnrichment$Term <- factor(goEnrichment$Term, levels=rev(goEnrichment$Term))

# Plot the results

library(ggplot2)
ggplot(goEnrichment, aes(x=Term, y=-log10(KS))) +
  stat_summary(geom = "bar", fun.y = mean, position = "dodge") +
  xlab("Biological process") +
  ylab("Enrichment") +
  ggtitle("Title") +
  scale_y_continuous(breaks = round(seq(0, max(-log10(goEnrichment$KS)), by = 2), 1)) +
  theme_bw(base_size=0.5) +
  theme(
    legend.position='none',
    legend.background=element_rect(),
    plot.title=element_text(angle=0, size=9, face="bold", vjust=1),
    axis.text.x=element_text(angle=0, size=9, face="bold", hjust=1.10),
    axis.text.y=element_text(angle=0, size=9, face="bold", vjust=0.5),
    axis.title=element_text(size=9, face="bold"),
    legend.key=element_blank(),   
    legend.key.size=unit(1, "cm"),
    legend.text=element_text(size=9), 
    title=element_text(size=9)) +
  guides(colour=guide_legend(override.aes=list(size=2.5))) +
  coord_flip()