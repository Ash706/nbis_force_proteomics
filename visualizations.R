
## Heatmap

```{r, fig8, fig.height=18, fig.width=12, echo=FALSE, message=FALSE}
library(gplots)
#basal.vs.lp.topgenes <- basal.vs.lp$ENTREZID[1:100]
#i <- which(v$genes$ENTREZID %in% basal.vs.lp.topgenes)
samples<-exprs(pre_processed) %>% colnames()
mycol <- colorpanel(1000,"blue","white","red")

png("../NBIS/Results/Normalized_NoGenotype/HeatMap_OXBaf_72_genes_noDand.png", units = "cm", height = 18, width = 12, res = 300)

heatmap.2(exprs(pre_processed)[OXBaf_genes,], scale="row",
          labRow=pre_processed@featureData@data[OXBaf_genes,"SYMBOL"],# labCol=as.character(Sum_Expe_filtered$group_name), 
          col=mycol, trace="none", density.info="none", cexRow = 0.15+ 0.5/log10(length(OXBaf_genes)),
          cexCol = 0.15 + 0.5/log10(length(samples)),offsetRow = 2,
          offsetCol = 0, key = TRUE,lhei=c(1, 10), lwid = c(1,3),
          keysize = 3,key.title = NA,
          key.xlab = NA,
          margin=c(4,4), adjRow = T,adjCol = T, dendrogram="none", labCol = pre_processed$treatmen,Colv=T
)

dev.off()
```

## Volcano plot
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
v1 <- EnhancedVolcano(results_DE_Baf$OxBafvsControl,
                      lab = results_DE_Baf$OxBafvsControl$SYMBOL,
                      x = 'logFC',
                      y = 'adj.P.Val',
                      pCutoff = 0.05,
                      ylim = c(0, max(-log10(results_DE_Baf$OxBafvsControl$adj.P.Val), na.rm=TRUE) + 1),
                      FCcutoff = 1,
                      #transcriptPointSize = 1.5,
                      #transcriptLabSize = 3, 
                      title = "OxBaf. Baf Control", subtitle = "", shadeBins = 4)


## Piechart

ggpie {ggpubr}
df <- data.frame(
    group = c("Male", "Female", "Child"),
    value = c(25, 25, 50))

head(df)
ggpie(df, "value", label = "group", color = c("red", "green", "blue"), fill = c("red", "green", "blue"))



