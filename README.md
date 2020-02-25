# Bioinformatics_tutorial DESeq2

#First we need some tools:

install.packages("htmltools")
library(htmltools)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

BiocManager::install("DESeq2")

library( "DESeq2" )

library(ggplot2)

countsName <- "https://raw.githubusercontent.com/BTomlinson/Bioinformatics_tutorial/master/AB5075_UTvT.csv"

download.file(countsName, destfile = "AB5075_UTvT.csv", method = "auto")

countData <- read.csv('AB5075_UTvT.csv', header = TRUE, sep = ",")

head(countData)

metaDataName <- "https://raw.githubusercontent.com/BTomlinson/Bioinformatics_tutorial/master/AB5075_UTvT_metadata.csv"

download.file(metaDataName, destfile = "AB5075_UTvT_metadata.csv", method = "auto")

metaData <- read.csv('AB5075_UTvT_metadata.csv', header = TRUE, sep = ",")

metaData

dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=metaData,
                              design=~dex, tidy = TRUE)

dds                

dds <- DESeq(dds)

?DESeq

res <- results(dds)
head(results(dds, tidy=TRUE))

summary(res)

res <- res[order(res$padj),]
head(res)

par(mfrow=c(2,3))
plotCounts(dds, gene="adeA", intgroup="dex")
plotCounts(dds, gene="ABUW_1318", intgroup="dex")


#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
