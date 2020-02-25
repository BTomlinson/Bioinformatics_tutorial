# Bioinformatics_tutorial DESeq2

#First we need some tools:

#Install htmltools
install.packages("htmltools")
library(htmltools)


#Install BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")


#Use BiocManager to install DESeq2
BiocManager::install("DESeq2")


#Call out DESeq2 and ggplot2
library("DESeq2")
library(ggplot2)

# Lets get our count data
#Save the URL of the counts table from github as ReadCountsURL
ReadCountsURL <- "https://raw.githubusercontent.com/BTomlinson/Bioinformatics_tutorial/master/AB_UTvT.csv"


#Download the ReadCountsURL file and save it as AB_UTvT.csv
download.file(ReadCountsURL, destfile = "AB_UTvT.csv", method = "auto")


#Read the downloaded file (AB_UTvT.csv) and save it as ReadCountsTable
ReadCountsTable <- read.csv('AB_UTvT.csv', header = TRUE, sep = ",")


#Lets look at the top of the table we just made
head(ReadCountsTable)


# Lets define our data using a metadatatable
#now we need to tell R what sample belongs to which experimental group. In this case, we have a control group (untreated, UT) and a group exposed to antibiotic (treated, T). This information is defined by the metaData table


#Save the URL of the metadata table from github as metaDataURL
metaDataURL <- "https://raw.githubusercontent.com/BTomlinson/Bioinformatics_tutorial/master/AB_UTvT_metadata.csv"


#Download the ReadCountsURL file and save it as AB_UTvT_metadata.csv
download.file(metaDataURL, destfile = "AB_UTvT_metadata.csv", method = "auto")


#Read the downloaded file (AB_UTvT_metadata.csv) and save it as metaDataTable
metaDataTable <- read.csv('AB_UTvT_metadata.csv', header = TRUE, sep = ",")


#Let's take a look at our table
metaDataTable
#The main thing to note about this file is each sample is matched to the appropriate experimental condition (untreated or treated)
#You are also able to define different experimental conditions or even GEO_ID tags in one column should you wish - for now we will simply focus on untreated vs. treated

# Now we can start DESeq
#Now we need to make the appropriate DESeqDataSet object. Heres how we can construct that from what we just did above:
DataSet <- DESeqDataSetFromMatrix(countData=ReadCountsTable,
                              colData=metaDataTable,
                              design=~condition, tidy = TRUE)
#Design specifies how the counts from each gene depend on our variables in the metadata.
#For this dataset, the factor we care about is our treatment status (condition)
#tidy=TRUE tells DESeq2 to output the results table with rownames as a first #column called 'row.
 
 
#OKAY, lets check out the object we just made
DataSet                
#This is a good way to make sure your object is constructed correctly. 
#As you can see, we made a DESeqDataSet. 
#dim: 3981 6 tells us there are 3981 rows and 6 columns
#rownames shows us our rows represent each gene name 
#colnames shows us each column represents each experimental sample


#If you would like some light reading on DESeq2 and all of the arguments...
?DESeq


#Excellent speed reading. Now that you are DESeq2 expert, we are ready to run the DESeq function on our DataSet!
DEDataSet <- DESeq(DataSet)
#estimateSizeFactors calculates the relative library depth of each sample 
#estimateDispersions estimates the dispersion of counts for each gene 
#nbinomWaldTest calculates the significance of coefficients in a Negative Binomial GLM using the size and dispersion outputs


#lets save these results as DESeqResults
DESeqResults <- results(DEDataSet)
#and take a quick peek at them
head(results(DEDataSet, tidy=TRUE))


#You can also get a quick summary of the results table, this is useful for identifying outliers too! But we have a much prettier way of doing this coming up...
summary(DESeqResults)


#Lets re-order our resuls to the most interesting, er I mean significant results are a the top. This sorts the results by adjusted p-value
DESeqResults <- DESeqResults[order(DESeqResults$padj),]
head(DESeqResults)


# Some useful visualizations
#we can use the plotCounts function to compare the normalized counts between untreated and treated conditions for a few genes of interest
par(mfrow=c(2,3))
plotCounts(DEDataSet, gene="adeA", intgroup="condition")
plotCounts(DEDataSet, gene="adeB", intgroup="condition")
plotCounts(DEDataSet, gene="adeC", intgroup="condition")
plotCounts(DEDataSet, gene="ABUW_1931", intgroup="condition")
plotCounts(DEDataSet, gene="ABUW_1318", intgroup="condition")


#To further decipher the data, I recommend blasting the genes with significantly altered expression to determine what their function is! Or, look into additional tools like PANTHER (http://www.pantherdb.org) or KEGG (https://www.kegg.jp)!


#reset par
par(mfrow=c(1,1))


#How about we make a basic volcano plot to view the differential expression spread?
with(DESeqResults, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))


#We can color these points to see where our significant results lie
#blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(DESeqResults, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(DESeqResults, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))


#You can do a Principal Component Analysis to see how your samples group by condition
#This brings out strong patterns from large and complex datasets and samples with similar expression profiles will cluster together
#For a great, simplified explanation on PCA - https://blog.bioturing.com/2018/06/14/principal-component-analysis-explained-simply/
#First we transform the raw count data, use the vst function to perform variance stabilizing transformation
VSData <- vst(DESeqResults, blind=FALSE)


#Then, use the DESEQ2 plotPCA function to plot Principal Component Analysis to see how these samples group by condition
plotPCA(VSData, intgroup="condition")
