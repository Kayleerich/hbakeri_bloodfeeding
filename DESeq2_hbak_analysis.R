library(DESeq2)
library(dplyr) 
library(stringr)
library(PCAtools)

setwd("C:\\path\\to\\files")

## select subset from data, if FALSE: all data, if TRUE: excludes day 6 data
d7_only <- TRUE

## set default variables
out_file <- "hbak_PRJEB15396.WBPS18_DESeq2_LightvDark_results.tsv"
desc <- "Days 6 and 7, condition Light vs. Dark"
mat_file <- "hbak_PRJEB15396.WBPS18_featurecounts.Rmatrix.tsv"
info_file <- "sample_table.tsv"
levels <- c("Dark", "Light")


## sample info table
coldata <- read.delim(info_file)
row.names(coldata) <- coldata$identifier
coldata[1] <- NULL

## count data from featureCounts 
cts <- read.delim(mat_file)
row.names(cts) <- cts$Geneid
cts$Geneid <- NULL

## change variables and dataframe content if subsetting data 
if(d7_only == TRUE){
  cts <- as.matrix(cts %>% select(contains("D7")))
  coldata <- coldata %>% filter(str_detect(day, "D7"))
  out_file <- "hbak_PRJEB15396.WBPS18_DESeq2_LightvDark_results_D7Only.tsv"
  desc <- "Day 7 only, condition Light vs. Dark"
}

## check if sample identifiers from count data match sample info
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))


## make DESeq2 object, set 'Dark' as reference level
dds <- DESeqDataSetFromMatrix(countData = round(cts),
                              colData = coldata,
                              design = ~ condition)
dds$condition <- factor(dds$condition, levels = levels)
dds <- DESeq(dds)
dds <- dds[rowSums(counts(dds)) > 0,]


## get results, change pval cut-off to 0.05 instead of (default) 0.1
res <- results(dds, alpha=0.05)
print(desc)
summary(res)
# Day 7 only, condition: Light vs. Dark
# out of 23631 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 389, 1.6%
# LFC < 0 (down)     : 162, 0.69%
# outliers [1]       : 645, 2.7%
# low counts [2]     : 5495, 23%
# (mean count < 10)

# Days 6 and 7, condition: Light vs. Dark
# out of 23631 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 389, 1.6%
# LFC < 0 (down)     : 162, 0.69%
# outliers [1]       : 645, 2.7%
# low counts [2]     : 5476, 23%
# (mean count < 10)

resOrdered <- res[order(res$padj),]
head(resOrdered[ order(resOrdered$padj, decreasing = FALSE), ], 75)
write.table(as.data.frame(resOrdered), sep = "\t",
            file=out_file)


## get only results with significant padj, top 75 up- and down-regulated, top 75 differentially expressed
resSig <- subset(resOrdered, padj < 0.05)
downReg <- head(resSig[ order(resSig$log2FoldChange), ], 75)
upReg <- head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ], 75)
topDiffExpr <- head(resSig[ order(abs(resSig$log2FoldChange), decreasing = TRUE), ], 75)


# Plot sample variance: transform read counts to log2 scale while minimizing difference for rows  
# with small counts and taking differences between library sizes of the samples into account
rld <- rlog(dds, blind=FALSE)
plotPCA(rld, intgroup = "condition") 
