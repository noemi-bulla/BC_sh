library("DESeq2")
library(ggplot2)
library(pheatmap)
library(ggfortify)
library(tidyverse)
library(openxlsx)
library(readxl)
library(limma)
library(sva)


be<- null.merged.gene_counts.SummarizedExperiment
counts <- be@assays@data@listData[["salmon.merged.gene_counts"]]
counts<- round(counts)#counts_paep<- counts[, grep("PAEP",colnames(counts), ignore.case= TRUE)]
sample<-be@colData@rownames
cond <- sub("_[^_]+$", "", sample)

condition <- factor(ifelse(cond %in% c("PT_shPAEP1", "PT_shPAEP2"), "PT_shPAEP", cond),
                    levels = c("PT_shSCR", "PT_shPAEP"))
condition_3 <- factor(cond, levels= c("PT_shSCR","PT_shPAEP1", "PT_shPAEP2")) 

colData <- data.frame(row.names = sample, condition = condition)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ condition)
dds

### change ensg name ###
gene_annot<- data.frame(
  ensg = be@elementMetadata@listData[["gene_id"]], 
  genename = be@elementMetadata@listData[["gene_name"]]
)

gene_annot_top_genes <- gene_annot[gene_annot$ensg %in% rownames(dds),]
gene_names_matched <- gene_annot_top_genes$genename[match(rownames(dds),
                                                          gene_annot_top_genes$ensg)]
rownames(dds) <- gene_names_matched

### 

smallestGroupSize <- 2
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

### DE analysis ###
dds <- DESeq(dds)
res <- results(dds, contrast=c("condition","PT_shPAEP","PT_shSCR"))
summary(res)

resLFC <- lfcShrink(dds, coef="condition_PT_shPAEP_vs_PT_shSCR", type="apeglm")
resOrdered <- resLFC[order(resLFC$padj),]
### visualization MA ###
plotMA(res, ylim=c(-2,2))
plotMA(resOrdered, ylim=c(-2,2)) #shrink apeglm
idx <- identify(res$baseMean, res$log2FoldChange) # if you want to find a gene in MA plot and recover id of gene

plotCounts(dds, gene="PAEP", intgroup="condition") #IMPORTANT , see condition with 3 groups


resSig <- subset(res, padj < 0.1)

rld <- rlog(dds, blind=FALSE)
resSig <- subset(rld@, padj < 0.1)
### PCA ###
plotPCA(rld, intgroup="condition")

### heat map visualization ###
top_genes <- res[order(abs(res$log2FoldChange), decreasing = TRUE), ]  # Sort by absolute log2FoldChange
select <- rownames(top_genes)[1:50]
df <- data.frame(condition = dds$condition)
rownames(df) <- colnames(assay(rld))
df$condition <- as.factor(df$condition)
pheatmap(assay(rld)[select,],
         cluster_rows=FALSE,
         show_rownames=TRUE,
         cluster_cols=TRUE,
         annotation_col=df,
         color = colorRampPalette(c("blue", "white", "red"))(50))

sampleDists <- dist(t(assay(rld))

