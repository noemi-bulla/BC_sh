library(edgeR)
library(ggplot2)
library(pheatmap)
library(ggfortify)
library(tidyverse)
library(openxlsx)
library(readxl)
library(limma)
library(sva)


be<- null.merged.gene_counts.SummarizedExperiment
counts <- be@assays@data@listData[["salmon.merged.gene_counts"]] #counts_paep<- counts[, grep("PAEP",colnames(counts), ignore.case= TRUE)]
sample<-be@colData@rownames
cond <- sub("_[^_]+$", "", sample)


condition <- factor(ifelse(cond %in% c("PT_shPAEP1", "PT_shPAEP2"), "PT_shPAEP", cond),
                    levels = c("PT_shSCR", "PT_shPAEP"))
condition_3 <- factor(cond, levels= c("PT_shSCR","PT_shPAEP1", "PT_shPAEP2")) 

y <- DGEList(counts = counts, group = condition) 

### change ensg name ###
gene_annot<- data.frame(
  ensg = be@elementMetadata@listData[["gene_id"]], 
  genename = be@elementMetadata@listData[["gene_name"]]
)

gene_annot_top_genes <- gene_annot[gene_annot$ensg %in% rownames(y$counts),]
gene_names_matched <- gene_annot_top_genes$genename[match(rownames(y$counts),
                                                          gene_annot_top_genes$ensg)]
rownames(y$counts) <- gene_names_matched

keep <- filterByExpr(y)
y<- y[keep, ,keep.lib.sizes=FALSE]
y<- calcNormFactors(y)


### DE edgeR ###
design<-model.matrix(~ 0 + condition)
y<- estimateDisp(y,design=design)

fit <- glmQLFit(y, design)

SCRvsPAEP <- makeContrasts(conditionPT_shPAEP - conditionPT_shSCR,
                           levels = design)

qlf<-glmQLFTest(fit, contrast=SCRvsPAEP)
t<-topTags(qlf, n=Inf)
tt<- t$table

ttt<-tt[tt$FDR<=0.1,]
zzz<-tt[tt$PValue<=0.05,]
ttt<- ttt %>% arrange(desc(logFC))
zzz<-zzz %>% arrange(desc(logFC))