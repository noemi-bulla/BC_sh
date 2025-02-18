library(tximport)
library(DESeq2)
library(readr)
library(tidyr)
library(dplyr)
library(sva)
library(edgeR)
library(tibble)
library(ggplot2)
library(ggfortify)

sample<- factor(c("PT_shSCR_1","PT_shSCR_2","PT_shSCR_3",
                                  "PT_shSCR_4","PT_shSCR_5","PT_shPAEP1_1",
                                  "PT_shPAEP1_2","PT_shPAEP1_3","PT_shPAEP1_4",
                                  "PT_shPAEP1_5","PT_shPAEP2_1","PT_shPAEP2_2",
                                  "PT_shPAEP2_3","PT_shPAEP2_4")) 

condition <- factor(ifelse(grepl("shPAEP1", sample), "PT_shPAEP1",
                           ifelse(grepl("shPAEP2", sample), "PT_shPAEP2",
                                  "PT_shSCR")),
                    levels = c("PT_shPAEP1", "PT_shPAEP2", "PT_shSCR"))

condition_2 <- factor(ifelse(grepl("shPAEP1", sample), "PT_shPAEP",
                             ifelse(grepl("shPAEP2", sample), "PT_shPAEP",
                                    "PT_shSCR")),
                      levels = c("PT_shPAEP", "PT_shSCR"))

cell_cycle<- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4",
               "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL",     
               "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2", "RPA2",   
               "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2",
               "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2",
               "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
               "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1",
               "E2F8", "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
               "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
               "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
               "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",  
               "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
               "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
               "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
               "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA")

sampletable <- data.frame(
  sample = sample,
  condition = condition  
)

files<- file.path(sample, "quant.sf")
names(files)<- sample
tx2gene<- read.delim("tx2gene.tsv", header=FALSE, stringsAsFactors = FALSE)
colnames(tx2gene) <- c("transcript_id", "gene_id", "gene_name")
txi<-tximport(files, type="salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)
dds <- DESeqDataSetFromTximport(txi, sampletable, ~condition)
counts<- dds@assays@data@listData[["counts"]]

gene_id_to_name <- tx2gene %>% distinct(gene_id, gene_name) %>% column_to_rownames("gene_id")
rownames(counts) <- gene_id_to_name[rownames(counts), "gene_name"]

# Check for NA values (some gene_ids may not have a corresponding gene_name)
counts <- counts[!is.na(rownames(counts)), ]

y <- DGEList(counts = counts, group = condition)
y$counts

keep <- filterByExpr(y)
y<- y[keep, ,keep.lib.sizes=FALSE]
y<- calcNormFactors(y)

### PCA ### 
cv_function <- function(x) {
  if (abs(mean(x)) <= 1e-8) return(0)  
  return(sd(x) / mean(x))
}
st_var <- function(x) {
  (x - mean(x)) / sd(x)
}


logCPM <- t(cpm(y, log=TRUE, prior.count = 2))
logCPM_t<- t(logCPM)

cv_values_cpm <- apply(logCPM, 2, cv_function)
cv_results_cpm <- data.frame(gene = colnames(logCPM), CV = cv_values_cpm)
top_cv_genes_cpm <- cv_results_cpm[order(-cv_results_cpm$CV), ][1:1000, ]  
pca_df_subset_cpm <- logCPM[, top_cv_genes_cpm$gene, drop=FALSE]  

st_var_values_cpm <- apply(pca_df_subset_cpm, 2, st_var)
corrected_logCPM_t<- t(corrected_logCPM)
pca<-prcomp(corrected_logCPM_t, scale. = TRUE)

autoplot(pca, data = data.frame(Group = condition), colour = "Group") + 
  ggtitle("PCA of RNA-seq Samples") +
  theme_minimal()

pc1_values <- pca$x[,1]
num_reads<-colSums(y$counts)

cell_cycle_reads <- logCPM_t[rownames(logCPM_t) %in% cell_cycle, ]
mean_cell_cycle_reads <- colMeans(cell_cycle_reads, na.rm=TRUE)

cor_pc1_cellcycle <- cor(pc1_values, mean_cell_cycle_reads, method="pearson")
cor_pc1_totalreads <- cor(pc1_values, num_reads, method="pearson")
batch <- factor(c(rep("PT_shSCR", 5), rep("PT_shPAEP1", 5), rep("PT_shPAEP2", 4)))
batch <- sample
covariates <- data.frame(total_reads=num_reads, cell_cycle=mean_cell_cycle_reads)
corrected_logCPM <- ComBat(dat=logCPM_t, batch=batch, mod=model.matrix(~ cell_cycle))
corrected_logCPM_num_reads <- ComBat(dat=logCPM_t, batch=, mod=model.matrix(~ total_reads, data=covariates))
subset_cc<- logCPM_t[rownames(logCPM_t) %in% cell_cycle, ]
mean_expression <- colMeans(subset_cc)

cv_values_cpm <- apply(logCPM, 2, cv_function)
cv_results_cpm <- data.frame(gene = colnames(logCPM), CV = cv_values_cpm)
top_cv_genes_cpm <- cv_results_cpm[order(-cv_results_cpm$CV), ][1:1000, ]  
pca_df_subset_cpm <- logCPM[, top_cv_genes_cpm$gene, drop=FALSE]  

st_var_values_cpm <- apply(pca_df_subset_cpm, 2, st_var)

column_variances <- apply(st_var_values_cpm, 2, var)
column_measn<- colMeans(st_var_values_cpm)

num_reads<- colSums(y$counts)

pca <- prcomp(st_var_values_cpm, scale = TRUE)
glimpse(pca)
pc1_scores<-pca$x[,1]


pairwise_correlation <- data.frame(
  PCA_Score = pc1_scores,
  cell_cycle_mean = mean_expression,
  Correlation = cor(pca$x[,1], mean_expression)
)

print(pairwise_correlation)

###combat ###
covariates <- data.frame(
  Tot_Reads = num_reads,      
  Mean_Expression = mean_expression
)
rownames(covariates) <- colnames(logCPM_t)
mod_matrix <- model.matrix(~ Tot_Reads + Mean_Expression, data = covariates)
combat_adjusted <- ComBat(dat = logCPM_t, batch = NULL, mod = modcombat)

pca_plot<-ggplot(pca, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +  
  ggtitle("PCA of RNA-seq Samples") +
  theme_bw()


### DE edgeR ###
design<-model.matrix(~ 0 + condition_2)
y<- estimateDisp(y,design=design)

fit <- glmQLFit(y, design)

SCRvsPAEP <- makeContrasts(condition_2PT_shPAEP - condition_2PT_shSCR,
                           levels = design)

qlf<-glmQLFTest(fit, contrast=SCRvsPAEP)
t<-topTags(qlf, n=Inf)
tt<- t$table

ttt<-tt[tt$FDR<=0.1,]
ttt<- ttt %>% arrange(desc(logFC))

### PAEP1 vs PAEP2 ### 
y_3 <- DGEList(counts = counts, group = condition)

### change ensg name ###
gene_annot<- data.frame(
  ensg = be@elementMetadata@listData[["gene_id"]], 
  genename = be@elementMetadata@listData[["gene_name"]]
)

gene_annot_top_genes <- gene_annot[gene_annot$ensg %in% rownames(y_3$counts),]
gene_names_matched <- gene_annot_top_genes$genename[match(rownames(y_3$counts),
                                                          gene_annot_top_genes$ensg)]
rownames(y_3$counts) <- gene_names_matched

keep <- filterByExpr(y_3)
y_3<- y_3[keep, ,keep.lib.sizes=FALSE]
y_3<- calcNormFactors(y_3)
design_3<-model.matrix(~ 0 + condition)
y_3<- estimateDisp(y,design=design_3)

fit <- glmQLFit(y_3, design_3)

PAEP1vsPAEP2 <- makeContrasts(conditionPT_shPAEP1 - conditionPT_shPAEP2,
                              levels = design_3)

qlf<-glmQLFTest(fit, contrast=PAEP1vsPAEP2)
x_3<-topTags(qlf, n=Inf)
xx_3<- x_3$table

xxx_3<-xx_3[xx_3$FDR<=0.1,]
xxx_3<-xxx_3 %>% arrange(desc(logFC))

### PAEP1 vs SCR ###
PAEP1vsSCR <- makeContrasts(conditionPT_shPAEP1 - conditionPT_shSCR,
                            levels = design_3)

qlf<-glmQLFTest(fit, contrast=PAEP1vsSCR)
t_3<-topTags(qlf, n=Inf)
tt_3<- t_3$table

ttt_3<-tt_3[tt_3$FDR<=0.1,]
ttt_3<-ttt_3 %>% arrange(desc(logFC))

### PAEP2vsSCR ###
PAEP2vsSCR <- makeContrasts(conditionPT_shPAEP2 - conditionPT_shSCR,
                            levels = design_3)

qlf<-glmQLFTest(fit, contrast=PAEP2vsSCR)
n_3<-topTags(qlf, n=Inf)
nn_3<- n_3$table

nnn_3<-nn_3[nn_3$FDR<=0.1,]
nnn_3<-nnn_3 %>% arrange(desc(logFC))

### Deseq2 ###
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds_filtered <- dds[keep,]
dds_final <- DESeq(dds_filtered)
rld <- rlog(dds_final, blind = TRUE)
pca_plot <- plotPCA(rld, intgroup = "condition")
