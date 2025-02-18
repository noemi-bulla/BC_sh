library(edgeR)
library(ggplot2)
library(pheatmap)
library(ggfortify)
library(tidyverse)
library(openxlsx)
library(readxl)
library(limma)
library(sva)
library(readr)



be<- null.merged.gene_counts.SummarizedExperiment
counts <- be@assays@data@listData[["salmon.merged.gene_counts"]] #counts_paep<- counts[, grep("PAEP",colnames(counts), ignore.case= TRUE)]
sample<-be@colData@rownames
cond <- sub("_[^_]+$", "", sample)

condition <- factor(ifelse(cond %in% c("PT_shPAEP1", "PT_shPAEP2"), "PT_shPAEP", cond),
                    levels = c("PT_shSCR", "PT_shPAEP"))
condition_3 <-factor(cond, levels= c("PT_shPAEP1", "PT_shPAEP2", "PT_shSCR"))
condition_all <- factor(sample, levels= c("PT_shSCR_1","PT_shSCR_2","PT_shSCR_3",
                                      "PT_shSCR_4","PT_shSCR_5","PT_shPAEP1_1",
                                      "PT_shPAEP1_2","PT_shPAEP1_3","PT_shPAEP1_4",
                                      "PT_shPAEP1_5","PT_shPAEP2_1","PT_shPAEP2_2",
                                      "PT_shPAEP2_3","PT_shPAEP2_4")) 

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

### PCA before DE ### 
logCPM <- cpm(y, log=TRUE)

pca <- prcomp(t(logCPM), scale = TRUE)
autoplot(pca, data = data.frame(Group = condition), colour = "Group") + 
  ggtitle("PCA of RNA-seq Samples") +
  theme_minimal()

### alternative PCA (CV) ###
pca_df<-t(y$counts)

cv_function <- function(x) {
  if (abs(mean(x)) <= 1e-8) return(0)  
  return(sd(x) / mean(x))
}

cv_values <- apply(pca_df, 2, cv_function)
cv_results <- data.frame(gene = colnames(pca_df), CV = cv_values)
top_cv_genes <- cv_results[order(-cv_results$CV), ]
top_2500_cv_genes <- top_cv_genes[1:2000, ]
pca_df_subset <- pca_df[, top_2500_cv_genes$gene, drop=FALSE]

st_var <- function(x) {
  (x - mean(x)) / sd(x)
}

st_var_values <- apply(pca_df_subset, 2, st_var)

column_variances <- apply(st_var_values, 2, var)
column_measn<- colMeans(st_var_values)

pca <- prcomp(st_var_values, center=TRUE, scale = FALSE)
autoplot(pca, data = data.frame(Group = condition), colour = "Group") + 
  ggtitle("PCA of RNA-seq Samples") +
  theme_minimal()

### hybrid approach ###
logCPM <- t(cpm(y, log=TRUE, prior.count = 2))
logCPM_t<- t(logCPM)

cv_values_cpm <- apply(logCPM, 2, cv_function)
cv_results_cpm <- data.frame(gene = colnames(logCPM), CV = cv_values_cpm)
top_cv_genes_cpm <- cv_results_cpm[order(-cv_results_cpm$CV), ][1:1000, ]  
pca_df_subset_cpm <- logCPM[, top_cv_genes_cpm$gene, drop=FALSE]  

st_var_values_cpm <- apply(pca_df_subset_cpm, 2, st_var)
st_var_values_cpm_t<- t(st_var_values_cpm)
column_variances <- apply(st_var_values_cpm, 2, var)
column_measn<- colMeans(st_var_values_cpm)

num_reads<- colSums(logCPM_t)

lm_fit <- lm(logCPM ~ num_reads)
residuals <- residuals(lm_fit)

# Now use the residuals for PCA
pca_input_adjusted <- residuals

# Perform PCA on the adjusted data
pca <- prcomp(pca_input_adjusted, scale = TRUE)

# Extract the PC1 scores and calculate the correlation
pc1_scores <- pca$x[, 1]
correlation_values <- cor(pc1_scores, num_reads, method = "pearson")

pca_data <- data.frame(
  PC1 = pc1_scores,
  PC2 = pca$x[, 2],  
  Group = condition_all,
  Correlation = round(correlation_values, 3) 
)

### PCA ###
pca_plot<-ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Correlation)) +
  geom_point(size = 3) + 
  geom_text(vjust=-2 , size = 4) + 
  ggtitle("PCA of RNA-seq Samples with Correlation (PC1)") +
  theme_bw()

ggsave(plot=pca_plot,"PCA_hybrid_5000_PC1_annot.png",dpi=300)


### check on PCA ###
explained_variance <- pca$sdev^2 / sum(pca$sdev^2)
cumulative_variance <- cumsum(explained_variance)
data.frame(PC = 1:length(explained_variance), 
           Explained_Variance = explained_variance, 
           Cumulative_Variance = cumulative_variance)

screeplot(pca, main = "Scree Plot", col = "blue", pch = 16)

library(GGally)  
pca_df <- data.frame(pca$x[, 1:5])  
ggpairs(pca_df)



pca_matrix <- pca$x[, 1:5]  # based on first 5PC PT_shPAEP2,2 and PT_shPAEP2_1 are the most variable. 
pheatmap(pca_matrix, cluster_cols = TRUE, 
         main = "Heatmap of First 5 PCs")



hc_samples <- hclust(dist(pca_matrix), method = "ward.D2")  
plot(hc_samples, main = "Hierarchical Clustering of First 5 PCs", 
     xlab = "Samples", sub = "", cex = 0.8)
     
## combat ##
batch <- factor(c(rep(1,5), rep(2,4), rep(3,5)))
logCPM_corrected <- ComBat(dat=logCPM, batch=batch, mod=NULL) #Useless because SCR,PAEP1,PAEP2 same batch

## removeBatchEffect ##

logCPM_corrected <- removeBatchEffect(logCPM_t, batch=batch)
logCPM_corrected <- ComBat(dat=logCPM_t, batch=batch, mod=NULL) 
logCPM_corrected_t<- t(logCPM_corrected)

### checks to understand PAEP2 variability ###
plotMDS(y,col=as.numeric(condition))
plotMDS(logCPM_t, col=as.numeric(condition))
logCPM_t<-t(logCPM)
st_var_values_cpm_t<-t(st_var_values_cpm)
log_CPM_paep2_t<-t(log_CPM_paep2)
boxplot(colSums(logCPM_t) ~ condition_3, main="PAEP Expression Across Conditions")
boxplot(colSums(st_var_values_cpm_t) ~ condition_3, main="PAEP Expression Across Conditions")

#
boxplot(colSums(st_var_values_cpm_t) ~ condition_3, main="PAEP Expression Across Conditions")
boxplot(colSums(logCPM_t) ~ condition_3, main="PAEP Expression Across Conditions")
boxplot(colSums(y$counts) ~ condition_3, main="Library Size Across Conditions")

#check
gene_to_remove <-rownames(st_var_values_cpm_t)
top_100<- gene_to_remove[1:100]
logCPM_t_filtered<- logCPM_t[!rownames(logCPM_t) %in% gene_to_remove,]
boxplot(logCPM_t[top_100, ] ~ condition_3, main="Expression of Removed Genes by Condition")

library(ggplot2)
library(reshape2) 
top_100_exp <- logCPM_t[gene_to_remove, ]

top_100_long <- melt(top_100_exp)
colnames(top_100_long) <- c("Gene", "Sample", "Expression")

# Add sample condition information
top_100_long$Condition <- condition_3[match(top_100_long$Sample, colnames(logCPM_t))]

# Boxplot of expression for top 100 genes by condition
ggplot(top_100_long, aes(x = Condition, y = Expression)) +
  geom_boxplot() +
  labs(title = "Expression of Top CV Genes by Condition",
       x = "Condition",
       y = "Expression Level") +
  theme_minimal()





logCPM_filtered<-t(logCPM_t_filtered)
cv_values_cpm <- apply(logCPM_filtered, 2, cv_function)
cv_results_cpm <- data.frame(gene = colnames(logCPM_filtered), CV = cv_values_cpm)
top_2500_cv_genes_cpm <- cv_results_cpm[order(-cv_results_cpm$CV), ][1:2500, ]  
pca_df_subset_cpm <- logCPM[, top_2500_cv_genes_cpm$gene, drop=FALSE]  

st_var_values_cpm <- apply(pca_df_subset_cpm, 2, st_var)
st_var_values_cpm_t<-t(st_var_values_cpm)
boxplot(colSums(st_var_values_cpm_t) ~ condition_3, main="PAEP Expression Across Conditions")
pca <- prcomp(st_var_values_cpm, scale = TRUE)
pca_plot<-autoplot(pca, data = data.frame(Group = condition_3), colour = "Group") + 
  ggtitle("PCA of RNA-seq Samples") +
  theme_bw()
ggsave(plot=pca_plot,"PCA_hybrid_5000.png",dpi=300)






barplot(logCPM_t, names.arg=colnames(logCPM_t), las=2, main="Library Sizes")

boxplot(logCPM_t[, condition_3 == "PT_shPAEP2"], main = "LogCPM Distribution for PAEP2")
boxplot(log_CPM_paep2_t, main = "LogCPM Distribution of PAEP2 Samples", col = "lightblue")
plotMDS(logCPM_t, col = as.numeric(factor(condition_3)))
boxplot(logCPM_t, main="Distribution of logCPM Values Across Samples")
plotMeanVar(y)
logCPM_PAEP2 <- logCPM[, condition_3 == "PT_shPAEP2", drop=FALSE]
var_logCPM <- apply(logCPM_PAEP2, 1, var)
boxplot(var_logCPM, var_voom, names=c("logCPM", "Voom"), main="Variance Comparison (PAEP2 Samples)")

### extract PC for GSEA ###
pca_loadings <- pca$rotation[, 1:2]  
loadings_df <- data.frame(
  gene = rownames(pca_loadings),
  PC1_loading = pca_loadings[, 1],
  PC2_loading = pca_loadings[, 2]
)
write.csv(loadings_df, "pca_loadings.csv", row.names = FALSE)

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
write.xlsx(ttt, "PT_shSCRvsPT_shPAEP_FDR.xlsx", rowNames=TRUE)
plotMeanVar(y, show.raw.vars=TRUE, show.tagwise.vars=TRUE, main="Mean-Variance Plot (logCPM)")

### check GRN of PAEP ###

grn_corr<- read_excel("top50_corr_tr_GRN_paep.xlsx", col_names = FALSE)
colnames(grn_corr) <- c("genes")
grn_chemor<-read_excel("top50_chemor_promet.xlsx", col_names = FALSE)
colnames(grn_chemor) <- c("genes")

grn_ensg_corr <- tt[rownames(tt) %in% grn_corr$genes,] %>%
  filter(FDR <= 0.1 | PValue < 0.05) %>% 
  arrange(desc(logFC))  

grn_ensg_chemor <- tt[rownames(tt) %in% grn_chemor$genes,] %>%
  filter(FDR <= 0.1 | PValue < 0.05) %>%  
  arrange(desc(logFC))  

write.xlsx(grn_ensg_corr, "DE_corr.xlsx", rowNames=TRUE)
write.xlsx(grn_ensg_chemor, "DE_chemor.xlsx", rowNames=TRUE)


# volcano plot #
volcano_plot <- function(zzz, title) {
  
  zzz$logFC <- as.numeric(as.character(zzz$logFC))
  zzz$threshold <- zzz$PValue < 0.05 & abs(zzz$logFC) > 1
  
  ggplot(zzz, aes(x = logFC, y = -log10(PValue), color = threshold)) + 
    geom_point(alpha = 0.8) + 
    scale_color_manual(values = c("black", "red")) + 
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white"),  
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "lightgrey", linewidth = 0.5), 
      panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.25)  
    ) + 
    ggtitle(title) + 
    xlab("Log Fold Change") + 
    ylab("-Log10 P-value") + 
    theme(legend.position = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.5) + 
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black", linewidth = 0.5)  
}

volcano_plot <- volcano_plot(zzz,title ='Differential Gene Expression')
ggsave(plot= volcano_plot, "volcano_plot_PAEP.png", dpi=300)



### heatmap ###

top_genes <- ttt #25 top DEGs
logCPM_t <- t(logCPM)
heat_data<- logCPM_t[rownames(logCPM_t) %in% rownames(top_genes), ]


color_palette<- colorRampPalette(c("blue", "white", "red"))(50)
heatmap_plot<- pheatmap(heat_data,
                        cluster_cols = TRUE,
                        scale="row",
                        annotation_col=data.frame(condition=condition, row.names = colnames(heat_data)),
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        color = color_palette
)
ggsave("heatmap_PAEP_corr.png", plot=heatmap_plot, width=12,height = 10,dpi=300)


### PAEP1 vs PAEP2 ### 
y_3 <- DGEList(counts = counts, group = condition_3)

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
design_3<-model.matrix(~ 0 + condition_3)
y_3<- estimateDisp(y,design=design_3)

fit <- glmQLFit(y_3, design_3)

PAEP1vsPAEP2 <- makeContrasts(condition_3PT_shPAEP1 - condition_3PT_shPAEP2,
                           levels = design_3)

qlf<-glmQLFTest(fit, contrast=PAEP1vsPAEP2)
x_3<-topTags(qlf, n=Inf)
xx_3<- x_3$table

xxx_3<-xx_3[xx_3$FDR<=0.1,]
zzz_3<-xx_3[xx_3$PValue<=0.05,]
xxx_3<-xxx_3 %>% arrange(desc(logFC))
write.xlsx(xxx_3, "PAEP1vsPAEP2_degs.xlsx", rowNames=TRUE)
plotMeanVar(y_3, show.raw.vars=TRUE, show.tagwise.vars=TRUE, main="Mean-Variance Plot (logCPM)")

### save for GSEA ###
write.csv(xxx_3, "PAEP1vsPAEP2_gsea.csv")

### heatmap PAEP1vsPAEP2 ###
logCPM <- t(cpm(y_3, log=TRUE, prior.count=2))
top_genes <- ttt_3[1:50,] #100 top DEGs
logCPM_t <- t(logCPM)
heat_data<- logCPM_t[rownames(logCPM_t) %in% rownames(top_genes), ]


color_palette<- colorRampPalette(c("blue", "white", "red"))(50)
heatmap_plot<- pheatmap(heat_data,
                        cluster_cols = TRUE,
                        scale="row",
                        annotation_col=data.frame(condition=condition_3, row.names = colnames(heat_data)),
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        color = color_palette
)

ggsave("heatmap_PAEP1vsPAEP2.png", plot=heatmap_plot, width=12,height = 15,dpi=300)


### PAEP1 vs SCR ###
PAEP1vsSCR <- makeContrasts(condition_3PT_shPAEP1 - condition_3PT_shSCR,
                              levels = design_3)

qlf<-glmQLFTest(fit, contrast=PAEP1vsSCR)
t_3<-topTags(qlf, n=Inf)
tt_3<- t_3$table

ttt_3<-tt_3[tt_3$FDR<=0.1,]
zzz_3<-tt_3[tt_3$PValue<=0.05,]
ttt_3<-ttt_3 %>% arrange(desc(logFC))
write.csv(ttt_3,"PAEP1vsSCR_gsea.csv")
write.xlsx(ttt_3, "PAEP1vsSCR_degs.xlsx", rowNames=TRUE)
### check gene_regulatory network PAEP1vsSCR ###
grn_corr<- read_excel("top50_corr_tr_GRN_paep.xlsx", col_names = FALSE)
colnames(grn_corr) <- c("genes")
grn_chemor<-read_excel("top50_chemor_promet.xlsx", col_names = FALSE)
colnames(grn_chemor) <- c("genes")

grn_ensg_corr <- ttt_3[rownames(ttt_3) %in% grn_corr$genes,] %>%
  filter(FDR <= 0.1 | PValue <= 0.1) %>% 
  arrange(desc(logFC))  

grn_ensg_chemor <- ttt_3[rownames(ttt_3) %in% grn_chemor$genes,] %>%
  filter(FDR <= 0.1 | PValue <= 0.1) %>% 
  arrange(desc(logFC))  

write.csv(grn_ensg_corr, "DE_corr_PAEP1vsSCR.csv")
write.csv(grn_ensg_chemor, "DE_chemor_PAEP1vsSCR.csv")

### heatmap PAEP1vsSCR ###
logCPM <- t(cpm(y_3, log=TRUE, prior.count=2))
top_genes <- ttt_3[1:50,] #50 top DEGs
logCPM_t <- t(logCPM)
heat_data<- logCPM_t[rownames(logCPM_t) %in% rownames(top_genes), ]


color_palette<- colorRampPalette(c("blue", "white", "red"))(50)
heatmap_plot<- pheatmap(heat_data,
                        cluster_cols = TRUE,
                        scale="row",
                        annotation_col=data.frame(condition=condition_3, row.names = colnames(heat_data)),
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        color = color_palette
)
ggsave("heatmap_PAEP1vsSCR.png", plot=heatmap_plot, width=12,height = 15,dpi=300)

### PAEP2vsSCR ###
PAEP2vsSCR <- makeContrasts(condition_3PT_shPAEP2 - condition_3PT_shSCR,
                            levels = design_3)

qlf<-glmQLFTest(fit, contrast=PAEP2vsSCR)
n_3<-topTags(qlf, n=Inf)
nn_3<- n_3$table

nnn_3<-nn_3[nn_3$FDR<=0.1,]
zzz_3<-tt_3[tt_3$PValue<=0.05,]
nnn_3<-nnn_3 %>% arrange(desc(logFC))
write.csv(nnn_3, "PAEP2vsSCR_gsea.csv")
write.xlsx(nnn_3, "PAEP2vsSCR_degs.xlsx", rowNames=TRUE)
### check gene_regulatory network PAEP2vsSCR ###
grn_corr<- read_excel("top50_corr_tr_GRN_paep.xlsx", col_names = FALSE)
colnames(grn_corr) <- c("genes")
grn_chemor<-read_excel("top50_chemor_promet.xlsx", col_names = FALSE)
colnames(grn_chemor) <- c("genes")

grn_ensg_corr <- nnn_3[rownames(nnn_3) %in% grn_corr$genes,] %>%
  filter(FDR <= 0.1 | PValue <= 0.05) %>% 
  arrange(desc(logFC))  

grn_ensg_chemor <- nnn_3[rownames(nnn_3) %in% grn_chemor$genes,] %>%
  filter(FDR <= 0.1 | PValue <= 0.05) %>%
  arrange(desc(logFC))  

write.csv(grn_ensg_corr, "DE_corr_PAEP2vsSCR.csv")
write.csv(grn_ensg_chemor, "DE_chemor_PAEP2vsSCR.csv")
write.xlsx(grn_ensg_corr, "DE_corr_PAEP2vsSCR.xlsx", rowNames=TRUE)
write.xlsx(grn_ensg_chemor, "DE_chemor_PAEP2vsSCR.xlsx",rowNames=TRUE)

### heatmap PAEP2vsSCR ###
logCPM <- t(cpm(y_3, log=TRUE, prior.count=2))
top_genes <- rbind(head(ttt_3,25), tail(ttt_3,25)) #50 top DEGs and 50 bottom DEGs
logCPM_t <- t(logCPM)
heat_data<- logCPM_t[rownames(logCPM_t) %in% rownames(top_genes), ]


color_palette<- colorRampPalette(c("blue", "white", "red"))(50)
heatmap_plot<- pheatmap(heat_data,
                        cluster_cols = TRUE,
                        scale="row",
                        annotation_col=data.frame(condition=condition_3, row.names = colnames(heat_data)),
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        color = color_palette
)
ggsave("heatmap_PAEP2vsSCR.png", plot=heatmap_plot, width=12,height = 15,dpi=300)


### compare DEGs in PAEP1 vs SCR with PAEP2 vs SCR ###
common_rows<- intersect(rownames(ttt_3),rownames(nnn_3))
ttt_subset<-ttt_3[common_rows,]
nnn_subset<-nnn_3[common_rows,]

colnames(ttt_subset)<- paste0(colnames(ttt_subset), "_1")
colnames(nnn_subset)<-paste0(colnames(nnn_subset),"_2")

df_new<- cbind(ttt_subset,nnn_subset)
write.csv(df_new, "DEG_intersection_paep12.csv")

### compare DEGS 
common<- intersect(rownames(df_new), rownames(ttt))

ttt_subset<-ttt_3[common,]
nnn_subset<-nnn_3[common,]

colnames(ttt_subset)<- paste0(colnames(ttt_subset), "_1")
colnames(nnn_subset)<-paste0(colnames(nnn_subset),"_2")

df_new<- cbind(ttt_subset,nnn_subset)

### Pearson correlation DEGs paep1vsscr & paep2vsscr ###
common_genes <- intersect(rownames(ttt_3), rownames(nnn_3))
logFC_paep1 <- ttt_3[common_genes, "logFC"]
logFC_paep2 <- nnn_3[common_genes, "logFC"]

plot(logFC_paep1, logFC_paep2, xlab="LogFC PAEP1 vs SCR", ylab="LogFC PAEP2 vs SCR", main="Correlation of LogFC")
abline(lm(logFC_paep2 ~ logFC_paep1), col="red")  # Regression line
cor(logFC_paep1, logFC_paep2)  

###  pca calculated above ###
pca
loadings_PC1 <- abs(pca$rotation[,1])
top_PC1_genes <- names(sort(loadings_PC1, decreasing=TRUE))  # Top 2500 genes

# Intersect with PAEP1 vs PAEP2 DEGs
common_genes <- intersect(top_PC1_genes, rownames(xxx_3))
print(length(common_genes))


### PCA deg paep1 vs paep2 
logCPM_t<- t(logCPM)
logCPM_subset <- logCPM_t[rownames(logCPM_t) %in% rownames(xxx_3), ]
pca_DEGs <- prcomp(t(logCPM_subset), scale=TRUE)
plot(pca_DEGs$x[,1:2], col=condition_3, pch=16, main="PCA on PAEP1 vs PAEP2 DEGs")

logCPM_filtered <- logCPM_t[!rownames(logCPM_t) %in% rownames(xxx_3), ]
pca_filtered <- prcomp(t(logCPM_filtered), scale=TRUE)
plot(pca_filtered$x[,1:2], col=condition_3, pch=16, main="PCA After Removing PAEP1 vs PAEP2 DEGs")


common_deg_logFC <- xxx_3[rownames(xxx_3) %in% common_genes, ]
hist(common_deg_logFC$logFC, main="Distribution of LogFC for Common DEGs", xlab="logFC")

expr_paep1 <- logCPM_t[, condition_3 == "PT_shPAEP1"]
expr_paep2 <- logCPM_t[, condition_3 == "PT_shPAEP2"]

# Select expression for common genes
expr_common_paep1 <- expr_paep1[common_genes, ]
expr_common_paep2 <- expr_paep2[common_genes, ]

# Calculate correlation
cor_expr <- cor(expr_common_paep1, expr_common_paep2)
print(cor_expr) 


### check intra variability of paep2 samples ###
expr_paep2 <- logCPM_t[, condition_3 == "PT_shPAEP2"]
dist_matrix <- dist(t(expr_paep2), method = "euclidean")
hclust_res <- hclust(dist_matrix, method = "ward.D2")
plot(hclust_res, labels = colnames(expr_paep2), main = "Hierarchical Clustering of PAEP2 Samples")


cor_paep2_samples <- cor(expr_paep2)
print(cor_paep2_samples)

plotMDS(expr_paep2)

### cv in paep2 ###
cv_values_paep2 <- apply(t(expr_paep2), 2, cv_function)
cv_results_paep2 <- data.frame(gene = colnames(t(expr_paep2)), CV = cv_values_paep2)
top_2500_cv_genes_paep2 <- cv_results_paep2[order(-cv_results_paep2$CV), ][1:2500, ]
pca_df_subset_paep2 <- logCPM[, !colnames(logCPM) %in% top_2500_cv_genes_paep2$gene, drop=FALSE]

st_var_values_paep2 <- apply(pca_df_subset_paep2, 2, st_var)

column_variances <- apply(st_var_values_paep2, 2,var)
column_measn<- colMeans(st_var_values_paep2)

pca <- prcomp(st_var_values_paep2, scale = TRUE)
pca_plot<-autoplot(pca, data = data.frame(Group = condition_3), colour = "Group") + 
  ggtitle("PCA of RNA-seq Samples") +
  theme_bw()
ggsave(plot=pca_plot,"PCA_hybrid_5000.png",dpi=300)

### cv in paep1 ###
expr_paep1 <- logCPM_t[, condition_3 == "PT_shPAEP1"]
cv_values_paep1 <- apply(t(expr_paep1), 2, cv_function)
cv_results_paep1 <- data.frame(gene = colnames(t(expr_paep1)), CV = cv_values_paep1)
top_2500_cv_genes_paep1 <- cv_results_paep1[order(-cv_results_paep1$CV), ][1:2500, ]



### cv in scr ###
expr_paep <- logCPM_t[, condition_3 == "PT_shSCR"]
cv_values_paep_scr <- apply(t(expr_paep), 2, cv_function)
cv_results_paep_scr <- data.frame(gene = colnames(t(expr_paep)), CV = cv_values_paep_scr)
top_2500_cv_genes_paep_scr <- cv_results_paep_scr[order(-cv_results_paep_scr$CV), ][1:2500, ]


### heat map intra expression paep2 ###
paep_diff<- c("PT_shPAEP2_1" ,"PT_shPAEP2_2")
expr_paep2 <- expr_paep2[,colnames(expr_paep2) %in% paep_diff]
cv_values_paep2 <- apply(t(expr_paep2), 2, cv_function)
tp_cv_genes <-names(sort(cv_values_paep2, decreasing =TRUE)[1:2500])
expr_paep2_tp_cv<- expr_paep2[tp_cv_genes,]
condition_paep2<- factor(c("PT_shPAEP2_1" ,"PT_shPAEP2_2"))
color_palette<- colorRampPalette(c("blue", "white", "red"))(50)
heatmap_plot<- pheatmap(expr_paep2_tp_cv,
                        cluster_cols = TRUE,
                        scale="row",
                        annotation_col=data.frame(condition=condition_paep2, row.names = colnames(expr_paep2_tp_cv)),
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        color = color_palette
)

ggsave(plot= heatmap_plot, "heatmap_intra_paep2.png", dpi=300)
### heatmap intra expression paep1 ###
cv_values_paep1 <- apply(t(expr_paep1), 2, cv_function)
tp_cv_genes <-names(sort(cv_values_paep1, decreasing =TRUE)[1:50])
expr_paep1_tp_cv<- expr_paep1[tp_cv_genes,]
condition_paep1<- factor(c("PT_shPAEP1_1" ,"PT_shPAEP1_2","PT_shPAEP1_3","PT_shPAEP1_4","PT_shPAEP1_5"))
color_palette<- colorRampPalette(c("blue", "white", "red"))(50)
heatmap_plot<- pheatmap(expr_paep1_tp_cv,
                        cluster_cols = TRUE,
                        scale="row",
                        annotation_col=data.frame(condition=condition_paep1, row.names = colnames(expr_paep1_tp_cv)),
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        color = color_palette
)
ggsave(plot= heatmap_plot, "heatmap_intra_paep1.png", dpi=300)
### heatmap intra expression scr ###
cv_values_paep_scr <- apply(t(expr_paep), 2, cv_function)
tp_cv_genes <-names(sort(cv_values_paep_scr, decreasing =TRUE)[1:50])
expr_paep_scr_tp_cv<- expr_paep[tp_cv_genes,]
condition_paep_scr<- factor(c("PT_shSCR_1" ,"PT_shSCR_2","PT_shSCR_3","PT_shSCR_4","PT_shSCR_5"))
color_palette<- colorRampPalette(c("blue", "white", "red"))(50)
heatmap_plot<- pheatmap(expr_paep_scr_tp_cv,
                        cluster_cols = TRUE,
                        scale="row",
                        annotation_col=data.frame(condition=condition_paep_scr, row.names = colnames(expr_paep_scr_tp_cv)),
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        color = color_palette
)
ggsave(plot= heatmap_plot, "heatmap_intra_scr.png", dpi=300)

### ###
subset_logCPM<- logCPM[ ,!colnames(logCPM) %in% tp_cv_genes, drop=FALSE]

st_var_values_no <- apply(subset_logCPM, 2, st_var)

column_variances <- apply(st_var_values_no, 2,var)
column_measn<- colMeans(st_var_values_no)

pca <- prcomp(st_var_values_no, scale = TRUE)
pca_plot<-autoplot(pca, data = data.frame(Group = condition_3), colour = "Group") + 
  ggtitle("PCA of RNA-seq Samples") +
  theme_bw()
ggsave(plot=pca_plot,"PCA_hybrid_5000.png",dpi=300)


### analysis redone without PAEP2_1 and PAEP2_2 ###

be<- null.merged.gene_counts.SummarizedExperiment
counts <- be@assays@data@listData[["salmon.merged.gene_counts"]]
counts_paep<- counts[, !grepl("PT_shPAEP2_1|PT_shPAEP2_2",colnames(counts), ignore.case= TRUE)]
sample<-be@colData@rownames
sample <- sample[!grepl("PT_shPAEP2_1|PT_shPAEP2_2", sample)]
cond <- sub("_[^_]+$", "", sample)

condition <- factor(ifelse(cond %in% c("PT_shPAEP1", "PT_shPAEP2"), "PT_shPAEP", cond),
                    levels = c("PT_shSCR", "PT_shPAEP"))
condition_3 <- factor(cond, levels= c("PT_shSCR","PT_shPAEP1", "PT_shPAEP2")) 

y <- DGEList(counts = counts_paep, group = condition) 

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


logCPM <- t(cpm(y, log=TRUE, prior.count=2))
logCPM_t<- t(logCPM)
cv_values_cpm <- apply(logCPM, 2, cv_function)
cv_results_cpm <- data.frame(gene = colnames(logCPM), CV = cv_values_cpm)
top_2500_cv_genes_cpm <- cv_results_cpm[order(-cv_results_cpm$CV), ][1:2500, ]  
pca_df_subset_cpm <- logCPM[, top_2500_cv_genes_cpm$gene, drop=FALSE]  

st_var_values_cpm <- apply(pca_df_subset_cpm, 2, st_var)

column_variances <- apply(st_var_values_cpm, 2, var)
column_measn<- colMeans(st_var_values_cpm)

pca <- prcomp(st_var_values_cpm, scale = TRUE)
pca_plot<-autoplot(pca, data = data.frame(Group = condition), colour = "Group") + 
  ggtitle("PCA of RNA-seq Samples") +
  theme_bw()
ggsave(plot=pca_plot,"PCA_hybrid_5000.png",dpi=300)

explained_variance <- pca$sdev^2 / sum(pca$sdev^2)
cumulative_variance <- cumsum(explained_variance)
data.frame(PC = 1:length(explained_variance), 
           Explained_Variance = explained_variance, 
           Cumulative_Variance = cumulative_variance)

screeplot(pca, main = "Scree Plot", col = "blue", pch = 16)

library(GGally)  
pca_df <- data.frame(pca$x[, 1:5])  # Extract first 5 PCs
ggpairs(pca_df)



pca_matrix <- pca$x[, 1:5]  # Extract first 5 PCs
pheatmap(pca_matrix, cluster_rows = TRUE, cluster_cols = TRUE, 
         main = "Heatmap of First 5 PCs")



hc_samples <- hclust(dist(pca_matrix), method = "ward.D2")  
plot(hc_samples, main = "Hierarchical Clustering of First 5 PCs", 
     xlab = "Samples", sub = "", cex = 0.8)


design<-model.matrix(~ 0 + condition_3)
PAEP1vsPAEP2 <- makeContrasts(condition_3PT_shPAEP1 - condition_3PT_shPAEP2,
                            levels = design)
y<- estimateDisp(y,design=design)
fit <- glmQLFit(y, design)
qlf<-glmQLFTest(fit, contrast=PAEP1vsPAEP2)
n_3<-topTags(qlf, n=Inf)
nn_3<- n_3$table

nnn_3<-nn_3[nn_3$FDR<=0.1,]
zzz_3<-tt_3[tt_3$PValue<=0.05,]
nnn_3<-nnn_3 %>% arrange(desc(logFC))




### heatmap 3x3 ###
library(readxl)
paep1vsscr<- read_excel("PAEP1vsSCR_degs.xlsx")
paep2vsscr<-read_excel("PAEP2vsSCR_degs.xlsx")
paep1vspaep2<-read_excel("PAEP1vsPAEP2_degs.xlsx")
paep1vsscr<-paep1vsscr[1:25,]
paep2vsscr<-paep2vsscr[1:25,]
paep1vspaep2<-paep1vspaep2[1:25,]
combined_top_genes <- unique(c(paep1vsscr$...1, 
                               paep2vsscr$...1, 
                              paep1vspaep2$...1))

logCPM <- t(cpm(y_3, log=TRUE, prior.count=2))
logCPM_t <- t(logCPM)  # Ensure you have the transposed logCPM matrix
heat_data <- logCPM_t[rownames(logCPM_t) %in% combined_top_genes, ]

# Create a color palette for the heatmap
color_palette <- colorRampPalette(c("blue", "white", "red"))(50)

# Create the heatmap
heatmap_plot <- pheatmap(
  heat_data,
  cluster_cols = TRUE,
  scale = "row",
  annotation_col = data.frame(condition = condition_3, row.names = colnames(heat_data)),
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = color_palette
)

# Save the plot as a PNG file
ggsave("heatmap_top_50_DEGs_across_conditions.png", plot = heatmap_plot, width = 12, height = 15, dpi = 300)

