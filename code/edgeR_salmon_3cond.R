library("readr")
library(edgeR)
library(tidyverse)
library("readxl")
library(ggplot2)

be<- null.merged.gene_counts.SummarizedExperiment
gene_counts<- be@assays@data@listData[["salmon.merged.gene_counts"]]
sample<-be@colData@rownames
condition_3 <- sub("_[^_]+$", "", sample)
condition_3 <- factor(condition_3, levels= c("PT_shSCR", "PT_shPAEP1", "PT_shPAEP2"))
design<-model.matrix(~ 0 + condition_3)


y<- DGEList(gene_counts,
            group=condition_3)

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

### estimate dispersion and fit model ###
y<- estimateDisp(y,design=design)
fit <- glmQLFit(y, design)

SCRvsPAEP1 <- makeContrasts(condition_3PT_shPAEP1 - condition_3PT_shSCR, levels = design)
SCRvsPAEP2<- makeContrasts(condition_3PT_shPAEP2 - condition_3PT_shSCR, levels = design)

qlf1<-glmQLFTest(fit, contrast=SCRvsPAEP1)
qlf2<-glmQLFTest(fit, contrast=SCRvsPAEP2)

tt1<-topTags(qlf1, n=Inf)
tt2<-topTags(qlf2, n=Inf)

tt1<- tt1$table
tt2<- tt2$table

ttt1<- tt1[tt1$FDR <= 0.1,] 
ttt2<-tt2[tt2$FDR <= 0.1,] 

ttt1<- ttt1 %>% arrange(desc(logFC))
ttt2<- ttt2 %>% arrange(desc(logFC))

write.csv(tt1, "PT_shSCRvsPT_shPAEP1.csv", row.names=TRUE)
write.csv(tt2, "PT_shSCRvsPT_shPAEP2.csv", row.names=TRUE)

### check GRN of PAEP ###
grn_corr<- read_excel("top50_corr_tr_GRN_paep.xlsx", col_names = FALSE)
colnames(grn_corr) <- c("genes")
grn_chemor<-read_excel("top50_chemor_promet.xlsx", col_names = FALSE)
colnames(grn_chemor) <- c("genes")

grn_ensg_corr_1 <- ttt1[rownames(ttt1) %in% grn_corr$genes,] %>%
  filter(FDR < 0.1 | PValue < 0.05) %>% 
  arrange(desc(logFC))  

grn_ensg_chemor_1 <- ttt1[rownames(ttt1) %in% grn_chemor$genes,] %>%
  filter(FDR < 0.1 | PValue < 0.05) %>%  
  arrange(desc(logFC))  

grn_ensg_corr_2 <- ttt2[rownames(ttt2) %in% grn_corr$genes,] %>%
  filter(FDR < 0.1 | PValue < 0.05) %>% 
  arrange(desc(logFC))  

grn_ensg_chemor_2 <- ttt2[rownames(ttt2) %in% grn_chemor$genes,] %>%
  filter(FDR < 0.1 | PValue < 0.05) %>%  
  arrange(desc(logFC))  

### plotting ### 


# volcano plot #
volcano_plot <- function(tt2, title) {
  
  tt2$logFC <- as.numeric(as.character(tt2$logFC))
  tt2$threshold <- tt2$FDR < 0.1 & abs(tt2$logFC) > 1
  ggplot(tt2, aes(x = logFC, y = -log10(PValue), color = threshold)) + 
    geom_point(alpha = 0.8) + 
    scale_color_manual(values = c("black", "red")) + 
    theme(
      legend.position = "none",
      panel.background = element_rect(fill = "white"),  # Set plot panel background to white
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "lightgrey", linewidth = 0.5),  # Keep major grid lines
      panel.grid.minor = element_line(color = "lightgrey", linewidth = 0.25)  # Keep minor grid lines
      # Set the entire plot background to white   # Remove minor grid lines
    ) + 
    ggtitle(title) + 
    xlab("Log Fold Change") + 
    ylab("-Log10 P-value") + 
    theme(legend.position = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black", linewidth = 0.5) + 
    geom_vline(xintercept = c(-1, 1), linetype = "dotted", color = "black", linewidth = 0.5)  
}

volcano_plot <- volcano_plot(tt2,title ='Differential Gene Expression')
ggsave(plot= volcano_plot, "volcano_plot_PAEP2.png", dpi=300)

# PCA # 

pca_df <- t(y$counts)

batch<- c(rep(1,5), rep(2,4), rep(3, 5))

cv_function <- function(x) {
  if (mean(x) == 0) return(0)  
  return(sd(x) / mean(x))
}
cv_values <- apply(pca_df, 2, cv_function)
cv_results <- data.frame(gene = colnames(pca_df), CV = cv_values)
top_cv_genes <- cv_results[order(-cv_results$CV), ]
top_2500_cv_genes <- top_cv_genes[1:2500, ]
pca_df_subset <- pca_df[, top_2500_cv_genes$gene, drop=FALSE]

pca_df_corrected <- ComBat(dat=t(pca_df_subset), batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
pca_df_corrected <- t(pca_df_corrected)

st_var <- function(x) {
  (x - mean(x)) / sd(x)
}
st_var_values <- apply(pca_df_corrected, 2, st_var)

column_variances <- apply(st_var_values, 2, var)
column_measn<- colMeans(st_var_values)

pca_res <- prcomp(st_var_values, center = TRUE, scale. = FALSE)

pca_df <- data.frame(
  sample = rownames(st_var_values),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  condition = condition_3
)

pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) + theme_light() +
  ggtitle("PCA of Gene Expression Data") +
  xlab("PC1") + ylab("PC2") +
  theme(legend.position = "right")

ggsave("PCA_PAEP_combat_3.png", plot = pca_plot, dpi = 300)

### heatmap ##

top_genes_1 <- ttt1
top_genes_2 <- ttt2 
topgenes_2 <-rbind(head(ttt2,25), tail(ttt2,25))
norm_counts <- cpm(y, log=TRUE)
heat_data<- norm_counts[rownames(norm_counts) %in% rownames(topgenes_2), ]


color_palette<- colorRampPalette(c("blue", "white", "red"))(50)
heatmap_plot<- pheatmap(heat_data,
                        cluster_cols = TRUE,
                        scale="row",
                        annotation_col=data.frame(condition=condition_3, row.names = colnames(heat_data)),
                        show_rownames = TRUE,
                        show_colnames = TRUE,
                        color = color_palette
)
ggsave("heatmap_PAEP_1.png", plot=heatmap_plot, width=12,height = 10,dpi=300)
ggsave("heatmap_PAEP_2.png", plot=heatmap_plot, width=12,height = 10,dpi=300)
