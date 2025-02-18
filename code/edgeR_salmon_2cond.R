library("readr")
library(edgeR)
library(tidyverse)
library("readxl")
library(ggfortify)
library(pheatmap)
library(openxlsx)
library(sva)

be<- null.merged.gene_counts.SummarizedExperiment
gene_counts<- be@assays@data@listData[["salmon.merged.gene_counts"]]
sample<-be@colData@rownames
condition <- sub("_[^_]+$", "", sample)
condition <- factor(ifelse(condition %in% c("PT_shPAEP1", "PT_shPAEP2"), "PT_shPAEP", condition),
                    levels = c("PT_shSCR", "PT_shPAEP"))
condition_3 <- factor(condition, levels= c("PT_shSCR", "PT_shPAEP1", "PT_shPAEP2"))
design<-model.matrix(~ 0 + condition)


m<- DGEList(gene_counts,
            group=condition)

gene_annot<- data.frame(
  ensg = be@elementMetadata@listData[["gene_id"]], 
  genename = be@elementMetadata@listData[["gene_name"]]
)

gene_annot_top_genes <- gene_annot[gene_annot$ensg %in% rownames(m$counts),]
gene_names_matched <- gene_annot_top_genes$genename[match(rownames(m$counts),
                                                          gene_annot_top_genes$ensg)]
rownames(m$counts) <- gene_names_matched

keep <- filterByExpr(m)
m<- m[keep, ,keep.lib.sizes=FALSE]
m<- calcNormFactors(m)

plotMDS(m, method="bcv", col=as.numeric(m$samples$group))
legend("bottomleft", as.character(unique(m$samples$group)), pch=20)
### estimate dispersion and fit model ###
m<- estimateDisp(m,design=design)
fit <- glmQLFit(m, design)

SCRvsPAEP <- makeContrasts(conditionPT_shPAEP - conditionPT_shSCR,
                           levels = design)

qlf<-glmQLFTest(fit, contrast=SCRvsPAEP)
tt<-topTags(qlf, n=Inf)
zzz<- tt$table

ttt<- zzz[zzz$FDR <= 0.1,] 
xxx<- zzz[zzz$PValue <= 0.05,]
xxx<- xxx %>% arrange(desc(logFC))
ttt<- ttt %>% arrange(desc(logFC))
zzz<- zzz %>% arrange(desc(logFC))
write.xlsx(zzz, "PT_shSCRvsPT_shPAEP_all.xlsx", rowNames=TRUE)
write.xlsx(ttt, "PT_shSCRvsPT_shPAEP_FDR.xlsx", rowNames=TRUE)
write.xlsx(xxx, "PT_shSCRvsPT_shPAEP_PValue.xlsx", rowNames=TRUE)

### check GRN of PAEP ###

grn_corr<- read_excel("top50_corr_tr_GRN_paep.xlsx", col_names = FALSE)
colnames(grn_corr) <- c("genes")
grn_chemor<-read_excel("top50_chemor_promet.xlsx", col_names = FALSE)
colnames(grn_chemor) <- c("genes")

grn_ensg_corr <- tt$table[rownames(tt$table) %in% grn_corr$genes,] %>%
  filter(FDR < 0.1 | PValue < 0.05) %>% 
  arrange(desc(logFC))  

grn_ensg_chemor <- tt$table[rownames(tt$table) %in% grn_chemor$genes,] %>%
  filter(FDR < 0.1 | PValue < 0.05) %>%  
  arrange(desc(logFC))  

write.xlsx(grn_ensg_corr, "DE_corr.xlsx", rowNames=TRUE)
write.xlsx(grn_ensg_chemor, "DE_chemor.xlsx", rowNames=TRUE)

### plotting ### 


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


### PCA ### 
normalized_counts <- m$counts / m$samples$norm.factors
pca_df <- t(normalized_counts)
write.csv(s,"raw_counts.csv", row.names = TRUE)

batch<-c(rep(1,9),rep(2,5))


cv_function <- function(x) {
  if (mean(x) == 0) return(0)  
  return(sd(x) / mean(x))
}
cv_values <- apply(pca_df, 2, cv_function)
cv_results <- data.frame(gene = colnames(pca_df), CV = cv_values)
top_cv_genes <- cv_results[order(-cv_results$CV), ]
top_2500_cv_genes <- top_cv_genes[1:2500, ]
pca_df_subset <- pca_df[, top_2500_cv_genes$gene, drop=FALSE]

num_reads <- rowSums(pca_df) 


num_reads_log<- log(num_reads)
mod <- model.matrix(~ num_reads) 
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

ggsave("PCA_PAEP_combatbatch2_num_reads.png", plot = pca_plot, dpi = 300)


### extract for GSEA ### 
pca_loadings <- pca_res$rotation[, 1:2]  
loadings_df <- data.frame(
  gene = rownames(pca_loadings),
  PC1_loading = pca_loadings[, 1],
  PC2_loading = pca_loadings[, 2]
)
write.csv(loadings_df, "pca_loadings.csv", row.names = FALSE)

loadings_df$abs_PC1 <- abs(loadings_df$PC1_loading)
sorted_loadings <- loadings_df[order(-loadings_df$abs_PC1), ]
write.csv(sorted_loadings, "ranked_genes_for_gsea.csv", row.names = FALSE)

condition_data <- data.frame(sample = rownames(st_var_values), condition = condition_3)
write.csv(condition_data, "sample_conditions.csv", row.names = FALSE)


### Heatmap ### 
xxx_topgenes<- rbind(head(xxx,25), tail(xxx,25))
top_genes <- ttt #25 top DEGs
norm_counts <- cpm(y, log=TRUE)
heat_data<- norm_counts[rownames(norm_counts) %in% rownames(top_genes), ]


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


### GO ###




