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
library(limma)
library(openxlsx)
library(readxl)

sample<- factor(c("PT_shSCR_1","PT_shSCR_2","PT_shSCR_3",
                                  "PT_shSCR_4","PT_shSCR_5","PT_shPAEP1_1",
                                  "PT_shPAEP1_2","PT_shPAEP1_3","PT_shPAEP1_4",
                                  "PT_shPAEP1_5","PT_shPAEP2_1","PT_shPAEP2_2",
                                  "PT_shPAEP2_3","PT_shPAEP2_4")) 

condition <- factor(ifelse(grepl("shPAEP1", sample), "PT_shPAEP1",
                           ifelse(grepl("shPAEP2", sample), "PT_shPAEP2",
                                  "PT_shSCR")),
                    levels  = c("PT_shSCR", "PT_shPAEP1", "PT_shPAEP2"))

condition_2 <- factor(ifelse(grepl("shPAEP1", sample), "PT_shPAEP",
                             ifelse(grepl("shPAEP2", sample), "PT_shPAEP",
                                    "PT_shSCR")),
                      levels = c("PT_shSCR", "PT_shPAEP"))

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
counts <- counts[!is.na(rownames(counts)), ]

y <- DGEList(counts = counts, group = condition)
y$counts

keep <- filterByExpr(y)
y<- y[keep, ,keep.lib.sizes=FALSE]
y<- calcNormFactors(y)
design<-model.matrix(~ 0 + condition_2)
y<- estimateDisp(y,design=design)

### PCA ### 
cv_function <- function(x) {
  if (abs(mean(x)) <= 1e-18) return(0)  
  return(sd(x) / mean(x))
}
st_var <- function(x) {
  (x - mean(x)) / sd(x)
}


logCPM <- t(cpm(y, log=TRUE, prior.count = 2))
logCPM_t<- t(logCPM)
num_reads<- colSums(y$counts)
cell_cycle_reads <- logCPM_t[rownames(logCPM_t) %in% cell_cycle, ]
mean_cell_cycle_reads <- colMeans(cell_cycle_reads, na.rm=TRUE)

### no correction ###
cv_values_cpm <- apply(logCPM, 2, cv_function)
cv_results_cpm <- data.frame(gene = colnames(logCPM), CV = cv_values_cpm)
top_cv_genes_cpm <- cv_results_cpm[order(-cv_results_cpm$CV), ][1:1000, ]  
pca_df_subset_cpm <- logCPM[, top_cv_genes_cpm$gene, drop=FALSE]  

st_var_values_cpm <- apply(pca_df_subset_cpm, 2, st_var)

pca<-prcomp(st_var_values_cpm, scale = TRUE)
pc1_values<-pca$x[,1]
cor_pc1_cellcycle <- cor(pc1_values, mean_cell_cycle_reads, method="pearson")
cor_pc1_totalreads <- cor(pc1_values, num_reads, method="pearson")
pca_df <- as.data.frame(pca$x)

pca_plot<-ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +  
  ggtitle("PCA of RNA-seq Samples : 5000") +
  theme_bw() 
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2),
           label = paste("PC1 vs Cell Cycle Correlation: ",
                         round(cor_pc1_cellcycle, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black") + 
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2) - 2,
           label = paste("PC1 vs Total Reads Correlation: ",
                         round(cor_pc1_totalreads, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black")

ggsave(filename = "/Users/ieo7295/Desktop/BC_sh/results/pca_plot/pca_5000_2.png", plot=pca_plot,width= 9, height= 6,dpi=300)


### regression num_reads
remove_num_reads_effect <- function(gene_expr) {
  model <- lm(gene_expr ~ num_reads)
  residuals(model)  
}

corrected_logCPM_num <- apply(logCPM, 2, remove_num_reads_effect)
cv_values_cpm <- apply(corrected_logCPM_num, 2, cv_function)
cv_results_cpm <- data.frame(gene = colnames(logCPM), CV = cv_values_cpm)
top_cv_genes_cpm <- cv_results_cpm[order(-cv_results_cpm$CV), ][1:1000, ]  
pca_df_subset_cpm <- logCPM[, top_cv_genes_cpm$gene, drop=FALSE]  

st_var_values_cpm <- apply(pca_df_subset_cpm, 2, st_var)

pca<-prcomp(st_var_values_cpm, scale. = TRUE)
pc1_values<-pca$x[,1]
cor_pc1_cellcycle <- cor(pc1_values, mean_cell_cycle_reads, method="pearson")
cor_pc1_totalreads <- cor(pc1_values, num_reads, method="pearson")
pca_df <- as.data.frame(pca$x)

pca_plot<-ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +  
  ggtitle("PCA of RNA-seq Samples : 4000") +
  theme_bw() +
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2),
           label = paste("PC1 vs Cell Cycle Correlation: ",
                         round(cor_pc1_cellcycle, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black") + 
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2) - 2,
           label = paste("PC1 vs Total Reads Correlation: ",
                         round(cor_pc1_totalreads, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black")

ggsave(filename = "/Users/ieo7295/Desktop/BC_sh/results/pca_plot/pca_5000num_reads.png", plot=pca_plot,width= 9, height= 6,dpi=300)

### regression cell_cycle ###
remove_cell_cycle_effect <- function(gene_expr) {
  model <- lm(gene_expr ~ mean_cell_cycle_reads)
  residuals(model)  
}

corrected_logCPM <- apply(logCPM, 2, remove_cell_cycle_effect)
cv_values_cpm <- apply(corrected_logCPM, 2, cv_function)
cv_results_cpm <- data.frame(gene = colnames(logCPM), CV = cv_values_cpm)
top_cv_genes_cpm <- cv_results_cpm[order(-cv_results_cpm$CV), ][1:1000, ]  
pca_df_subset_cpm <- logCPM[, top_cv_genes_cpm$gene, drop=FALSE]  

st_var_values_cpm <- apply(pca_df_subset_cpm, 2, st_var)

pca<-prcomp(st_var_values_cpm, scale. = TRUE)
pc1_values<-pca$x[,1]
cor_pc1_cellcycle <- cor(pc1_values, mean_cell_cycle_reads, method="pearson")
cor_pc1_totalreads <- cor(pc1_values, num_reads, method="pearson")
pca_df <- as.data.frame(pca$x)

pca_plot<-ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +  
  ggtitle("PCA of RNA-seq Samples : 5000") +
  theme_bw() +
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2),
           label = paste("PC1 vs Cell Cycle Correlation: ",
                         round(cor_pc1_cellcycle, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black") + 
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2) - 2,
           label = paste("PC1 vs Total Reads Correlation: ",
                         round(cor_pc1_totalreads, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black")

ggsave(filename = "/Users/ieo7295/Desktop/BC_sh/results/pca_plot/pca_5000_cell_cycle.png", plot=pca_plot,width= 9, height= 6,dpi=300)


### combat ###
num_reads_scaled <- scale(num_reads)
mod <- model.matrix(~ num_reads_scaled + mean_cell_cycle_reads)
batch <- factor(c(rep("PT_shSCR", 5), rep("PT_shPAEP1", 5), rep("PT_shPAEP2", 4)))
corrected_logCPM <- ComBat(dat = as.matrix(logCPM_t), batch = batch , mod = mod)

corrected_logCPM_t<-t(corrected_logCPM)
cv_values_cpm <- apply(corrected_logCPM_t, 2, cv_function)

cv_results_cpm <- data.frame(gene = colnames(logCPM), CV = cv_values_cpm)
top_cv_genes_cpm <- cv_results_cpm[order(-cv_results_cpm$CV), ][4000:5000, ]  
pca_df_subset_cpm <- logCPM[, top_cv_genes_cpm$gene, drop=FALSE]  

st_var_values_cpm <- apply(pca_df_subset_cpm, 2, st_var)


pca<-prcomp(st_var_values_cpm, scale. = TRUE)
pc1_values<-pca$x[,1]
cor_pc1_cellcycle <- cor(pc1_values, mean_cell_cycle_reads, method="pearson")
cor_pc1_totalreads <- cor(pc1_values, num_reads, method="pearson")
pca_df <- as.data.frame(pca$x)

pca_plot<-ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +  
  ggtitle("PCA of RNA-seq Samples : 5000") +
  theme_bw() +
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2),
           label = paste("PC1 vs Cell Cycle Correlation: ",
                         round(cor_pc1_cellcycle, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black") + 
  annotate("text", x = max(pca_df$PC1), y = max(pca_df$PC2) - 2,
           label = paste("PC1 vs Total Reads Correlation: ",
                         round(cor_pc1_totalreads, 2)),
           hjust = 1, vjust = 1, size = 5, color = "black")

ggsave(filename = "/Users/ieo7295/Desktop/BC_sh/results/pca_plot/pca_5000_combat.png", plot=pca_plot,width= 9, height= 6,dpi=300)


### DE edgeR ###

fit <- glmQLFit(y, design)

SCRvsPAEP <- makeContrasts(condition_2PT_shPAEP - condition_2PT_shSCR,
                           levels = design)

qlf<-glmQLFTest(fit, contrast=SCRvsPAEP)
t<-topTags(qlf, n=Inf)
tt<- t$table

ttt<-tt[tt$FDR<=0.1,]
ttt<- ttt %>% arrange(desc(logFC))
write.xlsx(ttt,"/Users/ieo7295/Desktop/BC_sh/results/pca_plot/Degs_paepvsscr.xlsx",rowNames = TRUE)

### PAEP1 vs PAEP2 ### 
y_3 <- DGEList(counts = counts, group = condition)


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
write.xlsx(xxx_3,"/Users/ieo7295/Desktop/BC_sh/results/pca_plot/Degs_paep1vspaep2.xlsx",rowNames = TRUE)
### PAEP1 vs SCR ###
PAEP1vsSCR <- makeContrasts(conditionPT_shPAEP1 - conditionPT_shSCR,
                            levels = design_3)

qlf<-glmQLFTest(fit, contrast=PAEP1vsSCR)
t_3<-topTags(qlf, n=Inf)
tt_3<- t_3$table

ttt_3<-tt_3[tt_3$FDR<=0.1,]
ttt_3<-ttt_3 %>% arrange(desc(logFC))
write.xlsx(ttt_3,"/Users/ieo7295/Desktop/BC_sh/results/pca_plot/Degs_paep1vscr.xlsx",rowNames = TRUE)
### PAEP2vsSCR ###
PAEP2vsSCR <- makeContrasts(conditionPT_shPAEP2 - conditionPT_shSCR,
                            levels = design_3)

qlf<-glmQLFTest(fit, contrast=PAEP2vsSCR)
n_3<-topTags(qlf, n=Inf)
nn_3<- n_3$table

nnn_3<-nn_3[nn_3$FDR<=0.1,]
nnn_3<-nnn_3 %>% arrange(desc(logFC))
write.xlsx(nnn_3,"/Users/ieo7295/Desktop/BC_sh/results/pca_plot/Degs_paep2vsscr.xlsx",rowNames = TRUE)

### heatmap top 50DEG of contrast sh1 vs scr, sh2 vs scr , sh1 vs sh2 ###
ttt_3$abs_logFC <- abs(ttt_3$logFC)
xxx_3$abs_logFC <- abs(xxx_3$logFC)
nnn_3$abs_logFC <- abs(nnn_3$logFC)


top_50_paep1scr <- ttt_3[order(-ttt_3$abs_logFC), ][1:50, ]
top_50_paep1paep2 <- xxx_3[order(-xxx_3$abs_logFC), ][1:50, ]
top_50_paep2scr <- nnn_3[order(-nnn_3$abs_logFC), ][1:50, ]

genes_paep1scr <- rownames(top_50_paep1scr)
genes_paep1paep2 <- rownames(top_50_paep1paep2)
genes_paep2scr <- rownames(top_50_paep2scr)


genes <-  unique(c(genes_paep1scr,genes_paep1paep2, genes_paep2scr))
logFC_matrix <- data.frame(Gene = genes)

logFC_matrix$sh1vsscr <- top_50_paep1scr$logFC[match(logFC_matrix$Gene, rownames(top_50_paep1scr))]
logFC_matrix$sh1vssh2 <- top_50_paep1paep2$logFC[match(logFC_matrix$Gene, rownames(top_50_paep1paep2))]
logFC_matrix$sh2vsscr <- top_50_paep2scr$logFC[match(logFC_matrix$Gene, rownames(top_50_paep2scr))]

# Replace NAs with 0 (for genes not in all contrasts)
logFC_matrix[is.na(logFC_matrix)] <- 0

# Convert to matrix for heatmap
rownames(logFC_matrix) <- logFC_matrix$Gene
logFC_matrix <- as.matrix(logFC_matrix[, -1]) 

# Heatmap plot
heatmap_contrast<- pheatmap(logFC_matrix, 
         cluster_cols = TRUE,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = TRUE,
         color = color_palette)



color_palette<- colorRampPalette(c("blue", "white", "red"))(50)
heatmap_plot <- pheatmap(expression_data,
                         cluster_cols = TRUE,
                         scale = "row",
                         annotation_col = data.frame(condition = condition, row.names = colnames(logCPM_t)),
                         show_rownames = TRUE,
                         show_colnames = TRUE,
                         color = color_palette)
ggsave("/Users/ieo7295/Desktop/BC_sh/results/pca_plot/heatmap_PAEP.png", plot=heatmap_contrast, width=12,height = 20,dpi=300)


### heatmap paep vs scr ###
top_genes<- ttt
logCPM_t <- t(logCPM)
heat_data<- logCPM_t[rownames(logCPM_t) %in% rownames(top_genes), ]



color_palette<- colorRampPalette(c("blue", "white", "red"))(50)
heatmap_plot <- pheatmap(heat_data,
                         cluster_cols = TRUE,
                         scale = "row",
                         annotation_col = data.frame(condition = condition, row.names = colnames(logCPM_t)),
                         show_rownames = TRUE,
                         show_colnames = TRUE,
                         color = color_palette)
ggsave("/Users/ieo7295/Desktop/BC_sh/results/pca_plot/heatmap_PAEP_vs_SCR.png", plot=heatmap_plot, width=12,height = 10,dpi=300)


### heatmap sh1 vs scr ###
top_genes<-nnn_3[1:50,]
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
ggsave("/Users/ieo7295/Desktop/BC_sh/results/pca_plot/heatmap_PAEP1_vs_scr.png", plot=heatmap_plot, width=12,height = 10,dpi=300)



### GRN all contrasts###
grn_corr<- read_excel("top50_corr_tr_GRN_paep.xlsx", col_names = FALSE)
colnames(grn_corr) <- c("genes")
grn_chemor<-read_excel("top50_chemor_promet.xlsx", col_names = FALSE)
colnames(grn_chemor) <- c("genes")

grn_ensg_corr <- tt_3[rownames(tt_3) %in% grn_corr$genes,] %>%
  filter(FDR <= 0.1 | PValue <= 0.1) %>% 
  arrange(desc(logFC))  

grn_ensg_chemor <- tt_3[rownames(tt_3) %in% grn_chemor$genes,] %>%
  filter(FDR <= 0.1 | PValue <= 0.1) %>% 
  arrange(desc(logFC)) 

write.xlsx(grn_ensg_corr ,"/Users/ieo7295/Desktop/BC_sh/results/pca_plot/corr_paep1_vs_scr.xlsx", rowNames=TRUE)
write.xlsx(grn_ensg_chemor,"/Users/ieo7295/Desktop/BC_sh/results/pca_plot/chemor_paep1_vs_scr.xlsx", rowNames=TRUE)

### volcano_plot ####
volcano_plot <- function(xx_3, title) {
  
  xx_3$logFC <- as.numeric(as.character(xx_3$logFC))
  xx_3$threshold <- xx_3$FDR < 0.1 & abs(xx_3$logFC) > 1
  
  ggplot(xx_3, aes(x = logFC, y = -log10(PValue), color = threshold)) + 
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

volcano_plot <- volcano_plot(xx_3,title ='Differential Gene Expression')
ggsave(plot= volcano_plot, "volcano_plot_PAEP1_vs_paep2.png", dpi=300)

### Deseq2 ###
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize

dds_filtered <- dds[keep,]
dds_final <- DESeq(dds_filtered)
rld <- rlog(dds_final, blind = TRUE)
combat_data <- ComBat(dat = assay(rld), batch = batch, mod = model)
assay(rld) <- combat_data
colData(dds)
pca_plot <- plotPCA(rld, intgroup = "condition")
ggsave(plot=pca_plot,"/Users/ieo7295/Desktop/BC_sh/results/pca_plot/deseq_pca_combat.png",dpi=300)

colData(dds)$num_reads_scaled <- num_reads_scaled
colData(dds)$mean_cell_cycle_reads <- mean_cell_cycle_reads
num_reads_scaled <- scale(num_reads)
model <- model.matrix(~ num_reads_scaled + mean_cell_cycle_reads)
batch <- factor(c(rep("PT_shSCR", 5), rep("PT_shPAEP1", 5), rep("PT_shPAEP2", 4)))
corrected_logCPM <- ComBat(dat = as.matrix(logCPM_t), batch = batch , mod = mod)



###heatmap deseq2 ###
library("RColorBrewer")
sampleDists <- dist(logCPM)

sampleDistMatrix <- as.matrix(sampleDists)

rownames(sampleDistMatrix) <- paste(rownames(logCPM))

colnames(sampleDistMatrix) <- paste(rownames(logCPM))

colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(50) 

clustering_plot <- pheatmap(sampleDistMatrix,
                            clustering_distance_rows = sampleDists,
                            clustering_distance_cols = sampleDists,
                            col = colors,
                            fontsize_col = 8,
                            fontsize_row = 8)

