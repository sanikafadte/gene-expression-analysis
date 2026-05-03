# Gene Expression Analysis (Tumor vs Normal)

# Load libraries
library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# Load dataset
gset <- getGEO("GSE44076", GSEMatrix = TRUE)

expr <- exprs(gset[[1]])
meta <- pData(gset[[1]])

# Filter Tumor vs Normal
keep <- meta$`sample type:ch1` %in% c("Tumor", "Normal")
expr2 <- expr[, keep]
meta2 <- meta[keep, ]

meta2$condition <- factor(meta2$`sample type:ch1`)

# Differential expression
design <- model.matrix(~ meta2$condition)
fit <- lmFit(expr2, design)
fit <- eBayes(fit)

results <- topTable(fit, coef = 2, number = Inf)

# Save results
write.csv(results, "DEG_results.csv")

# Volcano plot
results$direction <- ifelse(
  results$adj.P.Val < 0.05 & results$logFC > 1, "Up",
  ifelse(results$adj.P.Val < 0.05 & results$logFC < -1, "Down", "NS")
)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = direction), alpha = 0.6) +
  theme_minimal()

ggsave("volcano_plot_clean.png", width = 8, height = 6)

# Heatmap
top20 <- results[order(results$adj.P.Val), ][1:20, ]
mat <- expr2[rownames(top20), ]
mat <- t(scale(t(mat)))

pheatmap(mat,
         show_rownames = TRUE,
         show_colnames = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         filename = "heatmap_clean.png"
)

# GO enrichment
sig_genes <- results[results$adj.P.Val < 0.05, ]
genes <- sig_genes$GeneSymbol
genes <- genes[genes != "" & !is.na(genes)]

gene_ids <- bitr(genes,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db
)

ego <- enrichGO(
  gene = gene_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP"
)

png("GO_enrichment.png", width = 800, height = 600)
dotplot(ego, showCategory = 10)
dev.off()