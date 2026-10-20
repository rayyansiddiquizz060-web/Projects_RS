# Differential expression analysis: Flight vs Ground Control
# OSD-576 (RR-23), mouse tibialis anterior
# Input: filtered counts from notebooks/02_qc.ipynb

library(DESeq2)

# ---- Load data ----
dataframe  <- read.csv("data/processed/dataframe_filtered.csv", row.names = 1, check.names = FALSE)
samples <- read.csv("data/samples.csv")

dim(dataframe)                  # should be 16115 genes x 18 samples
table(samples$condition)     # should be 9 Flight, 9 Ground


# ---- Check sample order ----
# DESeq2 matches samples to groups by position, so the order must match exactly
stopifnot(all(colnames(dataframe) == samples$sample))

samples$condition <- factor(samples$condition, levels = c("Ground", "Flight"))

# ---- Set up DESeq2 ----
# Ground is the reference level, so positive log2FC = higher in Flight
rownames(samples) <- samples$sample


dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(dataframe),
  colData   = samples,
  design    = ~ condition
)
dds


# ---- Run DESeq2 ----
dds <- DESeq(dds)

saveRDS(dds, "data/processed/dds.rds")

res <- results(dds, contrast = c("condition", "Flight", "Ground"), alpha = 0.05)
summary(res)

# ---- Shrink fold changes ----
# Pulls noisy fold changes from low-count genes toward zero (used for plots and ranking)
resultsNames(dds)
res_shrunk <- lfcShrink(dds, coef = "condition_Flight_vs_Ground", type = "apeglm")

# ---- Save results ----
dir.create("results", showWarnings = FALSE)

out <- as.data.frame(res)
out$log2FC_shrunk <- res_shrunk$log2FoldChange
out$gene_id <- rownames(out)
out <- out[order(out$padj), ]
out <- out[, c("gene_id", setdiff(names(out), "gene_id"))] 
write.csv(out, "results/deseq2_flight_vs_ground.csv", row.names = FALSE)


# Variance-stabilized counts for the heatmap in Step 6
vsd <- vst(dds, blind = FALSE)
write.csv(assay(vsd), "results/vst_counts.csv")

# Record R and package versions for reproducibility
writeLines(capture.output(sessionInfo()), "results/sessionInfo_deseq2.txt")


# ---- Quick look at results ----
head(out[, c("gene_id", "baseMean", "log2FoldChange", "padj")], 10)
sum(out$padj < 0.05, na.rm = TRUE)