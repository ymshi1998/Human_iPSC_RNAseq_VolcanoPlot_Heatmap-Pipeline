library(DESeq2)
library(writexl)

# Read raw count file
counts <- read.table("GSE184419_Read_counts.txt",
                     header = TRUE,
                     sep = "\t",
                     stringsAsFactors = FALSE)

# Set gene_id (symbols) as rownames
rownames(counts) <- counts$Geneid
counts <- counts[ , -1]

# Construct group information
phenotypes <- c("GG", "GG", "GA", "GA")
colData <- data.frame(phenotypes = factor(phenotypes))
rownames(colData) <- colnames(counts)

# Construct DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ phenotypes
)

# run DESeq2
dds <- DESeq(dds)

# Prepare normalized counts for Heatmap
norm_counts <- counts(dds, normalized=TRUE)

# Extract DESeq results of GA vs GG
res <- results(dds, contrast = c("phenotypes", "GA", "GG"))
res$Gene_Symbol <- rownames(res)

# Save the result
write_xlsx(as.data.frame(res), "iPSC_DEG_GA_vs_GG.xlsx")
write.csv(norm_counts, "normalized_counts.csv")
