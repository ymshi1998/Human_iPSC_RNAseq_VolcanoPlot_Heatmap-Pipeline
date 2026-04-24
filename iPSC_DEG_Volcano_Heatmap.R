#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

# ---------- User inputs ----------
de_stats <- "iPSC_DEG_GA_vs_GG.xlsx"

out_heat_png <- "Human_iPSC_Basal_GA_vs_GG_Heatmap.png"
out_vol_png  <- "Human_iPSC_Basal_GA_vs_GG_Volcano.png"

padj_cut <- 0.05
lfc_cut <- 1

# ----------  Volcano ----------
message("Reading DE stats file: ", de_stats)
de <- read_excel(de_stats)

de2 <- de %>%
  mutate(
    Gene_symbol = as.character(Gene_Symbol),
    log2FC = as.numeric(log2FoldChange),
    p = as.numeric(pvalue),
    padj = as.numeric(padj),
    neglog10p = -log10(padj),
    is_sig = (padj < padj_cut) & (abs(log2FC) >= lfc_cut)
  )  

label_df <- de2 %>%
  filter(is_sig) %>%
  arrange(padj) %>%
  head(30)

vol <- ggplot(de2, aes(x = log2FC, y = neglog10p)) +
  geom_point(aes(color = is_sig), alpha = 0.8, size = 1.8) +
  scale_color_manual(values = c("FALSE" = "grey20", "TRUE" = "red")) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "grey50") +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = Gene_symbol),
    size = 3,
    max.overlaps = Inf,
  ) +
  labs(
    title = "Human iPSC-EC basal GA vs GG",
    x = "log2FC (GA vs GG)",
    y = "-log10(padj)",
    color = "Significance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face="bold")
  )
  # xlim(-7, 7) +        # FIXED X-AXIS RANGE
  # ylim(-1, 5)           # FIXED Y-AXIS RANGE

ggsave(out_vol_png, vol, width = 10, height = 9, dpi = 300)

message("Done.")

# ----------  Heatmap ----------
norm <- read.csv("normalized_counts.csv", row.names = 1)

# top 50
deg_sig <- de2[de2$is_sig == TRUE & !is.na(de2$is_sig), ]
deg_sig_sorted <- deg_sig[order(abs(deg_sig$log2FC), decreasing = TRUE), ]
top50 <- head(deg_sig_sorted$Gene_symbol, 50)

# subset expression matrix
mat50 <- norm[top50, ]

# row-wise scaling
mat_scaled <- t(scale(t(mat50)))

png(out_heat_png, width = 6000, height = 3000, res = 300)
p <- pheatmap(
  t(mat_scaled),
  fontsize_row = 13,
  fontsize_col = 13,
  cellwidth = 23,
  cellheight = 50
)
dev.off()
message("Done.")
