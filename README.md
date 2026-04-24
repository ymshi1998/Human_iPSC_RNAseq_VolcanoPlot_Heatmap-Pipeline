# Human_iPSC_RNAseq_VolcanoPlot_Heatmap-Pipeline
Stanford Internship: Human_iPSC_RNAseq_VolcanoPlot_Heatmap-Pipeline

## Introduction
This is a R language pipeline that firstly generates a DEG (Differentially expressed gene) list, then draw a volcano plot and a heatmap graph according to the generated DEG list.

### File Introduction
For file `Raw_Count_Generating_DEG.R`, this file takes `GSE184419_Read_counts.txt` as an raw gene count input. The data is downloaded from Omnibus website: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE184419

The output of this file is: a DEG list named "iPSC_DEG_GA_vs_GG.xlsx", and a normalized counts file "normalized_counts.csv" for heatmap drawing.

For file `iPSC_DEG_Volcano_Heatmap.R`, this is a following-up file which picks up the output of the previous running result, and draws the corresponding volcano plot and the heatmap. The volcano plot marks up 30 most important genes, whereas the heatmap draws top 50 most upregulated / downregulated genes.

### Cutoff Usage
We used log2FC cutoff as 1 and padj cutoff as 0.05.

## Running Instruction
The author uses RStudio to run the program. Firstly run `Raw_Count_Generating_DEG.R`, then run `iPSC_DEG_Volcano_Heatmap.R`. The two files are supposed to be put in the same folder.

## Attention
If you found you are missing any packages while running the program as a whole, you can always firstly install those packages either from code or from R console.
