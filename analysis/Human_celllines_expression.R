source("~/FGFR/Daniel/R/Nature_figures/Nature_figures/analysis/functions_FGFR_analysis.R")
library(stringr)
library(DEXSeq)


#################################################
# Cell line
#################################################
# isoform data - RSEM
data_path = "~/HPC/harris/Daniel_FGFR2/cell_lines/rsem/"
samp_label_cell = list.files(data_path)
samp_label_cell = str_split_fixed(samp_label_cell[grep("genes", samp_label_cell)], "_", 5)[,3]
FGFR2_iso_cell = rsem_iso_FGFR2(data_path, samp_label_cell)

# gene data - RSEM
FGFR2_gene_cell = rsem_gene_target("ENSG00000066468", data_path, samp_label_cell)
FGFR1_gene_cell = rsem_gene_target("ENSG00000077782", data_path, samp_label_cell)
FGFR3_gene_cell = rsem_gene_target("ENSG00000068078", data_path, samp_label_cell)

# exon data - Dexseq
data_path = "~/HPC/harris/Daniel_FGFR2/cell_lines/dexseq/"
FGFR2_exon_cell = dexseq_exon_FGFR2(data_path, samp_label_cell)
FGFR2_exon_cell$len = FGFR2_exon_cell$merge_len
FGFR2_exon_cell$merge_len = FGFR2_exon_cell$merge / FGFR2_exon_cell$len[,3]
rownames(FGFR2_exon_cell$merge_len) = rownames(FGFR2_exon_cell$merge_per)
rownames(FGFR2_exon_cell$RPKM) = rownames(FGFR2_exon_cell$merge_per)

# fusion data
# no fusion for FGFR1 or FGFR3
data_path = "~/HPC/harris/Daniel_FGFR2/cell_lines/star_fusion/"
FGFR2_fus_cell = fusion_target("FGFR2", data_path, samp_label_cell)
FGFR2_fus_cell = FGFR2_fus_cell[startsWith(rownames(FGFR2_fus_cell), "FGFR2"),]

# annotation
ind_sort = order(FGFR2_gene_cell$rsem_fpkm)
FGFR_gene_cell = rbind(log2(FGFR2_gene_cell$rsem_fpkm[ind_sort]+1),
                       log2(FGFR1_gene_cell$rsem_fpkm[ind_sort]+1),
                       log2(FGFR3_gene_cell$rsem_fpkm[ind_sort]+1))
rownames(FGFR_gene_cell) = c("FGFR2_Exp", "FGFR1_Exp", "FGFR3_Exp")

save.image("~/FGFR/Daniel/R/Nature_figures/data/Human_celllines/Human_celllines_expression.RData")