
library(dplyr)

#################################################
# CNV data
#################################################
# copy number data
data_path = "/DATA/projects/j.bhin/Daniel_FGFR2/cell_lines/cnv/"
samp_label_cell_cnv = list.files(data_path)
samp_label_cell_cnv = samp_label_cell_cnv[grepl("txt", samp_label_cell_cnv)]

ind_FGFR1 = 67287:67290
ind_FGFR2 = 83538:83545
ind_FGFR3 = 32070:32071
ind_FGFR4 = 49645:49647

cnv_FGFR = matrix(NA, 4, length(samp_label_cell_cnv))
for(i in 1:length(samp_label_cell_cnv)){
  t = read.table(paste(data_path, samp_label_cell_cnv[i], sep=""), header = T, sep = "\t", stringsAsFactors = F)
  cnv_FGFR[1,i] = (2^(mean(t$RATIO_CORRECTED[ind_FGFR1]))) *2
  cnv_FGFR[2,i] = (2^(median(t$RATIO_CORRECTED[ind_FGFR2]))) *2
  cnv_FGFR[3,i] = (2^(median(t$RATIO_CORRECTED[ind_FGFR3]))) *2
  cnv_FGFR[4,i] = (2^(median(t$RATIO_CORRECTED[ind_FGFR4]))) *2
}
rownames(cnv_FGFR) = c("FGFR1_CNV", "FGFR2_CNV", "FGFR3_CNV", "FGFR4_CNV")
colnames(cnv_FGFR) = str_split_fixed(samp_label_cell_cnv, "_", 6)[,3]

save.image("~/FGFR/Daniel/R/Nature_figures/data/Human_celllines/Human_celllines_CNV.RData")

