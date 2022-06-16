library(dplyr)

FGFR2_SV_PURPLE_ALL = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/HMF/HMF_data_purple.RData")[["FGFR2_PURPLE_ALL"]]

# Sample information
target_samp = unique(FGFR2_SV_PURPLE_ALL %>% pull(sample_name))

# sample information
samp_info_meta = read.table("~/FGFR/Daniel/R/Nature_figures/data/HMF/metadata.tsv", sep = "\t", stringsAsFactors = F, header = T, comment.char = "")
samp_info_target = samp_info_meta[match(target_samp, samp_info_meta$setName),]

purple_path = "/DATA/projects/j.bhin/Daniel_FGFR2/HMF/somatics/"

FGFR2_CNV = matrix(NA, length(target_samp), 2)
for(i in 1:length(target_samp)){
  t_path = paste(purple_path, target_samp[i], sep = "")
  t_files = list.files(t_path)
  # cnv
  t_df = read.table(paste(t_path, t_files[grepl("cnv.gene", t_files)], sep = "/"), header = T, stringsAsFactors = F, sep = "\t")
  FGFR2_CNV[i,1] = t_df %>% filter(gene == "FGFR2") %>% pull(maxCopyNumber)
  # ploidy
  t_df = read.table(paste(t_path, t_files[grepl("purple.purity.tsv", t_files)], sep = "/"), header = T, stringsAsFactors = F, sep = "\t")
  FGFR2_CNV[i,2] = FGFR2_CNV[i,1]/t_df$ploidy
}
rownames(FGFR2_CNV) = samp_info_target$sampleId
colnames(FGFR2_CNV) = c("FGFR2_CN", "FGFR2_CN_ploidy")

saveRDS(FGFR2_CNV, file = "~/FGFR/Daniel/R/Nature_figures/data/HMF/HMF_data_cnv.RDS")
