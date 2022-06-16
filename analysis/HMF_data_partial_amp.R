library(GenomicRanges)
library(dplyr)

##############################################################
# data processing
##############################################################
# HMF SV data
load("~/FGFR/Daniel/R/Nature_figures/data/HMF/HMF_data_purple.RData")

# loading hg19 genomic information
hg19_info = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/HMF/HMF_data_purple.RData")[["hg19_info"]]

# FGFR2 genomic range (uc021pzy)
ind = data.frame(hg19_info$trx) %>% filter(grepl("uc021pzy", tx_name)) %>% pull(tx_id)
t_exon = hg19_info$trx_exon[[ind]]
t_intron = rev(hg19_info$trx_intron[[ind]])
df1 = data.frame(data.frame(t_exon)[, c(1,2,3,5)], id = paste("Exon", 1:length(t_exon)))
df2 = rbind(data.frame(data.frame(t_intron)[, c(1,2,3,5)], id = paste("Intron", 1:length(t_intron))), NA)
FGFR2_info = (gdata::interleave(df1, df2))
FGFR2_info = GRanges(FGFR2_info[-nrow(FGFR2_info),])


purple_path = "/DATA/projects/j.bhin/Daniel_FGFR2/HMF/somatics/"
all_samples = list.files(purple_path)

FGFR2_CNV_brkpt = c()
for(i in 1:length(all_samples)){
  setwd(paste(purple_path, all_samples[i], sep = ""))
  t_files = list.files()
  t_cnv = read.table(t_files[grepl("purple.cnv.somatic.tsv", t_files)], sep = "\t", stringsAsFactors = F, header = T)
  
  t_start = t_cnv[,c(1,2,2)] %>% mutate(chromosome = paste("chr", chromosome, sep = ""))
  colnames(t_start) = colnames(t_start) = c("chromosome", "start", "end")
  t_start = GRanges(t_start)
  t_start_id = rep(NA, length(t_start))
  t_start_id[queryHits(findOverlaps(t_start, FGFR2_info))] = FGFR2_info$id[subjectHits(findOverlaps(t_start, FGFR2_info))]
  
  t_end = t_cnv[,c(1,3,3)] %>% mutate(chromosome = paste("chr", chromosome, sep = ""))
  colnames(t_end) = colnames(t_end) = c("chromosome", "start", "end")
  t_end = GRanges(t_end)
  t_end_id = rep(NA, length(t_end))
  t_end_id[queryHits(findOverlaps(t_end, FGFR2_info))] = FGFR2_info$id[subjectHits(findOverlaps(t_end, FGFR2_info))]
  
  FGFR2_CNV_brkpt = rbind(FGFR2_CNV_brkpt, 
                          t_cnv %>% 
                            mutate(start_id = t_start_id, end_id = t_end_id, sample_id = all_samples[i]) %>% 
                            filter(!is.na(start_id) | !is.na(end_id)))
  print(i)
}

t_samp = unique(FGFR2_CNV_brkpt %>% filter(start_id == "Intron 17" & is.na(end_id)) %>% pull(sample_id))

# for samples with entire E1-E17 CNV status
samp_with_brkpt_I17 = c()
for(i in 1:length(t_samp)){
  E1_E17 = FGFR2_CNV_brkpt %>% filter(sample_id == t_samp[i], start_id == "Intron 17" & is.na(end_id)) %>% pull(copyNumber)
  E18 = FGFR2_CNV_brkpt %>% filter(sample_id == t_samp[i], end_id == "Intron 17" & is.na(start_id)) %>% pull(copyNumber)
  if(length(E1_E17) == 1 & length(E18) == 1){
    samp_with_brkpt_I17 = rbind(samp_with_brkpt_I17, data.frame(E1_E17 = E1_E17, E18 = E18, CN_diff = E1_E17-E18, sample_id = t_samp[i]))
  }
}

save.image("~/FGFR/Daniel/R/Nature_figures/data/HMF/HMF_data_partial_amp.RData")
