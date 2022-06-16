library(dplyr)

samp_antisense = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/SB_transposon/SB_transposon_usage.RData")[["uniq_samp_antisense"]]
samp_sense = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/SB_transposon/SB_transposon_usage.RData")[["uniq_samp_sense"]]
samp_both = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/SB_transposon/SB_transposon_usage.RData")[["uniq_both"]]

junc_antisense = matrix(NA, length(samp_antisense), 2)
for(i in 1:length(samp_antisense)){
  t_sj = read.table(paste("/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/", samp_antisense[i], "/SJ.out.tab", sep = ""), sep = "\t", stringsAsFactors = F, header = F)
  t_gtf = read.table(paste("/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/", samp_antisense[i], "/reference/chr7_SB.gtf", sep = ""), sep = "\t", stringsAsFactors = F, header = F)
  t_gtf_sb = t_gtf %>% filter(grepl("transposon", V9), V3 == "exon")
  t_gtf_sb_E17 = t_gtf_sb %>% filter(grepl("ENSMUSE00001300180", V9))
  t_gtf_sb_tp = t_gtf_sb %>% filter(grepl("exon_id transposon", V9))
  t_gtf_sb_E18 = t_gtf_sb %>% filter(grepl("ENSMUSE00000708267", V9))
  
  junc_E17_E18 = as.numeric(as.matrix(t_sj %>% filter(V2 == t_gtf_sb_E18$V5+1, V3 == t_gtf_sb_E17$V4-1) %>% select(V6, V7)))
  junc_E17_E18 = sum(junc_E17_E18)

  junc_E17_TP = as.numeric(as.matrix(t_sj %>% filter(V2 %in% (t_gtf_sb_tp$V5+1), V3 %in% (t_gtf_sb_E17$V4-1)) %>% select(V6, V7)))
  junc_E17_TP = sum(junc_E17_TP)
  
  junc_antisense[i,1] = junc_E17_E18
  junc_antisense[i,2] = junc_E17_TP
}

junc_sense = matrix(NA, length(samp_sense), 2)
for(i in 1:length(samp_sense)){
  t_sj = read.table(paste("/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/", samp_sense[i], "/SJ.out.tab", sep = ""), sep = "\t", stringsAsFactors = F, header = F)
  t_gtf = read.table(paste("/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/", samp_sense[i], "/reference/chr7_SB.gtf", sep = ""), sep = "\t", stringsAsFactors = F, header = F)
  t_gtf_sb = t_gtf %>% filter(grepl("transposon", V9), V3 == "exon")
  t_gtf_sb_E17 = t_gtf_sb %>% filter(grepl("ENSMUSE00001300180", V9))
  t_gtf_sb_tp = t_gtf_sb %>% filter(grepl("exon_id transposon", V9))
  t_gtf_sb_E18 = t_gtf_sb %>% filter(grepl("ENSMUSE00000708267", V9))
  
  junc_E17_E18 = as.numeric(as.matrix(t_sj %>% filter(V2 == t_gtf_sb_E18$V5+1, V3 == t_gtf_sb_E17$V4-1) %>% select(V6, V7)))
  junc_E17_E18 = sum(junc_E17_E18)

  junc_E17_TP = as.numeric(as.matrix(t_sj %>% filter(V2 %in% (t_gtf_sb_tp$V5+1), V3 %in% (t_gtf_sb_E17$V4-1)) %>% select(V6, V7)))
  junc_E17_TP = sum(junc_E17_TP)
  
  junc_sense[i,1] = junc_E17_E18
  junc_sense[i,2] = junc_E17_TP
}

junc_both = matrix(NA, length(samp_both), 2)
# antisense, antisense, sense
for(i in 1:length(samp_both)){
  t_sj = read.table(paste("/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/", samp_both[i], "/SJ.out.tab", sep = ""), sep = "\t", stringsAsFactors = F, header = F)
  t_gtf = read.table(paste("/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/", samp_both[i], "/reference/chr7_SB.gtf", sep = ""), sep = "\t", stringsAsFactors = F, header = F)
  t_gtf_sb = t_gtf %>% filter(grepl("transposon", V9), V3 == "exon")
  t_gtf_sb_E17 = t_gtf_sb %>% filter(grepl("ENSMUSE00001300180", V9))
  t_gtf_sb_tp = t_gtf_sb %>% filter(grepl("exon_id transposon", V9))
  t_gtf_sb_E18 = t_gtf_sb %>% filter(grepl("ENSMUSE00000708267", V9))
  
  junc_E17_E18 = as.numeric(as.matrix(t_sj %>% filter(V2 == t_gtf_sb_E18$V5+1, V3 == t_gtf_sb_E17$V4-1) %>% select(V6, V7)))
  junc_E17_E18 = sum(junc_E17_E18)

  junc_E17_TP = as.numeric(as.matrix(t_sj %>% filter(V2 %in% (t_gtf_sb_tp$V5+1), V3 %in% (t_gtf_sb_E17$V4-1)) %>% select(V6, V7)))
  junc_E17_TP = sum(junc_E17_TP)
  
  junc_both[i,1] = junc_E17_E18
  junc_both[i,2] = junc_E17_TP
}

# sample 1 and 2 have junctions for only antisense insertion and sample 3 have junctions only for sense insertion
junc_antisense = rbind(junc_antisense, junc_both[1:2,])
rownames(junc_antisense) = c(samp_antisense, samp_both[1:2])

junc_sense = rbind(junc_sense, junc_both[3,])
rownames(junc_sense) = c(samp_sense, samp_both[3])


ratio_antisense = rbind(data.frame(ratio = junc_antisense[,1]/rowSums(junc_antisense), type = "E17-E18"),
                       data.frame(ratio = junc_antisense[,2]/rowSums(junc_antisense), type = "E17-TP"))


# randomized the sample order to mix sense and antisense for visualization of barplot
ratio_all = rbind(data.frame(ratio = c(junc_antisense[,1]/rowSums(junc_antisense), junc_sense[,1]/rowSums(junc_sense)), type = "E17-E18"),
                 data.frame(ratio = c(junc_antisense[,2]/rowSums(junc_antisense), junc_sense[,2]/rowSums(junc_sense)), type = "E17-TP"))
ratio_all$sample = rep(c(rownames(junc_antisense), rownames(junc_sense)), 2)

ratio_all = ratio_all %>% mutate(orientation = rep(NA, nrow(ratio_all)),
                                 orientation = replace(orientation, ratio_all$sample %in% rownames(junc_antisense), "antisense"),
                                 orientation = replace(orientation, ratio_all$sample %in% rownames(junc_sense), "sense"))
ratio_all$orientation = factor(ratio_all$orientation, levels = c("antisense", "sense"))

ratio_all = ratio_all[sample.int(nrow(ratio_all), nrow(ratio_all), replace = F),]
ratio_all = ratio_all %>% mutate(sample = factor(sample, levels = (ratio_all[ratio_all$type == "E17-TP",] %>% arrange(desc(ratio)) %>% pull(sample))))
ratio_all$annot = paste(ratio_all$type, ratio_all$orientation, sep = "_")

rm(i, t_gtf, t_gtf_sb, t_gtf_sb_tp, t_gtf_sb_E17, t_gtf_sb_E18, t_sj)
save.image("~/FGFR/Daniel/R/Nature_figures/data/SB_transposon/SB_transposon_junctions.RData")
