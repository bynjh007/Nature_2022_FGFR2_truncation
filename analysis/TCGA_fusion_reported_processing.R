library(GenomicRanges)
library(dplyr)

################################################
# FGFR2 genomic information
################################################
trx_info = R.utils::loadToEnv("~/FGFR/Daniel/R/R_files/FGFR2_RE_TCGA_20200804.RData")[["trx_info"]]
ind = data.frame(trx_info$trx) %>% filter(grepl("ENST00000358487", tx_name)) %>% pull(tx_id)
t_exon = trx_info$trx_exon[[ind]]
t_exon$id = paste("E", 1:length(t_exon), sep = "")
t_intron = rev(trx_info$trx_intron[[ind]])
t_intron$id = paste("I", 1:length(t_intron), sep = "")
FGFR2_info = sort(c(t_exon, t_intron))

################################################
# TCGA fusion genes from RNA-seq (previously reported from Cell Reports)
################################################
fusion = read.table("~/FGFR/Daniel/TCGA/fusion_RNAseq/TCGA_fusion.txt", header = T, sep = "\t", stringsAsFactors = F)
FGFR2_fusion = fusion %>% filter(grepl("FGFR2", Fusion))
FGFR2_fusion$sample_id = apply(stringr::str_split_fixed(FGFR2_fusion$Sample, "-", 7)[,1:4], 1, paste, collapse = "-")
FGFR2_fusion$submitter_id = apply(stringr::str_split_fixed(FGFR2_fusion$Sample, "-", 7)[,1:3], 1, paste, collapse = "-")
FGFR2_fusion$TN = stringr::str_split_fixed(FGFR2_fusion$Sample, "-", 7)[,4]

t = cbind(as.numeric(stringr::str_split_fixed(FGFR2_fusion$Breakpoint1, ":", 3)[,2]),
          as.numeric(stringr::str_split_fixed(FGFR2_fusion$Breakpoint2, ":", 3)[,2]))
t2 = t[,1]
t2[!startsWith(FGFR2_fusion$Fusion, "FGFR2")] = t[!startsWith(FGFR2_fusion$Fusion, "FGFR2"),2]

fusion_breakpoint = rep(NA, nrow(t))
fusion_breakpoint[t2<end(FGFR2_info[1])] = "E18"
fusion_breakpoint[t2==121480113|t2==start(FGFR2_info[3])|t2==end(FGFR2_info[1])] = "I17" # alternative exon 18
fusion_breakpoint[t2>start(FGFR2_info[3])] = "E1-E17" # 

FGFR2_fusion = FGFR2_fusion %>% mutate(breakpoint = fusion_breakpoint)
FGFR2_fusion$canonical = startsWith(FGFR2_fusion$Fusion, "FGFR2")
FGFR2_fusion$canonical[which(FGFR2_fusion$canonical == TRUE)] = "5'-fusion"
FGFR2_fusion$canonical[which(FGFR2_fusion$canonical == "FALSE")] = "3'-fusion"

FGFR2_fusion = FGFR2_fusion %>% 
  mutate(breakpoint = factor(breakpoint, levels = c("I17", "E18", "E1-E17")), 
         canonical = factor(canonical, levels = c("5'-fusion", "3'-fusion"))) %>% 
  arrange(canonical, breakpoint)


FGFR2_fusion$sample_id = substr(FGFR2_fusion$sample_id, 1, 15)
FGFR2_fusion = FGFR2_fusion[c(which(startsWith(FGFR2_fusion$Fusion, "FGFR2")), which(!startsWith(FGFR2_fusion$Fusion, "FGFR2"))),]

samples_fusion = read.table("~/FGFR/Daniel/TCGA/fusion_RNAseq/samples_used_fusion.txt", header = T, sep = "\t", stringsAsFactors = F)
colnames(samples_fusion) = c("Sample", "TCGA_code")
samples_fusion$submitter_id = apply(stringr::str_split_fixed(samples_fusion$Sample, "-", 5)[,1:3], 1, paste, collapse = "-")
samples_fusion$sample_id = apply(stringr::str_split_fixed(samples_fusion$Sample, "-", 5)[,1:4], 1, paste, collapse = "-")
samples_fusion$sample_id = substr(samples_fusion$sample_id, 1, 15)

save(FGFR2_fusion, samples_fusion, file = "~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_fusion_reported_processing.RData")
