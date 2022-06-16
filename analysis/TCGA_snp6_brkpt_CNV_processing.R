################################################
# retrieving ID
################################################
TCGAtranslateID = function(file_ids) {
  info = files() %>%
    filter( ~ file_id %in% file_ids) %>%
    select('cases.samples.submitter_id') %>%
    results_all()
  # The mess of code below is to extract TCGA barcodes
  # id_list will contain a list (one item for each file_id)
  # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
  id_list = lapply(info$cases,function(a) {
    a[[1]][[1]][[1]]})
  # so we can later expand to a data.frame of the right size
  barcodes_per_file = sapply(id_list,length)
  # And build the data.frame
  return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                    sample_id = unlist(id_list), stringsAsFactors = F))
}

library(dplyr)
library(tidyr)
manifest <- read.table("~/FGFR/Daniel/TCGA/snp6/gdc_manifest.2019-10-04.txt", header = T, stringsAsFactors = F)
manifest_masked = manifest %>% filter(grepl("nocnv", filename))

# Because it doesn't retreive IDs more than certain numbers, IDs might be partitoned.
library(GenomicDataCommons)
library(stringr)
file_uuids <- manifest_masked$id
snp6_TCGA_id = TCGAtranslateID(file_uuids)
snp6_TCGA_id = data.frame(snp6_TCGA_id, 
                          submitter_id = apply(str_split_fixed(snp6_TCGA_id$sample_id, "-", 4)[,1:3], 1, paste, collapse = "-"),
                          TN = str_split_fixed(snp6_TCGA_id$sample_id, "-", 4)[,4], stringsAsFactors = F)

detach("package:GenomicDataCommons", unload=TRUE)

# selecting the tumors
snp6_TCGA_id_tumor = snp6_TCGA_id %>% filter(grepl("01|02|03|04|05|06|07", TN))


################################################
# target data loading
################################################
seg_TCGA = lapply(snp6_TCGA_id_tumor$file_id, function(X){
  dir_samp = paste("~/FGFR/Daniel/TCGA/snp6/", X, sep = "")
  t = list.files(dir_samp)
  t_files = t[grep("txt", t)]
  read.table(paste(dir_samp, "/", t_files, sep = ""), sep = "\t", header = T, stringsAsFactors = F)
})


################################################
# truncation in FGFR2
################################################
trx_info = R.utils::loadToEnv("~/FGFR/Daniel/R/R_files/FGFR2_RE_TCGA_20200804.RData")[["trx_info"]]
ind = data.frame(trx_info$trx) %>% filter(grepl("ENST00000358487", tx_name)) %>% pull(tx_id)
t_exon = trx_info$trx_exon[[ind]]
t_exon$id = paste("E", 1:length(t_exon), sep = "")
t_intron = rev(trx_info$trx_intron[[ind]])
t_intron$id = paste("I", 1:length(t_intron), sep = "")
FGFR2_info = sort(c(t_exon, t_intron))

# take any samples with the breakpoint within FGFR2
FGFR2_RE_snp6 = c()
for(i in 1:length(seg_TCGA)){
  FGFR2_RE_snp6 = rbind(FGFR2_RE_snp6,
                        seg_TCGA[[i]] %>% filter(Chromosome == "10", End >= 121478341, End <= 121598084) %>%
                           mutate(sample_id = snp6_TCGA_id_tumor$sample_id[i], sample_index = i),
                        seg_TCGA[[i]] %>% filter(Chromosome == "10", Start >= 121478341, Start <= 121598084) %>%
                            mutate(sample_id = snp6_TCGA_id_tumor$sample_id[i], sample_index = i))
  print(i)
}
FGFR2_RE_snp6 = unique(FGFR2_RE_snp6)

t_start = (FGFR2_RE_snp6 %>% mutate(Chromosome = "chr10"))[,c(2,3,3)]
colnames(t_start) = c("chromosome", "start", "end")
t_start = GRanges(t_start)
t_start_id = rep(NA, length(t_start))
t_start_id[queryHits(findOverlaps(t_start, FGFR2_info))] = FGFR2_info$id[subjectHits(findOverlaps(t_start, FGFR2_info))]

t_end = (FGFR2_RE_snp6 %>% mutate(Chromosome = "chr10"))[,c(2,4,4)]
colnames(t_end) = c("chromosome", "start", "end")
t_end = GRanges(t_end)
t_end_id = rep(NA, length(t_end))
t_end_id[queryHits(findOverlaps(t_end, FGFR2_info))] = FGFR2_info$id[subjectHits(findOverlaps(t_end, FGFR2_info))]

FGFR2_RE_snp6$start_id = t_start_id
FGFR2_RE_snp6$end_id = t_end_id

FGFR2_RE_snp6$sample_id = substr(FGFR2_RE_snp6$sample_id, 1, 15)

t_samp = unique(FGFR2_RE_snp6 %>% filter(start_id == "I17" | end_id == "I17") %>% pull(sample_id))
FGFR2_RE_snp6_I17 = FGFR2_RE_snp6 %>% filter(sample_id %in% t_samp)

brkpt_I17_annot =c()
for(i in 1:length(t_samp)){
  t_re = FGFR2_RE_snp6 %>% filter(sample_id == t_samp[i])
  if(nrow(t_re)==2){
    if(is.na(t_re$start_id[1]) & is.na(t_re$end_id[2])){
      t_annot = data.frame(before_I17 = t_re$Segment_Mean[1], 
                           after_I17 = t_re$Segment_Mean[2], 
                           diff = t_re$Segment_Mean[2] - t_re$Segment_Mean[1]) %>%
        mutate(sample_id = t_samp[i])
      brkpt_I17_annot = rbind(brkpt_I17_annot, t_annot)
    } 
  }
}


################################################
# Copy number variation of FGFR2
################################################
FGFR2_alt = read.table("~/TCGA/cBioportal/alterations_across_samples.tsv", header = T, stringsAsFactors = F, sep = "\t")
FGFR2_cnv = read.table("~/TCGA/cBioportal/cna.txt", header = T, stringsAsFactors = F, sep = "\t")
FGFR2_cnv$FGFR2[which(FGFR2_cnv$SAMPLE_ID %in% FGFR2_alt$Sample.ID[which(FGFR2_alt$FGFR2..AMP == "no alteration")])] = "--"
FGFR2_cnv$FGFR2[which(FGFR2_cnv$FGFR2 == 2)] = "Amp"
FGFR2_cnv$FGFR2[which(FGFR2_cnv$FGFR2 == "-2")] = "Homdel"

save.image("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_snp6_brkpt_CNV_processing.RData")
