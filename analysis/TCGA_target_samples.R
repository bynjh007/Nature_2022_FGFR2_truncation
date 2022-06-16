# libraries
library(dplyr)

FGFR2_alt = read.table("~/FGFR/Daniel/R/Nature_figures/data/TCGA/cBioportal_alterations_across_samples.tsv", header = T, stringsAsFactors = F, sep = "\t")

########################################################
# collection of different types of variants
########################################################
# Mutations (hotspot, oncogenic, E18 truncation-frameshift, non-sense mutation)
load("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_mutation_processing.RData")
a = FGFR2_mut %>% filter(grepl("E18|Hotspot", Mutation_class)) %>% pull(SAMPLE_ID)

# low C1 expression (exon expression)
target_samp_C1_Med = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_exon_junction_processing.RData")[["target_samp_C1_Med"]]
exon_ratio_C1_Med = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_exon_junction_processing.RData")[["exon_ratio_C1_Med"]]
E17_junc_ratio = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_exon_junction_processing.RData")[["E17_junc_ratio"]]
samples_exon = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_exon_junction_processing.RData")[["samples_exon"]]
samples_junc = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_exon_junction_processing.RData")[["samples_junc"]]

# high C3 usage (junction reads from E17-C3)
target_samp_C3_junc = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_exon_junction_processing.RData")[["target_samp_C3_junc"]]

# high C4 usage (junction reads from E17-C4)
target_samp_C4_junc = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_exon_junction_processing.RData")[["target_samp_C4_junc"]]

# RSEM values
FGFR2_exp_rsem = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_exon_junction_processing.RData")[["FGFR2_exp_rsem"]]
FGFR2_exp_rsem = log2(FGFR2_exp_rsem+1)

# Copy number
FGFR2_cnv = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_snp6_brkpt_CNV_processing.RData")[["FGFR2_cnv"]]
b = FGFR2_cnv$SAMPLE_ID[which(FGFR2_cnv$FGFR2=="Amp" & FGFR2_cnv$SAMPLE_ID %in% names(FGFR2_exp_rsem))]

# breakpoints from SNP array
FGFR2_RE_snp6_I17 = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_snp6_brkpt_CNV_processing.RData")[["FGFR2_RE_snp6_I17"]]
brkpt_I17_annot = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_snp6_brkpt_CNV_processing.RData")[["brkpt_I17_annot"]]
c = unique(FGFR2_RE_snp6_I17$sample_id)

# partial amplification
FGFR2_partial_amp = brkpt_I17_annot %>% filter(after_I17>=0.3, diff>=0.3)

# Reported fusions
load("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_fusion_reported_processing.RData")
d = unique(FGFR2_fusion$sample_id)

####################################################################
# all samples with above variants
####################################################################
total_var_samples = unique(c(a, b, c, d, target_samp_C1_Med, target_samp_C3_junc, target_samp_C4_junc))

t_cnv = FGFR2_cnv %>% mutate(FGFR2 = ifelse(SAMPLE_ID %in% FGFR2_partial_amp$sample_id, "Partial amp", FGFR2))
df_samples = data.frame(TCGA_code_all = FGFR2_alt$Study.ID[match(total_var_samples, FGFR2_alt$Sample.ID)],
                        E18_trunc_mut = FGFR2_mut$Mutation_class[match(total_var_samples, FGFR2_mut$SAMPLE_ID)],
                        FGFR2_amp = t_cnv$FGFR2[match(total_var_samples, t_cnv$SAMPLE_ID)],
                        FGFR2_exp = FGFR2_exp_rsem[match(total_var_samples, names(FGFR2_exp_rsem))],
                        FGFR2_RE = FGFR2_fusion$canonical[match(total_var_samples, FGFR2_fusion$sample_id)],
                        FGFR2_RE_loc = FGFR2_fusion$breakpoint[match(total_var_samples, FGFR2_fusion$sample_id)],
                        SNP6_BP_I17 = total_var_samples %in% unique(unique(FGFR2_RE_snp6_I17$sample_id)),
                        Exon_C1_lack = total_var_samples %in% target_samp_C1_Med,
                        Exon_C1_ratio_exp = exon_ratio_C1_Med[match(total_var_samples, samples_exon$sample_id)],
                        Exon_C1_ratio_junc = E17_junc_ratio[2,match(total_var_samples, samples_junc$sample_id)],
                        Exon_C3_high = total_var_samples %in% target_samp_C3_junc,
                        Exon_C3_ratio_junc = E17_junc_ratio[4,match(total_var_samples, samples_junc$sample_id)], 
                        Exon_C4_high = total_var_samples %in% target_samp_C4_junc,
                        Exon_C4_ratio_junc = E17_junc_ratio[1,match(total_var_samples, samples_junc$sample_id)],
                        stringsAsFactors = F)

# One of the samples with lower C1 usage were not classified into certain category because of filtering criteria
# But that sample carried fusion, so it makes sense to categorize it to C1 lack
df_samples$Exon_C1_lack[which(df_samples$Exon_C1_ratio_exp<0.3963389)] = TRUE

df_samples = df_samples %>% mutate(E18_trunc_mut = replace(E18_trunc_mut, E18_trunc_mut == "Others", NA))


# reannotating for oncoplot
df_samples = 
  df_samples %>%
  mutate(TCGA_code_all = replace(TCGA_code_all, rownames(df_samples) == "TCGA-61-2095-02", "ov_tcga"),
         TCGA_code_all = stringr::str_replace(toupper(TCGA_code_all), "_TCGA", ""),
         FGFR2_amp = ifelse(FGFR2_amp == "Amp" | FGFR2_amp == "Partial amp" , FGFR2_amp, NA),
         FGFR2_RE = case_when(FGFR2_RE == "5'-fusion" ~ "3'-RE", FGFR2_RE == "3'-fusion" ~ "5'-RE"),
         SNP6_BP_I17 = ifelse(SNP6_BP_I17 == T, T, NA),
         Exon_C1_lack = ifelse(Exon_C1_lack == T, T, NA),
         Exon_C1_ratio_exp = ifelse(is.infinite(Exon_C1_ratio_exp) | is.na(Exon_C1_ratio_exp), 1, Exon_C1_ratio_exp),
         Exon_C3_high = ifelse(Exon_C3_high == T, T, NA),
         Exon_C3_ratio_junc = ifelse(is.infinite(Exon_C3_ratio_junc) | is.na(Exon_C3_ratio_junc), 0, Exon_C3_ratio_junc),
         Exon_C4_high = ifelse(Exon_C4_high == T, T, NA),
         Exon_C4_ratio_junc = ifelse(is.infinite(Exon_C4_ratio_junc) | is.na(Exon_C4_ratio_junc), 0, Exon_C4_ratio_junc),
         TCGA_code = replace(TCGA_code_all, TCGA_code_all %in%names(which(table(TCGA_code_all)<5)), "others"))

rownames(df_samples) = total_var_samples

#write.table(df_samples, "~/FGFR/Daniel/R/files/TCGA_target_samples_revision.txt", col.names = T, row.names = T, sep = "\t")


####################################################################
# post-processing of each type of variants for oncoplot
####################################################################

load("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_RE_RNAseq_processing.RData")
# all the breakpoint from STAR-Fusion is ref:in-strand type
t_right = stringr::str_split_fixed(FGFR2_star_fusion$right_annot, ":", 5)
t_right[,3] = paste("I", as.numeric(stringr::str_split_fixed(t_right[,3], " ", 2)[,2])-1, sep = "")
t_right = apply(t_right, 1, paste, collapse = ":")
a = FGFR2_star_fusion %>% mutate(Source = "star_fusion",
                                 RE_location = ifelse(startsWith(X.FusionName, "FGFR2"), "upstream", "downstream"),
                                 RE_type = ifelse(PROT_FUSION_TYPE == "INFRAME", "inframe", "out-of-frame"),
                                 FGFR2_brkpt = ifelse(grepl("chr10:121483698:-", LeftBreakpoint) | grepl("chr10:121480021:-|chr10:121480113:-", RightBreakpoint), "I17", "5' to E17"),
                                 C_trunc = LeftBreakpoint %in% "chr10:121483698:-",
                                 left_annot = gsub("Exon ", "I", left_annot),
                                 right_annot = t_right) %>%
  select(X.FusionName, LeftBreakpoint, RightBreakpoint, JunctionReadCount, FFPM, sample_id, RE_type, RE_location, 
         C_trunc, FGFR2_brkpt, left_annot, right_annot, Source)
colnames(a) = c("Fusion_name", "LeftBreakpoint", "RightBreakpoint", "Junc_count", "FFPM", "Sample_id", "RE_type", 
                "RE_location", "C_trunc", "FGFR2_brkpt", "left_annot", "right_annot", "Source")

# extract fusion information from RE_from_chimeric (** note that brkpt from RE_from_chimeric is one base different with star-fusion results **)
b = FGFR2_RE %>% mutate(Source = "RE_from_chimeric",
                        up.id = as.character(up.id), up.id = gsub("Exon ", "E", up.id), up.id = gsub("Intron ", "I", up.id),
                        up.id = ifelse(up.splice_type == "ref:in-strand", gsub("E", "I", up.id), up.id),
                        down.id = as.character(down.id), down.id = gsub("Exon ", "E", down.id), down.id = gsub("Intron ", "I", down.id),
                        down.id = ifelse(down.splice_type == "ref:in-strand", 
                                         paste("I", as.numeric(stringr::str_split_fixed(down.id, "E", 2)[,2])-1, sep = ""), down.id),
                        Fusion_name = paste(up.symbol, down.symbol, sep = "--"),
                        LeftBreakpoint = up.brkpt, RightBreakpoint = down.brkpt, Junc_count = count,
                        FFPM = FFPM, Sample_id = sample_id,
                        RE_location = ifelse(startsWith(Fusion_name, "FGFR2"), "upstream", "downstream"),
                        FGFR2_brkpt = ifelse(startsWith(Fusion_name, "FGFR2"), up.id, down.id),
                        FGFR2_brkpt = replace(FGFR2_brkpt, FGFR2_brkpt != "I17" & FGFR2_brkpt!="E18", "5' to E17"),
                        C_trunc_can = ((up.brkpt %in% "chr10:121483698:-") & (up.splice_type %in% "ref:in-strand") & (up.symbol %in% "FGFR2")),
                        C_trunc_noncan = ((FGFR2_brkpt %in% "I17") & (up.splice_type %in% "non_ref:in-strand") & (up.symbol %in% "FGFR2")),
                        C_trunc = ((FGFR2_brkpt %in% "I17") & grepl("in-strand", up.splice_type) & (up.symbol %in% "FGFR2")),
                        t_type = paste(up.splice_type, down.splice_type, sep = "//"),
                        RE_type = rep(NA, nrow(FGFR2_RE)),
                        RE_type = replace(RE_type, grepl("in-strand", t_type) & grepl("out-of-strand", t_type), "out-of-strand"),
                        RE_type = replace(RE_type, grepl("intergenic", t_type), "intergenic"),
                        RE_type = replace(RE_type, grepl("in-strand", up.splice_type) & grepl("in-strand", down.splice_type), "out-of-frame"),
                        left_annot = paste(up.symbol, up.txname, up.id, up.splice_type, sep = ":"),
                        right_annot = paste(down.symbol, down.txname, down.id, down.splice_type, sep = ":")) %>%
  select(Fusion_name, LeftBreakpoint, RightBreakpoint, Junc_count, FFPM, Sample_id, RE_type, RE_location, 
         C_trunc, FGFR2_brkpt, left_annot, right_annot, Source)

# combine the RE and Fusions
FGFR2_RE_all = rbind(a,b)

# remove out-of-strand fusion from upstream
FGFR2_RE_all = FGFR2_RE_all %>% filter(grepl("in-strand", left_annot)) %>%
  mutate(RE_type = factor(RE_type, levels = c("inframe", "out-of-frame", "intergenic", "out-of-strand"))) %>%
  filter(FFPM>=0.05)

# prioritizing the fusions
samp_fus = unique(FGFR2_RE_all$Sample_id)
FGFR2_RE_uniq = c()
for(i in 1:length(samp_fus)){
  t = FGFR2_RE_all %>% filter(Sample_id == samp_fus[i]) %>%
    arrange(desc(C_trunc), desc(RE_location), desc(FGFR2_brkpt), RE_type, desc(FFPM))
  FGFR2_RE_uniq = rbind(FGFR2_RE_uniq, t[1,])
}


#####################################################################################
# post-processing of each type of variants for oncoplot
# marking previously reported fusions to newly analyzed fusions/REs
#####################################################################################
# marking the known fusions to the re-analyzed fusions/RE
t = apply(FGFR2_fusion[, c("Breakpoint1", "Breakpoint2", "sample_id")], 1, paste, collapse = "_")
a = FGFR2_RE_uniq %>% 
  mutate(id = paste(LeftBreakpoint, RightBreakpoint, substr(Sample_id, 1, 15), sep = "_")) %>%
  mutate(known = id %in% t)
FGFR2_RE_uniq$known = a$known


#####################################################################################
# target samples for oncoplot
#####################################################################################
TCGA_target_samples = df_samples
TCGA_target_samples = data.frame(TCGA_target_samples, 
                                 FGFR2_RE_uniq[match(rownames(TCGA_target_samples), substr(FGFR2_RE_uniq$Sample_id, 1, 15)),])
TCGA_target_samples$Sample = rownames(TCGA_target_samples)

# inter- and intrachromosomal location
TCGA_target_samples = TCGA_target_samples %>% 
  mutate(FGFR2_chr_partner = ifelse(RE_location == "upstream", 
                                    stringr::str_split_fixed(RightBreakpoint, ":", 3)[,1], 
                                    stringr::str_split_fixed(LeftBreakpoint, ":", 3)[,1])) %>%
  mutate(FGFR2_chr_SV = ifelse(grepl("chr10", FGFR2_chr_partner), "intrachromosomal", "interchromosomal"),
         FGFR2_chr_SV = replace(FGFR2_chr_SV, is.na(FGFR2_chr_partner), NA))

# filtering out the C1 low samples without any detected fusion events
ind_filt_1 = TCGA_target_samples %>% mutate(index = 1:nrow(TCGA_target_samples)) %>% 
  filter(is.na(E18_trunc_mut), is.na(FGFR2_amp), is.na(FGFR2_RE), Exon_C1_lack == TRUE, 
         is.na(Exon_C3_high), is.na(Exon_C4_high), is.na(RE_type)) %>% pull(index)
# filtering out the SNP6 breakpoint without any detected fusion events
ind_filt_2 = TCGA_target_samples %>% mutate(index = 1:nrow(TCGA_target_samples)) %>% 
  filter(is.na(E18_trunc_mut), is.na(FGFR2_amp), is.na(FGFR2_RE), is.na(Exon_C1_lack), is.na(Exon_C3_high),
         SNP6_BP_I17 == TRUE, is.na(Exon_C4_high), is.na(RE_type)) %>% pull(index)
# filtering out the fusions without available RNA-seq but analyzed previously
ind_filt_3 = TCGA_target_samples %>% mutate(index = 1:nrow(TCGA_target_samples)) %>% 
  filter(!is.na(FGFR2_RE), is.na(RE_type)) %>% pull(index)

TCGA_target_samples = TCGA_target_samples[-c(ind_filt_1, ind_filt_2, ind_filt_3),]


#####################################################################################
# Gender information
#####################################################################################
clinicial = read.table("~/FGFR/Daniel/R/Nature_figures/data/TCGA/combined_study_clinical_data.tsv", 
                       header = T, sep = "\t", stringsAsFactors = F)

tumor_type = c("STAD", "BRCA", "UCEC", "BLCA", "HNSC", "LIHC", "OV", 
               "SARC", "CHOL", "LUSC", "LUAD","PAAD", "THCA", "GBM")

df_samples_oncoplot = TCGA_target_samples %>% mutate(TCGA_code = TCGA_code_all) %>%
  mutate(TCGA_code = replace(TCGA_code_all, !(TCGA_code_all %in% tumor_type), "others"))

df_samples_oncoplot = df_samples_oncoplot %>% mutate(Gender = clinicial$Sex[match(Sample, clinicial$Sample.ID)])
# all the missing samples from cBioportal is Female
df_samples_oncoplot$Gender[is.na(df_samples_oncoplot$Gender)] = "Female"
df_samples_oncoplot$Gender[df_samples_oncoplot$Gender == ""] = "Female"


#####################################################################################
# Gene expression
# download gene expression data for the target samples from GDC
#####################################################################################
library(TCGAbiolinks)
target_samples_revision = read.table("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_target_samples_revision.txt", header = T, sep = "\t", stringsAsFactors = F)
target_samples_revision = rownames(target_samples_revision)

FGFR2_alt = read.table("~/FGFR/Daniel/R/Nature_figures/data/TCGA/cBioportal_alterations_across_samples.tsv", header = T, stringsAsFactors = F, sep = "\t")
target_to_download = data.frame(sample_id = target_samples_revision, Study = FGFR2_alt$Study.ID[match(target_samples_revision, FGFR2_alt$Sample.ID)])
target_to_download$Study[target_to_download$sample_id == "TCGA-61-2095-02"] = "ov_tcga"
target_to_download$Study = paste("TCGA-", toupper(gsub("_tcga", "", target_to_download$Study)), sep = "")

projects <- unique(target_to_download$Study)
projects = c(projects, "TCGA-COAD", "TCGA-READ")
projects = projects[projects!="TCGA-COADREAD"]
# 
# dir.create("/DATA/projects/j.bhin/Daniel_FGFR2/TCGA/RNA/revision")
# setwd("/DATA/projects/j.bhin/Daniel_FGFR2/TCGA/RNA/revision")
query <- GDCquery(project = projects,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - FPKM",
                  barcode = target_to_download$sample_id)
match.file.cases <- getResults(query)
# GDCdownload(query, method = "api")

FGFR2_fpkm = c()
exp_samp = c()
for(i in projects){
  t_path = paste("/DATA/projects/j.bhin/Daniel_FGFR2/TCGA/RNA/revision/GDCdata/", i, 
                 "/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification", sep = "")
  t_files = list.files(t_path)
  for(j in t_files){
    t_gz = list.files(paste(t_path, j, sep = "/"))
    t_cnt = read.table(gzfile(paste(t_path, j, t_gz, sep = "/")), sep = "\t", stringsAsFactors = F, header = F)
    FGFR2_fpkm = cbind(FGFR2_fpkm, t_cnt[,2])
  }
  exp_samp = c(exp_samp, t_files)
}
FGFR2_fpkm = FGFR2_fpkm[grepl("ENSG00000066468", t_cnt$V1),]
names(FGFR2_fpkm) = substr(match.file.cases$sample.submitter_id[match(exp_samp, match.file.cases$id)], 1, 15)



df_samples_oncoplot$FGFR2_exp = log2(FGFR2_fpkm[match(df_samples_oncoplot$Sample, names(FGFR2_fpkm))]+1)
# all the samples with misssing FGFR2 expression (due to the lack of raw file) are the ones with FGFR2 amplification 
# filling the values with median expression of FGFR2 amplified samples
df_samples_oncoplot = df_samples_oncoplot %>% mutate(FGFR2_exp = replace(FGFR2_exp, is.na(FGFR2_exp), 
                                                                         median(FGFR2_exp[which(FGFR2_amp == "Amp")], na.rm = T)))

save.image("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_target_samples.RData")
