####################################################
# Sourcing functions to be used
####################################################
source("~/FGFR/Daniel/R/Nature_figures/functions/Functions_FM_FGFR2_RE_annotation.R")

####################################################
# loading data
####################################################
#Fusion list were cleaned from Oncoplot_2 where there are some mislabeled fusions amonng the samples with multiple fusions (10 samples).  
#For the samples with multiple fusions, representative fusion (FGFR2 is upstream of the fusion) was selected in each sample.  
#Also, FGFR2-BCCIP;DHX32 was removed because of inconsistency between matched gene and RE position
onco_data = read.table("~/FGFR/Daniel/R/Nature_figures/data/FM/FM_data/One patient per line_combined_Oncoplot_3.txt", 
                       header = T, sep = "\t", stringsAsFactors = F, quote = "", fill = T)

# combine all the minor tumor types into "others" group
#t = table(onco_data$TCGA_Oncoplot)
#onco_data$TCGA_Oncoplot[which(onco_data$TCGA_Oncoplot %in% names(t)[t<10])] = "Other"

# loading FGFR2 mutations
FGFR2_mut_all = read.table("~/FGFR/Daniel/R/Nature_figures/data/FM/FGFR2_mutation_all_MutationMapper.txt", header = T, sep = "\t", stringsAsFactors = F, quote = "", fill = T)
FGFR2_mut_hotspot = FGFR2_mut_all %>% mutate(Label = case_when(grepl("S252", Protein_Change) ~ "S252",
                                                           grepl("N549", Protein_Change) ~ "N549",
                                                           grepl("C382", Protein_Change) ~ "C382",
                                                           grepl("K659", Protein_Change) ~ "K659")) %>% 
                        filter(!is.na(Label), !grepl("fs", Protein_Change))

# FGFR2 RE
FGFR2_RE_in_strand = onco_data[,1:20] %>% filter(RE.type == "In-strand RE (Fusion)") %>% 
  mutate(Gene_1 = stringr::str_split_fixed(Rearrangement..RE., "-", 2)[,1],
         Gene_2 = stringr::str_split_fixed(Rearrangement..RE., "-", 2)[,2],
         Position.1_RE = round(Position.1_RE), Position.2_RE = round(Position.2_RE))
  
ind = match(FGFR2_RE_in_strand$Gene_1, hg19_info$trx_gene_AA$SYMBOL)
which(is.na(ind))
ind = match(FGFR2_RE_in_strand$Gene_2, hg19_info$trx_gene_AA$SYMBOL)
which(is.na(ind))
FGFR2_RE_in_strand$Gene_2[which(is.na(ind))] = c("LOC387723", "C11orf30", "LOC100132146") # LINC00959, EMSY, FAM240A

G_RE_up = FGFR2_RE_in_strand %>% 
  mutate(str = hg19_info$trx_gene_AA$STRAND[match(Gene_1, hg19_info$trx_gene_AA$SYMBOL)]) %>%
  mutate(str = replace(str, is.na(str), ".")) %>% 
  mutate(pos1_s = Position.1_RE, pos1_e = Position.1_RE) %>%
  select(Chr.1_RE, pos1_s, pos1_e, str)
colnames(G_RE_up) = c("seqnames", "tx_start", "tx_end", "strand")
G_RE_up = makeGRangesFromDataFrame(G_RE_up)

G_RE_down = FGFR2_RE_in_strand %>% 
  mutate(str = hg19_info$trx_gene_AA$STRAND[match(Gene_2, hg19_info$trx_gene_AA$SYMBOL)]) %>%
  mutate(str = replace(str, is.na(str), ".")) %>% 
  mutate(pos2_s = Position.2_RE, pos2_e = Position.2_RE) %>%
  select(Chr.2_RE, pos2_s, pos2_e, str)
colnames(G_RE_down) = c("seqnames", "tx_start", "tx_end", "strand")
G_RE_down = makeGRangesFromDataFrame(G_RE_down)

brkpt_annot_FM = c()
for(k in 1:nrow(FGFR2_RE_in_strand)){
  t_out = brk_annot_framework(G_up = G_RE_up[k], G_down = G_RE_down[k])
  brkpt_annot_FM = rbind(brkpt_annot_FM, data.frame(t_out %>% select(-c("up.AA_seq", "down.AA_seq")), effect = RE_type(t_out)))
  print(k)
}
brkpt_annot_FM = brkpt_annot_FM


# oncoplot table
onco_data_new = onco_data
# be careful not to use "replace" function which give rise to outcome
onco_data_new$Rearrangement..RE.[which(onco_data_new$RE.type == "In-strand RE (Fusion)")] = paste(brkpt_annot_FM$up.gene_symbol, brkpt_annot_FM$down.gene_symbol, sep = "-")
onco_data_new$RE.type[which(onco_data_new$RE.type == "In-strand RE (Fusion)")] = brkpt_annot_FM$effect

# FGFR2 amplification
FGFR2_amp = character(nrow(onco_data_new))
FGFR2_amp[!is.na(onco_data_new$FGFR2.Amp)] = "Full amplification"
FGFR2_amp[!is.na(onco_data_new$FGFR2.E1.17.Amp)] = "E1-17 amplification"

# FGFR2 truncating mutations
FGFR2_mut_trunc = onco_data_new[, c("E18.Truncating.Mutation", "E18.Splice.acceptor.Mutation")]
FGFR2_mut_trunc = FGFR2_mut_trunc %>% 
    mutate(E18.Truncating.Mutation = case_when(as.numeric(substr(E18.Truncating.Mutation, 2,4))<=783 ~ "E18-trunc-onco",
                                               as.numeric(substr(E18.Truncating.Mutation, 2,4))>783 ~ "E18-trunc-others"),
           E18.Splice.acceptor.Mutation = replace(E18.Splice.acceptor.Mutation, !is.na(E18.Splice.acceptor.Mutation), "E18-splice")) %>%
  as.matrix()
FGFR2_mut_trunc = apply(FGFR2_mut_trunc, 1, paste, collapse = ";")
FGFR2_mut_trunc = gsub("NA;NA", NA, FGFR2_mut_trunc)
FGFR2_mut_trunc = gsub(";NA", "", FGFR2_mut_trunc)
FGFR2_mut_trunc = gsub("NA;", "", FGFR2_mut_trunc)


# FGFR2 rearrangement
FGFR2_re_loc = rep(NA, nrow(onco_data_new))
FGFR2_re_loc[which(onco_data_new$Chr.1_RE == "chr10" & onco_data_new$Chr.2_RE == "chr10")] = "Intrachromosomal"
FGFR2_re_loc[which(onco_data_new$Chr.1_RE == "chr10" & onco_data_new$Chr.2_RE != "chr10")] = "Interchromosomal"
FGFR2_re_loc[which(onco_data_new$Chr.1_RE != "chr10" & onco_data_new$Chr.2_RE == "chr10")] = "Interchromosomal"

# FGFR2 truc location
onco_data_new$FGFR2.truncation.location[which(onco_data_new$FGFR2.truncation.location == "Intron 17_3' to C3")] = "I17"
onco_data_new$FGFR2.truncation.location[which(onco_data_new$FGFR2.truncation.location == "Intron 17_5' to C3")] = "I17"
onco_data_new$FGFR2.truncation.location[which(onco_data_new$FGFR2.truncation.location == "5' to Intron 17")] = "5' to E17"
onco_data_new$FGFR2.truncation.location[grepl("Exon 18", onco_data_new$FGFR2.truncation.location)] = "E18"
# onco_data$FGFR2.truncation.location[onco_data$FGFR2.truncation.location == "Exon 18_CDS"] = "E18_CDS"
# onco_data$FGFR2.truncation.location[onco_data$FGFR2.truncation.location == "Exon 18_3'-UTR"] = "E18_UTR"

# FGFR2 RE type
onco_data_new$RE.type[which(onco_data_new$RE.type == "Out-of-strand RE")] = "out-of-strand RE"
onco_data_new$RE.type[which(onco_data_new$RE.type == "RE without Partner")] = "intergenic RE"
onco_data_new$RE.type[which(onco_data_new$RE.type == "Internal RE")] = "re_internal"

t = rep(NA, nrow(onco_data_new))
t[grepl("5", onco_data_new$RE.location.to.FGFR2)] = "FGFR2 is downstream"
t[grepl("3", onco_data_new$RE.location.to.FGFR2)] = "FGFR2 is upstream"
onco_data_new$RE.location.to.FGFR2 = t

onco_data_new = data.frame(Specimen = onco_data_new$Specimen.name,
                           TCGA_type = onco_data_new$TCGA_Oncoplot,
                           FGFR2_mut_hot = rep(NA, nrow(onco_data_new)), 
                           FGFR2_mut_trunc = FGFR2_mut_trunc,
                           FGFR2_amp = FGFR2_amp,
                           RE_chr1 = onco_data_new$Chr.1_RE,
                           RE_pos1 = onco_data_new$Position.1_RE,
                           RE_chr2 = onco_data_new$Chr.2_RE,
                           RE_pos2 = onco_data_new$Position.2_RE,
                           RE_partner = onco_data_new$Rearrangement..RE.,
                           RE_FGFR2 = onco_data_new$RE.type,
                           RE_location = onco_data_new$RE.location.to.FGFR2,
                           FGFR2_brkpt = onco_data_new$FGFR2.truncation.location,
                           FGFR2_chr_SV = FGFR2_re_loc) %>%
    mutate(RE_FGFR2 = factor(RE_FGFR2, levels = c("in-frame", "in-strand (frame-unknown)", "intergenic RE", "out-of-strand RE", "re_internal")),
           FGFR2_brkpt = factor(FGFR2_brkpt, levels = c("I17", "E18", "5' to E17")),
           FGFR2_chr_SV = factor(FGFR2_chr_SV, levels = c("Intrachromosomal", "Interchromosomal")))


df_oncoplot = onco_data_new[, c("TCGA_type", "FGFR2_mut_hot", "FGFR2_mut_trunc", "FGFR2_amp", "RE_FGFR2", "RE_location", "FGFR2_brkpt", "FGFR2_chr_SV")] 
rownames(df_oncoplot) = onco_data$Specimen.name

######################################################################################################
# combining the tumors with only hotspot mutations
######################################################################################################
t_samp = setdiff(FGFR2_mut_hotspot$Sample_ID, onco_data$Specimen.name)
onco_data_mut = data.frame(Specimen.name = t_samp, TCGA_Oncoplot = FGFR2_mut_all$TCGA_code[match(t_samp, FGFR2_mut_all$Sample_ID)])

t_samp_mut_df = data.frame(onco_data_mut$TCGA_Oncoplot, matrix(NA, nrow(onco_data_mut), ncol(df_oncoplot)-1))
colnames(t_samp_mut_df) = colnames(df_oncoplot)

df_oncoplot = rbind(df_oncoplot, t_samp_mut_df)
rownames(df_oncoplot) = c(onco_data$Specimen.name, onco_data_mut$Specimen.name)

# correcting FGFR2 mutations
df_oncoplot = df_oncoplot %>% 
    mutat(FGFR2_mut_hot = replace(FGFR2_mut_hot, rownames(.) %in% FGFR2_mut_hotspot$Sample_ID, TRUE))

# TCGA tumor type
t = table(df_oncoplot$TCGA_type)
df_oncoplot$TCGA_type[which(df_oncoplot$TCGA_type %in% names(t)[t<15])] = "Other"


# filtering out the 5'-E17 rearrangement
# ******* for these samples, co-altering mutation list was requested *******
df_oncoplot_2 = df_oncoplot %>% filter(FGFR2_brkpt!="5' to E17" | is.na(FGFR2_brkpt), RE_location == "FGFR2 is upstream" | is.na(RE_location)) %>%
  mutate(FGFR2_brkpt = factor(FGFR2_brkpt, levels = c("I17", "E18"))) %>%
  mutate(Order = 1:nrow(.))


# target samples whose co-mutations are to be analysed.
write.table(df_oncoplot_2, file = "~/FGFR/Daniel/R/Nature_figures/data/FM/FM_tumors.txt", col.names = T, row.names = T, sep = "\t")

save.image("~/FGFR/Daniel/R/Nature_figures/data/FM/FM_FGFR2_vars.RData")
