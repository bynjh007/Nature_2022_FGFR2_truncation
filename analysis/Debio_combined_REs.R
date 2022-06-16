library(stringr)
library(dplyr)

################################################################
# FGFR rearrangement from Debio RNA-seq
################################################################
load("~/FGFR/Daniel/R/Nature_figures/data/Debio/Debio_FGFR2_RE_processing.RData")

################################################################
# FGFR rearrangement from CrownBio atabase
################################################################
gene_fusions = read.csv2("~/FGFR/Daniel/R/Nature_figures/data/Debio/CrownBio_gene fusion.csv", header = T, sep = ",",
                        stringsAsFactors = F)

################################################################
# Debio samples having FGFRi response
################################################################
samp_resp = xlsx::read.xlsx2("~/FGFR/Daniel/R/Nature_figures/data/Debio/Codes for PDX models Debio 1347_20211208.xlsx", sheetIndex = 1, header = T)
samp_info_RNA = read.table("~/FGFR/Daniel/R/Nature_figures/data/Debio/samples.csv", header = T, stringsAsFactors = F, sep = "\t")
# overlapping samples between drug response and RNA-seq
ovp_samp  = intersect(unique(substr(samp_info_RNA$SAMPLE_ID, 1, 6)), samp_resp$Model.name)

################################################################
# FGFR2 - Combine STAR fusion and Chimeric alignment output
################################################################
# all the breakpoint from STAR-Fusion is ref:in-strand type
t_right = str_split_fixed(FGFR2_star_fusion$right_annot, ":", 5)
t_right[,3] = paste("I", as.numeric(str_split_fixed(t_right[,3], " ", 2)[,2])-1, sep = "")
t_right = apply(t_right, 1, paste, collapse = ":")
a = FGFR2_star_fusion %>% mutate(Source = "star_fusion",
                                 RE_location = ifelse(startsWith(X.FusionName, "FGFR2"), "upstream", "downstream"),
                                 RE_type = ifelse(PROT_FUSION_TYPE == "INFRAME", "inframe", "out-of-frame"),
                                 FGFR2_brkpt = ifelse(grepl("Exon 17:ref:in-strand", left_annot), "I17", "out of I17"),
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
                                         paste("I", as.numeric(str_split_fixed(down.id, "E", 2)[,2])-1, sep = ""), down.id),
                        Fusion_name = paste(up.symbol, down.symbol, sep = "--"),
                        LeftBreakpoint = up.brkpt, RightBreakpoint = down.brkpt, Junc_count = count,
                        FFPM = FFPM, Sample_id = sample_id,
                        RE_location = ifelse(startsWith(Fusion_name, "FGFR2"), "upstream", "downstream"),
                        FGFR2_brkpt = ifelse(startsWith(Fusion_name, "FGFR2"), up.id, down.id),
                        FGFR2_brkpt = replace(FGFR2_brkpt, FGFR2_brkpt != "I17", "out of I17"),
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


################################################################
# FGFR3 - Combine STAR fusion and Chimeric alignment output
################################################################
# all the breakpoint from STAR-Fusion is ref:in-strand type
t_right = str_split_fixed(FGFR3_star_fusion$right_annot, ":", 5)
t_right[,3] = paste("I", as.numeric(str_split_fixed(t_right[,3], " ", 2)[,2])-1, sep = "")
t_right = apply(t_right, 1, paste, collapse = ":")
a = FGFR3_star_fusion %>% mutate(Source = "star_fusion",
                                 RE_location = ifelse(startsWith(X.FusionName, "FGFR3"), "upstream", "downstream"),
                                 RE_type = ifelse(PROT_FUSION_TYPE == "INFRAME", "inframe", "out-of-frame"),
                                 FGFR3_brkpt = ifelse(grepl("Exon 17:ref:in-strand", left_annot), "I17", "out of I17"),
                                 C_trunc = LeftBreakpoint %in% "chr4:1806934:+",
                                 left_annot = gsub("Exon ", "I", left_annot),
                                 right_annot = t_right) %>%
  select(X.FusionName, LeftBreakpoint, RightBreakpoint, JunctionReadCount, FFPM, sample_id, RE_type, RE_location, 
         C_trunc, FGFR3_brkpt, left_annot, right_annot, Source)
colnames(a) = c("Fusion_name", "LeftBreakpoint", "RightBreakpoint", "Junc_count", "FFPM", "Sample_id", "RE_type", 
                "RE_location", "C_trunc", "FGFR3_brkpt", "left_annot", "right_annot", "Source")

# extract fusion information from RE_from_chimeric (** note that brkpt from RE_from_chimeric is one base different with star-fusion results **)
b = FGFR3_RE %>% mutate(Source = "RE_from_chimeric",
                        up.id = as.character(up.id), up.id = gsub("Exon ", "E", up.id), up.id = gsub("Intron ", "I", up.id),
                        up.id = ifelse(up.splice_type == "ref:in-strand", gsub("E", "I", up.id), up.id),
                        down.id = as.character(down.id), down.id = gsub("Exon ", "E", down.id), down.id = gsub("Intron ", "I", down.id),
                        down.id = ifelse(down.splice_type == "ref:in-strand", 
                                         paste("I", as.numeric(str_split_fixed(down.id, "E", 2)[,2])-1, sep = ""), down.id), 
                        Fusion_name = paste(up.symbol, down.symbol, sep = "--"),
                        LeftBreakpoint = up.brkpt, RightBreakpoint = down.brkpt, Junc_count = count,
                        FFPM = FFPM, Sample_id = sample_id,
                        RE_location = ifelse(startsWith(Fusion_name, "FGFR3"), "upstream", "downstream"),
                        FGFR3_brkpt = ifelse(startsWith(Fusion_name, "FGFR3"), up.id, down.id),
                        FGFR3_brkpt = replace(FGFR3_brkpt, FGFR3_brkpt != "I17", "out of I17"),
                        C_trunc_can = ((up.brkpt %in% "chr4:1806934:+") & (up.splice_type %in% "ref:in-strand") & (up.symbol %in% "FGFR3")),
                        C_trunc_noncan = ((FGFR3_brkpt %in% "I17") & (up.splice_type %in% "non_ref:in-strand") & (up.symbol %in% "FGFR3")),
                        C_trunc = ((FGFR3_brkpt %in% "I17") & grepl("in-strand", up.splice_type) & (up.symbol %in% "FGFR3")),
                        t_type = paste(up.splice_type, down.splice_type, sep = "//"),
                        RE_type = rep(NA, nrow(FGFR3_RE)),
                        RE_type = replace(RE_type, grepl("in-strand", t_type) & grepl("out-of-strand", t_type), "out-of-strand"),
                        RE_type = replace(RE_type, grepl("intergenic", t_type), "intergenic"),
                        RE_type = replace(RE_type, grepl("in-strand", up.splice_type) & grepl("in-strand", down.splice_type), "out-of-frame"),
                        left_annot = paste(up.symbol, up.txname, up.id, up.splice_type, sep = ":"),
                        right_annot = paste(down.symbol, down.txname, down.id, down.splice_type, sep = ":")) %>%
  select(Fusion_name, LeftBreakpoint, RightBreakpoint, Junc_count, FFPM, Sample_id, RE_type, RE_location, 
         C_trunc, FGFR3_brkpt, left_annot, right_annot, Source)
                        
# combine the RE and Fusions
FGFR3_RE_all = rbind(a,b)

# remove out-of-strand fusion from upstream
FGFR3_RE_all = FGFR3_RE_all %>% filter(grepl("in-strand", left_annot)) %>%
  mutate(RE_type = factor(RE_type, levels = c("inframe", "out-of-frame", "intergenic", "out-of-strand"))) %>%
  filter(FFPM>=0.05)


################################################################
# Fusions from CrownBio
################################################################
ensembl_37 = biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh=37)
FGFR1_exon_coord_37 = biomaRt::getBM(attributes = c("ensembl_exon_id", "exon_chrom_start", 
                                        "exon_chrom_end", "strand", "transcript_biotype"), 
                         filters = "hgnc_symbol", values = "FGFR1", mart = ensembl_37)
FGFR1_exon_coord_37 = FGFR1_exon_coord_37 %>% filter(grepl("protein_coding|processed_transcript", transcript_biotype))

FGFR2_exon_coord_37 = biomaRt::getBM(attributes = c("ensembl_exon_id", "exon_chrom_start", 
                                        "exon_chrom_end", "strand", "transcript_biotype"), 
                         filters = "hgnc_symbol", values = "FGFR2", mart = ensembl_37)
FGFR2_exon_coord_37 = FGFR2_exon_coord_37 %>% filter(grepl("protein_coding|processed_transcript", transcript_biotype))

FGFR3_exon_coord_37 = biomaRt::getBM(attributes = c("ensembl_exon_id", "exon_chrom_start", 
                                        "exon_chrom_end", "strand", "transcript_biotype"), 
                         filters = "hgnc_symbol", values = "FGFR3", mart = ensembl_37)
FGFR3_exon_coord_37 = FGFR3_exon_coord_37 %>% filter(grepl("protein_coding|processed_transcript", transcript_biotype))

# functions
samp_fus = function(target){
  t_FGFR2 = gene_fusions %>% mutate(Sample_Name_2 = substr(Sample_Name, 1, 6)) %>% 
    filter(Sample_Name_2 %in% ovp_samp, Up_gene == target | Dw_gene == target)
}

RE_annot_FGFR_Crownbio = function(t_df, target, target_exon_coord, range_int17){
  
  target_str = target_exon_coord$strand[1]
  if(target_str == -1){
    exon_end = 2; exon_start = 3 
  } else {
    exon_end = 3; exon_start = 2 
  }
  
  # fusion location
  if(t_df$Up_gene == target){
    t_target_brkpt = t_df$Up_genome_pos
    t_partner = t_df[,c("Dw_genome_pos", "Dw_gene", "Dw_strand")]
    t_type = "upstream"
    t_cano_FGFR = t_target_brkpt %in% target_exon_coord[,exon_end]
  } else {
    t_target_brkpt = t_df$Dw_genome_pos
    t_partner = t_df[,c("Up_genome_pos", "Up_gene", "Up_strand")]
    t_type = "downstream"
    t_cano_FGFR = t_target_brkpt %in% target_exon_coord[,exon_start]
  }
  colnames(t_partner) = c("brkpt", "gene", "strand")
  
  # fusion location
  if(t_target_brkpt>=range_int17[1] & t_target_brkpt<=range_int17[2]){
    brk_I17 = "Yes"
  } else {
    brk_I17 = NA
  }
  
  if(target_str == -1){
    cano_I17_loc = range_int17[2]
  } else {
    cano_I17_loc = range_int17[1]
  }
  
  if(t_target_brkpt == cano_I17_loc & t_df$Up_gene == target){
    cano_I17_truc = "Yes"
  } else {
    cano_I17_truc = NA
  }
  
  # in-frame or out-of-frame
  partner_exon = biomaRt::getBM(attributes = c("ensembl_exon_id", "exon_chrom_start", 
                                               "exon_chrom_end", "strand", "transcript_biotype"),
                                filters = "hgnc_symbol", values = as.character(t_partner$gene), mart = ensembl_37)
  partner_exon = partner_exon %>% filter(grepl("protein_coding|processed_transcript", transcript_biotype))
  
  if(t_df$Up_gene == target){ # target is upstream of a fusion
    if(all(partner_exon$strand == 1)){
      t_cano_parner = t_partner$brkpt %in% partner_exon$exon_chrom_start
    } else if(all(partner_exon$strand == -1)){
      t_cano_parner = t_partner$brkpt %in% partner_exon$exon_chrom_end
    } else {
      stop("not consistent across the exon's strand")
    }
  } else { # target is downstream of a fusion
    if(all(partner_exon$strand == 1)){
      t_cano_parner = t_partner$brkpt %in% partner_exon$exon_chrom_end
    } else if(all(partner_exon$strand == -1)){
      t_cano_parner = t_partner$brkpt %in% partner_exon$exon_chrom_start
    } else {
      stop("not consistent across the exon's strand")
    }
  }
  
  t_cano = c(t_cano_FGFR, t_cano_parner)
  if(sum(t_cano) == 2){
    cano = "inframe"
  } else {
    cano = "out-of-frame"
  }
  
  fus_name = paste(t_df$Up_gene, "--", t_df$Dw_gene, sep = "")
  fus_info = c(fus_name, cano, gsub("-", "", t_df$Sample_Name), t_type, brk_I17, cano_I17_truc)
  n_spanning = suppressWarnings(sum(as.numeric(t_df[,12:15]), na.rm = T))
  out = list(fus_info = fus_info, n_span = n_spanning)
  return(out)
}

# Function for RE refinement
RE_final = function(FGFR_name){
    Crown_FGFR = samp_fus(FGFR_name)
    t_samp = unique(Crown_FGFR$Sample_Name_2)

    if(FGFR_name == "FGFR1"){
        FGFR_exon_coord_37 = FGFR1_exon_coord_37; range_int17 = c(38271322, 38271436)
    } else if(FGFR_name == "FGFR2"){
        FGFR_exon_coord_37 = FGFR2_exon_coord_37; range_int17 = c(123239535, 123243212)
    } else if(FGFR_name == "FGFR3"){
        FGFR_exon_coord_37 = FGFR3_exon_coord_37; range_int17 = c(1808661, 1808843)
    } else {
        break;
    }

    FGFR_Crownbio = c()
    for(i in 1:length(t_samp)){
        t_fgfr = Crown_FGFR %>% filter(Sample_Name_2 == t_samp[i])
        t_annot = c()
        for(j in 1:nrow(t_fgfr)){
            t = RE_annot_FGFR_Crownbio(t_fgfr[j,], FGFR_name, FGFR_exon_coord_37, range_int17)
            t_annot = rbind(t_annot, t$fus_info)
        }
        t_annot = as.data.frame(t_annot)
        colnames(t_annot) = c("X.FusionName", "PROT_FUSION_TYPE", "sample_id", "type", "I17", "canonical_I17")

        # selecting the representative one
        t_annot = t_annot %>% arrange(desc(type), I17, PROT_FUSION_TYPE)
        FGFR_Crownbio = rbind(FGFR_Crownbio, t_annot[1,])
    }
    return(FGFR_Crownbio)
}

FGFR1_Crownbio = RE_final("FGFR1")
FGFR2_Crownbio = RE_final("FGFR2")
FGFR3_Crownbio = RE_final("FGFR3")


################################################################
# Combining REs from CrownBio and our analysis
################################################################
# combining FGFR1 REs
t_all = FGFR1_Crownbio %>% select(X.FusionName, PROT_FUSION_TYPE, sample_id, type, I17)
t_all$sample_id = substr(t_all$sample_id, 1, 6)

t_samp = as.character(unique(t_all$sample_id))
FGFR1_RE_comb = c()
for(i in t_samp){
  t_df = t_all %>% filter(sample_id == i) %>% arrange(desc(type), I17, PROT_FUSION_TYPE)
  FGFR1_RE_comb = rbind(FGFR1_RE_comb, t_df[1,])
}

# combining FGFR2 REs
t_star = FGFR2_RE_all %>% arrange(desc(C_trunc), desc(RE_location), FGFR2_brkpt, RE_type) %>% 
  mutate(I17 = ifelse(FGFR2_brkpt == "I17", "Yes", NA)) %>%
  mutate(Sample_id = substr(Sample_id, 1, 6)) %>%
  select(Fusion_name, RE_type, Sample_id, RE_location, I17)
t_samp = unique(t_star$Sample_id)
FGFR2_RE_rep = c()
for(i in 1:length(t_samp)){
  t_star_2 = t_star %>% filter(Sample_id == t_samp[i])
  FGFR2_RE_rep = rbind(FGFR2_RE_rep, t_star_2[1,])
}
colnames(FGFR2_RE_rep) = c("X.FusionName", "PROT_FUSION_TYPE", "sample_id", "type", "I17")

t_all = rbind(FGFR2_Crownbio %>% select(X.FusionName, PROT_FUSION_TYPE, sample_id, type, I17),
              FGFR2_RE_rep)
t_all$sample_id = substr(t_all$sample_id, 1, 6)

t_samp = as.character(unique(t_all$sample_id))
FGFR2_RE_comb = c()
for(i in t_samp){
  t_df = t_all %>% filter(sample_id == i) %>% arrange(desc(type), I17, PROT_FUSION_TYPE)
  FGFR2_RE_comb = rbind(FGFR2_RE_comb, t_df[1,])
}

# combining FGFR3 REs
t_star = FGFR3_RE_all %>% arrange(desc(C_trunc), desc(RE_location), FGFR3_brkpt, RE_type) %>% 
  mutate(I17 = ifelse(FGFR3_brkpt == "I17", "Yes", NA)) %>%
  mutate(Sample_id = substr(Sample_id, 1, 6)) %>%
  select(Fusion_name, RE_type, Sample_id, RE_location, I17)
t_samp = unique(t_star$Sample_id)
FGFR3_RE_rep = c()
for(i in 1:length(t_samp)){
  t_star_2 = t_star %>% filter(Sample_id == t_samp[i])
  FGFR3_RE_rep = rbind(FGFR3_RE_rep, t_star_2[1,])
}
colnames(FGFR3_RE_rep) = c("X.FusionName", "PROT_FUSION_TYPE", "sample_id", "type", "I17")

t_all = rbind(FGFR3_Crownbio %>% select(X.FusionName, PROT_FUSION_TYPE, sample_id, type, I17),
              FGFR3_RE_rep)
t_all$sample_id = substr(t_all$sample_id, 1, 6)

t_samp = as.character(unique(t_all$sample_id))
FGFR3_RE_comb = c()
for(i in t_samp){
  t_df = t_all %>% filter(sample_id == i) %>% arrange(desc(type), I17, PROT_FUSION_TYPE)
  FGFR3_RE_comb = rbind(FGFR3_RE_comb, t_df[1,])
}

save.image("~/FGFR/Daniel/R/Nature_figures/data/Debio/Debio_combined_REs.RData")

