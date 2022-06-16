#############################################
# Libraries & functions
#############################################
library(rstatix)
library(stringr)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(patchwork)

split_cell_names = function(X){str_split_fixed(X, "_", 2)[,1]}

#############################################
# Drug response data (GDSC and PharmacoDB)
#############################################
# AZD4547 response for CCLE cell lines obtained from PharmacoDB
cell_models_PDB = read.table("~/FGFR/Daniel/R/nvim_R/data/PharmacoDB_cell_annot.csv", sep = ",", header = T, stringsAsFactors = F)

drug_resp_PDB = read.table("~/FGFR/Daniel/R/nvim_R/data/PharmacoDB_drug_response.csv", sep = ",", header = T, stringsAsFactors = F) %>%
    mutate(AUC = 1 - AAC)
drug_resp_PDB = data.frame(drug_resp_PDB,
                           cell_models_PDB[match(drug_resp_PDB$cell.name, cell_models_PDB$unique.cellid),
                                           c("unique.cellid", "CCLE.cellid", "CTRPv2.cellid", "GDSC1000.cellid")]) %>% 
                arrange(desc(AUC))

AZD4547_PDB = drug_resp_PDB %>% filter(drugname == "AZD4547")
AZD4547_PDB = AZD4547_PDB[match(unique(AZD4547_PDB$cell.name), AZD4547_PDB$cell.name),] %>% arrange(AUC)

# GDSC data loading
# for PD173074, GDSC data was used because of notorious inconsistency between CCLE ID and PharmacoDB ID for PD173074
PD173074_GDSC = read.table("~/FGFR/Daniel/R/nvim_R/data/PD173074_GDSC2.csv", sep = ",", header = T, stringsAsFactors = F)
PD173074_GDSC = PD173074_GDSC %>% mutate(label = gsub("-", "", Cell.line))


#############################################
# CCLE data loading
#############################################
# CCLE cell line
cell_models_CCLE = read.csv2("~/FGFR/Daniel/R/Nature_figures/data/CCLE_GDSC_PharmacoDB/sample_info.csv", 
                             stringsAsFactors = F, header = T, sep = ",")
cell_models_CCLE$Cell_name = cell_models_CCLE$stripped_cell_line_name

# CCLE FGFR2 exon usage
FGFR2_exon_usage_CCLE = t(read.table("~/FGFR/Daniel/R/Nature_figures/data/CCLE_GDSC_PharmacoDB/FGFR2_ExonUsageRatio.txt", header = F, stringsAsFactors = F, sep = "\t"))
colnames(FGFR2_exon_usage_CCLE) = FGFR2_exon_usage_CCLE[1,]
FGFR2_exon_usage_CCLE = FGFR2_exon_usage_CCLE[-c(1,2),]
rownames(FGFR2_exon_usage_CCLE) = FGFR2_exon_usage_CCLE[,1]
FGFR2_exon_usage_CCLE = FGFR2_exon_usage_CCLE[,-1]
class(FGFR2_exon_usage_CCLE) = "numeric"
FGFR2_exon_usage_CCLE = FGFR2_exon_usage_CCLE[,order(colnames(FGFR2_exon_usage_CCLE))]

# CNV data
cnv_CCLE = read.table("~/FGFR/Daniel/R/Nature_figures/data/CCLE_GDSC_PharmacoDB/CCLE_copynumber_byGene_2013-12-03.txt", header = T, sep = "\t", stringsAsFactors = F)
rownames(cnv_CCLE) = cnv_CCLE[,2]; cnv_CCLE = as.matrix(cnv_CCLE[, -c(1:5)])
class(cnv_CCLE) = "numeric"

# SV data
sv_CCLE = read.table("~/FGFR/Daniel/R/Nature_figures/data/CCLE_GDSC_PharmacoDB/CCLE_SVs.txt", header = T, sep = "\t", stringsAsFactors = F)
sv_FGFR2_CCLE = sv_CCLE[unique(c(grep("FGFR2", sv_CCLE$gene1), grep("FGFR2", sv_CCLE$gene2))),]
sv_FGFR1_CCLE = sv_CCLE[unique(c(grep("FGFR1", sv_CCLE$gene1), grep("FGFR1", sv_CCLE$gene2))),]
sv_FGFR3_CCLE = sv_CCLE[unique(c(grep("FGFR3", sv_CCLE$gene1), grep("FGFR3", sv_CCLE$gene2))),]
sv_FGFR4_CCLE = sv_CCLE[unique(c(grep("FGFR4", sv_CCLE$gene1), grep("FGFR4", sv_CCLE$gene2))),]

# RNA-seq data
t = read.csv2("~/FGFR/Daniel/R/Nature_figures/data/CCLE_GDSC_PharmacoDB/FGFR2_CCLE_RNAseq_genes_rpkm_20180929.txt", header = T, sep = "\t", stringsAsFactors = F, check.names = F)
exp_FGFR2_CCLE = as.numeric(t(t[,-c(1,2)]))
names(exp_FGFR2_CCLE) = colnames(t)[-c(1,2)]
t = read.csv2("~/FGFR/Daniel/R/Nature_figures/data/CCLE_GDSC_PharmacoDB/FGFR1_CCLE_RNAseq_genes_rpkm_20180929.txt", header = F, sep = "\t", stringsAsFactors = F, check.names = F)
exp_FGFR1_CCLE = as.numeric(t(t[,-c(1,2)]))
names(exp_FGFR1_CCLE) = names(exp_FGFR2_CCLE)
t = read.csv2("~/FGFR/Daniel/R/Nature_figures/data/CCLE_GDSC_PharmacoDB/FGFR3_CCLE_RNAseq_genes_rpkm_20180929.txt", header = F, sep = "\t", stringsAsFactors = F, check.names = F)
exp_FGFR3_CCLE = as.numeric(t(t[,-c(1,2)]))
names(exp_FGFR3_CCLE) = names(exp_FGFR2_CCLE)
t = read.csv2("~/FGFR/Daniel/R/Nature_figures/data/CCLE_GDSC_PharmacoDB/FGFR4_CCLE_RNAseq_genes_rpkm_20180929.txt", header = F, sep = "\t", stringsAsFactors = F, check.names = F)
exp_FGFR4_CCLE = as.numeric(t(t[,-c(1,2)]))
names(exp_FGFR4_CCLE) = names(exp_FGFR2_CCLE)

gmean_norm = function(x){(x-exp(mean(log(x))))/sd(x)} 
exp_FGFR_comp = colSums(rbind(gmean_norm(exp_FGFR1_CCLE+1), 
                              gmean_norm(exp_FGFR2_CCLE+1),
                              gmean_norm(exp_FGFR3_CCLE+1),
                              gmean_norm(exp_FGFR4_CCLE+1)))
names(exp_FGFR_comp) = names(exp_FGFR2_CCLE)

# Mutation
mut_CCLE = read.table("~/FGFR/Daniel/R/Nature_figures/data/CCLE_GDSC_PharmacoDB/CCLE_DepMap_18q3_maf_20180718.txt", header = T, sep = "\t", stringsAsFactors = F)
FGFR2_mut_CCLE = mut_CCLE %>% filter(Hugo_Symbol=="FGFR2" & Variant_Classification != "Silent")
FGFR1_mut_CCLE = mut_CCLE %>% filter(Hugo_Symbol=="FGFR1" & Variant_Classification != "Silent")
FGFR3_mut_CCLE = mut_CCLE %>% filter(Hugo_Symbol=="FGFR3" & Variant_Classification != "Silent")
FGFR4_mut_CCLE = mut_CCLE %>% filter(Hugo_Symbol=="FGFR4" & Variant_Classification != "Silent")
FGFR2_mut_CCLE_hotspot = FGFR2_mut_CCLE %>% filter(grepl("S252|C382|N549|K659", Protein_Change))
FGFR3_mut_CCLE_hotspot = FGFR3_mut_CCLE %>% filter(grepl("R248|S249|Y373|K650", Protein_Change))

# Fusion
fusions_CCLE = read.table("~/FGFR/Daniel/R/Nature_figures/data/CCLE_GDSC_PharmacoDB/CCLE_Fusions_20181130.txt", header = T, stringsAsFactors = F, sep = "\t", quote = "")
fusions_unfilt_CCLE = read.table("~/FGFR/Daniel/R/Nature_figures/data/CCLE_GDSC_PharmacoDB/CCLE_Fusions_unfiltered_20181130.txt", header = T, stringsAsFactors = F, sep = "\t", quote = "")


#############################################
# Fusion annotation
#############################################
fusion_annotation = function(FGFR){
    if(FGFR == "FGFR2"){
        I17_loc = c(123239535, 123243212); exp_FGFR_CCLE = exp_FGFR2_CCLE
    } else if(FGFR == "FGFR1"){
        I17_loc = c(38271322, 38271436); exp_FGFR_CCLE = exp_FGFR1_CCLE
    } else if(FGFR == "FGFR3"){
        I17_loc = c(1808661, 1808843); exp_FGFR_CCLE = exp_FGFR3_CCLE
    } else {
        break;
    }

    FGFR_fusion_CCLE = df_fusion %>% 
        mutate(Expression = log2(exp_FGFR_CCLE[match(X.sample, names(exp_FGFR_CCLE))]+1)) %>%
        filter(LeftGene == FGFR | RightGene == FGFR, sumFFPM>=0.1, SpanningFragCount>=5, Expression >=1,
           !grepl("chrM", LeftBreakpoint), !grepl("chrM", RightBreakpoint)) %>%
        mutate(type = rep(NA, nrow(.)),
               type = replace(type, startsWith(X.FusionName, paste(FGFR, "--", sep = "")), "upstream"),
               type = replace(type, !startsWith(X.FusionName, paste(FGFR, "--", sep = "")), "downstream"))
    
    # annotating the fusions
    t = rep(NA, nrow(FGFR_fusion_CCLE))
    ind_l = which(FGFR_fusion_CCLE$type == "upstream")
    ind_r = which(FGFR_fusion_CCLE$type == "downstream")
    t[ind_l] = as.numeric(str_split_fixed(FGFR_fusion_CCLE$LeftBreakpoint[ind_l], ":", 3)[,2])
    t[ind_r] = as.numeric(str_split_fixed(FGFR_fusion_CCLE$RightBreakpoint[ind_r], ":", 3)[,2])
    t[t>=I17_loc[1] & t<=I17_loc[2]] = "I17"
    t[t!="I17"] = "5'-E17"
    FGFR_fusion_CCLE = FGFR_fusion_CCLE %>% mutate(locus = t)
    
    return(FGFR_fusion_CCLE)
}


df_fusion = fusions_unfilt_CCLE %>% mutate(cellid = split_cell_names(clnam))
# FGFR1
FGFR1_fusion_CCLE = fusion_annotation("FGFR1")
FGFR2_fusion_CCLE = fusion_annotation("FGFR2")
FGFR3_fusion_CCLE = fusion_annotation("FGFR3")


##########################################################################################
# representative fusions for each cell line & manual curation to check in- or out-of-frame
##########################################################################################
t_uniq = unique(FGFR2_fusion_CCLE$X.sample)
FGFR2_fusion_rep = c()
for(i in t_uniq){
  t = FGFR2_fusion_CCLE %>% filter(X.sample == i) %>% arrange(desc(type), desc(locus), desc(SpliceType))
  FGFR2_fusion_rep = rbind(FGFR2_fusion_rep, t[1,])
}
t = c("SNU16_STOMACH" = "out-of-frame", "NCIH2066_LUNG" = "inframe", "NCIH716_LARGE_INTESTINE" = "inframe",
      "KATOIII_STOMACH" = "out-of-frame", "HCC1187_BREAST" = "inframe")
FGFR2_fusion_rep = FGFR2_fusion_rep %>% 
    mutate(partner = t[match(X.sample, names(t))])

t_uniq = unique(FGFR1_fusion_CCLE$X.sample)
FGFR1_fusion_rep = c()
for(i in t_uniq){
  t = FGFR1_fusion_CCLE %>% filter(X.sample == i) %>% arrange(desc(type), desc(locus), desc(SpliceType))
  FGFR1_fusion_rep = rbind(FGFR1_fusion_rep, t[1,])
}
t = c("KG1_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE" = "inframe", "SNU761_LIVER" = "out-of-frame")
FGFR1_fusion_rep = FGFR1_fusion_rep %>% 
    mutate(partner = t[match(X.sample, names(t))])

t_uniq = unique(FGFR3_fusion_CCLE$X.sample)
FGFR3_fusion_rep = c()
for(i in t_uniq){
  t = FGFR3_fusion_CCLE %>% filter(X.sample == i) %>% arrange(desc(type), desc(locus), desc(SpliceType))
  FGFR3_fusion_rep = rbind(FGFR3_fusion_rep, t[1,])
}
t = c("RT4_URINARY_TRACT" = "out-of-frame", "SW780_URINARY_TRACT" = "inframe", 
      "RT112_URINARY_TRACT" = "inframe", "SKNMC_BONE" = "inframe")
FGFR3_fusion_rep = FGFR3_fusion_rep %>% 
    mutate(partner = t[match(X.sample, names(t))])


################################################
# FGFR2 C3 usage
################################################
# chr10_123239535_123237848 : C1
# chr10_123239627_123239565 : C2
# chr10_123239999_123239740 : unknown (intronic)
# chr10_123241691_123241349 : C3
# chr10_123243317_123243212 : E17

# C3 alternative exon usage
ind_C1 = grep("chr10_123239535_123237848", colnames(FGFR2_exon_usage_CCLE))[2]
ind_C2 = grep("chr10_123239627_123239565", colnames(FGFR2_exon_usage_CCLE))[2]
ind_C3 = grep("chr10_123241691_123241349", colnames(FGFR2_exon_usage_CCLE))[2]
df_sp = data.frame(CCLE_ID = rownames(FGFR2_exon_usage_CCLE),
                   C1 = FGFR2_exon_usage_CCLE[,ind_C1],
                   C2 = FGFR2_exon_usage_CCLE[,ind_C2],
                   C3 = FGFR2_exon_usage_CCLE[,ind_C3]) %>%
  mutate(FGFR2_exp = exp_FGFR2_CCLE[match(CCLE_ID, names(exp_FGFR2_CCLE))]) %>%
  mutate(AUC = AZD4547_PDB$AUC[match(split_cell_names(CCLE_ID), AZD4547_PDB$CTRPv2.cellid)]) %>%
  filter(!is.na(AUC))

# C3 high sample (log2 (for normality) --> Z-score --> one-side P-value (0.01))
# most of the cell lines lowly express FGFR2 (median 0.3 RPKM), so select the samples with top 30% (2.14 RPKM)
t = FGFR2_exon_usage_CCLE[,ind_C3]
#t[t==0] = NA
t[is.na(t)] = 0
t2 = cbind(exp_FGFR2_CCLE[match(rownames(FGFR2_exon_usage_CCLE), names(exp_FGFR2_CCLE))], 
           log2(exp_FGFR2_CCLE[match(rownames(FGFR2_exon_usage_CCLE), names(exp_FGFR2_CCLE))]+1),
           t, scale(log2(t)), pnorm(-scale(log2(t))))
C3_high_ccle = names(which(t2[,1] >= quantile(exp_FGFR2_CCLE, probs = c(0.7)) & t2[,3]>0.5))


################################################
# FGF/FGFR amplification
################################################
# FGFR1/2/3 amplification
ind_fgfr = c(which(rownames(cnv_CCLE) == "FGFR1"),
             which(rownames(cnv_CCLE) == "FGFR2"),
             which(rownames(cnv_CCLE) == "FGFR3"),
             which(rownames(cnv_CCLE) == "FGFR4"))
FGFR_cnv = t(cnv_CCLE[ind_fgfr,])
FGFR_cnv = FGFR_cnv >= 2

FGFR1_amp =  names(which(FGFR_cnv[,1] == T))
FGFR2_amp =  names(which(FGFR_cnv[,2] == T))
FGFR3_amp =  names(which(FGFR_cnv[,3] == T))
FGFR4_amp =  names(which(FGFR_cnv[,4] == T))

# FGF3/4/19 amplification
ind_fgf = c(which(rownames(cnv_CCLE) == "FGF3"),
        which(rownames(cnv_CCLE) == "FGF4"),
        which(rownames(cnv_CCLE) == "FGF19"))
FGF_cnv = t(cnv_CCLE[ind_fgf,])
FGF_cnv = FGF_cnv >= 2

FGF_amp =  names(which(rowSums(FGF_cnv*1)==3))


#############################################
# mapping drug response data to each cell line
#############################################
# AZD4547
AZD4547_PDB_annot = AZD4547_PDB %>% 
    mutate(FGFR1_RE = apply(FGFR1_fusion_rep[match(CTRPv2.cellid, split_cell_names(FGFR1_fusion_rep$X.sample)), 
                                c("type", "locus", "partner")], 1, paste, collapse = ":"),
           FGFR1_RE = replace(FGFR1_RE, FGFR1_RE == "NA:NA:NA", "none"),
           FGFR1_RE = ifelse(CTRPv2.cellid %in% split_cell_names(names(exp_FGFR_comp)), FGFR1_RE, NA),
           FGFR2_RE = apply(FGFR2_fusion_rep[match(CTRPv2.cellid, split_cell_names(FGFR2_fusion_rep$X.sample)), 
                                c("type", "locus", "partner")], 1, paste, collapse = ":"),
           FGFR2_RE = replace(FGFR2_RE, FGFR2_RE == "NA:NA:NA", "none"),
           FGFR2_RE = ifelse(CTRPv2.cellid %in% split_cell_names(names(exp_FGFR_comp)), FGFR2_RE, NA),
           FGFR3_RE = apply(FGFR3_fusion_rep[match(CTRPv2.cellid, split_cell_names(FGFR3_fusion_rep$X.sample)), 
                                c("type", "locus", "partner")], 1, paste, collapse = ":"),
           FGFR3_RE = replace(FGFR3_RE, FGFR3_RE == "NA:NA:NA", "none"),
           FGFR3_RE = ifelse(CTRPv2.cellid %in% split_cell_names(names(exp_FGFR_comp)), FGFR3_RE, NA)) %>%
    mutate(E18_C3_usage = CTRPv2.cellid %in% split_cell_names(C3_high_ccle)) %>%
    mutate(FGF3_4_9_AMP = ifelse(CTRPv2.cellid %in% split_cell_names(FGF_amp), "FGF3/4/19_AMP", "no AMP"),
           FGF3_4_9_AMP = ifelse(CTRPv2.cellid %in% split_cell_names(colnames(cnv_CCLE)), FGF3_4_9_AMP, NA),
           FGF3_4_9_AMP = factor(FGF3_4_9_AMP, levels = c("no AMP", "FGF3/4/19_AMP"))) %>%
    mutate(FGFR_AMP = rep("no AMP", nrow(AZD4547_PDB)),
           FGFR_AMP = replace(FGFR_AMP, CTRPv2.cellid %in% split_cell_names(FGFR1_amp), "FGFR1"),
           FGFR_AMP = replace(FGFR_AMP, CTRPv2.cellid %in% split_cell_names(FGFR2_amp), "FGFR2"),
           FGFR_AMP = replace(FGFR_AMP, CTRPv2.cellid %in% split_cell_names(FGFR3_amp), "FGFR3"),
           FGFR_AMP = replace(FGFR_AMP, CTRPv2.cellid %in% split_cell_names(FGFR4_amp), "FGFR4"),
           FGFR_AMP = ifelse(CTRPv2.cellid %in% split_cell_names(colnames(cnv_CCLE)), FGFR_AMP, NA),
           FGFR_AMP = factor(FGFR_AMP, levels = c("no AMP", "FGFR1", "FGFR2", "FGFR3", "FGFR4"))) %>%
    mutate(FGFR1_exp = log2(exp_FGFR1_CCLE[match(CTRPv2.cellid, split_cell_names(names(exp_FGFR1_CCLE)))]+1),
           FGFR2_exp = log2(exp_FGFR2_CCLE[match(CTRPv2.cellid, split_cell_names(names(exp_FGFR2_CCLE)))]+1),
           FGFR3_exp = log2(exp_FGFR3_CCLE[match(CTRPv2.cellid, split_cell_names(names(exp_FGFR3_CCLE)))]+1),
           FGFR4_exp = log2(exp_FGFR4_CCLE[match(CTRPv2.cellid, split_cell_names(names(exp_FGFR4_CCLE)))]+1),
           FGFR_exp_comp = exp_FGFR_comp[match(CTRPv2.cellid, split_cell_names(names(exp_FGFR4_CCLE)))]) %>%
    mutate(FGFR2_exon.C3 = FGFR2_exon_usage_CCLE[match(CTRPv2.cellid, split_cell_names(rownames(FGFR2_exon_usage_CCLE))), ind_C3]) %>%
    mutate(FGFR_hotspot_mut = rep("none", nrow(.)),
           FGFR_hotspot_mut = replace(FGFR_hotspot_mut, CTRPv2.cellid %in% split_cell_names(FGFR2_mut_CCLE_hotspot$Tumor_Sample_Barcode), "FGFR2"),
           FGFR_hotspot_mut = replace(FGFR_hotspot_mut, CTRPv2.cellid %in% split_cell_names(FGFR3_mut_CCLE_hotspot$Tumor_Sample_Barcode), "FGFR3"),
           FGFR_hotspot_mut = ifelse(CTRPv2.cellid %in% split_cell_names(unique(mut_CCLE$Tumor_Sample_Barcode)), FGFR_hotspot_mut, NA),
           FGFR_hotspot_mut = factor(FGFR_hotspot_mut, levels = c("none", "FGFR2", "FGFR3")))

# PD173074
PD173074_GDSC_annot = PD173074_GDSC %>%
    mutate(FGFR1_RE = apply(FGFR1_fusion_rep[match(label, split_cell_names(FGFR1_fusion_rep$X.sample)), 
                                c("type", "locus", "partner")], 1, paste, collapse = ":"),
           FGFR1_RE = replace(FGFR1_RE, FGFR1_RE == "NA:NA:NA", "none"),
           FGFR1_RE = ifelse(label %in% split_cell_names(names(exp_FGFR_comp)), FGFR1_RE, NA),
           FGFR2_RE = apply(FGFR2_fusion_rep[match(label, split_cell_names(FGFR2_fusion_rep$X.sample)), 
                                c("type", "locus", "partner")], 1, paste, collapse = ":"),
           FGFR2_RE = replace(FGFR2_RE, FGFR2_RE == "NA:NA:NA", "none"),
           FGFR2_RE = ifelse(label %in% split_cell_names(names(exp_FGFR_comp)), FGFR2_RE, NA),
           FGFR3_RE = apply(FGFR3_fusion_rep[match(label, split_cell_names(FGFR3_fusion_rep$X.sample)), 
                                c("type", "locus", "partner")], 1, paste, collapse = ":"),
           FGFR3_RE = replace(FGFR3_RE, FGFR3_RE == "NA:NA:NA", "none"),
           FGFR3_RE = ifelse(label %in% split_cell_names(names(exp_FGFR_comp)), FGFR3_RE, NA)) %>% 
    mutate(E18_C3_usage = label %in% split_cell_names(C3_high_ccle)) %>%
    mutate(FGF3_4_9_AMP = ifelse(label %in% split_cell_names(FGF_amp), "FGF3/4/19_AMP", "no AMP"),
           FGF3_4_9_AMP = ifelse(label %in% split_cell_names(colnames(cnv_CCLE)), FGF3_4_9_AMP, NA),
           FGF3_4_9_AMP = factor(FGF3_4_9_AMP, levels = c("no AMP", "FGF3/4/19_AMP"))) %>%
    mutate(FGFR_AMP = rep("no AMP", nrow(.)),
           FGFR_AMP = replace(FGFR_AMP, label %in% split_cell_names(FGFR1_amp), "FGFR1"),
           FGFR_AMP = replace(FGFR_AMP, label %in% split_cell_names(FGFR2_amp), "FGFR2"),
           FGFR_AMP = replace(FGFR_AMP, label %in% split_cell_names(FGFR3_amp), "FGFR3"),
           FGFR_AMP = replace(FGFR_AMP, label %in% split_cell_names(FGFR4_amp), "FGFR4"),
           FGFR_AMP = ifelse(label %in% split_cell_names(colnames(cnv_CCLE)), FGFR_AMP, NA),
           FGFR_AMP = factor(FGFR_AMP, levels = c("no AMP", "FGFR1", "FGFR2", "FGFR3", "FGFR4"))) %>%
    mutate(FGFR1_exp = log2(exp_FGFR1_CCLE[match(label, split_cell_names(names(exp_FGFR1_CCLE)))]+1),
           FGFR2_exp = log2(exp_FGFR2_CCLE[match(label, split_cell_names(names(exp_FGFR2_CCLE)))]+1),
           FGFR3_exp = log2(exp_FGFR3_CCLE[match(label, split_cell_names(names(exp_FGFR3_CCLE)))]+1),
           FGFR4_exp = log2(exp_FGFR4_CCLE[match(label, split_cell_names(names(exp_FGFR4_CCLE)))]+1),
           FGFR_exp_comp = exp_FGFR_comp[match(label, split_cell_names(names(exp_FGFR4_CCLE)))]) %>%
    mutate(FGFR2_exon.C3 = FGFR2_exon_usage_CCLE[match(label, split_cell_names(rownames(FGFR2_exon_usage_CCLE))), ind_C3]) %>%
    mutate(FGFR_hotspot_mut = rep("none", nrow(.)),
           FGFR_hotspot_mut = replace(FGFR_hotspot_mut, label %in% split_cell_names(FGFR2_mut_CCLE_hotspot$Tumor_Sample_Barcode), "FGFR2"),
           FGFR_hotspot_mut = replace(FGFR_hotspot_mut, label %in% split_cell_names(FGFR3_mut_CCLE_hotspot$Tumor_Sample_Barcode), "FGFR3"),
           FGFR_hotspot_mut = ifelse(label %in% split_cell_names(unique(mut_CCLE$Tumor_Sample_Barcode)), FGFR_hotspot_mut, NA),
           FGFR_hotspot_mut = factor(FGFR_hotspot_mut, levels = c("none", "FGFR2", "FGFR3")))

save.image("~/FGFR/Daniel/R/Nature_figures/data/CCLE_GDSC_PharmacoDB/CCLE_processing.RData")


