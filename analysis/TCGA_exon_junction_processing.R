################################################################################
# Mining the samples where C-terminal truncation of FGFR2 might occurs
# Exon-level quantification files were obtained from GDC Legacy Archive
# Filtering out the samples with low expression of FGFR2 gene : log2(RSEM)<5 [gene level] & median(RPKM)>=1 [exon level]

################################################################################
# gene expression profile
################################################################################
# zscore from cBioportal
FGFR2_exp_zscore = read.table("~/FGFR/Daniel/R/Nature_figures/data/TCGA/mRNA Expression z-Scores (RNA Seq V2 RSEM).txt", header = T, stringsAsFactors = F, sep = "\t")
# gene expression - normalized rsem read counts from firehose
FGFR2_exp_rsem = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/FGFR2_gene_expression_TCGA.RData")[["FGFR2_data"]]
names(FGFR2_exp_rsem) = substring(text = gsub("[.]", "-", names(FGFR2_exp_rsem)), first = 1, last = 15)
# remove redundant samples by taking the samples with high exon expression
FGFR2_exp_rsem = sort(FGFR2_exp_rsem, decreasing = T)
t = unique(names(FGFR2_exp_rsem))
FGFR2_exp_rsem = FGFR2_exp_rsem[match(t, names(FGFR2_exp_rsem))]
# shuffling the samples for heatmap visualization
set.seed(123)
FGFR2_exp_rsem = FGFR2_exp_rsem[sample(length(FGFR2_exp_rsem))]
FGFR2_exp_log = log2(FGFR2_exp_rsem+1)


################################################################################
# Exon-level quantification (RPKM) from GDC Archeive regacy
################################################################################
samples_exon = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/FGFR2_exon_expression_TCGA_Legacy.RData")[["samples_exon"]]
samples_exon$sample_id = substr(samples_exon$cases, 1, 15)
exon_list = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/FGFR2_exon_expression_TCGA_Legacy.RData")[["exon_list"]]
FGFR2_exon_rpkm = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/FGFR2_exon_expression_TCGA_Legacy.RData")[["FGFR2_exon"]]
# sorting according to the rpkm value
ind = order(apply(FGFR2_exon_rpkm, 2, median), decreasing = T)
FGFR2_exon_rpkm = FGFR2_exon_rpkm[,ind]
samples_exon = samples_exon[ind,]
# selecting the redundant samples based on exon expression
t = unique(samples_exon$sample_id)
ind = match(t, samples_exon$sample_id)
FGFR2_exon_rpkm = FGFR2_exon_rpkm[, ind]
samples_exon = samples_exon[ind,]


########################################################
# C1 exon usage based on exon-level expression
########################################################
#Identification of the samples with low expression of C1  
#Genomic rearrangement (at intron 17) or alternative splicing (C3 usage) could be the reason of the low expression of C1 (truncation)  
#No difference of relative C1 expression between normal vs tumor  
#Some of the tumors showed lower rel.C1 than minimun of those of normal

# selecting the samples expressing E17 (because E17 is a kinase domain, it should be expressed for oncogenic activity)
ind_N = samples_exon$sample_id[which(stringr::str_split_fixed(samples_exon$sample_id, "-", 4)[,4] == "11" & FGFR2_exon_rpkm[3,]>0)]
ind_T = samples_exon$sample_id[which(stringr::str_split_fixed(samples_exon$sample_id, "-", 4)[,4] != "11" & FGFR2_exon_rpkm[3,]>0)]

# filtering the samples with low expression of FGFR2 gene <5 (log2 scale)
exp_samp_gene = names(which(FGFR2_exp_log>=5))
# filtering the samples with low expression of FGFR2 exon <1 (log2 scale)
ind_exons_filt = c(2, 12, 13, 20) # non-canonical exons / alternatively expressed exons --> minor exons
t_exon = FGFR2_exon_rpkm[-ind_exons_filt,]
exp_samp_exon = samples_exon$sample_id[which(apply(t_exon, 2, median) >=1)]
# intersection of the samples from gene-based and exon-based filtering
exp_samp = intersect(exp_samp_gene, exp_samp_exon)

ind_N = which(samples_exon$sample_id %in% intersect(ind_N, exp_samp))
ind_T = which(samples_exon$sample_id %in% intersect(ind_T, exp_samp))

# exon expresion ratio between C1 and median of E1-E17
exon_ratio_C1_Med = t_exon[1,] / apply(t_exon[-1,], 2, median)
names(exon_ratio_C1_Med) = samples_exon$sample_id
#boxplot(exon_ratio_C1_Med[ind_N], exon_ratio_C1_Med[ind_T], names = c("Normal", "Tumor"), ylab = "C1/(E1-E17)")
#boxplot(exon_ratio_C1_Med[ind_N], exon_ratio_C1_Med[ind_T], names = c("Normal", "Tumor"), ylab = "C1/(E1-E17)", ylim = c(0,3))
#length(which(exon_ratio_C1_Med[ind_T] < min(exon_ratio_C1_Med[ind_N])))

# samples with C1 expression is lower than minimum of normal samples
target_samp_C1_Med = names(which(exon_ratio_C1_Med[ind_T] < min(exon_ratio_C1_Med[ind_N])))


########################################################
# Exon usage based on junction reads
########################################################
#Junction list (1: E17-C4, 2: E17-C1, 3: E17-C2, 4: E17-C3)  
#Like exon expression-based C3 usage, junction-based C3 usage also showed clear downregulation in tumor comparing to normal  
#should do the same thing for second type junction files

junc_cnt = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/FGFR2_junc_quantification_TCGA.RData")[["cnt_merge_v2"]]
samples_junc = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/TCGA/FGFR2_junc_quantification_TCGA.RData")[["samples_junc"]]
samples_junc$sample_id = substr(samples_junc$cases, 1, 15)

E17_junc_cnt = junc_cnt[1:4,]
E17_junc_ratio = t(t(E17_junc_cnt) / colSums(E17_junc_cnt))

# normal vs tumor
t = stringr::str_split_fixed(samples_junc$cases, "-", 7)
tumor_solid = samples_junc[!grepl("11", t[,4]),]
normal_solid = samples_junc[grepl("11", t[,4]),]

tumor_E17_junc_cnt = E17_junc_cnt[,which(colnames(E17_junc_cnt) %in% tumor_solid$id)]
normal_E17_junc_cnt = E17_junc_cnt[,which(colnames(E17_junc_cnt) %in% normal_solid$id)]

tumor_E17_junc_ratio = t(t(tumor_E17_junc_cnt) / colSums(tumor_E17_junc_cnt))
normal_E17_junc_ratio = t(t(normal_E17_junc_cnt) / colSums(normal_E17_junc_cnt))

colnames(tumor_E17_junc_cnt) = samples_junc$cases[match(colnames(tumor_E17_junc_cnt), samples_junc$id)]
colnames(tumor_E17_junc_ratio) = samples_junc$cases[match(colnames(tumor_E17_junc_ratio), samples_junc$id)]

# sorting based on C3 ratio
tumor_E17_junc_cnt = tumor_E17_junc_cnt[,order(tumor_E17_junc_ratio[2,])]
tumor_E17_junc_ratio = tumor_E17_junc_ratio[,order(tumor_E17_junc_ratio[2,])]

# E17 junction reads larger than 5
ind_t1 = which(colSums(tumor_E17_junc_cnt) >= 5 & tumor_E17_junc_cnt[1,]>2)
ind_t2 = which(colSums(tumor_E17_junc_cnt) >= 5 & tumor_E17_junc_cnt[2,]>2)
ind_t3 = which(colSums(tumor_E17_junc_cnt) >= 5 & tumor_E17_junc_cnt[3,]>2)
ind_t4 = which(colSums(tumor_E17_junc_cnt) >= 5 & tumor_E17_junc_cnt[4,]>2)
ind_n = colSums(normal_E17_junc_cnt) >= 5

# boxplot(100*tumor_E17_junc_ratio[1, ind_t1], 100*normal_E17_junc_ratio[1, ind_n], ylab = "C4 usage (%)", names = c("Tumor", "Normal"), main = ("C4-usage (juntion-reads ratio)"))
# boxplot(tumor_E17_junc_cnt[1, ind_t1], normal_E17_junc_cnt[1, ind_n], ylab = "C4 usage (read counts)", names = c("Tumor", "Normal"), main = ("C4-usage (juntion-reads counts)"))
# 
# boxplot(100*tumor_E17_junc_ratio[2, ind_t2], 100*normal_E17_junc_ratio[2, ind_n], ylab = "C1 usage (%)", names = c("Tumor", "Normal"), main = ("C1-usage (juntion-reads ratio)"))
# boxplot(tumor_E17_junc_cnt[2, ind_t2], normal_E17_junc_cnt[2, ind_n], ylab = "C1 usage (read counts)", names = c("Tumor", "Normal"), main = ("C1-usage (juntion-reads counts)"))
# 
# boxplot(100*tumor_E17_junc_ratio[3, ind_t3], 100*normal_E17_junc_ratio[3, ind_n], ylab = "C2 usage (%)", names = c("Tumor", "Normal"), main = ("C2-usage (juntion-reads ratio)"))
# boxplot(tumor_E17_junc_cnt[3, ind_t3], normal_E17_junc_cnt[3, ind_n], ylab = "C2 usage (read counts)", names = c("Tumor", "Normal"), main = ("C2-usage (juntion-reads counts)"))
# 
# boxplot(100*tumor_E17_junc_ratio[4, ind_t4], 100*normal_E17_junc_ratio[4, ind_n], ylab = "C3 usage (%)", names = c("Tumor", "Normal"), main = ("C3-usage (juntion-reads ratio)"))
# boxplot(tumor_E17_junc_cnt[4, ind_t4], normal_E17_junc_cnt[4, ind_n], ylab = "C3 usage (read counts)", names = c("Tumor", "Normal"), main = ("C3-usage (juntion-reads counts)"))

# Target samples with high C3 usage (larger than maximum of normal)
t = names(which(tumor_E17_junc_ratio[4, ind_t4] > max(normal_E17_junc_ratio[4, ind_n])))
target_samp_C3_junc = unique(substr(t, 1, 15))
# Target samples with high C4 usage (larger than maximum of normal)
t = names(which(tumor_E17_junc_ratio[1, ind_t1] > max(normal_E17_junc_ratio[1, ind_n])))
target_samp_C4_junc = unique(substr(t, 1, 15))

save.image("~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_exon_junction_processing.RData")

