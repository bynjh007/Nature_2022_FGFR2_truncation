library(dplyr)

# PURPLE and meta data loading
FGFR2_SV_PURPLE_ALL = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/HMF/HMF_data_purple.RData")[["FGFR2_PURPLE_ALL"]]
FGFR2_CNV = readRDS("~/FGFR/Daniel/R/Nature_figures/data/HMF/HMF_data_cnv.RDS")
FGFR2_SV_PURPLE_RNA = R.utils::loadToEnv("~/FGFR/Daniel/R/Nature_figures/data/HMF/HMF_data_SV_RNA.RData")[["FGFR2_PURPLE_ALL"]]

target_samp = unique(FGFR2_SV_PURPLE_ALL %>% pull(sample_name))
samp_info_meta = read.table("~/FGFR/Daniel/R/Nature_figures/data/HMF/metadata.tsv", sep = "\t", stringsAsFactors = F, header = T, comment.char = "")
samp_info_target = samp_info_meta[match(target_samp, samp_info_meta$setName),]

# Redundant samples with multiple biopysies 
FGFR2_SV_PURPLE_ALL$hmfPatientId = samp_info_meta$hmfPatientId[match(FGFR2_SV_PURPLE_ALL$sample_name, samp_info_meta$setName)]

redun_samp = unique(FGFR2_SV_PURPLE_ALL %>% select(hmfPatientId, sample_name)) %>% 
    group_by(hmfPatientId) %>% summarise(n = n()) %>% filter(n>1) %>% pull(hmfPatientId)
redun_samp_C_trunc = unique(FGFR2_SV_PURPLE_ALL %>% filter(hmfPatientId %in% redun_samp, C_trunc == T) %>% pull(hmfPatientId))

# Remove redundant samples
samp_nonred_no_C_trunc = FGFR2_SV_PURPLE_ALL$sample_id[match(setdiff(redun_samp, redun_samp_C_trunc), FGFR2_SV_PURPLE_ALL$hmfPatientId)]
samp_nonred_C_trunc = c(47, 80, 82, 35)
samp_others = FGFR2_SV_PURPLE_ALL %>% filter(hmfPatientId %in% setdiff(unique(FGFR2_SV_PURPLE_ALL$hmfPatientId), redun_samp)) %>% pull(sample_id)
FGFR2_SV_PURPLE_ALL_noredunant = FGFR2_SV_PURPLE_ALL %>% filter(sample_id %in% c(samp_others, samp_nonred_C_trunc, samp_nonred_no_C_trunc))

# Discarding SVs that have unidentified partner or upstream is out-of-strand FGFR2
FGFR2_SV_PURPLE_TARGET = FGFR2_SV_PURPLE_ALL_noredunant %>% 
  filter(RE_type != "SGL", # removing SGL breakend
         !(up.RE_category == "out-of-strand" & RE_type != "internal")) # removing the SVs with upstream is out-of-strand orientation

FGFR2_SV_PURPLE_TARGET = FGFR2_SV_PURPLE_TARGET %>% 
  mutate(RE_FGFR2 = as.character(RE_type)) %>%
  mutate(Tumor_type = samp_info_target$primaryTumorLocation[match(sample_name, samp_info_target$setName)],
         Tumor_type = factor(Tumor_type, levels = names(sort(table(Tumor_type), decreasing = T)))) %>%
  mutate(FGFR2_cnv = FGFR2_CNV[match(sample_name, samp_info_target$setName),1]) %>%
  mutate(RE_FGFR2 = replace(RE_FGFR2, grepl("out-of-frame", RE_FGFR2), "Frame unknown"),
         RE_FGFR2 = replace(RE_FGFR2, RE_FGFR2 == "in-frame", "In-frame fusion"),
         RE_FGFR2 = replace(RE_FGFR2, RE_FGFR2 == "intergenic", "Intergenic space"),
         RE_FGFR2 = replace(RE_FGFR2, RE_FGFR2 == "out-of-strand", "Out-of-strand"),
         RE_FGFR2 = replace(RE_FGFR2, RE_FGFR2 == "internal", "Internal")) %>%
  mutate(FGFR2_SV_loc = ifelse(up.gene_symbol == "FGFR2", up.loc_in_trx, down.loc_in_trx),
         FGFR2_brkpt = rep(NA, nrow(FGFR2_SV_PURPLE_TARGET)),
         FGFR2_brkpt = replace(FGFR2_brkpt, stringr::str_split_fixed(FGFR2_SV_loc, ":", 6)[,6] == "Intron 17", "I17"),
         FGFR2_brkpt = replace(FGFR2_brkpt, stringr::str_split_fixed(FGFR2_SV_loc, ":", 6)[,6] == "Exon 18", "E18"),
         FGFR2_brkpt = replace(FGFR2_brkpt, is.na(FGFR2_brkpt), "5' to E17"),
         FGFR2_brkpt = factor(FGFR2_brkpt, levels = c("I17", "E18", "5' to E17")),
         parter_SV_chr = ifelse(up.gene_symbol == "FGFR2", as.character(down.seqnames), as.character(up.seqnames)),
         FGFR2_chr_SV = ifelse(parter_SV_chr == "chr10", "Intrachromosomal", "Interchromosomal"),
         RE_location = ifelse(FGFR2_is_upstream == TRUE, "FGFR2 is upstream", "FGFR2 is downstream"),
         fusion_id = paste(up.gene_symbol, down.gene_symbol, sep = "_")) %>%
  select(-RE_type)

a = do.call("paste", c(FGFR2_SV_PURPLE_TARGET[, c("up.seqnames", "up.start", "up.loc_in_trx", "up.orientation", 
                                                  "down.seqnames", "down.start", "down.loc_in_trx", "down.orientation",
                                                  "sample_name")], sep = "_"))
b = do.call("paste", c(FGFR2_SV_PURPLE_RNA[, c("up.seqnames", "up.start", "up.loc_in_trx", "up.orientation", 
                                               "down.seqnames", "down.start", "down.loc_in_trx", "down.orientation",
                                               "sample_name")], sep = "_"))
FGFR2_SV_PURPLE_TARGET = data.frame(FGFR2_SV_PURPLE_TARGET, 
                                    FGFR2_SV_PURPLE_RNA[match(a,b), c("STAR_Fusion", "n_RNA", "STAR_chim")])
FGFR2_SV_PURPLE_TARGET = FGFR2_SV_PURPLE_TARGET %>% mutate(n_RNA = replace(n_RNA, RE_FGFR2 == "Internal", NA))

# This BP is located in the last intron of the longest isoform in RNLS (hg19) but that isoform is excluded in hg38, therefore this RE should be re-classified as intergenic RE
FGFR2_SV_PURPLE_TARGET$RE_FGFR2[which(FGFR2_SV_PURPLE_TARGET$STAR_chim == "FGFR2_Intron 17_chr10:121483697_NA_intergenic::PTEN--RNLS_chr10:88183205")] = "Intergenic space"
rm(a,b)


#############################################################
# copy numbers and RE heterogeneity
#############################################################
samp_id = unique(FGFR2_SV_PURPLE_TARGET %>% pull(sample_id))
trunc_samp_id = unique(FGFR2_SV_PURPLE_TARGET %>% filter(C_trunc == T) %>% pull(sample_id))

# one SV can have multiple chain ID --> unfold
df_unfolded = c()
for(i in 1:nrow(FGFR2_SV_PURPLE_TARGET)){
  t = strsplit(FGFR2_SV_PURPLE_TARGET$Cluster_Chain[i], "/")[[1]]
  df_unfolded = rbind(df_unfolded, 
                      FGFR2_SV_PURPLE_TARGET[i,] %>% dplyr::slice(rep(1:n(), each = length(t))) %>% mutate(Cluster_Chain_2 = t))
}

CN_samp = matrix(NA, length(samp_id), 3)
df_CN_samp = c()
for(i in 1:length(samp_id)){

    ####################################
    # chain links information
    ####################################
    sample_name = unique(FGFR2_SV_PURPLE_TARGET %>% filter(sample_id == samp_id[i]) %>% pull(sample_name))
    t = list.files(paste("/DATA/projects/j.bhin/Daniel_FGFR2/HMF/linx", sample_name, sep = "/"))
    links_data = read.table(paste("/DATA/projects/j.bhin/Daniel_FGFR2/HMF/linx", sample_name, 
                                  t[grep("linx.links.tsv", t)], sep = "/"), sep = "\t", stringsAsFactors = F, header = T) %>%
        mutate(Cluster_Chain = paste(clusterId, chainId, sep = "_"))
 
    # unfolded SVs-chains
    samp_SV = df_unfolded %>% filter(sample_id == samp_id[i])
    if(samp_id[i] %in% trunc_samp_id){
        samp_SV_trunc = samp_SV %>% filter(C_trunc == T)

        ####################################
        # chains with truncation SVs
        ####################################
        chains_trunc = unlist(strsplit(unique(samp_SV_trunc %>% pull(Cluster_Chain)), "/"))
        chains_trunc = chains_trunc[!is.na(chains_trunc)]
          
        # chains with truncation SVs do not exist in linx output --> Ploidy information from PURPLE output
        unmatched_links = chains_trunc[which(is.na(match(chains_trunc, links_data$Cluster_Chain)))]
        unmatched_links_trunc_PL = unlist(lapply(unmatched_links, function(x){mean(unique(samp_SV_trunc %>% 
                                                                                            filter(Cluster_Chain_2 %in% x)) %>% pull(Ploidy))}))
        # chains with truncation SVs exist in linx output --> Ploidy information from LINX_link
        matched_links = setdiff(chains_trunc, unmatched_links)
        matched_links_trunc_PL = links_data$ploidy[match(matched_links, links_data$Cluster_Chain)]
          
        # Ploidy information for C-terminal truncated chains
        df_trunc = data.frame(cluster_chain = c(unmatched_links, matched_links), 
                        ploidy = c(unmatched_links_trunc_PL, matched_links_trunc_PL),
                        group = "truncated", sample = samp_id[i]) %>% 
            mutate(SVID = samp_SV$SVID[match(cluster_chain, samp_SV$Cluster_Chain_2)])
    } else {
        chains_trunc = c()
        unmatched_links_trunc_PL = c()
        matched_links_trunc_PL = c()
        df_trunc = c()
    }

    ####################################
    # chains without truncation SVs
    ####################################
    chains_notrunc = setdiff(unlist(strsplit(unique(samp_SV %>% pull(Cluster_Chain)), "/")), chains_trunc)
    chains_notrunc = chains_notrunc[!is.na(chains_notrunc)]

    # chains with truncation SVs do not exist in linx output --> Ploidy information from PURPLE output
    unmatched_links = chains_notrunc[which(is.na(match(chains_notrunc, links_data$Cluster_Chain)))]
    unmatched_links_notrunc_PL = unlist(lapply(unmatched_links, function(x){mean(unique(samp_SV %>% filter(Cluster_Chain_2 %in% x)) %>% pull(Ploidy))}))
    # chains with truncation SVs exist in linx output --> Ploidy information from LINX_link
    matched_links = setdiff(chains_notrunc, unmatched_links)
    matched_links_notrunc_PL = links_data$ploidy[match(matched_links, links_data$Cluster_Chain)]

    # Ploidy information for chains without C-terminal truncations
    if(length(chains_notrunc)!=0){
    df_notrunc = data.frame(cluster_chain = c(unmatched_links, matched_links), 
                            ploidy = c(unmatched_links_notrunc_PL, matched_links_notrunc_PL),
                            group = "non_truncated", sample = samp_id[i]) %>%
      mutate(SVID = samp_SV$SVID[match(cluster_chain, samp_SV$Cluster_Chain_2)])

    } else {
    df_notrunc = NULL
    }

    # ratio between truncs and nontruncs chains
    notrunc_ploidy = sum(c(unmatched_links_notrunc_PL, matched_links_notrunc_PL), na.rm = T)
    trunc_ploidy = sum(c(unmatched_links_trunc_PL, matched_links_trunc_PL), na.rm = T)
    trunc_ratio = trunc_ploidy/(notrunc_ploidy + trunc_ploidy)
    CN_samp[i,] = c(notrunc_ploidy, trunc_ploidy, trunc_ratio)

    # list of truncated and non-truncated chains
    t_df = rbind(df_trunc, df_notrunc)
    t_uniq_svid = unique(t_df$SVID)
    t_df = t_df %>% mutate(chains_label = 1:nrow(.), SVID_label = match(SVID, t_uniq_svid))
    df_CN_samp = rbind(df_CN_samp, t_df)

}
colnames(CN_samp) = c("nontrunc_PL", "trunc_PL", "trunc_ratio")

# annotation
df_CN_samp = df_CN_samp %>% 
    mutate(sample_name = samp_info_target$sampleId[match(FGFR2_SV_PURPLE_TARGET$sample_name[match(sample, FGFR2_SV_PURPLE_TARGET$sample_id)], 
                                                         samp_info_target$setName)]) %>%
    mutate(ID = paste(sample, cluster_chain, SVID, sep = "_")) %>%
    mutate(RE_FGFR2 = as.character(df_unfolded$RE_FGFR2[match(ID, paste(df_unfolded$sample_id, df_unfolded$Cluster_Chain_2, df_unfolded$SVID, sep = "_"))])) %>%
    mutate(RE_loc = df_unfolded$FGFR2_brkpt[match(ID, paste(df_unfolded$sample_id, df_unfolded$Cluster_Chain_2, df_unfolded$SVID, sep = "_"))])

# all the I17/E18 REs are truncated REs : 
df_CN_samp = df_CN_samp %>%
    mutate(RE_FGFR2 = ifelse(group == "non_truncated", "5'-E17", RE_FGFR2))

rm(i, t_df, t, t_uniq_svid, trunc_ploidy, trunc_samp_id, trunc_ratio, unmatched_links, unmatched_links_notrunc_PL, unmatched_links_trunc_PL,
   chains_notrunc, chains_trunc, chains_label, matched_links, matched_links_notrunc_PL, matched_links_trunc_PL)

save.image("~/FGFR/Daniel/R/Nature_figures/data/HMF/HMF_data_samples_filtering.RData")
