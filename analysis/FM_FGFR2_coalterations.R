library(dplyr)

#####################################################################################
# data loading 
#####################################################################################
# samples to be analyzed
FMI_tumors = read.table("~/FGFR/Daniel/R/Nature_figures/data/FM/FM_tumors.txt", header = T, sep = "\t", stringsAsFactors = F)

# loading FGFR2 mutations
FGFR2_mut_all = read.table("~/FGFR/Daniel/R/Nature_figures/data/FM/FGFR2_mutation_all_MutationMapper.txt", header = T, sep = "\t", stringsAsFactors = F, quote = "", fill = T)


# loading co-mutations based on the "FMI_tumors"
FGFR2_mut_concurrent_pancancer = xlsx::read.xlsx2("~/FGFR/Daniel/R/Nature_figures/data/FM/FM FGFR2 Coalts.xlsx", sheetIndex = "All FMI tumors", header = T) %>% filter(variantType != "RE")
FGFR2_mut_concurrent_breast = xlsx::read.xlsx2("~/FGFR/Daniel/R/Nature_figures/data/FM/FM FGFR2 Coalts.xlsx", sheetIndex = "Breast", header = T) %>% filter(variantType != "RE")

####################################################################################
# all tumors
#####################################################################################
t_samp = rownames(FMI_tumors) 
t_mut_all = unique(names(sort(table(FGFR2_mut_concurrent_pancancer$gene), decreasing = T)))
t_mat_all = matrix(NA, length(t_samp), 30)
for(i in 1:length(t_samp)){
  for(j in 1:30){
    t_df = FGFR2_mut_concurrent_pancancer %>% filter(deidentifiedSpecimenName == t_samp[i], gene == t_mut_all[j])
    t_cn = t_df %>% filter(variantType == "CN") %>% pull(amplificationDeletion_CN)
    t_sv = t_df %>% filter(variantType == "SV") %>% pull(proteinEffect_SV)
    t_mat_all[i,j] = paste(c(t_cn, t_sv), collapse = ", ")
  }
  print(i)
}
colnames(t_mat_all) = t_mut_all
t_mat_all[t_mat_all ==""] = NA

FMI_tumors_all = data.frame(Sample_ID = rownames(FMI_tumors), FMI_tumors, t_mat_all)

#####################################################################################
# breast tumors
#####################################################################################
t_samp = rownames((FMI_tumors %>% filter(TCGA_type == "BRCA"))) 
t_mut_BRCA = unique(names(sort(table(FGFR2_mut_concurrent_breast$gene), decreasing = T)))
t_mat_BRCA = matrix(NA, length(t_samp), 30)
for(i in 1:length(t_samp)){
  for(j in 1:30){
    t_df = FGFR2_mut_concurrent_breast %>% filter(deidentifiedSpecimenName == t_samp[i], gene == t_mut_BRCA[j])
    t_cn = t_df %>% filter(variantType == "CN") %>% pull(amplificationDeletion_CN)
    t_sv = t_df %>% filter(variantType == "SV") %>% pull(proteinEffect_SV)
    t_mat_BRCA[i,j] = paste(c(t_cn, t_sv), collapse = ", ")
  }
  print(i)
}
colnames(t_mat_BRCA) = t_mut_BRCA
t_mat_BRCA[t_mat_BRCA ==""] = NA

FMI_tumors_BRCA = data.frame(Sample_ID = t_samp, FMI_tumors %>% filter(TCGA_type == "BRCA"), t_mat_BRCA)

#####################################################################################
# mutation re-labeling
#####################################################################################
# Mutation in the other genes
mut_sum_func = function(X){
  
  mut_type = function(Y){
    # splice site
    ss = grepl("splice", Y)
    # truncation (nonsense, frameshift)
    tr = grepl("[*]", Y)
    # non-frameshift mutation by indel
    nonfs = grepl("ins|del", Y) & !grepl("splice", Y) & !grepl("deletion", Y)
    # promoter mut
    pr = grepl("promoter", Y)
    # amplification
    amp = grepl("amplification", Y)
    # deletion
    del = grepl("deletion", Y)
    tmp_type = rowSums(rbind(ss, tr, nonfs, pr, amp, del))
    return(tmp_type)
  }
  
  # multiple variations
  if(grepl("[,]", X)){
    n = stringr::str_count(X, "[,]") + 1
    tmp = stringr::str_split_fixed(X, "[,]", n)
    mut_sum = mut_type(tmp) # mutation type except for missense
    
    # add missense mutation
    mut_sum = c(n-sum(mut_sum), mut_sum)
    mut_sum[mut_sum>=1] = 1
    
    # unique variations  
  } else {
    mut_sum = mut_type(X) # mutation type except for missense
    
    # add missense mutation
    if(sum(mut_sum)==0){
      mut_sum = c(1, mut_sum)
    } else {
      mut_sum = c(0, mut_sum)
    }
  }
  
  names(mut_sum) = c("missense", "splice_site", "truncation", "nonframeshift_indel", "promoter", "amplification", "deletion")
  out_code_sum = paste(names(mut_sum[mut_sum!=0]), collapse = ";")
  
  # take representative mutation in case the samples with multiple mutations (Truncation >=Missense mutation)
  mut_sum_rep = c(mut_sum[1], mut_sum[3], sum(mut_sum[c(2,4)]), mut_sum[5], mut_sum[6:7])
  mut_sum_rep[mut_sum_rep>=1] = 1
  names(mut_sum_rep) = c("missense", "truncation", "splice_nonfs",  "promoter", "amplification", "deletion")
  if(mut_sum_rep[1]==1 & mut_sum_rep[2]==1){ # Truncation > Missense mutation
    mut_sum_rep[1] =0
  }
  if(mut_sum_rep[2]==1 & mut_sum_rep[3]==1){ # Truncation > splice & nonfs
    mut_sum_rep[3] =0
  }
  if(mut_sum_rep[1]==1 & mut_sum_rep[3]==1){ # splice & nonfs > Missense mutation
    mut_sum_rep[1] =0
  }
  out_code_sum_rep = paste(names(mut_sum_rep[mut_sum_rep!=0]), collapse = ";")
  
  out_code = list(sum = out_code_sum_rep, all = out_code_sum)
  return(out_code)
}

# for all tumors
mut_sum_mat_all = matrix(NA, nrow(t_mat_all), ncol(t_mat_all))
mut_all_mat_all = matrix(NA, nrow(t_mat_all), ncol(t_mat_all))
for(i in 1:nrow(t_mat_all)){
  for(j in 1:ncol(t_mat_all)){
    if(!is.na(t_mat_all[i,j])){
      t = mut_sum_func(t_mat_all[i,j])
      mut_sum_mat_all[i,j] = t$sum
      mut_all_mat_all[i,j] = t$all
    }
  }
}
colnames(mut_sum_mat_all) = colnames(t_mat_all)
colnames(mut_all_mat_all) = colnames(t_mat_all)
rownames(mut_sum_mat_all) = FMI_tumors_all$Sample_ID
rownames(mut_all_mat_all) = FMI_tumors_all$Sample_ID
mut_binary_mat_all = (!is.na(mut_sum_mat_all)) * 1

# for breast tumors
mut_sum_mat_BRCA = matrix(NA, nrow(t_mat_BRCA), ncol(t_mat_BRCA))
mut_all_mat_BRCA = matrix(NA, nrow(t_mat_BRCA), ncol(t_mat_BRCA))
for(i in 1:nrow(t_mat_BRCA)){
  for(j in 1:ncol(t_mat_BRCA)){
    if(!is.na(t_mat_BRCA[i,j])){
      t = mut_sum_func(t_mat_BRCA[i,j])
      mut_sum_mat_BRCA[i,j] = t$sum
      mut_all_mat_BRCA[i,j] = t$all
    }
  }
}
colnames(mut_sum_mat_BRCA) = colnames(t_mat_BRCA)
colnames(mut_all_mat_BRCA) = colnames(t_mat_BRCA)
rownames(mut_sum_mat_BRCA) = FMI_tumors_BRCA$Sample_ID
rownames(mut_all_mat_BRCA) = FMI_tumors_BRCA$Sample_ID
mut_binary_mat_BRCA = (!is.na(mut_sum_mat_BRCA)) * 1

rm(i, j, t, t_cn, t_df, t_mat_all, t_mat_BRCA, t_mut_all, t_mut_BRCA, t_sv, t_samp)

save.image("~/FGFR/Daniel/R/Nature_figures/data/FM/FM_FGFR2_coalterations.RData")
#
