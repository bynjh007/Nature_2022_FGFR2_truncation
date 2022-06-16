library(TCGAbiolinks)
library(dplyr)
library(DT)
library(stringr)

projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]

#######################################
# Illumina Hiseq - using TCGAbiolinks
#######################################
for(proj in projects){
  print(proj)
  query <- GDCquery(project = proj,
                    data.category = "Gene expression",
                    data.type = "Exon junction quantification", 
                    legacy = TRUE,
                    platform = "Illumina HiSeq")
  tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20, 
                       directory = "~/FGFR/Daniel/TCGA/junction_quantification/Illumina_Hiseq/"),
           error = function(e) GDCdownload(query, method = "client",
                                           directory = "~/FGFR/Daniel/TCGA/junction_quantification/Illumina_Hiseq/"))
}

#######################################
# Illumina GA
#######################################
# manifest file was obtained from GDC legacy database
# ~/tools/gdc-client download -m gdc_manifest.2020-02-16_GA.txt


# Both Hiseq and GA data could be downloaded by not setting platform option,,, but I missed it. So I additionally downloaded the GA files.
#query <- GDCquery(project = proj,
                  #data.category = "Gene expression",
                  #data.type = "Exon junction quantification",
                  #legacy = TRUE,
                  #experimental.strategy = "RNA-Seq")

#######################################
# file - sample mapping
#######################################
# sample information
for(i in 1:length(projects)){
  proj = projects[i]
  query <- GDCquery(project = proj,
                    data.category = "Gene expression",
                    data.type = "Exon junction quantification",
                    legacy = TRUE,
                    experimental.strategy = "RNA-Seq")
  match.file.cases <- getResults(query)
  match.file.cases$project <- proj
  if(i == 1){
    match.file.cases.all = match.file.cases
  } else {
    x = match.file.cases[,match(colnames(match.file.cases.all), colnames(match.file.cases))]
    match.file.cases.all = rbind(match.file.cases.all, x)
  }
  print(i)
}

t = c()
for(i in 1:nrow(match.file.cases.all)){
  t = c(t, paste(match.file.cases.all$tags[[i]], collapse = "_"))
}
match.file.cases.all$tags = t
write.table(match.file.cases.all, file = "~/FGFR/Daniel/TCGA/junction_quantification/sample_information.txt", col.names = T, row.names = F, quote = F, sep = "\t")


#######################################
# version 2 samples
#######################################
samples_junc = data.frame(match.file.cases.all[grep("v2", match.file.cases.all$tags),], stringsAsFactors = F)
samples_junc$submitter_id = apply(str_split_fixed(samples_junc$cases, "[-]", 5)[,1:3], 1, paste, collapse = "-")


#######################################
# Illumina Hiseq data loading
#######################################
path_default = "~/FGFR/Daniel/TCGA/junction_quantification/Illumina_Hiseq"
path_folder = list.files(path_default)
path_folder = path_folder[grep("TCGA", path_folder)]

cnt_new = c()
junc_new = c()
for (i in 1:33){
  path_tumor = paste(path_default, path_folder[i], "legacy", "Gene_expression", "Exon_junction_quantification", sep = "/")
  path_files = list.files(path_tumor)
  
  for(j in 1:length(path_files)){
    path_target = paste(path_tumor, path_files[j], sep = "/")
    file_junc = list.files(path_target)
    tab_junc = read.table(paste(path_target, file_junc, sep = "/"), header = T, stringsAsFactors = F)
    
    # extract FGFR2 regions
    t = str_split_fixed(tab_junc$junction, ":|\\-|\\,|\\+", 8)
    t_loc = cbind(as.numeric(t[,2]), as.numeric(t[,6]))
    t_loc = cbind(apply(t_loc, 1, min), apply(t_loc, 1, max))
    ind = which(t[,1] == "chr10" & t_loc[,1] >= 123237844 & t_loc[,2] <= 123357972)
    # target region read counts
    tab_junc = as.matrix(tab_junc$raw_counts[ind])
    colnames(tab_junc) = path_files[j]
    # target region table
    df = cbind(t[ind,1], t_loc[ind,1], t_loc[ind,2])
    t_reg = apply(df, 1, paste, collapse = "_")
    rownames(tab_junc) = t_reg
    # there are different types of junction sections
    
    if(i == 1 & j ==1){
      junc_info = t_reg
      cnt = tab_junc
    } else {
      if(all(junc_info == t_reg)){
        cnt = cbind(cnt, tab_junc)
      } else {
        tab_junc = as.numeric(tab_junc)
        names(tab_junc) = t_reg
        cnt_new = c(cnt_new, tab_junc)
        junc_new = c(junc_new, rep(path_files[j], length(tab_junc)))
      }
    }
    print(paste(j, "in", i, "th tumors"))
  }
}

# summarize new type junction information
# Two additional juncion information
n = unique(junc_new)
for(i in 1:length(n)){
  ind = which(junc_new ==n[i])
  t_junc = names(cnt_new[ind])
  if(i == 1){
    junc_info_extra = cbind(t_junc, 1)
  } else {
    t_ind = as.numeric(unique(junc_info_extra[,2]))
    t_presence = numeric(length(t_ind))
    for(j in 1:length(t_ind)){
      t_ref = junc_info_extra[which(as.numeric(junc_info_extra[,2]) == t_ind[j]),1]
      if(all(t_ref == t_junc)){
        t_presence[j] = 1
      }
    }
    if(sum(t_presence)==0){
      junc_info_extra = rbind(junc_info_extra, 
                              cbind(t_junc, as.numeric(junc_info_extra[nrow(junc_info_extra),2])+1))
    }
  }
  print(i)
}

# make junction count matrix for extra junction information
n = unique(junc_new)
cnt_2 = c()
cnt_3 = c()
for(i in 1:length(n)){
  ind = which(junc_new == n[i])
  t_junc = names(cnt_new[ind])
  t_cnt = as.matrix(cnt_new[ind])
  colnames(t_cnt) = n[i]
  if(all(t_junc == junc_info_extra[junc_info_extra[,2] == "1", 1])){
    cnt_2 = cbind(cnt_2, t_cnt)
  } else if(all(t_junc == junc_info_extra[junc_info_extra[,2] == "2", 1])){
    cnt_3 = cbind(cnt_3, t_cnt)
  } else {
    print("Neither 2 nor 3")
    exit()
  }
}

# combining all type of junctions
cnt_all_Hiseq = list(cnt, cnt_2, cnt_3)
cnt_all_Hiseq_v2 = list(cnt[,colnames(cnt) %in% samples_junc$file_id],
                        cnt_2[,colnames(cnt_2) %in% samples_junc$file_id],
                        cnt_3[,colnames(cnt_3) %in% samples_junc$file_id])


######################################
# Illumina GA data loading
######################################
path_default = "~/FGFR/Daniel/TCGA/junction_quantification/Illumina_GA"
path_folder = list.files(path_default)
path_folder = path_folder[!grepl("txt", path_folder)]

cnt_new_GA = c()
junc_new_GA = c()
for (i in 1:length(path_folder)){
  path_tumor = paste(path_default, path_folder[i], sep = "/")
  path_files = list.files(path_tumor)
  path_files = path_files[grep("quantification.txt", path_files)]
  
  tab_junc = read.table(paste(path_tumor, path_files, sep = "/"), header = T, stringsAsFactors = F)
    
  # extract FGFR2 regions
  t = str_split_fixed(tab_junc$junction, ":|\\-|\\,|\\+", 8)
  t_loc = cbind(as.numeric(t[,2]), as.numeric(t[,6]))
  t_loc = cbind(apply(t_loc, 1, min), apply(t_loc, 1, max))
  ind = which(t[,1] == "chr10" & t_loc[,1] >= 123237844 & t_loc[,2] <= 123357972)
  # target region read counts
  tab_junc = as.matrix(tab_junc$raw_counts[ind])
  colnames(tab_junc) = path_folder[i]
  # target region table
  df = cbind(t[ind,1], t_loc[ind,1], t_loc[ind,2])
  t_reg = apply(df, 1, paste, collapse = "_")
  rownames(tab_junc) = t_reg
  # there are different types of junction sections
    
  if(i == 1){
    junc_info_GA = t_reg
    cnt_GA = tab_junc
  } else {
    if(all(junc_info_GA == t_reg)){
      cnt_GA = cbind(cnt_GA, tab_junc)
    } else {
      tab_junc = as.numeric(tab_junc)
      names(tab_junc) = t_reg
      cnt_new_GA = c(cnt_new_GA, tab_junc)
      junc_new_GA = c(junc_new_GA, rep(path_folder[i], length(tab_junc)))
    }
  }
  print(i)
}

# summarize new type junction information
# Two additional juncion information
n = unique(junc_new_GA)
for(i in 1:length(n)){
  ind = which(junc_new_GA ==n[i])
  t_junc = names(cnt_new_GA[ind])
  if(i == 1){
    junc_info_extra_GA = cbind(t_junc, 1)
  } else {
    t_ind = as.numeric(unique(junc_info_extra_GA[,2]))
    t_presence = numeric(length(t_ind))
    for(j in 1:length(t_ind)){
      t_ref = junc_info_extra_GA[which(as.numeric(junc_info_extra_GA[,2]) == t_ind[j]),1]
      if(all(t_ref == t_junc)){
        t_presence[j] = 1
      }
    }
    if(sum(t_presence)==0){
      junc_info_extra_GA = rbind(junc_info_extra_GA, 
                              cbind(t_junc, as.numeric(junc_info_extra_GA[nrow(junc_info_extra_GA),2])+1))
    }
  }
  print(i)
}

# make junction count matrix for extra junction information
n = unique(junc_new_GA)
cnt_2_GA = c()
cnt_3_GA = c()
for(i in 1:length(n)){
  ind = which(junc_new_GA == n[i])
  t_junc = names(cnt_new_GA[ind])
  t_cnt = as.matrix(cnt_new_GA[ind])
  colnames(t_cnt) = n[i]
  if(all(t_junc == junc_info_extra_GA[junc_info_extra_GA[,2] == "1", 1])){
    cnt_2_GA = cbind(cnt_2_GA, t_cnt)
  } else if(all(t_junc == junc_info_extra_GA[junc_info_extra_GA[,2] == "2", 1])){
    cnt_3_GA = cbind(cnt_3_GA, t_cnt)
  } else {
    print("Neither 2 nor 3")
    exit()
  }
}

# combining all type of junctions
cnt_all_GA = list(cnt_GA, cnt_2_GA, cnt_3_GA)
cnt_all_GA_v2 = list(cnt_GA[,colnames(cnt_GA) %in% samples_junc$file_id],
                     cnt_2_GA[,colnames(cnt_2_GA) %in% samples_junc$file_id],
                     cnt_3_GA[,colnames(cnt_3_GA) %in% samples_junc$file_id])


######################################
# combining Hiseq and GA data
######################################
# all the version 2 data have the same junction list --> merge Hiseq-v2 and GA-v2 data
all(rownames(cnt_all_Hiseq_v2[[1]]) == rownames(cnt_all_GA_v2[[1]]))
cnt_merge_v2 = cbind(cnt_all_Hiseq_v2[[1]], cnt_all_GA_v2[[1]])
# matching samples_junc order and cnt file
cnt_merge_v2 = cnt_merge_v2[,match(samples_junc$id, colnames(cnt_merge_v2))]
all(colnames(cnt_merge_v2) == samples_junc$id)


save.image("~/FGFR/Daniel/R/Nature_figures/data/TCGA/FGFR2_junc_quantification_TCGA.RData")

