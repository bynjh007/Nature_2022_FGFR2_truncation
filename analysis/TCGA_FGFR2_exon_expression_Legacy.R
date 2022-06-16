library(TCGAbiolinks)
library(dplyr)
library(DT)
library(stringr)

################################################
# data download
################################################
# data was downloaded by gdc-clients with the manifest file ("~/FGFR/Daniel/TCGA/exon_quantification/gdc_manifest.2019-10-04_exon.txt")

################################################
# sample annotation
################################################
projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]
for(i in 1:length(projects)){
  proj = projects[i]
  query <- GDCquery(project = proj,
                    data.category = "Gene expression",
                    data.type = "Exon quantification",
                    legacy = TRUE,
                    experimental.strategy = "RNA-Seq")
  
  match.file.cases <- getResults(query)
  match.file.cases$project <- proj
  if(i == 1){
    match.file.cases.all = match.file.cases[, c("version", "cases", "tags", "platform", "file_name", "file_id", "experimental_strategy",
                                                "id", "data_type", "project", "tissue.definition")]
    t_lab = colnames(match.file.cases)
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
write.table(match.file.cases.all, file = "~/FGFR/Daniel/TCGA/exon_quantification/sample_information.txt", col.names = T, row.names = F, quote = F, sep = "\t")

################################################
# data loading
################################################
# Using verson2 file (most of the samples with v1 are overlapped with those of v2)
# all the version 2 data have the same FGFR2 exon information #chr10:123237845-123239627:- to chr10:123357476-123357972:-
samples_exon = data.frame(match.file.cases.all[grep("v2", match.file.cases.all$tags),], stringsAsFactors = F)
samples_exon$submitter_id = apply(str_split_fixed(samples_exon$cases, "[-]", 5)[,1:3], 1, paste, collapse = "-")

FGFR2_exon = matrix(NA, 21, nrow(samples_exon))
for(i in 1:nrow(FGFR2_exon)){
  t = list.files(paste("~/FGFR/Daniel/TCGA/exon_quantification/", samples_exon$file_id[i], sep = ""))
  t = t[grep("quantification.txt", t)]
  t_exon = read.table(paste("~/FGFR/Daniel/TCGA/exon_quantification/", samples_exon$file_id[i], "/", t, sep = ""), sep = "\t", stringsAsFactors = F, header = T)
  
  ind = which(t_exon$exon=="chr10:123237845-123239627:-") : which(t_exon$exon=="chr10:123357476-123357972:-")
  FGFR2_exon[,i] = t_exon$RPKM[ind]
  if(i == 1){
    rownames(FGFR2_exon) = t_exon$exon[ind]
  } else {
    if(rownames(FGFR2_exon) == t_exon$exon[ind]){
      print(i)
    } else{
      break
    }
  }
}

rm(t, t_exon, ind, i, folders, files)

# matching to "ENST00000358487"
exon_list = data.frame(exon_list, exon_num_ENST00000358487 = c(18, 17.5, 17:8, 7.5, 7:2, 1.5, 1),
                       stringsAsFactors = F)
colnames(FGFR2_exon) = samples_exon$id

save.image("~/FGFR/Daniel/R/Nature_figures/data/TCGA/FGFR2_exon_expression_TCGA.RData")


