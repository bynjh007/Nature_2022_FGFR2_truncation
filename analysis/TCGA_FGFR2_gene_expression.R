## data loading
# Exclude FFEP data
library(stringr)
folders=list.files("~/TCGA/data/")
folders=folders[folders!="FPPP"]
samples=list()
FGFR2_data=list()
for(i in  1:length(folders)){
  tmp_folder=paste("~/TCGA/data/", folders[i], "/20160128", sep = "")
  setwd(tmp_folder)
  tmp_files=list.files(tmp_folder)
  x=tmp_files[grepl("RSEM_genes_normalized", tmp_files) & grepl("Level_3[.]", tmp_files) & !grepl("md5", tmp_files)]
  untar(x[1])
  setwd(str_split_fixed(string = x, pattern = "[.]tar.gz", 2)[,1][1])
  t_list=list.files()
  t_data=read.table(t_list[grep("RSEM", t_list)], header = T, sep = "\t", stringsAsFactors = F)
  t_data=t_data[-1,]; rownames(t_data)=t_data[,1]; t_data=t_data[,-1]
  
  # normalized counts
  ind=grep("FGFR2", rownames(t_data))
  t_fgfr2=t_data[ind,]
  t_fgfr2=as.numeric(t_fgfr2)
  names(t_fgfr2)=colnames(t_data)
  FGFR2_data[[i]]=t_fgfr2
  t_type=cbind(folders[i], str_split_fixed(names(t_fgfr2), "[.]", 5)[,4], names(t_fgfr2))
  samples[[i]]=t_type
  
  unlink(paste("~/TCGA/data/", folders[i], "/20160128/", str_split_fixed(string = x, pattern = "[.]tar.gz", 2)[,1][1], sep = ""), recursive = T)
  rm(t_data)
  
  print(i)
}
names(FGFR2_data)=folders

FGFR2_data=unlist(FGFR2_data)
samples=do.call(rbind, samples)
colnames(samples)=c("Type", "TN", "sample")
names(FGFR2_data)=samples[,3]

rm(tmp_folder, tmp_files, x, t_list, ind, t_fgfr2, t_type, i)

save.image("~/FGFR/Daniel/R/Nature_figures/data/TCGA/FGFR2_gene_expression_TCGA.RData")
