################################# 
# RSEM isoform data
################################# 
rsem_iso_FGFR2 = function(data_path, samp_label){
  # data_path: path for rsem (ex, "~/HPC/harris/Daniel_FGFR2/cell_lines/rsem/")
  
  FGFR2_ISO = read.table("~/HPC/turing/FGFR/Daniel/ENSEMBL_FGFR2_ISO.txt", 
                         header = T, sep = "\t", stringsAsFactors = F, quote = "")
  folders = list.files(data_path)
  files = folders[grep("isoforms", folders)]
  
  files = paste(data_path, files, sep = "")
  rsem_isoforms = do.call("cbind", sapply(X = files, FUN = read.table, header = T, sep = "\t", stringsAsFactors = F,
                                          colClasses = c(rep("NULL", 6), "numeric", "NULL")))
  colnames(rsem_isoforms) = samp_label
  isoforms = read.table(files[1], header = T, sep = "\t", stringsAsFactors = F)[,1:4]
  
  rsem_FGFR2 = rsem_isoforms[match(FGFR2_ISO$ENSEMBL_ID, isoforms$transcript_id),]
  isoforms_FGFR2 = isoforms[match(FGFR2_ISO$ENSEMBL_ID, isoforms$transcript_id),]
  rownames(rsem_FGFR2) = FGFR2_ISO$ENSEMBL_ID
  return(list(rsem_fpkm = rsem_FGFR2, iso_info = isoforms_FGFR2))
}


################################# 
# RSEM gene data
################################# 
rsem_gene_target = function(target_ENSEMBL_ID, data_path, samp_label){
  # data_path: path for rsem (ex, "~/HPC/harris/Daniel_FGFR2/cell_lines/rsem/")
  
  folders = list.files(data_path)
  files = folders[grep("genes", folders)]
  
  files = paste(data_path, files, sep = "")
  rsem_genes = do.call("cbind", sapply(X = files, FUN = read.table, header = T, sep = "\t", stringsAsFactors = F,
                                       colClasses = c(rep("NULL", 6), "numeric")))
  
  colnames(rsem_genes) = samp_label
  genes = read.table(files[1], header = T, sep = "\t", stringsAsFactors = F)[,1:4]
  rsem_target = rsem_genes[which(genes$gene_id == target_ENSEMBL_ID),]
  genes_target = genes[which(genes$gene_id == target_ENSEMBL_ID),]
  
  return(list(rsem_fpkm = rsem_target, gene_info = genes_target))
}



################################# 
# DEXseq exon data
################################# 
dexseq_exon_FGFR2 = function(data_path, samp_label){
  # data_path: path for rsem (ex, "~/HPC/harris/Daniel_FGFR2/cell_lines/dexseq/")
  
  library(DESeq2)
  # FGFR2 exon
  flat_fgfr2 = read.table("~/HPC/turing/Resources/human_GRCh38/DEXSeq/FGFR2_dexseq.gff", 
                          header = F, sep = "\t", stringsAsFactors = F)
  
  files = paste(data_path, list.files(data_path), sep = "")
  sampleTable = data.frame(sampleNames = samp_label, fileName = list.files(data_path), condition = samp_label)
  dxd = DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = data_path, design = ~condition)
  dxd = estimateSizeFactors(dxd)
  norm_exon <- counts(dxd, normalized=TRUE)
  colnames(norm_exon) = samp_label
  
  norm_exon_fgfr2 = norm_exon[grep("ENSG00000066468", rownames(norm_exon)),]
  # remove exon part 12 which is not annotated in UCSC
  norm_exon_fgfr2 = norm_exon_fgfr2[!grepl(":012", rownames(norm_exon_fgfr2)),]
  
  # merge the overlapping exon parts
  ind = unique(flat_fgfr2$V6)
  fgfr2_merge = c()
  fgfr2_len = matrix(NA, length(ind), 3)
  for(i in 1:length(ind)){
    t = norm_exon_fgfr2[which(flat_fgfr2[,6] == ind[i]),]
    t_flat = flat_fgfr2[which(flat_fgfr2[,6] == ind[i]),]
    
    if (is.vector(t)){
      fgfr2_merge = cbind(fgfr2_merge, t)
      fgfr2_len[i,1:2] = unlist(c(t_flat[4], t_flat[5]))
      fgfr2_len[i,3] = as.numeric(t_flat[5] - t_flat[4])
      
    } else {
      fgfr2_merge = cbind(fgfr2_merge, colSums(t))
      fgfr2_len[i,1:2] = c(t_flat[1,4], t_flat[nrow(t_flat),5])
      fgfr2_len[i,3] =as.numeric(t_flat[nrow(t_flat),5] - t_flat[1,4])
    }
  }
  fgfr2_merge = t(fgfr2_merge)
  fgfr2_merge_len = fgfr2_merge / fgfr2_len[,3]
  fgfr2_RPKM = t(t(fgfr2_merge*(10^9)) / (colSums(norm_exon))) / fgfr2_len[,3]
  fgfr2_merge_per = t(t(fgfr2_merge_len) / colSums(fgfr2_merge_len))
  rownames(fgfr2_merge_per) = sprintf("exon_part:%d (chr10:%d-%d)", 1:length(unique(flat_fgfr2$V6)), fgfr2_len[,1], fgfr2_len[,2])
  return(list(merge = fgfr2_merge, merge_len = fgfr2_len, merge_per = fgfr2_merge_per, RPKM = fgfr2_RPKM))
}


################################# 
# STAR-Fusion
################################# 
fusion_target = function(target, data_path, samp_label){
  # data_path: path for rsem (ex, "~/HPC/harris/Daniel_FGFR2/cell_lines/star_fusion/")
  
  # folders = list.files(data_path)
  # files = sapply(folders, FUN = function(X) {paste(data_path, X, "/FusionInspector-validate/finspector.FusionInspector.fusions.abridged.tsv.annotated.coding_effect", sep ="")})
  # #files = sapply(folders, FUN = function(X) {paste(data_path, X, "/star-fusion.fusion_predictions.abridged.tsv", sep ="")})
  # 
  # 
  # FGFR_fus_FI = c()
  # for(i in 1:length(files)){
  #   fus_list = read.table(file = files[i], sep = "\t", header = T, stringsAsFactors = F, quote = "", comment.char = "")
  #   if((sum(grepl(pattern = "FGFR", fus_list[,1]))) !=0){
  #     FGFR_fus_FI = rbind(FGFR_fus_FI, data.frame(fus_list[grep("FGFR", fus_list$X.FusionName),], sample = samp_label[i]))
  #   }
  # }
  # 
  # 
  
  folders = list.files(data_path)
  FGFR_fus_FI = c()
  for(i in 1:length(folders)){
    tmp_folder = paste(data_path, folders[i], sep = "")
    if("FusionInspector-validate" %in% list.files(tmp_folder)){
      tmp_path = paste(data_path, folders[i], '/FusionInspector-validate/finspector.FusionInspector.fusions.abridged.tsv.annotated.coding_effect', sep = "")
      fus_list = read.table(tmp_path, sep = "\t", header = T, stringsAsFactors = F, quote = "", comment.char = "")
      FGFR_fus_FI = rbind(FGFR_fus_FI, data.frame(fus_list, sample = samp_label[i]))
    }
  }
  FGFR_fus_FI = FGFR_fus_FI[grep(target, FGFR_fus_FI$X.FusionName),]
  
  # Filtering FFPM
  FGFR_fus_filter = FGFR_fus_FI[which(FGFR_fus_FI$FFPM>=0.1),]
  
  fus_uniq = unique(FGFR_fus_filter$X.FusionName)
  fus_mat = matrix(0, length(fus_uniq), length(samp_label))
  for(i in 1:length(fus_uniq)){
    t = FGFR_fus_filter[which(FGFR_fus_filter$X.FusionName == fus_uniq[i]),]
    t_samp = unique(t$sample)
    for(j in 1:length(t_samp)){
      fus_mat[i, samp_label==t_samp[j]] = sum(t$FFPM[t$sample == t_samp[j]])
    }
  }
  colnames(fus_mat) = samp_label
  rownames(fus_mat) = fus_uniq
  
  return(fus_mat)
}





################################# 
# Heatmap generation
################################# 
# draw_heatmap = function(FGFR2_gene, FGFR2_iso, FGFR2_fus, annot){
#   library(ComplexHeatmap)
#   ha = HeatmapAnnotation(Gene_counts = norm_gene_fgfr2_mc[ind_sort],
#                          col = list(Gene_counts = colorRamp2(c(-3,0,3), c("blue", "white", "red"))))
#   
#   
#   ord = names(order(gene_mat, decreasing = F))
#   FGFR2_gene = gene_mat[match(ord, names(gene_mat))]
#   FGFR2_iso = iso_mat[, match(ord, colnames(iso_mat))]
#   FGFR2_fus = fus_mat[, match(ord, colnames(fus_mat))]
# 
# 
# }

