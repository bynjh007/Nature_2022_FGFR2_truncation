# In the Debio Server
path_folder = "path/bam" # location of bam directory from the pipeline ("path" should be changed)
output = "path/FGFR2_output.RData" # output file after extracting FGFR2 information ("path" should be changed)

#######################################################
# library loading
#######################################################
library(dplyr)
library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes.gr <- genes(txdb)
symbols.eg <- select(Homo.sapiens,keys=elementMetadata(genes.gr)$gene_id,columns="SYMBOL",keytype="GENEID")
elementMetadata(genes.gr)$symbol <- symbols.eg$SYMBOL
exons.gr = exonicParts(txdb)
exons.gr_FGFR2 = data.frame(subset(exons.gr, any(gene_id %in% "2263")))

########################################################
# junction count function
########################################################
junction_count = function(chim.out, target, grHs, min_num){
  # range - character vector (chr, start, end, strand)
  # grHs - reference gene annotation
  
  tmp_data = read.table(chim.out, header = F, stringsAsFactors = F, sep = "\t")
  colnames(tmp_data)[1:9] = c("chr_donorA", "brkpt_donorA", "strand_donorA", 
                         "chr_acceptorB", "brkpt_acceptorB", "strand_acceptorB",
                         "junction_type", "repeat_left_lenA", "repeat_right_lenB")
  
  #############################
  # spanning reads analysis
  #############################
  # donor
  target_donor_junc = tmp_data %>% filter(chr_donorA == target[1], brkpt_donorA >= as.numeric(target[2]), brkpt_donorA <= as.numeric(target[3]), strand_donorA == target[4])
  target_donor_junc_summary = target_donor_junc %>% group_by(chr_donorA, brkpt_donorA, strand_donorA, chr_acceptorB, brkpt_acceptorB, strand_acceptorB, junction_type) %>% summarise(n = n())
  # select the junction reads with spanning reads >=5
  target_donor_junc_summary = target_donor_junc_summary %>% filter(n>=min_num, junction_type != -1)
  
  # acceptor
  target_acceptor_junc = tmp_data %>% filter(chr_acceptorB == target[1], brkpt_acceptorB >= as.numeric(target[2]), brkpt_acceptorB <= as.numeric(target[3]), strand_acceptorB == target[4])
  target_acceptor_junc_summary = target_acceptor_junc %>% group_by(chr_donorA, brkpt_donorA, strand_donorA, chr_acceptorB, brkpt_acceptorB, strand_acceptorB, junction_type) %>% summarise(n = n())
  # select the junction with spanning reads >=5
  target_acceptor_junc_summary = target_acceptor_junc_summary %>% filter(n>=min_num, junction_type != -1)
  
  # merge the two results
  target_merge = data.frame(rbind(target_donor_junc_summary, target_acceptor_junc_summary))
  target_merge = target_merge %>% 
    mutate(brkpt_donorA = replace(brkpt_donorA, strand_donorA == "-", brkpt_donorA[strand_donorA == "-"]+1)) %>%
    mutate(brkpt_donorA = replace(brkpt_donorA, strand_donorA == "+", brkpt_donorA[strand_donorA == "+"]-1)) %>%
    mutate(brkpt_acceptorB = replace(brkpt_acceptorB, strand_acceptorB == "-", brkpt_acceptorB[strand_acceptorB == "-"]-1)) %>%
    mutate(brkpt_acceptorB = replace(brkpt_acceptorB, strand_acceptorB == "+", brkpt_acceptorB[strand_acceptorB == "+"]+1))
  
  # encompassing reads
  target_encomp = tmp_data %>% filter(junction_type == -1) 
  
  #####################
  # gene annotation
  #####################
  if(nrow(target_merge)!=0){
    #KnownGene1/2
    g1 <- NA; g2 <- NA
    
    #creating objects
    fusionList <- c()
    cat("\n")
    for(i in 1:nrow(target_merge)){
      cat(".")
      print(i)
      
      strand1 <- as.character(target_merge$strand_donorA[i])
      strand2 <- as.character(target_merge$strand_acceptorB[i])
      #spanning reads		 
      seed1 <- target_merge$n[i]
      #junction type
      junc <- target_merge$junction_type[i]
      #pos.tmp <- strsplit(as.character(target_merge$chrs.fusion[i]), "-")
      
      #defining start1
      #coverage1.tmp <- strsplit(as.character(target_merge$coverage1[i]), " ")
      start1 <- target_merge$brkpt_donorA[i] - 30
      end1 <- target_merge$brkpt_donorA[i]
      chr1 <- as.character(target_merge$chr_donorA[i])
      #defining start2
      #coverage2.tmp <- strsplit(as.character(target_merge$coverage2[i]), " ")
      start2 <- target_merge$brkpt_acceptorB[i]
      end2 <- target_merge$brkpt_acceptorB[i] + 30
      chr2 <- as.character(target_merge$chr_acceptorB[i])
      #defining gene1
      grG1 <-  GRanges(seqnames = chr1, ranges = IRanges(start = start1, end= end1))
      grG2 <-  GRanges(seqnames = chr2, ranges = IRanges(start = start2, end= end2))
      tmpG1 <- findOverlaps(grG1, grHs, type = "any", select = "first", ignore.strand = T)
      if(!is.na(tmpG1)){
        g1 <- elementMetadata(grHs[tmpG1])$symbol	            	
      }else{
        g1 <- paste(seqnames(grG1), paste(start(grG1),end(grG1), sep="-"),sep=":")
      }
      tmpG2 <- findOverlaps(grG2, grHs, type = "any", select = "first", ignore.strand = T)
      if(!is.na(tmpG2)){
        g2 <- elementMetadata(grHs[tmpG2])$symbol	            	
      }else{
        g2 <- paste(seqnames(grG2), paste(start(grG2),end(grG2), sep="-"),sep=":")
      }
      
      
      gr1 <- GRanges(seqnames = chr1,
                     ranges = IRanges(start = start1, end= end1),
                     strand = strand1,
                     KnownGene = as.character(g1))
      gr2 <- GRanges(seqnames = chr2,
                     ranges = IRanges(start = start2, end= end2),
                     strand = strand2,
                     KnownGene = as.character(g2))
      grl <- data.frame(GRangesList("gene1" = gr1, "gene2" = gr2))
      
      
      # encompassing reads
      if(target_merge$strand_donorA[i] == "-"){
        t = target_encomp %>% filter(chr_donorA == target_merge$chr_donorA[i],
                                     brkpt_donorA >= target_merge$brkpt_donorA[i],
                                     brkpt_donorA <= target_merge$brkpt_donorA[i]+50)
        if(target_merge$strand_acceptorB[i] == "-"){
          t = t %>% filter(chr_acceptorB == target_merge$chr_acceptorB[i],
                           brkpt_acceptorB <= target_merge$brkpt_acceptorB[i],
                           brkpt_acceptorB >= target_merge$brkpt_acceptorB[i]-50)
        } else { # target_merge$strand_acceptorB[i] == "+"
          t = t %>% filter(chr_acceptorB == target_merge$chr_acceptorB[i],
                           brkpt_acceptorB >= target_merge$brkpt_acceptorB[i],
                           brkpt_acceptorB <= target_merge$brkpt_acceptorB[i]+50)
        }
        
      } else {
        # target_merge$strand_donorA[i] == "+"
        t = target_encomp %>% filter(chr_donorA == target_merge$chr_donorA[i],
                                     brkpt_donorA <= target_merge$brkpt_donorA[i],
                                     brkpt_donorA >= target_merge$brkpt_donorA[i]-50)
        if(target_merge$strand_acceptorB[i] == "+"){
          t = t %>% filter(chr_acceptorB == target_merge$chr_acceptorB[i],
                           brkpt_acceptorB >= target_merge$brkpt_acceptorB[i],
                           brkpt_acceptorB <= target_merge$brkpt_acceptorB[i]+50)
        } else {  # target_merge$strand_acceptorB[i] == "-"
          t = t %>% filter(chr_acceptorB == target_merge$chr_acceptorB[i],
                           brkpt_acceptorB <= target_merge$brkpt_acceptorB[i],
                           brkpt_acceptorB >= target_merge$brkpt_acceptorB[i]-50)
        }
      }
      
      
      df = data.frame(Fusion_gene = paste(g1,g2, sep="->"),
                      gene_donor = grl$KnownGene[1], target_merge[i, c(1:3)],
                      gene_acceptor = grl$KnownGene[2], target_merge[i, c(4:6)],
                      junction_type = junc,
                      n_spanning = target_merge$n[i], n_encompassing = nrow(t))
      
      fusionList <- rbind(fusionList, df)
    }
    cat("\n")
  } else {
    fusionList = NA
  }
  return(fusionList)
}

########################################################
# Chimera alignment
########################################################

path_file = list.files(path_folder)
path_file = path_file[!grepl("graft|summary|done|_db", path_file)]

FGFR2_fusion = c()
for(i in 1:length(path_file)){
  setwd(paste(path_folder, "/", path_file[i], sep = ""))
  options(warn=0)
  t_fusion = junction_count(chim.out = "Chimeric.out.junction", target = c("chr10", "121478334", "121598458", "-"), 
                            grHs = genes.gr, min_num = 1)
  if(!is.na(t_fusion)[1]){
    t_fusion = t_fusion[!grepl("chrM", t_fusion$Fusion_gene),]
    t_fusion = t_fusion[!grepl("FGFR2->FGFR2", t_fusion$Fusion_gene),]
    if(nrow(t_fusion)!=0){
      FGFR2_fusion = rbind(FGFR2_fusion, data.frame(t_fusion, sample = path_file[i]))
    }
  }
  print(i)
}


########################################################
# junction analysis
########################################################
E17_junc_list = c()
for (x in path_file){
  tmp_data = read.table(paste(path_folder, "/", x, "/SJ.out.tab", sep = ""), header = F, stringsAsFactors = F, sep = "\t")
  t_E17_junc = tmp_data[which(tmp_data[,1] == "chr10" & tmp_data[,3] == 121483697),1:3]
  E17_junc_list = rbind(E17_junc_list, t_E17_junc)
  print(x)
}
E17_junc_list = unique(E17_junc_list)
E17_junc_list = E17_junc_list[order(E17_junc_list[,2]),]

E17_junc_cnt = matrix(NA, nrow(E17_junc_list), length(path_file))
for(i in 1:length(path_file)){
  tmp_data = read.table(paste(path_folder, "/", path_file[i], "/SJ.out.tab", sep = ""), header = F, stringsAsFactors = F, sep = "\t")
  tmp_coord = apply(tmp_data[,1:3], 1, paste, collapse = "_")
  ind = match(apply(E17_junc_list, 1, paste, collapse = "_"), tmp_coord)
  E17_junc_cnt[,i] = tmp_data[ind,7]
  print(i)
}
E17_junc_cnt[is.na(E17_junc_cnt)] = 0
E17_junc_list = cbind.data.frame(E17_junc_list, rowSums(E17_junc_cnt))
colnames(E17_junc_list) = c("chr", "loc_1", "loc_2", "sum_cnt")

# All the above data was analyzed in Debio server
save.image("~/FGFR/Daniel/R/Nature_figures/data/Debio/Debio_FGFR2_junction_processing.RData")

