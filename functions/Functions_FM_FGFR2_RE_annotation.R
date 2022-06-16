####################################################
# library
####################################################
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(org.Hs.eg.db)
library(RColorBrewer)
library(GenomicFeatures)
library(Biostrings)
library(dplyr)
options(stringsAsFactors = F)


##############################################################################
# Function 
##############################################################################

# Function for generating transcriptome information from genome annotation database
generating_genome_info = function(txdb, id_sym_tab, genome){
  # id_sym_tab: col1-TXNAME, col2-SYMBOL
  
  # list of transcript in the DB
  trx = transcripts(txdb)
  gene_exon = exonsBy(txdb, by = "gene")
  
  # intron and exon information
  trx_intron = intronsByTranscript(txdb)
  trx_exon = exonsBy(txdb)
  
  # TXID is same with index of the "trx", "trx_gene"
  trx_gene = AnnotationDbi::select(txdb, key = as.character(trx$tx_id), keytype = "TXID", columns = c("TXNAME", "GENEID"))
  exon_info = AnnotationDbi::select(txdb, key = as.character(trx$tx_id), keytype = "TXID", columns = c("EXONID", "EXONNAME", "TXID",  "TXNAME", "GENEID"))
  
  # adding SYMBOL and STRAND information
  trx_gene$SYMBOL = id_sym_tab[match(trx_gene$TXNAME, id_sym_tab[,1]), 2]
  trx_gene$STRAND = as.character(strand(trx))
  
  # protein coding region in each transcript (GRange object)
  G_cds = cdsBy(txdb, by="tx")
  
  # list of non- scaffold and mitochondrial chromosome
  ind_canonical = which(!(grepl("_|M", as.character(seqnames(trx)))))
  G_cds = G_cds[which(as.numeric(names(G_cds)) %in% ind_canonical)]
  
  # DNA and AA sequence
  cds_seqs = extractTranscriptSeqs(genome, G_cds)
  cds_seqs_aa = suppressWarnings(Biostrings::translate(cds_seqs))
  
  # filtering out the non-canonical transcript (stop codon in the middle of protein sequence which is caused by non-canonical translation start site)
  ind_canonical = which(!(grepl("\\*", substr(as.character(cds_seqs_aa), 1, nchar(as.character(cds_seqs_aa))-1))))
  G_cds = G_cds[ind_canonical]
  cds_seqs = cds_seqs[ind_canonical]
  cds_seqs_aa = cds_seqs_aa[ind_canonical]
  
  # adding the length of the AA and DNA of each transcript and sorting based on 1) AA and 2) DNA length
  trx_gene_AA = trx_gene %>% mutate(aa_length = width(cds_seqs_aa)[match(TXID, names(cds_seqs_aa))],
                                    trx_length = width(trx)) %>% 
    arrange(desc(aa_length), desc(trx_length))
  
  # selecting the major transcript for each gene based on the length of AA (longest transcript)
  # trx with non-canonical protein sequence is also annotated as NA
  total_gene_list = unique(trx_gene_AA$SYMBOL)
  trx_gene_AA_major = trx_gene_AA[match(total_gene_list, trx_gene_AA$SYMBOL),] %>% arrange(TXID)
  trx_gene_AA_major = trx_gene_AA_major %>% arrange(TXID)
  
  # GRange object for the major transcript
  trx_major = trx[trx_gene_AA_major$TXID]
  
  # re-ordering trx table based on the TXID
  trx_gene_AA = trx_gene_AA %>% arrange(TXID)
  
  trx_info = list(trx = trx, gene_exon = gene_exon, 
                  exon_info = exon_info, trx_intron = trx_intron, trx_exon = trx_exon, 
                  trx_major = trx_major, trx_gene_AA = trx_gene_AA, 
                  trx_gene_AA_major = trx_gene_AA_major, 
                  G_cds = G_cds, cds_seqs = cds_seqs, cds_seqs_aa = cds_seqs_aa)
  
  return(trx_info)
  
}

# Function for annotating breakpoints in FGFR2 and its partners
brk_annot = function(ref_trx, index, target_region, up_down){
  target_gene_trx = ref_trx[index,]
  target_gene_trx = target_gene_trx %>% mutate(Protein_coding = ifelse(!is.na(aa_length), TRUE, FALSE))
  
  if(nrow(target_gene_trx)>1){
    stop("more than 1")
  }
  
  ################################################## 
  # build the GRange object for the target transcript
  ##################################################
  ind = as.numeric(target_gene_trx$TXID)
  t_exon = hg19_info$trx_exon[[ind]]
  t_intron = hg19_info$trx_intron[[ind]]
  t_str = unique(as.character(strand(t_exon)))
  
  if(length(t_intron)!=0){
    if(t_str == "-"){
      t_intron = rev(t_intron)
    }
    df1 = data.frame(data.frame(t_exon)[, c(1,2,3,5)], id = paste("Exon", 1:length(t_exon)))
    df2 = rbind(data.frame(data.frame(t_intron)[, c(1,2,3,5)], id = paste("Intron", 1:length(t_intron))), NA)
    df_comb = (gdata::interleave(df1, df2))
    df_comb = df_comb[-nrow(df_comb),]
    G_ref = makeGRangesFromDataFrame(df_comb)
    G_ref$id = df_comb$id
  } else {
    df1 = data.frame(data.frame(t_exon)[, c(1,2,3,5)], id = paste("Exon", 1:length(t_exon)))
    G_ref = makeGRangesFromDataFrame(df1)
    G_ref$id = df1$id
  }
  
  #######################################################################################
  # Find the location of the breakpoint and downstream exons within the target transcript
  #######################################################################################
  # Compare the Ranges between the target transcript and target region (breakpoint)
  temp_ovp = findOverlaps(G_ref, target_region)
  loc_in_trx = G_ref[queryHits(temp_ovp)]
  ind = queryHits(temp_ovp)
  
  if(target_gene_trx$Protein_coding == T){
    G_cds_target = hg19_info$G_cds[[which(names(hg19_info$G_cds) == target_gene_trx$TXID)]]
    target_cds_seq_AA = as.character(hg19_info$cds_seqs_aa[names(hg19_info$cds_seqs_aa) == target_gene_trx$TXID][[1]])
    
    if(up_down == "up"){
      # breakpoint is in intronic region
      cds_target_region = G_cds_target[queryHits(findOverlaps(query = G_cds_target, subject = G_ref[1:ind]))]
      if(grepl("Exon", loc_in_trx$id)){
        if(length(queryHits(findOverlaps(query = G_cds_target, subject = target_region)))!=0){ # breakpoint in UTR region in the exon covering both UTR and coding region
          if(as.character(strand(G_ref))[1] == "+"){
            end(ranges(cds_target_region[length(cds_target_region)])) = end(ranges(target_region))
          } else {
            start(ranges(cds_target_region[length(cds_target_region)])) = start(ranges(target_region))
          }
        } else {
          cds_target_region = NULL
        }
      }
    } else { # downstream of the breakpoint
      cds_target_region = G_cds_target[queryHits(findOverlaps(query = G_cds_target, subject = G_ref[ind:length(G_ref)]))]
      if(grepl("Exon", loc_in_trx$id)){
        if(length(queryHits(findOverlaps(query = G_cds_target, subject = target_region)))!=0){ # breakpoint in UTR region in the exon covering both UTR and coding region
          if(as.character(strand(G_ref))[1] == "+"){
            start(ranges(cds_target_region[1])) = start(ranges(target_region))
          } else {
            end(ranges(cds_target_region[1])) = end(ranges(target_region))
          }
        } else {
          cds_target_region = NULL
        }
      }
    }
    
    # breakpoint at UTR region
    if(length(cds_target_region)!=0){
      # AA sequence before (if target is upstream) or after (if target is downstream) the breakpoint
      cds_target_region_DNA = do.call(xscat, getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, cds_target_region))
      cds_target_region_DNA_res = nchar(as.character(cds_target_region_DNA)) %% 3
      
      options(warn = -1)
      if(cds_target_region_DNA_res == 0 | up_down == "up"){
        cds_target_region_AA = translate(cds_target_region_DNA)
      } else {# if downstream of a breakpoint
        cds_target_region_AA = translate(cds_target_region_DNA[-c(1:cds_target_region_DNA_res)])
      }
      options(warn = 0)
      
      if(grepl(cds_target_region_AA, target_cds_seq_AA, fixed = T)){
        annot_target_region = data.frame("base_left" = cds_target_region_DNA_res, "AA_seq" = as.character(cds_target_region_AA), "AA_length" = nchar(cds_target_region_AA))
      } else {
        annot_target_region = data.frame("base_left" = cds_target_region_DNA_res, "AA_seq" = "not_consistent with reference", "AA_length" = nchar(cds_target_region_AA))
      }
      
    } else {
      annot_target_region = data.frame("base_left" = NA, "AA_seq" = "UTR", "AA_length" = NA)
    }
    
  } else { # breakpoint at non-coding gene
    annot_target_region = data.frame("base_left" = NA, "AA_seq" = "non-coding RNA", "AA_length" = NA)
  }
  
  # output
  if(as.character(loc_in_trx$id) != "NA"){
    t_loc = apply(data.frame(loc_in_trx), 1, paste, collapse = ":")
  } else {
    t_loc = "out_of_trx"
  }
  
  target_annot = data.frame(target_region)[,1:2] %>% mutate(gene_id = target_gene_trx$GENEID,
                                                            gene_symbol = target_gene_trx$SYMBOL,
                                                            trx_name = target_gene_trx$TXNAME,
                                                            aa_length = target_gene_trx$aa_length,
                                                            trx_length = target_gene_trx$trx_length,
                                                            loc_in_trx = t_loc,
                                                            base_left = annot_target_region$base_left,
                                                            AA_seq = annot_target_region$AA_seq,
                                                            AA_length = annot_target_region$AA_length)
  
  return(target_annot)
}
    
# Function for FGFR2 RE annotations
brk_annot_framework = function(G_up, G_down){
  # G_up: GRanges object for upstream breakpoint (seqnames, ranges, strand)
  # G_down: GRanges object for downstream breakpoint (seqnames, ranges, strand)
 
  # find the overlapping transcript
  index_ovp = list(up = queryHits(findOverlaps(hg19_info$trx_major, G_up)),
                   down = queryHits(findOverlaps(hg19_info$trx_major, G_down)))
  
  G_target = list(up = G_up, down = G_down)
  # if there is no overlap for the representative transcript
  # 1 is upstream and 2 is downstream
  target_annot_all = c()
  for(i in 1:2){
    
    up_down = ifelse(i == 1, "up", "down")
    
    # # if there is no overlap for the representative transcript, find the target region in the all the transcripts
    if(length(index_ovp[[i]]) == 0){
      index = queryHits(findOverlaps(hg19_info$trx, G_target[[i]]))
      if(length(index) !=0){
        index = hg19_info$trx_gene_AA[index,] %>% arrange(desc(aa_length), desc(trx_length)) %>% pull(TXID)
        target_annot = brk_annot(ref_trx = hg19_info$trx_gene_AA, index = index[1], target_region = G_target[[i]], up_down = up_down)
      } else {
        target_annot = data.frame(G_target) %>% mutate(gene_id = hg19_info$target_gene_trx$GENEID,
                                                       gene_symbol = target_gene_trx$SYMBOL,
                                                       trx_name = target_gene_trx$TXNAME,
                                                       loc_in_trx = t_loc,
                                                       down_exon = as.character(down_exon$id[1]))
      }
    
    } else {
      index = index_ovp[[i]][(hg19_info$trx_gene_AA_major[index_ovp[[i]],] %>% mutate(ind = 1:length(index_ovp[[i]])) %>% arrange(desc(aa_length), desc(trx_length)) %>% pull(ind))]
      target_annot = brk_annot(ref_trx = hg19_info$trx_gene_AA_major, index = index[1], target_region = G_target[[i]], up_down = up_down)
    }
    
    target_annot_all = rbind(target_annot_all, target_annot)
  }
  target_annot_out = data.frame(up = target_annot_all[1,], down = target_annot_all[2,])
  return(target_annot_out)
}


# function to define in-frame or out-of-frame from purple annotation
RE_type = function(SV_annot){
  # in-frame: both up and downstream is inframe
    base_left_sum = SV_annot$up.base_left + SV_annot$down.base_left
  # up and down are FGFR2
  if((SV_annot$up.gene_symbol %in% "FGFR2") & (SV_annot$down.gene_symbol %in% "FGFR2")){
    RE = "internal"
    # one of the partner is intergenic
  } else if(is.na(SV_annot$up.gene_symbol) | is.na(SV_annot$down.gene_symbol)){
    RE = "intergenic"
  } else if((base_left_sum %% 3) == 0 & !is.na(base_left_sum)){
    
    RE = "in-frame"
  } else {
    RE = "in-strand (frame-unknown)"
  }
  return(RE)
}


####################################################
# UCSC hg19 ID mapping
####################################################
hg_id_19_ENST = read.table("~/FGFR/Daniel/R/Nature_figures/resouces/UCSC.hg.19_ID_to_SYMBOL.txt", header = T, stringsAsFactors = F, sep = "\t", comment.char = "")
hg_id_19 = hg_id_19_ENST[,c(1,4)]
hg19_info = generating_genome_info(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                   id_sym_tab = hg_id_19, genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens)



