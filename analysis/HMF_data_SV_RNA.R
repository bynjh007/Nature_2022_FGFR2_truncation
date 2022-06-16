library(GenomicFeatures)
library(dplyr)
options(stringsAsFactors = F)

hg_id_19_ENST = read.table("~/FGFR/Daniel/R/Nature_figures/resources/UCSC.hg.19_ID_to_SYMBOL.txt", header = T, stringsAsFactors = F, sep = "\t", comment.char = "")
hg_id_19 = hg_id_19_ENST[,c(1,4)]

hg_id_38_ENST = read.table("~/FGFR/Daniel/R/Nature_figures/resources/GRCh38_genecode_v32_trx.txt", header = T, stringsAsFactors = F, sep = "\t")
hg38_ensembl = biomaRt::useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = "100")
hg38_entrezid = biomaRt::getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), mart = hg38_ensembl)
hg_id_38_ENST$GENEID = hg38_entrezid$entrezgene_id[match(stringr::str_split_fixed(hg_id_38_ENST$ENSG_ID, "[.]", 2)[,1],
                                                   hg38_entrezid$ensembl_gene_id)]
hg_id_38 = hg_id_38_ENST[, c(2,4)]

# hg19 and hg38 - importing the gene / transcript information
##############################################################################
# function for generating required information from genome annotation database
##############################################################################

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


##############################################################
# generation of information required for breakpoint annotation
##############################################################
hg19_info = generating_genome_info(txdb = TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                   id_sym_tab = hg_id_19, genome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens)
hg38_info = generating_genome_info(txdb = AnnotationDbi::loadDb(file = '~/FGFR/Daniel/R/Nature_figures/resources/txdb.GRCh38_genecode_v32_CTAT_lib.sqlite'), 
                                   id_sym_tab = hg_id_38, genome = BSgenome.Hsapiens.UCSC.hg38::Hsapiens)

##############################################################
# Functions
##############################################################

# function for breakpoint annotation
brk_annot_framework = function(G_up, up_orient, G_down, down_orient, genome_info, BSgenome){
  ################################################## 
  # find the overlapping transcript
  ##################################################
  index_ovp = list(up = queryHits(suppressWarnings(findOverlaps(genome_info$trx_major, G_up))),
                   down = queryHits(suppressWarnings(findOverlaps(genome_info$trx_major, G_down))))
  
  G_target = list(up = G_up, down = G_down)
  
  # if there is no overlap for the representative transcript
  # 1 is upstream and 2 is downstream
  target_annot_all = c()
  for(i in 1:2){
    
    if(i == 1){
      target_orient = up_orient
      up_down = "up"
    } else {
      target_orient = down_orient
      up_down = "down"
    }
    
    if(length(index_ovp[[i]]) == 0){
      # find the target region in the all the transcripts
      
      index = queryHits(suppressWarnings(findOverlaps(genome_info$trx, G_target[[i]])))
      if(length(index) !=0){
        index = genome_info$trx_gene_AA[index,] %>% arrange(desc(aa_length), desc(trx_length)) %>% pull(TXID)
        target_annot = brk_annot(ref_trx = genome_info$trx_gene_AA, index = index[1], target_region = G_target[[i]], target_orient = target_orient, 
                                 up_down = up_down, genome_info = genome_info, BSgenome = BSgenome)
        
        # no overlapping and SGL
      } else if(as.character(seqnames(G_target[[i]])) == "chr-1"){
        target_annot = data.frame(G_target[[i]])[,1:2] %>% mutate(gene_id = NA,
                                                                  gene_symbol = NA,
                                                                  trx_name = NA,
                                                                  aa_length = NA,
                                                                  trx_length = NA,
                                                                  loc_in_trx = NA,
                                                                  down_exon = NA,
                                                                  down_exon_frame = "no partner",
                                                                  RE_category = "SGL",
                                                                  orientation = target_orient)
        
      } else {
        target_annot = data.frame(G_target[[i]])[,1:2] %>% mutate(gene_id = NA,
                                                                  gene_symbol = NA,
                                                                  trx_name = NA,
                                                                  aa_length = NA,
                                                                  trx_length = NA,
                                                                  loc_in_trx = NA,
                                                                  down_exon = NA,
                                                                  down_exon_frame = "no partner",
                                                                  RE_category = "intergenic",
                                                                  orientation = target_orient)
      }
      
    } else {
      index = index_ovp[[i]][(genome_info$trx_gene_AA_major[index_ovp[[i]],] %>% mutate(ind = 1:length(index_ovp[[i]])) %>% arrange(desc(aa_length), desc(trx_length)) %>% pull(ind))]
      target_annot = brk_annot(ref_trx = genome_info$trx_gene_AA_major, index = index[1], target_region = G_target[[i]], target_orient = target_orient, up_down = up_down,
                               genome_info = genome_info, BSgenome = BSgenome)
    }
    
    target_annot_all = rbind(target_annot_all, target_annot)
  }
  
  # if both up and downstream are out-of-strand --> this is acually in-strand fusion
  if(sum(target_annot_all$RE_category == "out-of-strand")==2){
    target_annot_all = target_annot_all[c(2,1),]
    target_annot_all$RE_category = "in-strand"
  }
  
  return(target_annot_all)
}



# function to define the breakpoint
brk_annot = function(ref_trx, index, target_region, target_orient, up_down, genome_info, BSgenome){
  
  target_gene_trx = ref_trx[index,]
  target_gene_trx = target_gene_trx %>% mutate(Protein_coding = ifelse(!is.na(aa_length), TRUE, FALSE))
  
  if(nrow(target_gene_trx)>1){
    stop("more than one gene in the brkpt")
  }
  
  # target region is upstream of the RE
  if(up_down == "up"){
    
    if((target_gene_trx$STRAND == "+" & target_orient == 1) | (target_gene_trx$STRAND == "-" & target_orient == (-1))){
      re_categ = "in-strand"
    } else {
      re_categ = "out-of-strand"
    }
    # target region is downstream of the RE
  } else if (up_down == "down") {
    if((target_gene_trx$STRAND == "+" & target_orient == (-1)) | (target_gene_trx$STRAND == "-" & target_orient == 1)){
      re_categ = "in-strand"
    } else {
      re_categ = "out-of-strand"
    }
  } else {
    stop("'up_down' parameter should be either 'up' or 'down'")
  }
  
  
  ################################################## 
  # build the GRange object for the target transcript
  ##################################################
  ind = as.numeric(target_gene_trx$TXID)
  t_exon = genome_info$trx_exon[[ind]]
  t_intron = genome_info$trx_intron[[ind]]
  t_str = unique(as.character(strand(t_exon)))
  
  if(t_str == "-"){
    t_intron = rev(t_intron)
  }
  
  df1 = data.frame(data.frame(t_exon)[, c(1,2,3,5)], id = paste("Exon", 1:length(t_exon)))
  df2 = rbind(data.frame(data.frame(t_intron)[, c(1,2,3,5)], id = paste("Intron", 1:length(t_intron))), NA)
  df_comb = (gdata::interleave(df1, df2))
  df_comb = df_comb[-nrow(df_comb),]
  G_ref = makeGRangesFromDataFrame(df_comb)
  G_ref$id = df_comb$id
  
  #######################################################################################
  # Find the location of the breakpoint and downstream exons within the target transcript
  #######################################################################################
  # Compare the Ranges between the target transcript and target region (breakpoint)
  temp_ovp = findOverlaps(G_ref, target_region)
  loc_in_trx = G_ref[queryHits(temp_ovp)]
  
  ind = queryHits(temp_ovp)
  
  # following exon from the breakpoint
  if(ind == length(G_ref)){
    down_exon = G_ref[length(G_ref)]
  } else {
    down_exon = G_ref[(ind+1):length(G_ref)]
  }
  down_exon = down_exon[grepl("Exon", down_exon$id)]
  
  #######################################################################################
  # check whether the corresponding trx is protein-coding or not
  #######################################################################################
  if(target_gene_trx$Protein_coding == T){
    # protein coding region in the transcript
    G_cds_target = genome_info$G_cds[[which(names(genome_info$G_cds) == target_gene_trx$TXID)]]
    # Protein sequence for the transcript
    target_cds_seq_AA = as.character(genome_info$cds_seqs_aa[names(genome_info$cds_seqs_aa) == target_gene_trx$TXID][[1]])
    
    #######################################################################################
    # check the following exon is in-frame or out-of-frame or unknown
    #######################################################################################
    # the nearby exon is completely overlapped with exon of full CDS on the transcript (exon in the middle of transcript)
    if(sum(G_cds_target == down_exon[1]) == 1){
      # check whether the protein sequence of the neigbor exon (previous exon for upstream and next exon for downstream) is in-frame or not
      # this can be identified by the aa sequence of the next exon of the breakpoint
      t_seq = as.character(suppressWarnings(Biostrings::translate(Biostrings::getSeq(BSgenome, G_cds_target[which(G_cds_target == down_exon[1])])))[[1]])
      # target exon is inframe --> fully overlapped with protein sequence of the transcript
      if(grepl(t_seq, target_cds_seq_AA, fixed = T)){
        down_exon_frame = "inframe:CDS_full"
      } else {
        down_exon_frame = "out_of_frame:CDS_full"
      }
      
      # the following exon is partly overlapped with CDS of the transcript (exons with translation start (5') or end (3'))
    } else if(length(queryHits(findOverlaps(down_exon[1], G_cds_target))) == 1){
      
      # if the following exon include translation start site -> the following exon should include the first CDS exon
      if(length(queryHits(findOverlaps(down_exon[1], G_cds_target[1])))==1){
        down_exon_frame = "unknown:CDS_start"
        
        # if the following exon include translation end site -> the following exon should include the last CDS exon
      } else if(length(queryHits(findOverlaps(down_exon[1], G_cds_target[length(G_cds_target)])))==1){
        # protein sequence of the last CDS exon
        t_seq = as.character(Biostrings::translate(Biostrings::getSeq(BSgenome, G_cds_target[length(G_cds_target)]))[[1]])
        # check whether the last CDS exon is in-frame or not
        if(grepl(t_seq, target_cds_seq_AA, fixed = T)){
          down_exon_frame = "inframe:CDS_end"
        } else {
          down_exon_frame = "out_of_frame:CDS_end"
        }
        
      } else {
        stop("unexpected REs with partial CDS")
      }
      
      # the following exon is completely UTR
    } else if(length(queryHits(findOverlaps(down_exon[1], G_cds_target))) == 0){
      
      # if there is remaining CDS exons (5'-UTR)
      if(length(queryHits(findOverlaps(down_exon[1], G_cds_target)))!=0){
        down_exon_frame = "unknown:UTR_5"
        # if there is no remaining CDS exons (3'-UTR)
      } else {
        down_exon_frame = "unknown:UTR_3"
      }
      
    } else {
      stop("unexpected RE scenario")
    }
    
    # if target transcript is non-coding  
  } else {
    down_exon_frame = "unknown:noncoding trx"
  }
  down_exon = as.character(down_exon$id[1])
  
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
                                                            down_exon = down_exon,
                                                            down_exon_frame = down_exon_frame,
                                                            RE_category = re_categ,
                                                            orientation = target_orient)
  return(target_annot)
}


# annotation for the vcf file
get_orient = function(purple_vcf){
  if(endsWith(purple_vcf$V5, "[")){
    orient = c(1, -1)
  } else if (startsWith(purple_vcf$V5, "]")){
    orient = c(-1, 1)
  } else if (startsWith(purple_vcf$V5, "[")){
    orient = c(-1, -1)
  } else if (endsWith(purple_vcf$V5, "]")){
    orient = c(1, 1)
    # SGL
  } else {
    if(endsWith(purple_vcf$V5, ".")){
      orient = c(1, 0)
    } else {
      orient = c(-1, 0)
    }
  }
  return(orient)
}

trx_annot = function(ref_trx_id, target_region){
  t_enst_id = hg_id_19_ENST %>% filter(X.hg19.knownGene.name %in% ref_trx_id) %>% pull(ENST_ID)
  index = grep(t_enst_id, hg38_info$trx_gene_AA$TXNAME)
  
  target_gene_trx = hg38_info$trx_gene_AA[index,]
  target_gene_trx = target_gene_trx %>% mutate(Protein_coding = ifelse(!is.na(aa_length), TRUE, FALSE))
  
  ################################################## 
  # build the GRange object for the target transcript
  ##################################################
  ind = as.numeric(target_gene_trx$TXID)
  t_exon = hg38_info$trx_exon[[ind]]
  t_intron = hg38_info$trx_intron[[ind]]
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
  temp_ovp = findOverlaps(target_region, G_ref, ignore.strand = T)
  ind_map = queryHits(temp_ovp)
  ind_ref = subjectHits(temp_ovp)
  
  target_region_mapped = target_region[ind_map]
  
  loc_in_trx_map = G_ref[ind_ref]
  loc_in_trx_map$brkpt = paste(seqnames(target_region_mapped), ranges(target_region_mapped), sep = ":")
  loc_in_trx_map$junction_type = target_region_mapped$junction_type
  
  # check whether the breakpoint is mapped to end of exon (so, actual breakpoint is in somewhere of intronic region)
  junc_site = cbind(paste(seqnames(loc_in_trx_map), start(loc_in_trx_map), sep = ":"),
                    paste(seqnames(loc_in_trx_map), end(loc_in_trx_map), sep = ":"))
  a = as.numeric(junc_site[,1] == loc_in_trx_map$brkpt & grepl("Intron", loc_in_trx_map$id)) 
  b = as.numeric(junc_site[,2] == loc_in_trx_map$brkpt & grepl("Intron", loc_in_trx_map$id))
  
  loc_in_trx_map$canonical_junction = as.logical(a+b)
  
  # unmapped
  ind_unmap = setdiff(1:length(target_region), ind_map)
  if(length(ind_unmap)!=0){
    loc_in_trx_unmap = makeGRangesFromDataFrame(data.frame(seqnames = rep(unique(as.character(seqnames(G_ref))), length(ind_unmap)), 
                                                           tx_start = rep(0, length(ind_unmap)), 
                                                           tx_end = rep(0, length(ind_unmap)), 
                                                           strand = rep("*", length(ind_unmap))))
    loc_in_trx_unmap$id = NA
    loc_in_trx_unmap$brkpt = do.call("paste", c(data.frame(target_region)[ind_unmap,1:2], sep = ":"))
    loc_in_trx_unmap$junction_type = target_region$junction_type[ind_unmap]
    loc_in_trx_unmap$canonical_junction = FALSE
    
    # merging
    loc_in_trx = c(loc_in_trx_map, loc_in_trx_unmap)
    loc_in_trx = loc_in_trx[order(c(ind_map, ind_unmap))]
  } else {
    loc_in_trx = loc_in_trx_map
  }
  return(loc_in_trx)
}


neighbor_genes = function(target_ranges, genome_info){
  t_prot = which(!is.na(genome_info$trx_gene_AA_major$aa_length))
  t_genome = data.frame(genome_info$trx_major[t_prot], stringsAsFactors = F) %>% 
    mutate(gene = genome_info$trx_gene_AA_major$SYMBOL[t_prot]) %>% arrange(seqnames, start)
  
  t_tar_left = t_genome %>%
    filter(seqnames == as.character(seqnames(target_ranges)), end < start(target_ranges))
  t_tar_right = t_genome %>%
    filter(seqnames == as.character(seqnames(target_ranges)), start > start(target_ranges))
  
  output = paste(tail(t_tar_left$gene, 1), head(t_tar_right$gene, 1), sep = "--")
  return(output)
}


intergenic_annot = function(ref_pos, target_pos){
  
  ##################################
  # neighbor genes for the reference (hg19)
  ##################################
  
  ref_pos = data.frame(ref_pos$seqnames, ref_pos$brkpos, ref_pos$brkpos, ".")
  colnames(ref_pos) = c("seqnames", "tx_start", "tx_end", "strand")
  ref_pos = makeGRangesFromDataFrame(ref_pos)
  
  neighbor_genes_ref = neighbor_genes(target_ranges = ref_pos, genome_info = hg19_info)
  
  ##################################
  # neighbor genes for the reference (hg38)
  ##################################
  neighbor_genes_tar = rep(NA, length(target_pos))
  for(i in 1:length(target_pos)){
    ind = queryHits(findOverlaps(hg38_info$trx, target_pos[i], ignore.strand = T))
    # target breakpoint not overlapped with the known transcripts
    if(length(ind)==0){
      
      neighbor_genes_tar[i] = neighbor_genes(target_ranges = target_pos[i], genome_info = hg38_info)
      
      # target breakpoint overlapped with the known transcripts
    } else {
      neighbor_genes_tar[i] = (hg38_info$trx_gene_AA[ind,] %>% 
                                 mutate(canonical = !grepl("[.]", SYMBOL)) %>% 
                                 arrange(desc(canonical)) %>% pull(SYMBOL))[1]
    }
    print(i)
  }
  
  ##################################
  # comparison
  ##################################
  target_pos$neighbors = neighbor_genes_tar
  target_pos$ref = neighbor_genes_ref
  
  t_ovp_1 = as.numeric(stringr::str_split_fixed(target_pos$neighbors, "--", 2)[,1] %in% unique(stringr::str_split_fixed(target_pos$ref, "--", 2)))
  t_ovp_2 = as.numeric(stringr::str_split_fixed(target_pos$neighbors, "--", 2)[,2] %in% unique(stringr::str_split_fixed(target_pos$ref, "--", 2)))
  t_ovp = t_ovp_1 + t_ovp_2
  target_pos$overlap = t_ovp
  
  return(target_pos)
}


############################################# 
# target breakpoints from purple (WGS) 
#############################################
# we don't know whether reads from the STAR have the same orientation
# thus, both orientation should be considered
# e.g) upstream-downstream, downstream-upstream from STAR

comp_brk = function(t_annot, dir, target_region, up_down){
  t_gene = c(t_annot$up.trx_name, t_annot$down.trx_name)
  t_brk = list(t_annot[, c("up.seqnames", "up.start", "up.loc_in_trx", "up.down_exon_frame", "up.RE_category", "up.gene_symbol")],
               t_annot[, c("down.seqnames", "down.start", "down.loc_in_trx", "down.down_exon_frame", "down.RE_category", "down.gene_symbol")])
  colnames(t_brk[[1]]) = c("seqnames", "brkpos", "loc_in_trx", "down_exon_frame", "RE_category", "symbol")
  colnames(t_brk[[2]]) = c("seqnames", "brkpos", "loc_in_trx", "down_exon_frame", "RE_category", "symbol")
  
  # Forward direction
  if(dir == "F" & up_down == "up"){
    t_target = target_region[[1]]
  } else if (dir == "F" & up_down == "down"){
    t_target = target_region[[2]]
  } else if(dir == "R" & up_down == "up"){
    t_target = target_region[[2]]
  } else if(dir == "R" & up_down == "down"){
    t_target = target_region[[1]]
  } else {
    stop("dir should be either F/R and up_down should be either up/down")
  }
  
  # up or down
  ind = ifelse(up_down == "up", 1, 2)
  
  # target region is located on the gene
  if(!is.na(t_gene[ind]) & t_brk[[ind]]$RE_category == "in-strand" & !grepl("UTR_3", t_brk[[ind]]$down_exon_frame)){
    t_brk_pos = unlist(strsplit(as.character(t_brk[[ind]]$loc_in_trx), ":"))[6]
    t_loc = trx_annot(ref_trx_id = t_gene[ind], target_region = t_target)
    # breakpoint is in the intronic region --> reads to canonical junction in the intron
    if(grepl("Intron", as.character(t_brk[[ind]]$loc_in_trx))){
      target_ind = which(t_loc$id == t_brk_pos & t_loc$canonical_junction == TRUE)
    # breakpoint is in the exonic region --> reads to mapped in any of the region in the gene
    } else {
      target_ind = which(!is.na(t_loc$id))
    }
  } else if(!is.na(t_gene[ind]) & t_brk[[ind]]$RE_category == "out-of-strand"){
    t_brk_pos = unlist(strsplit(as.character(t_brk[[ind]]$loc_in_trx), ":"))[6]
    t_loc = trx_annot(ref_trx_id = t_gene[ind], target_region = t_target)
    target_ind = which(!is.na(t_loc$id))
  } else if(up_down == "down" & grepl("UTR_3", t_brk[[ind]]$down_exon_frame)){
    ref_pos = t_brk[[ind]][1:2]
    colnames(ref_pos) = c("seqnames", "brkpos")
    t_loc = intergenic_annot(ref_pos = ref_pos, target_pos = t_target)
    target_ind = grep(t_brk[[ind]]$symbol, t_loc$neighbors)
  # target region is located on the intergenic region
  } else if(as.character(t_brk[[ind]]$seqnames) != "chr-1"){
    ref_pos = t_brk[[ind]][1:2]
    colnames(ref_pos) = c("seqnames", "brkpos")
    t_loc = intergenic_annot(ref_pos = ref_pos, target_pos = t_target)
    if(grepl("--FGFR2", t_loc$ref[1])){
      target_ind = grep("--FGFR2", t_loc$neighbors)
    } else {
      target_ind = which(t_loc$overlap >0 & t_loc$neighbors!="FGFR2")
    }
    
    # SGL type
  } else {
    target_ind = NA
  }
  
  return(target_ind)
}

create_Gref = function(genome_info, TXID){
  ind = as.numeric(TXID)
  t_exon = genome_info$trx_exon[[ind]]
  t_intron = genome_info$trx_intron[[ind]]
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
  return(G_ref)
}

########################################################################
# RNA-seq and WGS data loading
########################################################################
load("~/FGFR/Daniel/R/R_files/FGFR2_RE_HMF_v3.RData")

list_linx = list.files("/DATA/projects/j.bhin/Daniel_FGFR2/HMF/linx/")
temp = c()
for(i in 1:length(list_linx)){
  t = list.files(paste("/DATA/projects/j.bhin/Daniel_FGFR2/HMF/linx/", list_linx[i], sep = ""))
  t = t[grep("T", t)][1]
  temp = c(temp, stringr::str_split_fixed(t, "[.]", 2)[1])
}
list_linx = data.frame(file = list_linx, sample = temp, stringsAsFactors = F)

list_RNA = list.files("/DATA/projects/j.bhin/Daniel_FGFR2/HMF/RNAseq/bam/")
list_RNA = data.frame(file = list_RNA, sample = stringr::str_split_fixed(list_RNA, "_", 3)[,1], stringsAsFactors = F)


########################################################################
# FGFR2 SVs from Linx/purple/STAR
########################################################################

linx_path = "/DATA/projects/j.bhin/Daniel_FGFR2/HMF/linx"
purple_path = "/DATA/projects/j.bhin/Daniel_FGFR2/HMF/somatics/"
rna_path = "/DATA/projects/j.bhin/Daniel_FGFR2/HMF/RNAseq/bam/"

target_samp = intersect(list_linx$sample, list_RNA$sample)
target_in_linx = list_linx[match(target_samp, list_linx$sample),]
target_in_RNA = list_RNA[match(target_samp, list_RNA$sample),]

FGFR2_SV = vector(mode = "list", length(target_samp))
for(i in 1:length(target_samp)){
  ######################
  # linx
  ######################
  # fusions defined by linx
  t = list.files(paste(linx_path, target_in_linx$file[i], sep = "/"))
  linx_fus = read.table(paste(linx_path, target_in_linx$file[i], t[grep("linx.fusion.tsv", t)], sep = "/"), sep = "\t", stringsAsFactors = F, header = T)
  linx_fus_FGFR2 = linx_fus %>% filter(grepl("FGFR2", Name))
  
  # SVs defined by linx
  linx_SV = read.table(paste(linx_path, target_in_linx$file[i], t[grep("vis_sv_data", t)], sep = "/"), sep = "\t", stringsAsFactors = F, header = T)
  linx_SV_FGFR2 = unique(rbind(linx_SV %>% filter(ChrStart == 10 , PosStart >= 123237844, PosStart <= 123357972),
                               linx_SV %>% filter(ChrEnd == 10 , PosEnd >= 123237844, PosEnd <= 123357972)))
  if(nrow(linx_SV_FGFR2)!=0){
    # annotation
    linx_SV_annot = c()
    for(j in 1:nrow(linx_SV_FGFR2)){
      G_1 = data.frame(linx_SV_FGFR2[j, c("ChrStart", "PosStart", "PosStart")], ".") %>%
        mutate(Orient = linx_SV_FGFR2$OrientStart[j], ChrStart = paste("chr", ChrStart, sep = ""))
      colnames(G_1) = c("seqnames", "tx_start", "tx_end", "strand", "Orient")
      G_1 = makeGRangesFromDataFrame(G_1)
      
      G_2 = data.frame(linx_SV_FGFR2[j, c("ChrEnd", "PosEnd", "PosEnd")], ".") %>%
        mutate(Orient = linx_SV_FGFR2$OrientEnd[j], ChrEnd = paste("chr", ChrEnd, sep = ""))
      colnames(G_2) = c("seqnames", "tx_start", "tx_end", "strand", "Orient")
      G_2 = makeGRangesFromDataFrame(G_2)
      
      t_out = brk_annot_framework(G_up = G_1, up_orient = linx_SV_FGFR2$OrientStart[j],
                                  G_down = G_2, down_orient = linx_SV_FGFR2$OrientEnd[j],
                                  BSgenome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                                  genome_info = hg19_info)
      t_out = data.frame(up = t_out[1,], down = t_out[2,])
      linx_SV_annot = rbind(linx_SV_annot, t_out)
    }
  } else {
    linx_SV_annot = data.frame(NULL)
  }
  
  
  ########################################
  # RE evidence from RNA-seq (STAR and STAR-Fusion)
  ########################################
  star_fus = FGFR2_star_fusion %>% filter(sample_id == target_in_RNA$file[i])
  star_re = FGFR2_RE %>% filter(sample_id == target_in_RNA$file[i])
  
  
  ########################################
  # purple
  ########################################
  t = list.files(paste(purple_path, target_in_linx$file[i], sep = "/"))
  purple_SV = read.table(gzfile(paste(purple_path, target_in_linx$file[i], t[grepl("purple.sv.ann", t) & !grepl("tbi", t)], sep = "/")), sep = "\t", stringsAsFactors = F, header = F)
  purple_SV_FGFR2 = purple_SV %>% filter(V1 == 10, V2 >=123237844, V2 <=123357972)
  
  if(nrow(purple_SV_FGFR2)!=0){
    t_partner = stringr::str_split_fixed(purple_SV_FGFR2$V5, "\\]|\\[|\\:",4)[,2:3]
    if(nrow(purple_SV_FGFR2)==1){
      t_partner = data.frame(t(t_partner), stringsAsFactors = F)
    } else {
      t_partner = data.frame(t_partner, stringsAsFactors = F)
    }
    t_partner[t_partner[,1] == "", 1] = "-1"
    t_partner[t_partner[,2] == "", 2] = 0
    colnames(t_partner) = c("V1", "V2")
    
    # annotation
    purple_SV_annot = c()
    for(j in 1:nrow(purple_SV_FGFR2)){
      orient = get_orient(purple_SV_FGFR2[j,])
      G_1 = data.frame(purple_SV_FGFR2[j, c("V1", "V2", "V2")], ".") %>%
        mutate(Orient = orient[1], V1 = paste("chr", V1, sep = ""))
      colnames(G_1) = c("seqnames", "tx_start", "tx_end", "strand", "Orient")
      G_1 = makeGRangesFromDataFrame(G_1)
      
      G_2 = data.frame(t_partner[j,1:2], t_partner[j,2],  ".") %>%
        mutate(Orient = orient[2], V1 = paste("chr", V1, sep = ""))
      colnames(G_2) = c("seqnames", "tx_start", "tx_end", "strand", "Orient")
      G_2 = makeGRangesFromDataFrame(G_2)
      
      t_out = suppressWarnings(brk_annot_framework(G_up = G_1, up_orient =  orient[1],
                                                   G_down = G_2, down_orient =  orient[2],
                                                   BSgenome = BSgenome.Hsapiens.UCSC.hg19::Hsapiens,
                                                   genome_info = hg19_info))
      t_out = data.frame(up = t_out[1,], down = t_out[2,])
      purple_SV_annot = rbind(purple_SV_annot, t_out)
    }
    # remove redundant SVs
    purple_SV_annot = purple_SV_annot %>% arrange(up.RE_category)
    t_sv = t(apply(purple_SV_annot %>% select(up.start, down.start), 1, sort))
    t_id = apply(t_sv, 1, paste, collapse = "-")
    purple_SV_annot = purple_SV_annot[match(unique(t_id), t_id),]
    
    # check whether the RE was detected by LINX and STAR-Fusion
    t_la = do.call("paste", c(purple_SV_annot[, c("up.gene_symbol", "down.gene_symbol")], sep = "_"))
    t_linx_lab = c(); t_sf_lab = c()
    for(j in 1:length(t_la)){
      t_la_2 = t_la[j]
      
      if(!("NA" %in% (unlist(strsplit(t_la_2, "_"))))){
        # LINX-Fusion
        if(!is.na(match(t_la_2, linx_fus_FGFR2$Name))){
          t_linx_lab_2 = do.call("paste", c(linx_fus[match(t_la_2, linx_fus$Name), c("Name", "FusedExonUp", "FusedExonDown")], sep = "_"))
        } else {
          t_linx_lab_2 = NA
        }
        
        # STAR-Fusion
        la_l = unlist(strsplit(t_la_2, "_"))[1]
        la_r = unlist(strsplit(t_la_2, "_"))[2]
        la_l = (hg_id_38_ENST %>% filter(GENEID == (hg19_info$trx_gene_AA_major %>% filter(SYMBOL == la_l) %>% pull(GENEID))) %>% pull(GENE_SYMBOL))[1]
        la_r = (hg_id_38_ENST %>% filter(GENEID == (hg19_info$trx_gene_AA_major %>% filter(SYMBOL == la_r) %>% pull(GENEID))) %>% pull(GENE_SYMBOL))[1]
        t_la_2 = paste(la_l, la_r, sep = "--")
        if(!is.na(match(t_la_2, star_fus$X.FusionName))){
          t_sf_lab_2 = paste((star_fus[match(t_la_2, star_fus$X.FusionName),] %>% 
                              mutate(lab = paste(left_annot, right_annot, sep = "_")) %>% pull(lab)), collapse = "//")
        } else {
          t_sf_lab_2 = NA
        }
      } else {
        t_linx_lab_2 = NA
        t_sf_lab_2 = NA
      }
      
      t_linx_lab = c(t_linx_lab, t_linx_lab_2)
      t_sf_lab = c(t_sf_lab, t_sf_lab_2)
    }
    purple_SV_annot$LINX_Fusion = t_linx_lab
    purple_SV_annot$STAR_Fusion = t_sf_lab
    
    
    # RNA level evidence from STAR
    t_star = read.table(paste(rna_path, target_in_RNA$file[i], "/Chimeric.out.junction", sep = ""), sep = "\t", stringsAsFactors = F, header = T)
    t_star_fgfr2 = unique(rbind(t_star %>% filter(chr_donorA == "chr10", brkpt_donorA >= 121478332, brkpt_donorA <= 121598458),
                                t_star %>% filter(chr_acceptorB == "chr10", brkpt_acceptorB >= 121478332,brkpt_acceptorB <= 121598458)))
    
    # too many reads in 33 sample
    if(i == 33){
      t_star_fgfr2 = t_star_fgfr2 %>% filter(brkpt_donorA == 121483697 | brkpt_acceptorB == 121483697)
    }
    
    t_star_fgfr2 = t_star_fgfr2 %>% filter(chr_donorA %in% paste("chr", c(1:22, "X"), sep = ""), chr_acceptorB %in% paste("chr", c(1:22, "X"), sep = ""))
    target_region_up = t_star_fgfr2[, c("chr_donorA", "brkpt_donorA", "brkpt_donorA", "strand_donorA")]
    target_region_down = t_star_fgfr2[, c("chr_acceptorB", "brkpt_acceptorB", "brkpt_acceptorB", "strand_acceptorB")]
    colnames(target_region_up) = c("seqnames", "tx_start", "tx_end", "strand")
    colnames(target_region_down) = c("seqnames", "tx_start", "tx_end", "strand")
    
    target_region_up = makeGRangesFromDataFrame(target_region_up)
    target_region_down = makeGRangesFromDataFrame(target_region_down)
    target_region_up$junction_type = t_star_fgfr2$junction_type
    target_region_down$junction_type = t_star_fgfr2$junction_type
    target_region = list(target_region_up, target_region_down)
    
    # too many SVs in 33 sample
    if(i == 33){
      purple_SV_annot = purple_SV_annot %>% filter(grepl("Intron 17", up.loc_in_trx))
      ind = t_star_fgfr2 %>% mutate(index= 1:nrow(t_star_fgfr2)) %>% 
        filter((chr_donorA == "chr10" & brkpt_donorA == 121483697) | (chr_acceptorB == "chr10" & brkpt_acceptorB == 121483697)) %>% pull(index)
      target_region = list(target_region_up[ind], target_region_down[ind])
    }
    
    n_RNA = c(); annot = c()
    for(k in 1:nrow(purple_SV_annot)){
      t_annot = purple_SV_annot[k,]
      
      # forward
      target_ind_u = comp_brk(t_annot = t_annot, dir = "F", target_region = target_region, up_down = "up")
      target_ind_d = comp_brk(t_annot = t_annot, dir = "F", target_region = target_region, up_down = "down")
      if(NA %in% c(target_ind_u, target_ind_d)){
        target_ind_F = setdiff(c(target_ind_u, target_ind_d), NA)
      } else {
        target_ind_F = intersect(target_ind_u, target_ind_d)
      }
      
      # reverse
      target_ind_u = comp_brk(t_annot = t_annot, dir = "R", target_region = target_region, up_down = "up")
      target_ind_d = comp_brk(t_annot = t_annot, dir = "R", target_region = target_region, up_down = "down")
      if(NA %in% c(target_ind_u, target_ind_d)){
        target_ind_R = setdiff(c(target_ind_u, target_ind_d), NA)
      } else {
        target_ind_R = intersect(target_ind_u, target_ind_d)
      }
      
      target_ind = c(target_ind_F, target_ind_R)
      
      # if multiple coordinates of the RNA-seq reads are potentially matched to WGS-SVs, selecting the ones with canonical splicing acceptor & donors
      target_ind = t_star_fgfr2[target_ind,] %>% mutate(target_ind = target_ind) %>% filter(junction_type>0) %>% pull(target_ind)
      
      t_tab = do.call("paste", c(t_star_fgfr2[target_ind, 1:7], sep = ":"))
      t_tab_sum = sort(table(t_tab), decreasing = T)
      t_target_ind = target_ind[match(names(t_tab_sum)[1], t_tab)]
      
      if(length(target_ind)!=0){
        
        t_tab = do.call("paste", c(t_star_fgfr2[target_ind, 1:7], sep = ":"))
        t_tab_sum = sort(table(t_tab), decreasing = T)
        t_target_ind = target_ind[match(names(t_tab_sum)[1], t_tab)]
        
        # annotating WGS-driven intergenic RE
        t_purple_gene = as.character(t_annot[, c("up.gene_symbol", "down.gene_symbol")])
        t_star_up = hg38_info$trx_gene_AA$SYMBOL[queryHits(findOverlaps(hg38_info$trx, target_region_up[t_target_ind], ignore.strand = T))]
        t_star_up = names(sort(table(t_star_up), decreasing = T))
        t_star_up = c(t_star_up[!(grepl("[.]", t_star_up))], t_star_up[(grepl("[.]", t_star_up))])[1]
        t_star_up = ifelse(is.null(t_star_up), NA, t_star_up)
        t_star_down = hg38_info$trx_gene_AA$SYMBOL[queryHits(findOverlaps(hg38_info$trx, target_region_down[t_target_ind], ignore.strand = T))]
        t_star_down = names(sort(table(t_star_down), decreasing = T))
        t_star_down = c(t_star_down[!(grepl("[.]", t_star_down))], t_star_down[(grepl("[.]", t_star_down))])[1]
        t_star_down = ifelse(is.null(t_star_down), NA, t_star_down)
        
        t_star_gene = c(t_star_up, t_star_down)
        
        # if the supporting reads from STAR has the same orientation with purple
        if(which(t_purple_gene == "FGFR2") == which(t_star_gene == "FGFR2")){
          t_ind = c(1,2)
        } else {
          t_ind = c(2,1)
        }
        t_trx_id_u = hg19_info$trx_gene_AA_major %>% filter(SYMBOL %in% t_star_gene[t_ind[1]]) %>% pull(TXNAME)
        if(length(t_trx_id_u)==0){
          t_trx_id_u = hg19_info$trx_gene_AA_major %>% filter(GENEID == (hg_id_38_ENST %>% filter(GENE_SYMBOL == t_star_gene[t_ind[1]]) %>% pull(GENEID))[1]) %>% pull(TXNAME)
        }
        t_trx_id_d = hg19_info$trx_gene_AA_major %>% filter(SYMBOL %in% t_star_gene[t_ind[2]]) %>% pull(TXNAME)
        if(length(t_trx_id_d)==0){
          t_trx_id_d = hg19_info$trx_gene_AA_major %>% filter(GENEID == (hg_id_38_ENST %>% filter(GENE_SYMBOL == t_star_gene[t_ind[2]]) %>% pull(GENEID))[1]) %>% pull(TXNAME)
        }
        
        # check whether the breakpoint is in the gene or intergenic region
        if(length(t_trx_id_u)>0){
          t_loc_u = trx_annot(ref_trx_id = t_trx_id_u, target_region = unique(target_region[[t_ind[1]]][t_target_ind]))
        } else {
          t_brkpt = target_region[[t_ind[1]]][t_target_ind]
          t_loc_u = data.frame(id = paste("intergenic", neighbor_genes(target_ranges = unique(t_brkpt), genome_info = hg38_info), sep = "::"),
                               brkpt = paste(as.character(seqnames(t_brkpt)), start(t_brkpt), sep = ":"))
        }
        
        if(length(t_trx_id_d)>0){
          t_loc_d = trx_annot(ref_trx_id = t_trx_id_d, target_region = unique(target_region[[t_ind[2]]][t_target_ind]))
        } else {
          t_brkpt = target_region[[t_ind[2]]][t_target_ind]
          t_loc_d = data.frame(id = paste("intergenic", neighbor_genes(target_ranges = unique(target_region[[t_ind[2]]][t_target_ind]), genome_info = hg38_info), sep = "::"),
                               brkpt = paste(as.character(seqnames(t_brkpt)), start(t_brkpt), sep = ":"))
        }
        
        t_sv_annot = paste(t_star_gene[t_ind[1]], t_loc_u$id, t_loc_u$brkpt, t_star_gene[t_ind[2]], t_loc_d$id, t_loc_d$brkpt, sep = "_")
        
      } else if(nrow(linx_fus_FGFR2)!=0){
        # make artificial t_annot file for linx fusion
        t_gene = stringr::str_split_fixed(linx_fus_FGFR2$Name, "_", 2)
        t_cnt = c(); t_fus= c()
        for(z in 1:nrow(linx_fus_FGFR2)){
          
          up_down = ifelse(startsWith(linx_fus_FGFR2$Name, "FGFR2"), "up", "down")
          trx_fgfr2 = hg38_info$trx_gene_AA_major %>% filter(SYMBOL == "FGFR2") %>% pull(TXID)
          trx_fgfr2_ref = create_Gref(genome_info = hg38_info, trx_fgfr2)
          
          # FGFR2 is upstream
          if(up_down == "up"){
            trx_fgfr2_ref_pos = trx_fgfr2_ref[which(trx_fgfr2_ref$id == paste("Intron", linx_fus_FGFR2$FusedExonUp[z]))]
            trx_fgfr2_ref_pos_map = end(trx_fgfr2_ref_pos)
            trx_partner = data.frame(hg38_info$trx_intron[hg38_info$trx_gene_AA_major %>% filter(SYMBOL == t_gene[z,2]) %>% pull(TXID)])
            if(nrow(trx_partner)==0){
              t_gene_2 = (hg_id_38_ENST %>% filter(GENEID == hg19_info$trx_gene_AA_major %>% filter(SYMBOL == t_gene[z,2]) %>% pull(GENEID)) %>% pull(GENE_SYMBOL))[1]
              trx_partner = data.frame(hg38_info$trx_intron[hg38_info$trx_gene_AA_major %>% filter(SYMBOL == t_gene_2) %>% pull(TXID)])
            }
          } else {
            trx_fgfr2_ref_pos = trx_fgfr2_ref[which(trx_fgfr2_ref$id == paste("Intron", linx_fus_FGFR2$FusedExonDown[z]))]
            trx_fgfr2_ref_pos_map = start(trx_fgfr2_ref_pos)
            trx_partner = data.frame(hg38_info$trx_intron[hg38_info$trx_gene_AA_major %>% filter(SYMBOL == t_gene[z,1]) %>% pull(TXID)])
            if(nrow(trx_partner)==0){
              t_gene_2 = (hg_id_38_ENST %>% filter(GENEID == hg19_info$trx_gene_AA_major %>% filter(SYMBOL == t_gene[z,1]) %>% pull(GENEID)) %>% pull(GENE_SYMBOL))[1]
              trx_partner = data.frame(hg38_info$trx_intron[hg38_info$trx_gene_AA_major %>% filter(SYMBOL == t_gene_2) %>% pull(TXID)])
            }
          }
          # forward
          t_1 = t_star_fgfr2 %>% filter(chr_donorA == as.character(seqnames(trx_fgfr2_ref_pos)), brkpt_donorA == trx_fgfr2_ref_pos_map)
          t_1_id = do.call("paste", c(t_1[,4:5], sep = "_"))
          # reverse
          t_2 = t_star_fgfr2 %>% filter(chr_acceptorB == as.character(seqnames(trx_fgfr2_ref_pos)), brkpt_acceptorB == trx_fgfr2_ref_pos_map)
          t_2_id = do.call("paste", c(t_2[,1:2], sep = "_"))
          
          t_id_sum = table(rbind(t_1_id, t_2_id))
          t_id_sum = c(t_id_sum[names(t_id_sum) %in% do.call("paste", c(trx_partner[, c(3,4)], sep = "_"))],  
                       t_id_sum[names(t_id_sum) %in% do.call("paste", c(trx_partner[, c(3,5)], sep = "_"))])
          t_cnt = c(t_cnt, t_id_sum)
          t_fus = c(t_fus, paste(linx_fus_FGFR2[k, c("Name", "FusedExonUp", "FusedExonDown")], collapse = "_"))
        }
        target_ind = c(which(do.call("paste", c(t_star_fgfr2[, c(1,2)], sep = "_")) %in% names(t_cnt)), 
                       which(do.call("paste", c(t_star_fgfr2[, c(4,5)], sep = "_")) %in% names(t_cnt)))
        t_sv_annot = paste(t_fus, collapse = "//", paste(names(t_cnt), collapse = "//"), sep = " & ") 
        
      } else {
        t_sv_annot = NA 
      }
      
      # remove the optical replicates
      t_RNA = unique(t_star_fgfr2[target_ind, -c(10)])
      n_RNA = c(n_RNA, nrow(t_RNA))
      
      annot = c(annot, t_sv_annot)
      
    }
    
    purple_SV_annot$n_RNA = n_RNA
    purple_SV_annot$STAR_chim = annot
    
  } else {
    purple_SV_annot = data.frame(NULL)
  }
  
  t_out = list(linx_fus = linx_fus_FGFR2, linx_SV = linx_SV_FGFR2, linx_SV_annot = linx_SV_annot, 
               purple_SV = purple_SV_FGFR2, purple_SV_annot = purple_SV_annot,
               STAR_Fusion = star_fus, STAR_RE = star_re)
  
  FGFR2_SV[[i]] = t_out
  
  print(i)
  
}
names(FGFR2_SV) = target_samp

# 
FGFR2_SV_n = c()
for(i in 1:length(FGFR2_SV)){
  t = unlist(lapply(FGFR2_SV[[i]], nrow))[c(1,3,5,6,7)]
  FGFR2_SV_n = rbind(FGFR2_SV_n, t)
}
rownames(FGFR2_SV_n) = target_samp


RE_type = function(purple_annot){
  # up and down are FGFR2
  if((purple_annot$up.gene_symbol %in% "FGFR2") & (purple_annot$down.gene_symbol %in% "FGFR2")){
    RE = "internal"
    # one of the partner is intergenic
  } else if(is.na(purple_annot$up.gene_symbol) | is.na(purple_annot$down.gene_symbol)){
    if(purple_annot$up.RE_category == "SGL" | purple_annot$down.RE_category == "SGL"){
      RE = "SGL"
    } else {
      RE = "intergenic"
    }
    # in-strand
  } else if((purple_annot$up.RE_category %in% "in-strand") & (purple_annot$down.RE_category %in% "in-strand")){
    # in-frame: both up and downstream is inframe
    if(grepl("inframe", purple_annot$up.down_exon_frame) & grepl("inframe", purple_annot$down.down_exon_frame)){
      if(grepl("Exon", purple_annot$up.loc_in_trx) | grepl("Exon", purple_annot$down.loc_in_trx)){
        RE = "out-of-frame (brkpt is in exon)"
      } else {
        RE = "in-frame"
      }
    } else {
      RE = "out-of-frame"
    }
  } else if((purple_annot$up.RE_category %in% "out-of-strand") | (purple_annot$down.RE_category %in% "out-of-strand")){
    RE = "out-of-strand"
  } else {
    RE = "unknown"
  }
  return(RE)
}

################################################################################
# combined results from PURPLE
################################################################################
ind_tar = which(rowSums(FGFR2_SV_n[,1:3]) >0)
mat_linx_fusion = rep(NA, length(ind_tar))
mat_purple_SV = matrix(NA, length(ind_tar), 4)
mat_star_fusion = matrix(NA, length(ind_tar), 4)
mat_star_RE = matrix(NA, length(ind_tar), 2)
FGFR2_PURPLE_ALL = c()

for(i in 1:length(ind_tar)){
  
  # Linx fusion
  t_linx = FGFR2_SV[[ind_tar[i]]]$linx_fus %>% arrange(desc(FusedExonUp)) %>% select(Name, FusedExonUp, FusedExonDown)
  t_linx$Name = gsub("TENC1", "TNS2", t_linx$Name)
  
  mat_linx_fusion[i] = paste(do.call("paste", c(t_linx, sep = "_")), collapse = " / ")
  
  t_RE = c()
  for(j in 1:nrow(FGFR2_SV[[ind_tar[i]]]$purple_SV_annot)){
    t_RE = c(t_RE, RE_type(FGFR2_SV[[ind_tar[i]]]$purple_SV_annot[j,]))
  }
  
  # Purple SV
  t_purple = FGFR2_SV[[ind_tar[i]]]$purple_SV_annot %>% 
    mutate(C_trunc = (up.down_exon == "Exon 18" & up.RE_category == "in-strand" & up.gene_symbol == "FGFR2")) %>%
    mutate(Exon_up = ifelse(up.RE_category == "in-strand", 
                            as.numeric(stringr::str_split_fixed(up.down_exon, " ", 2)[,2])-1,
                            as.numeric(stringr::str_split_fixed(up.down_exon, " ", 2)[,2])),
           Exon_down = ifelse(down.RE_category == "in-strand", 
                              as.numeric(stringr::str_split_fixed(down.down_exon, " ", 2)[,2]),
                              as.numeric(stringr::str_split_fixed(down.down_exon, " ", 2)[,2])-1)) %>%
    mutate(RE_type = factor(t_RE, levels = c("in-frame", "out-of-frame (brkpt is in exon)", "out-of-frame", "intergenic", "out-of-strand", "internal", "SGL")), 
           FGFR2_is_upstream = (up.gene_symbol == "FGFR2")) %>%
    mutate(sample_name = target_in_linx$file[ind_tar[i]], sample_id = ind_tar[i]) %>%
    arrange(desc(C_trunc), desc(FGFR2_is_upstream), RE_type)
  
  FGFR2_PURPLE_ALL = rbind(FGFR2_PURPLE_ALL, t_purple)
  
  mat_purple_SV[i,1] = do.call("paste", c(t_purple[, c("up.gene_symbol", "down.gene_symbol", "Exon_up", "Exon_down")], sep = "_"))[1]
  mat_purple_SV[i,2] = as.character(t_purple$RE_type[1])
  mat_purple_SV[i,3] = paste(do.call("paste", c(t_purple[, c("up.gene_symbol", "down.gene_symbol", "Exon_up", "Exon_down")], sep = "_")),
                             collapse = " / ")
  mat_purple_SV[i,4] = paste(t_purple$RE_type, collapse = " / ")
  
  # STAR-Fusion
  t_star_fusion = FGFR2_SV[[ind_tar[i]]]$STAR_Fusion %>% 
    mutate(C_trunc = (left_annot == "FGFR2:ENST00000457416.6:Exon 17:ref:in-strand"),
           FGFR2_is_upstream = startsWith(X.FusionName, "FGFR2")) %>%
    arrange(desc(C_trunc), desc(FGFR2_is_upstream))
  
  mat_star_fusion[i,1] = do.call("paste", c(t_star_fusion[, c("left_annot", "right_annot")], sep = "_"))[1]
  mat_star_fusion[i,2] = t_star_fusion$PROT_FUSION_TYPE[1]
  mat_star_fusion[i,3] = paste(do.call("paste", c(t_star_fusion[, c("left_annot", "right_annot")], sep = "_")),
                               collapse = " / ")
  mat_star_fusion[i,4] = paste(t_star_fusion$PROT_FUSION_TYPE, collapse = " / ")
  
  # STAR-RE
  t_star_re = FGFR2_SV[[ind_tar[i]]]$STAR_RE %>%
    mutate(C_trunc = (up.symbol == "FGFR2" & up.id == "Exon 18" & up.splice_type == "ref:in-strand"),
           FGFR2_is_upstream = (up.symbol == "FGFR2")) %>%
    arrange(desc(C_trunc), desc(FGFR2_is_upstream), count)
  
  mat_star_RE[i,1] = do.call("paste", c(t_star_re[, c("up.symbol", "down.symbol", "up.id", "down.id")], sep = "_"))[1]
  mat_star_RE[i,2] = paste(do.call("paste", c(t_star_re[, c("up.symbol", "down.symbol", "up.id", "down.id")], sep = "_")),
                           collapse = " / ")
  
}

mat_SV = cbind(mat_purple_SV, mat_linx_fusion, mat_star_fusion[,3], mat_star_RE[,2])
colnames(mat_SV) = c("major_SV", "major_SV_type", "all_SV", "all_SV_type", "LINX", "STAR_Fusion", "STAR_RE")
rownames(mat_SV) = ind_tar


################################################################################
# combined results from LINX
################################################################################
FGFR2_LINX_ALL = c()
for(i in 1:length(ind_tar)){
  t_linx = FGFR2_SV[[ind_tar[i]]]$linx_fus
  if(is.null(t_linx)){
    t_linx = data.frame(NULL)
  }
  
  if(nrow(t_linx)!=0){
    t_linx = t_linx %>% arrange(desc(FusedExonUp)) %>% select(Name, FusedExonUp, FusedExonDown)
    t_linx = do.call("paste", c(t_linx, sep = "_"))
    
    t_purple = FGFR2_SV[[ind_tar[i]]]$purple_SV_annot
    t_star_fus = FGFR2_SV[[ind_tar[i]]]$STAR_Fusion
    
    t_df = c()
    for(j in 1:length(t_linx)){
      # comparing with PURPLE
      if(nrow(t_purple)!=0){
        t_lab = t_linx[j]
        t_df_2 = data.frame(t_lab, t_purple[match(t_linx[j], t_purple$LINX_Fusion), c(8:11, 20:23)])
      } else {
        t_df_2 = data.frame(t_linx[j], rep(NA, 8))
      }

      # comparing with STAR-Fusion
      if(nrow(t_star_fus)!=0){
        t_lab = paste(stringr::str_split_fixed(t_linx[j], "_", 4)[, 1:2], collapse = "_")
        t_lab = gsub("TENC1", "TNS2", t_lab)
        t_df_3 = t_star_fus[match(t_lab, gsub("--", "_", t_star_fus$X.FusionName)), c("left_annot", "right_annot")]
      } else {
        t_df_3 = data.frame(matrix(NA, ncol = 2, nrow = 1))
      }
      t_df_comb = data.frame(t_df_2, t_df_3)
      colnames(t_df_comb) = c("LINX_Fusion", "up.loc_in_trx", "up.down_exon", "up.down_exon_frame", "up.RE_category", 
                           "down.loc_in_trx", "down.down_exon", "down.down_exon_frame", "down.RE_category", 
                           "STAR_Fusion.Left_annot", "STAR_Fusion.Right_annot")
      
      t_df_comb = data.frame(sample_id = ind_tar[[i]], t_df_comb)
      t_df = rbind(t_df, t_df_comb)
    }
    
    FGFR2_LINX_ALL = rbind(FGFR2_LINX_ALL, t_df)
    
  }
}

save(hg19_info, hg38_info, FGFR2_SV, FGFR2_SV_n, FGFR2_PURPLE_ALL, FGFR2_LINX_ALL, list_linx, list_RNA, target_in_linx, target_in_RNA, file = "~/FGFR/Daniel/R/R_files/HMF_data_SV_RNA_v2.RData")

