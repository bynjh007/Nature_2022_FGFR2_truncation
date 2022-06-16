library(BSgenome.Mmusculus.UCSC.mm10)
library(seqinr)
library(dplyr)

# chr7 fasta file
chr7_seq = read.fasta("~/Resources/genome/Chr7_GRCm38.fa", seqtype = "DNA")
name = unlist(strsplit(attr(chr7_seq[[1]], "Annot"), ":"))
chr7_seq = paste(as.character(chr7_seq[[1]]), collapse = "")

# chr7 gtf file
chr7_gtf = read.table("~/Resources/genome/Chr7_GRCm38.gtf", header = F, stringsAsFactors = F, sep = "\t")
chr7_gtf = chr7_gtf %>% filter(!(grepl("ENSMUST00000124096", V9))) # remove the longest Fgfr2 isoform
chr7_gtf = chr7_gtf %>% mutate(V5 = replace(V5, V3 == "gene" & grepl("Fgfr2", V9), max(as.numeric(chr7_gtf %>% filter(grepl("Fgfr2", V9), V3 != "gene") %>% pull(V5)))))

# chr7 fasta file for antisense insertion
sb_seq_anti = "ttcagccgatgatgaaattgccgcactggttgttagcaacgtagccggtatgtgaaagatggattcgcgggaatttagtggatcccccgggctgcaggaattcgatctgaagcctatagagtacgagccatagataaaataaaagattttatttagtctccagaaaaaggggggaatgaaagaccccacctgtaggtttggcaagctagcttaagtaacgccattttgcaaggcatggaaaatacataactgagaatagagaagttcagatcaaggttaggaacagagagacagcagaatatgggccaaacaggatatctgtggtaagcagttcctgccccggctcagggccaagaacagatggtccccagatgcggtcccgccctcagcagtttctagagaaccatcagatgtttccagggtgccccaaggacctgaaaatgaccctgtgccttatttgaactaaccaatcagttcgcttctcgcttctgttcgcgcgcttctgctccccgagctcaataaaagagcccacaacccctcactcggcgcgccagtcctccgatagactgcgtcgcccatcaagcttgctactagcaccagaacgcccgcgaggatctctcaggtaataaagagcgccaaggctggctgcaagcggagcctctgagagcctctgagggccagggctactgcacccttggtcctcaacgctggggtcttcagaactagaatgctgggggtggggtggggattcggttccctattccatcgcgcgttaagatacattgatgagtttggacaaaccacaactagaatgcagtgaaaaaaatgctttatttgtgaaatttgtgatgctattgctttatttgtaaccattataagctgcaataaacaagttaacaacaacaattgcattcattttatgtttcaggttcagggggaggtgtgggaggttttttaaagcaagtaaaacctctacaaatgtggtatggctgattatgatcagttatctagatccggtggatcccgggcccgcggtaccgtcgactgcagaattcgatgatcctctagagtccagatctgcgatctgcgttcttcttctttggttttcgggacctgggac"
sb_seq_sense = "gtcccaggtcccgaaaaccaaagaagaagaacgcagatcgcagatctggactctagaggatcatcgaattctgcagtcgacggtaccgcgggcccgggatccaccggatctagataactgatcataatcagccataccacatttgtagaggttttacttgctttaaaaaacctcccacacctccccctgaacctgaaacataaaatgaatgcaattgttgttgttaacttgtttattgcagcttataatggttacaaataaagcaatagcatcacaaatttcacaaataaagcatttttttcactgcattctagttgtggtttgtccaaactcatcaatgtatcttaacgcgcgatggaatagggaaccgaatccccaccccacccccagcattctagttctgaagaccccagcgttgaggaccaagggtgcagtagccctggccctcagaggctctcagaggctccgcttgcagccagccttggcgctctttattacctgagagatcctcgcgggcgttctggtgctagtagcaagcttgatgggcgacgcagtctatcggaggactggcgcgccgagtgaggggttgtgggctcttttattgagctcggggagcagaagcgcgcgaacagaagcgagaagcgaactgattggttagttcaaataaggcacagggtcattttcaggtccttggggcaccctggaaacatctgatggttctctagaaactgctgagggcgggaccgcatctggggaccatctgttcttggccctgagccggggcaggaactgcttaccacagatatcctgtttggcccatattctgctgtctctctgttcctaaccttgatctgaacttctctattctcagttatgtattttccatgccttgcaaaatggcgttacttaagctagcttgccaaacctacaggtggggtctttcattcccccctttttctggagactaaataaaatcttttattttatctatggctcgtactctataggcttcagatcgaattcctgcagcccgggggatccactaaattcccgcgaatccatctttcacataccggctacgttgctaacaaccagtgcggcaatttcatcatcggctgaa"
######################################################
# Fgfr2 exon information
######################################################
Fgfr2_exon_info = read.table("~/FGFR/Daniel/sb/NM_010207.2.txt", header = F, stringsAsFactors = F, sep = "\t")
colnames(Fgfr2_exon_info) = c("chr", "start", "end", "unit")
Fgfr2_trx_info = c()
for(i in 1:(nrow(Fgfr2_exon_info)-1)){
  int_s = Fgfr2_exon_info[i+1,3]+1
  int_e = Fgfr2_exon_info[i,2]-1
  Fgfr2_trx_info = rbind(Fgfr2_trx_info, Fgfr2_exon_info[i,], 
                         data.frame(chr = "chr7", start = int_s, end = int_e, unit = paste("Intron", i)))
}
Fgfr2_trx_info = rbind(Fgfr2_trx_info, tail(Fgfr2_exon_info,1))
Fgfr2_trx_info$width = Fgfr2_trx_info$end - Fgfr2_trx_info$start


######################################################
# insertion data
######################################################
options(stringsAsFactors = F)
Fgfr2_insertion = xlsx::read.xlsx2("/DATA/projects/j.bhin/Daniel_FGFR2/sb/Fgfr2 Transposon Insertions.xlsx", sheetIndex = 1, header = T)
Fgfr2_insertion = Fgfr2_insertion %>% filter(gene_name == "Fgfr2")
Fgfr2_insertion$position = as.numeric(as.character(Fgfr2_insertion$position))
Fgfr2_insertion$sample = gsub("13SKA008-Leg-right", "13SKA008-leg-right", Fgfr2_insertion$sample)

Fgfr2_insertion_RNAseq = xlsx::read.xlsx2("/DATA/projects/j.bhin/Daniel_FGFR2/sb/Fgfr2 Transposon Insertions.xlsx", sheetIndex = 5, header = T)
Fgfr2_insertion_RNAseq$Sample = gsub("13SKA008-Leg-right", "13SKA008-leg-right", Fgfr2_insertion_RNAseq$Sample)

######################################################
# antisense insertion
######################################################
samp_antisense_I17 = Fgfr2_insertion %>% filter(gene_orientation == "antisense", position >= Fgfr2_trx_info$start[34], position <= Fgfr2_trx_info$end[34])
samp_sense_I17 = Fgfr2_insertion %>% filter(gene_orientation == "sense", position >= Fgfr2_trx_info$start[34], position <= Fgfr2_trx_info$end[34])

uniq_samp_antisense = unique(as.character(samp_antisense_I17$sample))
uniq_samp_sense = unique(as.character(samp_sense_I17$sample))

uniq_both = intersect(uniq_samp_antisense, uniq_samp_sense)
uniq_samp_antisense = setdiff(uniq_samp_antisense, uniq_both)
uniq_samp_sense = setdiff(uniq_samp_sense, uniq_both)

dir.create("/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/", showWarnings = FALSE)
for(x in uniq_samp_antisense){
  t_path = paste("/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/", x, sep = "")
  dir.create(t_path, showWarnings = FALSE)
  dir.create(paste(t_path, "/reference", sep = ""), showWarnings = FALSE)
  
  # fasta file
  t_name = name
  pos_ins = sort(samp_antisense_I17 %>% filter(sample %in% x) %>% pull(position))
  chr7_seq_sb = chr7_seq
  n_gap = 0
  for(j in 1:length(pos_ins)){
    chr7_seq_sb = paste(substr(chr7_seq_sb, 1, pos_ins[j] + n_gap), sb_seq_anti, substr(chr7_seq_sb, pos_ins[j] + n_gap + 1, nchar(chr7_seq_sb)), sep = "")
    n_gap = nchar(sb_seq_anti) + n_gap
  }
  
  t_name[6] = nchar(chr7_seq_sb)
  t_name = paste(t_name, collapse = ":")
  t_name = gsub(">", "", t_name)
  
  write.fasta(sequences = chr7_seq_sb, names = t_name, as.string = T, file.out =  paste(t_path, "/reference/chr7_SB.fasta", sep = ""))
  
  # gtf file
  ind_start = which(chr7_gtf$V4 > max(pos_ins))
  ind_end = which(chr7_gtf$V5 > max(pos_ins))
  chr7_gtf_SB = chr7_gtf
  chr7_gtf_SB$V4[ind_start] = chr7_gtf_SB$V4[ind_start] + n_gap
  chr7_gtf_SB$V5[ind_end] = chr7_gtf_SB$V5[ind_end] + n_gap
  
  # add SB location in the gtf
  transposon_gtf = chr7_gtf_SB[grep("ENSMUST00000122054", chr7_gtf_SB$V9),]
  transposon_gtf$V9 = gsub("ENSMUST00000122054", "SB_transposon", transposon_gtf$V9)
  E17_ind = grep("exon_number 17", transposon_gtf$V9)[2]
  transposon_gtf_E1_E17 = transposon_gtf[1:E17_ind,]
  transposon_gtf_E18 = transposon_gtf[(E17_ind+1):nrow(transposon_gtf),]
  transposon_gtf_E18$V9 = gsub("exon_number 18", paste("exon_number", 18+length(pos_ins)), transposon_gtf_E18$V9)
  
  t_loc = stringr::str_locate_all(string = chr7_seq_sb, pattern = sb_seq_anti)
  transposon_gtf_sb = c()
  for(j in 1:length(pos_ins)){
    transposon_gtf_sb = rbind(transposon_gtf_sb, 
                        data.frame(V1 = "7", V2 = "protein_coding", V3 = "exon", V4 = t_loc[[1]][j,1], V5 = t_loc[[1]][j,2],
                                   V6 = ".", V7 = "-", "V8" = ".", 
                                   V9 = "gene_id ENSMUSG00000030849; transcript_id SB_transposon; exon_number XX; gene_name Fgfr2; gene_source ensembl_havana; gene_biotype protein_coding; transcript_name Fgfr2-SB; transcript_source ensembl_havana; tag CCDS; ccds_id CCDS52413; exon_id transposon;"),
                        data.frame(V1 = "7", V2 = "protein_coding", V3 = "CDS", V4 = t_loc[[1]][j,1], V5 = t_loc[[1]][j,2],
                                   V6 = ".", V7 = "-", "V8" = ".", 
                                   V9 = "gene_id ENSMUSG00000030849; transcript_id SB_transposon; exon_number XX; gene_name Fgfr2; gene_source ensembl_havana; gene_biotype protein_coding; transcript_name Fgfr2-SB; transcript_source ensembl_havana; tag CCDS; ccds_id CCDS52413; exon_id transposon;"))
    t_exon_number = 17+j
    transposon_gtf_sb$V9 = gsub("exon_number XX", paste("exon_number", t_exon_number), transposon_gtf_sb$V9)
  }
  
  transposon_gtf = rbind(transposon_gtf_E1_E17, transposon_gtf_sb, transposon_gtf_E18)
  chr7_gtf_SB = rbind(chr7_gtf_SB, transposon_gtf)
  
  write.table(chr7_gtf_SB, paste(t_path, "/reference/chr7_SB.gtf", sep = ""), col.names = F, row.names = F, sep = "\t", quote = F)

}

write.table(uniq_samp_antisense, "/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/samp_antisense.txt", col.names = F, row.names = F, quote = F)

conda activate /DATA/projects/j.bhin/Daniel_FGFR2/pipeline/.snakemake/conda/00312989 # same STAR version
input="/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/samp_antisense.txt"
while read line; do
fastq=$(ls /DATA/projects/ilc-sb-screen/analysis/data/interim/sb/rnaseq/fastq/downloaded/ | grep $line | grep -v md5)
STAR --runThreadN 24 --runMode genomeGenerate --genomeSAindexNbases 12 --sjdbOverhang 50 --sjdbGTFfile /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference/chr7_SB.gtf --genomeDir /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference --genomeFastaFiles /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference/chr7_SB.fasta
STAR --runThreadN 12 --genomeDir /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference --readFilesIn /DATA/projects/ilc-sb-screen/analysis/data/interim/sb/rnaseq/fastq/downloaded/$fastq --outFileNamePrefix /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/ --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c
samtools index /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/Aligned.sortedByCoord.out.bam
samtools faidx /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference/chr7_SB.fasta
echo $line
done < "$input"


######################################################
# sense insertion
######################################################
for(x in uniq_samp_sense){
  t_path = paste("/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/", x, sep = "")
  dir.create(t_path, showWarnings = FALSE)
  dir.create(paste(t_path, "/reference", sep = ""), showWarnings = FALSE)
  
  # fasta file
  t_name = name
  pos_ins = sort(samp_sense_I17 %>% filter(sample %in% x) %>% pull(position))
  chr7_seq_sb = chr7_seq
  n_gap = 0
  for(j in 1:length(pos_ins)){
    chr7_seq_sb = paste(substr(chr7_seq_sb, 1, pos_ins[j] + n_gap), sb_seq_sense, substr(chr7_seq_sb, pos_ins[j] + n_gap + 1, nchar(chr7_seq_sb)), sep = "")
    n_gap = nchar(sb_seq_sense) + n_gap
  }
  
  t_name[6] = nchar(chr7_seq_sb)
  t_name = paste(t_name, collapse = ":")
  t_name = gsub(">", "", t_name)
  
  write.fasta(sequences = chr7_seq_sb, names = t_name, as.string = T, file.out =  paste(t_path, "/reference/chr7_SB.fasta", sep = ""))
  
  # gtf file
  ind_start = which(chr7_gtf$V4 > max(pos_ins))
  ind_end = which(chr7_gtf$V5 > max(pos_ins))
  chr7_gtf_SB = chr7_gtf
  chr7_gtf_SB$V4[ind_start] = chr7_gtf_SB$V4[ind_start] + n_gap
  chr7_gtf_SB$V5[ind_end] = chr7_gtf_SB$V5[ind_end] + n_gap
  
  # add SB location in the gtf
  transposon_gtf = chr7_gtf_SB[grep("ENSMUST00000122054", chr7_gtf_SB$V9),]
  transposon_gtf$V9 = gsub("ENSMUST00000122054", "SB_transposon", transposon_gtf$V9)
  E17_ind = grep("exon_number 17", transposon_gtf$V9)[2]
  transposon_gtf_E1_E17 = transposon_gtf[1:E17_ind,]
  transposon_gtf_E18 = transposon_gtf[(E17_ind+1):nrow(transposon_gtf),]
  transposon_gtf_E18$V9 = gsub("exon_number 18", paste("exon_number", 18+length(pos_ins)), transposon_gtf_E18$V9)
  
  t_loc = stringr::str_locate_all(string = chr7_seq_sb, pattern = sb_seq_sense)
  transposon_gtf_sb = c()
  for(j in 1:length(pos_ins)){
    transposon_gtf_sb = rbind(transposon_gtf_sb, 
                              data.frame(V1 = "7", V2 = "protein_coding", V3 = "exon", V4 = t_loc[[1]][j,1], V5 = t_loc[[1]][j,2],
                                         V6 = ".", V7 = "-", "V8" = ".", 
                                         V9 = "gene_id ENSMUSG00000030849; transcript_id SB_transposon; exon_number XX; gene_name Fgfr2; gene_source ensembl_havana; gene_biotype protein_coding; transcript_name Fgfr2-SB; transcript_source ensembl_havana; tag CCDS; ccds_id CCDS52413; exon_id transposon;"),
                              data.frame(V1 = "7", V2 = "protein_coding", V3 = "CDS", V4 = t_loc[[1]][j,1], V5 = t_loc[[1]][j,2],
                                         V6 = ".", V7 = "-", "V8" = ".", 
                                         V9 = "gene_id ENSMUSG00000030849; transcript_id SB_transposon; exon_number XX; gene_name Fgfr2; gene_source ensembl_havana; gene_biotype protein_coding; transcript_name Fgfr2-SB; transcript_source ensembl_havana; tag CCDS; ccds_id CCDS52413; exon_id transposon;"))
    t_exon_number = 17+j
    transposon_gtf_sb$V9 = gsub("exon_number XX", paste("exon_number", t_exon_number), transposon_gtf_sb$V9)
  }
  
  transposon_gtf = rbind(transposon_gtf_E1_E17, transposon_gtf_sb, transposon_gtf_E18)
  chr7_gtf_SB = rbind(chr7_gtf_SB, transposon_gtf)
  
  write.table(chr7_gtf_SB, paste(t_path, "/reference/chr7_SB.gtf", sep = ""), col.names = F, row.names = F, sep = "\t", quote = F)
  
}

write.table(uniq_samp_sense, "/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/samp_sense.txt", col.names = F, row.names = F, quote = F)

conda activate /DATA/projects/j.bhin/Daniel_FGFR2/pipeline/.snakemake/conda/00312989 # same STAR version
input="/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/samp_sense.txt"
while read line; do
fastq=$(ls /DATA/projects/ilc-sb-screen/analysis/data/interim/sb/rnaseq/fastq/downloaded/ | grep $line | grep -v md5)
STAR --runThreadN 24 --runMode genomeGenerate --genomeSAindexNbases 12 --sjdbOverhang 50 --sjdbGTFfile /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference/chr7_SB.gtf --genomeDir /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference --genomeFastaFiles /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference/chr7_SB.fasta
STAR --runThreadN 12 --genomeDir /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference --readFilesIn /DATA/projects/ilc-sb-screen/analysis/data/interim/sb/rnaseq/fastq/downloaded/$fastq --outFileNamePrefix /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/ --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c
samtools index /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/Aligned.sortedByCoord.out.bam
samtools faidx /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference/chr7_SB.fasta
echo $line
done < "$input"

######################################################
# antisense + sense insertion
######################################################

for(x in uniq_both){
  t_path = paste("/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/", x, sep = "")
  dir.create(t_path, showWarnings = FALSE)
  dir.create(paste(t_path, "/reference", sep = ""), showWarnings = FALSE)
  
  # fasta file
  t_name = name
  pos_ins = rbind(samp_antisense_I17 %>% filter(sample %in% x),
                  samp_sense_I17 %>% filter(sample %in% x))
  pos_ins = pos_ins %>% arrange(position)
  chr7_seq_sb = chr7_seq
  n_gap = 0
  for(j in 1:nrow(pos_ins)){
    t_sb_seq = ifelse(pos_ins$gene_orientation[j] == "antisense", sb_seq_anti, sb_seq_sense)
    chr7_seq_sb = paste(substr(chr7_seq_sb, 1, pos_ins$position[j] + n_gap), t_sb_seq, substr(chr7_seq_sb, pos_ins$position[j] + n_gap + 1, nchar(chr7_seq_sb)), sep = "")
    n_gap = nchar(t_sb_seq) + n_gap
  }
  
  t_name[6] = nchar(chr7_seq_sb)
  t_name = paste(t_name, collapse = ":")
  t_name = gsub(">", "", t_name)
  
  write.fasta(sequences = chr7_seq_sb, names = t_name, as.string = T, file.out =  paste(t_path, "/reference/chr7_SB.fasta", sep = ""))
  
  # gtf file
  ind_start = which(chr7_gtf$V4 > max(pos_ins$position))
  ind_end = which(chr7_gtf$V5 > max(pos_ins$position))
  chr7_gtf_SB = chr7_gtf
  chr7_gtf_SB$V4[ind_start] = chr7_gtf_SB$V4[ind_start] + n_gap
  chr7_gtf_SB$V5[ind_end] = chr7_gtf_SB$V5[ind_end] + n_gap
  
  # add SB location in the gtf
  transposon_gtf = chr7_gtf_SB[grep("ENSMUST00000122054", chr7_gtf_SB$V9),]
  transposon_gtf$V9 = gsub("ENSMUST00000122054", "SB_transposon", transposon_gtf$V9)
  E17_ind = grep("exon_number 17", transposon_gtf$V9)[2]
  transposon_gtf_E1_E17 = transposon_gtf[1:E17_ind,]
  transposon_gtf_E18 = transposon_gtf[(E17_ind+1):nrow(transposon_gtf),]
  transposon_gtf_E18$V9 = gsub("exon_number 18", paste("exon_number", 18+length(pos_ins)), transposon_gtf_E18$V9)
  
  t_loc_anti = stringr::str_locate_all(string = chr7_seq_sb, pattern = sb_seq_anti)
  t_loc_sense = stringr::str_locate_all(string = chr7_seq_sb, pattern = sb_seq_sense)
  t_loc = rbind(data.frame(t_loc_anti[[1]], type = "antisense"),
                data.frame(t_loc_sense[[1]], type = "sense"))
  t_loc = t_loc %>% arrange(start)
  transposon_gtf_sb = c()
  for(j in 1:nrow(t_loc)){
    transposon_gtf_sb = rbind(transposon_gtf_sb, 
                              data.frame(V1 = "7", V2 = "protein_coding", V3 = "exon", V4 = t_loc[j,1], V5 = t_loc[j,2],
                                         V6 = ".", V7 = "-", "V8" = ".", 
                                         V9 = "gene_id ENSMUSG00000030849; transcript_id SB_transposon; exon_number XX; gene_name Fgfr2; gene_source ensembl_havana; gene_biotype protein_coding; transcript_name Fgfr2-SB; transcript_source ensembl_havana; tag CCDS; ccds_id CCDS52413; exon_id transposon;"),
                              data.frame(V1 = "7", V2 = "protein_coding", V3 = "CDS", V4 = t_loc[j,1], V5 = t_loc[j,2],
                                         V6 = ".", V7 = "-", "V8" = ".", 
                                         V9 = "gene_id ENSMUSG00000030849; transcript_id SB_transposon; exon_number XX; gene_name Fgfr2; gene_source ensembl_havana; gene_biotype protein_coding; transcript_name Fgfr2-SB; transcript_source ensembl_havana; tag CCDS; ccds_id CCDS52413; exon_id transposon;"))
    t_exon_number = 17+j
    transposon_gtf_sb$V9 = gsub("exon_number XX", paste("exon_number", t_exon_number), transposon_gtf_sb$V9)
    transposon_gtf_sb$V9 = gsub("exon_id transposon", paste("exon_id transposon", t_loc$type[j], sep = "_"), transposon_gtf_sb$V9)
  }
  
  transposon_gtf = rbind(transposon_gtf_E1_E17, transposon_gtf_sb, transposon_gtf_E18)
  chr7_gtf_SB = rbind(chr7_gtf_SB, transposon_gtf)
  
  write.table(chr7_gtf_SB, paste(t_path, "/reference/chr7_SB.gtf", sep = ""), col.names = F, row.names = F, sep = "\t", quote = F)
  
}


write.table(uniq_both, "/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/samp_both.txt", col.names = F, row.names = F, quote = F)

conda activate /DATA/projects/j.bhin/Daniel_FGFR2/pipeline/.snakemake/conda/00312989 # same STAR version
input="/DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/samp_both.txt"
while read line; do
fastq=$(ls /DATA/projects/ilc-sb-screen/analysis/data/interim/sb/rnaseq/fastq/downloaded/ | grep $line | grep -v md5)
STAR --runThreadN 24 --runMode genomeGenerate --genomeSAindexNbases 12 --sjdbOverhang 50 --sjdbGTFfile /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference/chr7_SB.gtf --genomeDir /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference --genomeFastaFiles /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference/chr7_SB.fasta
STAR --runThreadN 12 --genomeDir /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference --readFilesIn /DATA/projects/ilc-sb-screen/analysis/data/interim/sb/rnaseq/fastq/downloaded/$fastq --outFileNamePrefix /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/ --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c
samtools index /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/Aligned.sortedByCoord.out.bam
samtools faidx /DATA/projects/j.bhin/Daniel_FGFR2/sb/sashimi_realign/$line/reference/chr7_SB.fasta
echo $line
done < "$input"

save.image("~/FGFR/Daniel/R/Nature_figures/data/SB_transposon/SB_transposon_usage.RData")

