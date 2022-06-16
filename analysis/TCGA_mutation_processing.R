library(dplyr)

################################################################################
# mutation
################################################################################
FGFR2_alt = read.table("~/FGFR/Daniel/R/Nature_figures/data/TCGA/cBioportal_alterations_across_samples.tsv", header = T, stringsAsFactors = F, sep = "\t")
FGFR2_mut = read.table("~/FGFR/Daniel/R/Nature_figures/data/TCGA/cBioportal/mutations.txt", header = T, stringsAsFactors = F, sep = "\t")
FGFR2_mut$FGFR2[which(FGFR2_mut$SAMPLE_ID %in% FGFR2_alt$Sample.ID[which(FGFR2_alt$FGFR2..MUT == "no alteration")])] = "--"
FGFR2_mut_info = read.table("~/FGFR/Daniel/R/Nature_figures/data/TCGA/cBioportal_mutation_info.tsv", header = T, stringsAsFactors = F, sep = "\t", fill = T, comment.char = "")

# samples with multiple types of mutations
t = FGFR2_mut_info %>%
  mutate(Mutation.Type = replace(Mutation.Type, stringr::str_detect(Mutation.Type, "Frame"), "Frame_Shift_Mut")) %>%
  mutate(Mutation.Type = replace(Mutation.Type, stringr::str_detect(Mutation.Type, "Splice"), "Splice_Mut")) %>%
  mutate(Mutation.Type = gsub("Mutation", "Mut", Mutation.Type)) %>%
  mutate(putative_driver = ifelse(Annotation != "OncoKB: Unknown, level NA;CIViC: NA;MyCancerGenome: not present;CancerHotspot: no;3DHotspot: no" &
                                    Annotation != "OncoKB: Likely Neutral, level NA;CIViC: NA;MyCancerGenome: not present;CancerHotspot: no;3DHotspot: no",
                                  "Hotspot_Oncogenic_Mut", NA)) %>%
  mutate(Mut_loc = ifelse(Exon == "18/18", "E18", "E1-E17")) %>% 
  mutate(Mut_class = case_when(Mut_loc == "E18" & (Mutation.Type == "Nonsense_Mut" | Mutation.Type == "Frame_Shift_Mut") ~ "E18_truncation",
                               Mut_loc == "E1-E17" & putative_driver == "Hotspot_Oncogenic_Mut" ~ "Hotspot/Oncogenic"),
         Mut_class = ifelse(is.na(Mut_class), "Others", Mut_class)) %>%
  select(Mutation.Type, putative_driver, Mut_loc, Mut_class, Sample.ID) %>% 
  distinct()

t_lab = rep(NA, nrow(FGFR2_mut))
for(i in 1:nrow(FGFR2_mut)){
  t_mut_info_1 = FGFR2_mut[i,]
  t_mut_info_2 = t %>% filter(Sample.ID == FGFR2_mut$SAMPLE_ID[i]) %>% arrange(Mut_class)
  t_mut_info_2_sig = t_mut_info_2 %>% filter(Mut_class != "Others")
  if(nrow(t_mut_info_2)>=2){
    if(nrow(t_mut_info_2_sig)>=1){
      t_lab[i] = paste(t_mut_info_2_sig$Mut_class, collapse = " & ")
    } else {
      t_lab[i] = "Others"
    }
  } else if(nrow(t_mut_info_2)==1){
    t_lab[i] = t_mut_info_2$Mut_class
  } else {
    t_lab[i] = NA
  }
  print(i)
}

FGFR2_mut$Mutation_class = t_lab

save(FGFR2_mut, FGFR2_mut_info, file = "~/FGFR/Daniel/R/Nature_figures/data/TCGA/TCGA_mutation_processing.RData")
