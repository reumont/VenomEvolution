library(stringr)
library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(gggenes)
library(lamisc)
library(tibble)
library(gridExtra)
options(scipen = 999)


setwd("/Users/marianavelasque/Dropbox/Pesquisa/Post-doc/Ivan_research/200722_venom_transcriptome/201215_ortholog_synteny/210122_melittin_synteny")


blast_orthologs<- read_table2("210120_melittin_flank/blast_results_flank_filtered.blast", 
                              col_names = FALSE) %>%
   select(other_spp = X1, apis_mellifera =X2 ) 

blast_orthologs<- blast_orthologs %>%
   separate( other_spp, into = c("gene_desc", "other_spp", "duplicate"), sep = "_") %>%
   select("gene_desc", "other_spp", apis_mellifera)

blast_orthologs$other_spp <- paste(blast_orthologs$gene_desc, "_", blast_orthologs$other_spp)

blast_orthologs$other_spp<- sub(" ", "", blast_orthologs$other_spp)
blast_orthologs$other_spp<- sub(" ", "", blast_orthologs$other_spp)
blast_orthologs$gene_desc = NULL

blast_orthologs<- blast_orthologs %>%
   separate( apis_mellifera, into = c("gene_desc", "apis_mellifera", "duplicate"), sep = "_") %>%
   select("gene_desc", "other_spp", apis_mellifera)

blast_orthologs$apis_mellifera <- paste(blast_orthologs$gene_desc, "_", blast_orthologs$apis_mellifera)

blast_orthologs$apis_mellifera<- sub(" ", "", blast_orthologs$apis_mellifera)
blast_orthologs$apis_mellifera<- sub(" ", "", blast_orthologs$apis_mellifera)
blast_orthologs$gene_desc = NULL

#flanking_genes<- read_table2("apis_flanking_genes_filtered_mcscan", 
#                             col_names = FALSE)
#blast_orthologs<- blast_orthologs %>%
#   filter(apis_mellifera %in% flanking_genes$X1)

##get the species that have an ortholog
species<- read_table2("201231_mcscan_files/melittin.collinearity_filtered_species", 
                      col_names = FALSE)

species<- unique(c(species$X3))

prepare_bedfile<- function(bed_file){
   bed_file<- bed_file %>%
      filter(X8 == "CDS")%>%
      select(scafolds = X1, start = X2, end = X3, gene = X4, direction = X6, X8) %>%
      mutate(gene = str_replace(gene, "cds-", "")) %>%
      filter(gene %in% blast_orthologs$other_spp )
   
   bed_file<- bed_file %>%
      left_join(blast_orthologs, by = c("gene"="other_spp")) %>%
      select(!gene) %>%
      select( scafolds, start, end, direction,  gene = "apis_mellifera") 
   
   genes_insterest<- c("NP_001011607.1")
   
   scafolds<- bed_file %>% 
      filter(gene %in% genes_insterest) %>% 
      select(scafolds)
   
   scafolds_filt<- unique(c(scafolds$scafolds))
   
   tally_bed<- bed_file %>%     
      distinct(scafolds,  gene, .keep_all = TRUE) %>%
      group_by(scafolds) %>%
     dplyr:: summarise(n = n()) %>%
      filter(n>2)

   bed_file <- bed_file %>%
      distinct(scafolds, gene, .keep_all = TRUE) %>%
      filter(scafolds %in% tally_bed$scafolds) %>%
      distinct(scafolds,  gene, .keep_all = TRUE) %>%
      rank_in_group2(group_var = scafolds,
                     arrange_var = start) %>%
      mutate(start = rank*10, 
             end = start+10) %>%
      select(!rank)
   
}


apis_mellifera<- read_table2("GCF_003254395.2_Amel_HAv3._1_merged.bed",
                             col_names = F,comment = "#") %>%
   prepare_bedfile()%>%
   mutate(species = c("Apis mellifera"))  

apis_florea<- read_table2("GCF_000184785.3_Aflo_1.1_cds_filtered_merged.bed",
                          col_names = F,comment = "#") %>%
   prepare_bedfile() %>%
   mutate(species = c("Apis florea"))  


Solenopsis_invicta<- read_table2("GCF_000188075.2_Si_gnH_cds_filtered_merged.bed",
                                 col_names = F,comment = "#")  %>%
   prepare_bedfile() %>%
   mutate(species = c("Solenopsis invicta")) 


Megachile_rotundata<- read_table2("GCF_000220905.1_MROT_cds_filtered_merged.bed",
                                  col_names = F,comment = "#") %>%
   prepare_bedfile() %>%
   mutate(species = c("Megachile rotundata"))  


Habropoda_laboriosa<- read_table2("GCF_001263275.1_ASM126327v1_cds_filtered_merged.bed",
                                  col_names = F,comment = "#")  %>%
   prepare_bedfile() %>%
   mutate(species = c("Habropoda laboriosa")) 


Polistes_canadensis<- read_table2("GCF_001263275.1_ASM126327v1_cds_filtered_merged.bed",
                                  col_names = F,comment = "#") %>%
   prepare_bedfile() %>%
   mutate(species = c("Polistes canadensis"))  

Polistes_dominula<- read_table2("GCF_001465965.1_Pdom_cds_filtered_merged.bed",
                                col_names = F,comment = "#")  %>%
   prepare_bedfile() %>%
   mutate(species = c("Polistes dominula"))  



## select color for data set 
library(RColorBrewer)

color_genes <- as.data.frame(matrix(ncol=length(unique(apis_mellifera$gene)), nrow=1)) #create data table 
names(color_genes)=unique(apis_mellifera$gene)

color_genes[1,]= brewer.pal(length(unique(apis_mellifera$gene)), "Dark2")




apis_synteny<- rbind( apis_mellifera, Solenopsis_invicta,  Polistes_dominula,
                      Megachile_rotundata, Habropoda_laboriosa, 
                      Polistes_canadensis)

apis_synteny %>%
   filter(species == "Apis mellifera") %>%
   ggplot( aes(xmin = start, xmax = end, y = scafolds, fill = gene, color = gene)) +
   geom_gene_arrow() +
   #xlim(0, apis_synteny$size[apis_synteny$species == "Apis mellifera"][1]) +
   labs(y ='\n\nApis mellifera\n\n\n') +
   scale_color_manual(values=color_genes,
                      name = "Genes") +
   scale_fill_manual(values=color_genes,
                     name = "Genes") +
   theme_genes()+
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12,face="bold"))

ggsave("../../Figures_genomic_locations/Melittin_plot_filtered/apis_melifera_melittin.png", height=8, width=8, dpi=800)


apis_synteny %>%
   dplyr::filter(species == "Habropoda laboriosa") %>%
   ggplot( aes(xmin = start, xmax = end, y = scafolds, fill = gene, color = gene)) +
   geom_gene_arrow() +
   #xlim(0, apis_synteny$size[apis_synteny$species == "Apis mellifera"][1]) +
   labs(y ='\n\nHabropoda laboriosa\n\n\n') +
   scale_color_manual(values=color_genes,
                      name = "Genes") +
   scale_fill_manual(values=color_genes,
                     name = "Genes") +
   theme_genes()+
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12,face="bold"))

ggsave("../../Figures_genomic_locations/Melittin_plot_filtered/Habropoda_laboriosa_melittin.png", height=8, width=8, dpi=800)

apis_synteny %>%
   filter(species == "Megachile rotundata") %>%
   ggplot( aes(xmin = start, xmax = end, y = scafolds, fill = gene, color = gene)) +
   geom_gene_arrow() +
   #xlim(0, apis_synteny$size[apis_synteny$species == "Apis mellifera"][1]) +
   labs(y ='\n\nMegachile rotundata\n\n\n') +
   scale_color_manual(values=color_genes,
                      name = "Genes") +
   scale_fill_manual(values=color_genes,
                     name = "Genes") +
   theme_genes()+
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12,face="bold"))
ggsave("../../Figures_genomic_locations/Melittin_plot_filtered/Megachile_rotundata_melittin.png", height=8, width=8, dpi=800)


apis_synteny %>%
   filter(species == "Polistes canadensis") %>%
   ggplot( aes(xmin = start, xmax = end, y = scafolds, fill = gene, color = gene)) +
   geom_gene_arrow() +
   #xlim(0, apis_synteny$size[apis_synteny$species == "Apis mellifera"][1]) +
   labs(y ='\n\nPolistes canadensis\n\n\n') +
   scale_color_manual(values=color_genes,
                      name = "Genes") +
   scale_fill_manual(values=color_genes,
                     name = "Genes") +
   theme_genes()+
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12,face="bold"))
ggsave("../../Figures_genomic_locations/Melittin_plot_filtered/Polistes_canadensis_melittin.png", height=8, width=8, dpi=800)

apis_synteny %>%
   filter(species == "Polistes dominula") %>%
   ggplot( aes(xmin = start, xmax = end, y = scafolds, fill = gene, color = gene)) +
   geom_gene_arrow() +
   #xlim(0, apis_synteny$size[apis_synteny$species == "Apis mellifera"][1]) +
   labs(y ='\n\nOPolistes dominula\n\n\n') +
   scale_color_manual(values=color_genes,
                      name = "Genes") +
   scale_fill_manual(values=color_genes,
                     name = "Genes") +
   theme_genes()+
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12,face="bold"))
ggsave("../../Figures_genomic_locations/Melittin_plot_filtered/Polistes_dominula_melittin.png", height=8, width=8, dpi=800)

apis_synteny %>%
   filter(species == "Solenopsis invicta") %>%
   ggplot( aes(xmin = start, xmax = end, y = scafolds, fill = gene, color = gene)) +
   geom_gene_arrow() +
   #xlim(0, apis_synteny$size[apis_synteny$species == "Apis mellifera"][1]) +
   labs(y ='\n\nSolenopsis invicta\n\n\n') +
   scale_color_manual(values=color_genes,
                      name = "Genes") +
   scale_fill_manual(values=color_genes,
                     name = "Genes") +
   theme_genes()+
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12,face="bold"))
ggsave("../../Figures_genomic_locations/Melittin_plot_filtered/Solenopsis_invicta_melittin.png", height=8, width=8, dpi=800)


prepare_bedfile<- function(bed_file){
   bed_file<- bed_file %>%
      filter(X8 == "CDS")%>%
      select(scafolds = X1, start = X2, end = X3, gene = X4, direction = X6, X8) %>%
      mutate(gene = str_replace(gene, "cds-", "")) %>%
      filter(gene %in% blast_orthologs$other_spp )
   
   bed_file<- bed_file %>%
      left_join(blast_orthologs, by = c("gene"="other_spp")) %>%
      select(!gene) %>%
      select( scafolds, start, end, direction,  gene = "apis_mellifera") 
   
   genes_insterest<- c("NP_001011607.1")
   
   scafolds<- bed_file %>% 
      filter(gene %in% genes_insterest) %>% 
      select(scafolds)
   
   scafolds_filt<- unique(c(scafolds$scafolds))
   
   bed_file <- bed_file %>%
      distinct(scafolds, gene, .keep_all = TRUE) %>%
      distinct(scafolds,  gene, .keep_all = TRUE) %>%
      rank_in_group2(group_var = scafolds,
                     arrange_var = start) %>%
      mutate(start = rank*10, 
             end = start+10) %>%
      select(!rank)
   
}


apis_mellifera<- read_table2("GCF_003254395.2_Amel_HAv3._1_merged.bed",
                             col_names = F,comment = "#") %>%
   prepare_bedfile()%>%
   mutate(species = c("Apis mellifera"))  

apis_florea<- read_table2("GCF_000184785.3_Aflo_1.1_cds_filtered_merged.bed",
                          col_names = F,comment = "#") %>%
   prepare_bedfile() %>%
   mutate(species = c("Apis florea"))  


Solenopsis_invicta<- read_table2("GCF_000188075.2_Si_gnH_cds_filtered_merged.bed",
                                 col_names = F,comment = "#")  %>%
   prepare_bedfile() %>%
   mutate(species = c("Solenopsis invicta")) 


Megachile_rotundata<- read_table2("GCF_000220905.1_MROT_cds_filtered_merged.bed",
                                  col_names = F,comment = "#") %>%
   prepare_bedfile() %>%
   mutate(species = c("Megachile rotundata"))  


Habropoda_laboriosa<- read_table2("GCF_001263275.1_ASM126327v1_cds_filtered_merged.bed",
                                  col_names = F,comment = "#")  %>%
   prepare_bedfile() %>%
   mutate(species = c("Habropoda laboriosa")) 


Polistes_canadensis<- read_table2("GCF_001263275.1_ASM126327v1_cds_filtered_merged.bed",
                                  col_names = F,comment = "#") %>%
   prepare_bedfile() %>%
   mutate(species = c("Polistes canadensis"))  

Polistes_dominula<- read_table2("GCF_001465965.1_Pdom_cds_filtered_merged.bed",
                                col_names = F,comment = "#")  %>%
   prepare_bedfile() %>%
   mutate(species = c("Polistes dominula"))  



## select color for data set 
library(RColorBrewer)

color_genes <- as.data.frame(matrix(ncol=length(unique(apis_mellifera$gene)), nrow=1)) #create data table 
names(color_genes)=unique(apis_mellifera$gene)

color_genes[1,]= brewer.pal(length(unique(apis_mellifera$gene)), "Dark2")




apis_synteny<- rbind( apis_mellifera, Solenopsis_invicta,  Polistes_dominula,
                      Megachile_rotundata, Habropoda_laboriosa, 
                      Polistes_canadensis)

apis_synteny %>%
   filter(species == "Apis mellifera") %>%
   ggplot( aes(xmin = start, xmax = end, y = scafolds, fill = gene, color = gene)) +
   geom_gene_arrow() +
   #xlim(0, apis_synteny$size[apis_synteny$species == "Apis mellifera"][1]) +
   labs(y ='\n\nApis mellifera\n\n\n') +
   scale_color_manual(values=color_genes,
                      name = "Genes") +
   scale_fill_manual(values=color_genes,
                     name = "Genes") +
   theme_genes()+
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12,face="bold"))

ggsave("../../Figures_genomic_locations/Melittin_plot_unfiltered/apis_melifera_melittin.png", height=8, width=8, dpi=800)


apis_synteny %>%
   dplyr::filter(species == "Habropoda laboriosa") %>%
   ggplot( aes(xmin = start, xmax = end, y = scafolds, fill = gene, color = gene)) +
   geom_gene_arrow() +
   #xlim(0, apis_synteny$size[apis_synteny$species == "Apis mellifera"][1]) +
   labs(y ='\n\nHabropoda laboriosa\n\n\n') +
   scale_color_manual(values=color_genes,
                      name = "Genes") +
   scale_fill_manual(values=color_genes,
                     name = "Genes") +
   theme_genes()+
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12,face="bold"))

ggsave("../../Figures_genomic_locations/Melittin_plot_unfiltered/Habropoda_laboriosa_melittin.png", height=8, width=8, dpi=800)

apis_synteny %>%
   filter(species == "Megachile rotundata") %>%
   ggplot( aes(xmin = start, xmax = end, y = scafolds, fill = gene, color = gene)) +
   geom_gene_arrow() +
   #xlim(0, apis_synteny$size[apis_synteny$species == "Apis mellifera"][1]) +
   labs(y ='\n\nMegachile rotundata\n\n\n') +
   scale_color_manual(values=color_genes,
                      name = "Genes") +
   scale_fill_manual(values=color_genes,
                     name = "Genes") +
   theme_genes()+
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12,face="bold"))
ggsave("../../Figures_genomic_locations/Melittin_plot_unfiltered/Megachile_rotundata_melittin.png", height=8, width=8, dpi=800)


apis_synteny %>%
   filter(species == "Polistes canadensis") %>%
   ggplot( aes(xmin = start, xmax = end, y = scafolds, fill = gene, color = gene)) +
   geom_gene_arrow() +
   #xlim(0, apis_synteny$size[apis_synteny$species == "Apis mellifera"][1]) +
   labs(y ='\n\nPolistes canadensis\n\n\n') +
   scale_color_manual(values=color_genes,
                      name = "Genes") +
   scale_fill_manual(values=color_genes,
                     name = "Genes") +
   theme_genes()+
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12,face="bold"))
ggsave("../../Figures_genomic_locations/Melittin_plot_unfiltered/Polistes_canadensis_melittin.png", height=8, width=8, dpi=800)

apis_synteny %>%
   filter(species == "Polistes dominula") %>%
   ggplot( aes(xmin = start, xmax = end, y = scafolds, fill = gene, color = gene)) +
   geom_gene_arrow() +
   #xlim(0, apis_synteny$size[apis_synteny$species == "Apis mellifera"][1]) +
   labs(y ='\n\nOPolistes dominula\n\n\n') +
   scale_color_manual(values=color_genes,
                      name = "Genes") +
   scale_fill_manual(values=color_genes,
                     name = "Genes") +
   theme_genes()+
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12,face="bold"))
ggsave("../../Figures_genomic_locations/Melittin_plot_unfiltered/Polistes_dominula_melittin.png", height=8, width=8, dpi=800)

apis_synteny %>%
   filter(species == "Solenopsis invicta") %>%
   ggplot( aes(xmin = start, xmax = end, y = scafolds, fill = gene, color = gene)) +
   geom_gene_arrow() +
   #xlim(0, apis_synteny$size[apis_synteny$species == "Apis mellifera"][1]) +
   labs(y ='\n\nSolenopsis invicta\n\n\n') +
   scale_color_manual(values=color_genes,
                      name = "Genes") +
   scale_fill_manual(values=color_genes,
                     name = "Genes") +
   theme_genes()+
   theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.text=element_text(size=12),
         axis.title=element_text(size=12,face="bold"))
ggsave("../../Figures_genomic_locations/Melittin_plot_unfiltered/Solenopsis_invicta_melittin.png", height=8, width=8, dpi=800)




