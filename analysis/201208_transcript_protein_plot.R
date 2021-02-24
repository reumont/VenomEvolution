setwd("~/Dropbox/Pesquisa/Post-doc/Ivan_research/200722_venom_transcriptome/scripts")

library(ggthemes)
library(tidyverse)
library(limma)
library(ggplot2)
library(readr)
library(NOISeq)
library(VennDiagram)
library(data.table)
library(magrittr)
library(kableExtra)
library(ggpubr)
library(RColorBrewer)
library(stringi)
library(dplyr)
library(viridis)
library(circlize)
library(png)
options("scipen"=100, "digits"=4)
Xylocopa_genes<- read_csv("../figures/Xylocopa_proteome_figure.csv") %>%
   group_by(protein_id) %>%
   summarise(Xylocopa_TPM = sum(Xylocopa_TPM, na.rm = TRUE))


Apis_genes<- read_csv("../figures/Apis_proteome_figure.csv") %>%
   group_by(protein_id) %>%
   summarise(Apis_TPM = sum(Apis_TPM, na.rm = TRUE))
Halictus_genes<- read_csv("../figures/Halictus_proteome_figure.csv")%>%
   group_by(protein_id) %>%
   summarise(Halictu_TPM = sum(Halictus_TPM, na.rm = TRUE))

##Common core


all_species_tr_plot<- Apis_genes  %>%
   full_join(Halictus_genes, by = "protein_id") %>%
   full_join( Xylocopa_genes, by = "protein_id") 


all_spp<- all_species_tr_plot %>%
   drop_na() %>%
   set_colnames(c("protein_id", "Apis", "Halictus", "Xylocopa")) %>%
   gather(species,TPM, Apis:Xylocopa, na.rm = T) %>% 
   ggplot(aes(x = protein_id, y = TPM,fill=species)) +
   geom_bar(position="stack", stat="identity") +
   labs(y= "Venom gland protein composition", x = "Proteins") +
   coord_flip()  + 
   theme_bw() +
   scale_fill_manual(values=c("Apis"     = "#078285", 
                              "Halictus" = "#F1962C",
                              "Xylocopa" = "#4A0933"),
                     name = "") +
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         plot.caption = element_text(hjust = 0),
         text = element_text(size=15)) 
all_spp
ggsave(filename = "../201208_figures/common_core_raw_tpm", width = 7, height = 6, dpi =300, units = "in", device='svg', limitsize = FALSE)
   
all_species_tr_plot %>%
   drop_na() %>%
   set_colnames(c("protein_id", "Apis", "Halictus", "Xylocopa")) %>%
   gather(species,TPM, Apis:Xylocopa, na.rm = T) %>% 
   ggplot(aes(x = protein_id, y = TPM,fill=species)) +
   geom_bar(position = "fill",stat = "identity") +
   labs(y= "Venom gland TPM", x = "Proteins") +
   coord_flip()  + 
   theme_bw() +
   scale_fill_manual(values=c("Apis"     = "#078285", 
                              "Halictus" = "#F1962C",
                              "Xylocopa" = "#4A0933"),
                     name = "") +
   scale_y_continuous(labels = scales::percent_format())  +
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         plot.caption = element_text(hjust = 0),
         text = element_text(size=15)) 

ggsave(filename = "../201208_figures/common_core_perc", width = 7, height = 6, dpi =300, units = "in", device='svg', limitsize = FALSE)

all_species_tr_plot %>%
   set_colnames(c("protein_id", "Apis", "Halictus", "Xylocopa")) %>%
   select("protein_id", "Xylocopa") %>%
   drop_na() %>%
   ggplot(aes(x = protein_id, y = Xylocopa)) +
   geom_col(color = "#4A0933", fill = "#4A0933") +
   labs(y= "Xylocopa venom gland TPM", x = "Proteins") +
   coord_flip()  + 
   theme_bw() +
   theme(legend.position = "none", 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         plot.caption = element_text(hjust = 0),
         text = element_text(size=15)) 

ggsave(filename = "../201208_figures/Xylocopa_tpm", width = 7, height = 6, dpi =300, units = "in", device='svg', limitsize = FALSE)

all_species_tr_plot %>%
   set_colnames(c("protein_id", "Apis", "Halictus", "Xylocopa")) %>%
   select("protein_id", "Apis") %>%
   drop_na() %>%
   ggplot(aes(x = protein_id, y = Apis)) +
   geom_col(color ="#078285", fill = "#078285") +
   labs(y= "Apis venom gland TPM", x = "Proteins") +
   coord_flip()  + 
   theme_bw() +
   theme(legend.position = "none", 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         plot.caption = element_text(hjust = 0),
         text = element_text(size=15)) 

ggsave(filename = "../201208_figures/Apis_tpm", width = 7, height = 6, dpi =300, units = "in", device='svg', limitsize = FALSE)
all_species_tr_plot %>%
   set_colnames(c("protein_id", "Apis", "Halictus", "Xylocopa")) %>%
   select("protein_id", "Halictus") %>%
   drop_na() %>%
   ggplot(aes(x = protein_id, y = Halictus)) +
   geom_col(color ="#F1962C", fill = "#F1962C") +
   labs(y= "Halictus venom gland TPM", x = "Proteins") +
   coord_flip()  + 
   theme_bw() +
   scale_fill_manual(values=c("Apis"     = "#078285", 
                              "Halictus" = "#F1962C",
                              "Xylocopa" = "#4A0933"),
                     name = "") +
   theme(legend.position = "none", 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         plot.caption = element_text(hjust = 0),
         text = element_text(size=15)) 

ggsave(filename = "../201208_figures/Halictus_tpm", width = 7, height = 6, dpi =300, units = "in", device='svg', limitsize = FALSE)

library(laviz)
col_colors <- c("#7e6ba4",
                "#6daf4b",
                "#9a47be",
                "#c99048",
                "#5ea9a1",
                "#be4b57",
                "#5b572d")

all_species_tr_plot %>%
   drop_na() %>%
   set_colnames(c("Proteins", "Apis", "Halictus", "Xylocopa")) %>%
   gather(species,TPM, Apis:Xylocopa, na.rm = T) %>% 
   ggplot(aes(x = species, y = TPM,fill= Proteins)) +
   geom_bar(position = "fill",stat = "identity") +
   labs(y= "Proportion of components in venom gland", x = "Species") +
   theme_bw() +
   #scale_fill_manual(values = col_colors)   + 
   scale_fill_viridis(discrete = TRUE, option = "D") +
   #scale_fill_viridis_d(begin = 0, end = 1) +
   #scale_fill_brewer(palette = "Set1", direction = -2) +
   scale_y_continuous(labels = scales::percent_format())  +
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         plot.caption = element_text(hjust = 0),
         text = element_text(size=15)) + 
   coord_flip()

ggsave(filename = "../201208_figures/common_core_perc_spe", width = 7, height = 6, dpi =300, units = "in", device='svg', limitsize = FALSE)

nb.cols <- 22
mycolors <- colorRampPalette(brewer.pal(8, "Accent"))(nb.cols)

all_species_tr_plot %>%
  
   set_colnames(c("Proteins", "Apis", "Halictus", "Xylocopa")) %>%
   gather(species,TPM, Apis:Xylocopa, na.rm = T) %>% 
   ggplot(aes(x = species, y = TPM,fill= Proteins)) +
   geom_bar(position = "fill",stat = "identity") +
   labs(y= "Proportion of components in venom gland ", x = "Species") +
   theme_bw() +
   scale_fill_manual(values = mycolors, name = "") +
   #scale_fill_manual(values = col_colors)   + 
   #scale_fill_viridis(discrete = TRUE, option = "D",name = "") +
   #scale_fill_viridis_d(begin = 0, end = 1) +
   #scale_fill_brewer(palette = "Set1", direction = -2) +
   scale_y_continuous(labels = scales::percent_format())  +
   theme(panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         plot.caption = element_text(hjust = 0),
         text = element_text(size=12),
         legend.position="bottom",
         legend.justification="left",
        legend.text = element_text(color = "black", size = 9),
        plot.margin = unit(c(1, 1, 1, 1), "cm")) +
   coord_flip() +
   guides(colour = guide_legend(override.aes = list(size=2)))
   

ggsave(filename = "../201208_figures/perc_spe_all", width = 11, height = 6, dpi =300, units = "in", device='svg', limitsize = FALSE)

#all_species_tr_plot$percent<-all_species_tr_plot$Xylocopa_TPM/ sum(all_species_tr_plot$Xylocopa_TPM)*100 


a<- all_species_tr_plot %>%
   set_colnames(c("protein_id", "Apis", "Halictus", "Xylocopa")) %>%
   select("protein_id", "Xylocopa") %>%
   drop_na() %>%
   mutate(Percentage=round(Xylocopa/sum(Xylocopa)*100,4)) %>%
   ggplot(aes(x = protein_id, y = Percentage)) +
   geom_col(color = "#4A0933", fill = "#4A0933") +
   labs(y= "", x = "") +
   scale_y_continuous(limits = c(0,100)) +
   coord_flip()  + 
   theme_bw() +
   theme(legend.position = "none", 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         plot.caption = element_text(hjust = 0),
         text = element_text(size=15)) 


b<-all_species_tr_plot %>%
   set_colnames(c("protein_id", "Apis", "Halictus", "Xylocopa")) %>%
   select("protein_id", "Apis") %>%
   drop_na() %>%
   mutate(Percentage=round(Apis/sum(Apis)*100,4)) %>%
   ggplot(aes(x = protein_id, y = Percentage))  +
   geom_col(color ="#078285", fill = "#078285") +
   labs(y= "", x = "") +
   scale_y_continuous(limits = c(0,100)) +
   coord_flip()  + 
   theme_bw() +
   theme(legend.position = "none", 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         plot.caption = element_text(hjust = 0),
         text = element_text(size=15)) 

c<-all_species_tr_plot %>%
   set_colnames(c("protein_id", "Apis", "Halictus", "Xylocopa")) %>%
   select("protein_id", "Halictus") %>%
   drop_na() %>%
   mutate(Percentage=round(Halictus/sum(Halictus)*100,4)) %>%
   ggplot(aes(x = protein_id, y = Percentage))  +
   geom_col(color ="#F1962C", fill = "#F1962C") +
   labs(y= "Proportion of venom components (%)", x = "") +
   scale_y_continuous(limits = c(0,100)) +
  
    coord_flip()  + 
   theme_bw() +
   theme(legend.position = "none", 
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         plot.caption = element_text(hjust = 0),
         text = element_text(size=15)) 

legend<- get_legend(all_spp)

ggarrange(ggarrange(a,b,c, ncol =1, heights = c(0.85, 1, 0.8)),
          legend, ncol=2,  widths =  c(1, 0.15))
ggsave(filename = "../201208_figures/percentage_all_spp", width = 10, height = 8, dpi =300, units = "in", device='svg', limitsize = FALSE)



# Chord diagram

#generate and adjacency matrix list
library(circlize)


all_species_proteins <- all_species_tr_plot %>%
   set_colnames(c("Proteins", "Apis", "Halictus", "Xylocopa")) %>%
   data.frame() %>%
   gather(key, value, Apis:Xylocopa) %>%
   drop_na()  %>%
   group_by(key) %>%
   mutate(percent = value / sum(value) * 100) %>%
   select(Proteins, key, value = percent) %>%
   arrange(desc(value))


circos.clear()
# Selecting colour palette
mycolor <- viridis(length(unique(all_species_proteins$Proteins)),alpha = 1, begin = 0, end = 1, option = "D")
mycolor <- sample(mycolor, length(unique(all_species_proteins$Proteins)))


grid.col = NULL # just create the variable

grid.col[all_species_proteins$Proteins] = mycolor
tiff('../201208_figures/Chord_diagram1.png', units="in", width=8, height=8, res=400, compression = 'lzw')

circos.clear() 
circos.par(start.degree = 100, points.overflow.warning = T, gap.degree = 4)
par(cex = 2.5, mar = c(0, 0, 0, 0))
g<-{chordDiagram(all_species_proteins,grid.col = grid.col,row.col = 1:3, annotationTrack = "grid", annotationTrackHeight = c(0.08, 0.08), big.gap = 20, directional = TRUE, diffHeight = 0.06, transparency = 0.3, preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(all_species_proteins))))))
   # we go back to the first track and customize sector labels
   circos.track(track.index = 1, panel.fun = function(x, y) {
      circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", cex = 0.3, niceFacing = TRUE, adj = c(0,0.4)) }, bg.border = NA, track.height = 10)
}
dev.off()


all_species_filtered <- all_species_tr_plot %>%
   set_colnames(c("protein_id", "Apis", "Halictus", "Xylocopa")) %>%
   gather(species,TPM, Apis:Xylocopa, na.rm = T) %>%
   drop_na()

venn.diagram(
   x = list(    
      all_species_filtered %>% filter(species=="Apis" ) %>% pull(protein_id), 
      all_species_filtered %>% filter(species=="Halictus" ) %>% pull(protein_id), 
      all_species_filtered %>% filter(species=="Xylocopa" ) %>% pull(protein_id)), 
   category.names = c("Apis (17)" , "Halictus (11)", "Xylocopa (14)"),
   filename = '../201208_figures/venom_venn_filt.png',
   output = TRUE ,
   imagetype="png" ,
   height = 500 , 
   width = 500 , 
   resolution = 400,
   compression = "lzw",
   lwd = 1,
   col=c( '#078285', '#F1962C', '#4A0933' ),
   fill = c(alpha('#078285',0.3), alpha('#F1962C',0.3), alpha('#4A0933',0.3)),
   cex = 0.4,
   fontfamily = "sans",
   cat.cex = 0.4,
   cat.pos = c(9, -11, -1),
   cat.dist = c(0.4, 0.33, 0.3),
   cat.default.pos = "text",
   cat.fontfamily = "sans",
   cat.col = c('#078285', '#F1962C', '#4A0933'))




