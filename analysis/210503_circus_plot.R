library(circlize)
library(RColorBrewer)
library(viridis)
library(Cairo)
library(VennDiagram)
library(tidyverse)
library(limma)
library(cowplot)
library(circlize)
library(scales)
library(svglite)
library(venn)
library(ggplot2)
library(gridExtra)
library(forcats)
library(randomcoloR)
library(readr)
library(NOISeq)
library(data.table)
library(ggforce)
library(magrittr)
library(stringi)
library(kableExtra)
library(ggpubr)
library(tximport)
library(stringi)
library(dplyr)
#options(scipen=999)
## negative function
'%ni%' <- Negate('%in%')


##Function to clean up data

clean_transcripts <- function(filter_data) {
  data <-  filter_data %>%
    dplyr::select("Transcript_ID" = "Transcript ID", "tpm_value"= "Expression value"  ,
                  "Blast_toxprot" = "Blast-P toxprot (sseqid)",
                  "Blast_Apocrita" = "Blast-P UniProt Apocrita (sseqid)",
                  "Blast_UniProt" = "Blast-P UniProt All (sseqid)", "Physiological_association"= "Physiological association" , 
                  "Venom_protein_family"= "Venom protein family") %>%
    mutate(Venom_protein_family=replace(Venom_protein_family, Venom_protein_family=="Uncharacterized protein", NA)) %>%
    mutate(Venom_protein_family=replace(Venom_protein_family, 
                                        Venom_protein_family=="Mellitin", "Melittin")) %>%
    mutate(Venom_protein_family=replace(Venom_protein_family, 
                                        Venom_protein_family=="Peroxiredoxin4", "Peroxiredoxin 4")) %>%
    mutate(Venom_protein_family=replace(Venom_protein_family, 
                                        Venom_protein_family=="Phospholipase  A2", "Phospholipase A2")) %>%
    mutate(Venom_protein_family=replace(Venom_protein_family, 
                                        Venom_protein_family=="Venom Serine Protease", "Venom serine protease")) %>%
    mutate(Venom_protein_family=replace(Venom_protein_family, 
                                        Venom_protein_family=="Major royal jelly protein", "Major royal jelly protein")) %>%
    mutate(Venom_protein_family=replace(Venom_protein_family, 
                                        Venom_protein_family=="Major Royal Jelly", "Major royal jelly protein")) %>%
    
    mutate(Venom_protein_family=replace(Venom_protein_family, 
                                        Venom_protein_family=="Dipeptidy peptidase 4", "Dipeptidyl peptidase 4")) %>%
    group_by(Venom_protein_family) %>% 
    summarise(tpm_value = sum(tpm_value)) 
  
}

##load Xylococopa data set
Xylocopa_data <- read_delim("../data/SupplementaryTable3_XylocopaOnlyVenomComponents.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE, col_names = T) %>% 
  clean_transcripts() %>%
  dplyr::select(Venom_protein_family, "Xylocopa_TPM" = "tpm_value")

Apis_data <- read_delim("../data/SupplementaryTableS1_ApisONlyVenomComponents.txt", 
                        "\t", escape_double = FALSE, trim_ws = TRUE, col_names = T) %>% 
  clean_transcripts() %>%
  dplyr::select(Venom_protein_family, Apis_TPM = tpm_value)

Halcitus_data <- read_delim("../data/SupplementaryTableS2_Halcitus_OnlyVenomCOmponents.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE, col_names = T) %>% 
  clean_transcripts() %>%
  dplyr::select(Venom_protein_family, Halcitus_TPM = tpm_value)


## create new data set merging all species tables

hymen_plots <- Apis_data %>%
  full_join(Halcitus_data, by = "Venom_protein_family") %>%
  full_join(Xylocopa_data, by = "Venom_protein_family") %>%
  set_colnames(c("Proteins", "Apis", "Halictus", "Xylocopa")) %>%
  arrange(desc(Halictus))


all_species <- hymen_plots %>%
  set_colnames(c("protein_id", "Apis", "Halictus", "Xylocopa")) %>%
  gather(species,TPM, Apis:Xylocopa, na.rm = T) %>%
  arrange(desc(TPM))

## Group data based on species (turn a long table into a wide one)
all_species_proteins <- hymen_plots %>%
  set_colnames(c("Proteins", "Apis", "Halictus", "Xylocopa")) %>%
  data.frame() %>%
  gather(key, value, Apis:Xylocopa) %>%
  drop_na() %>%
  group_by(key) %>%
  mutate(percent = value / sum(value) * 100) %>%
  dplyr::select(Proteins, key, value = percent) %>%
  arrange(desc(value))

##dplyr::select color scheme for the plot
n= length(hymen_plots$Proteins)

mycolor <- distinctColorPalette(n)

## Inicialise the chord diagram

circos.clear()
grid.col = NULL # just create the variable

grid.col[hymen_plots$Proteins] = mycolor
pdf('../results/Chord_diagram_all.pdf', width=9, height=9,  useDingbats=FALSE)

circos.clear() 
circos.par(start.degree = 170, points.overflow.warning = T, gap.degree = 3)
par(cex = 2.5, mar = c(0, 0, 0, 0))
g<-{chordDiagram(all_species_proteins,grid.col = grid.col,row.col = 1:3, annotationTrack = "grid", annotationTrackHeight = c(0.08, 0.08), big.gap = 2, directional = TRUE, diffHeight = 0.003, transparency = 0.3, preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(all_species_proteins))))))
  # we go back to the first track and customize sector labels
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, facing = "clockwise", cex = 0.38, niceFacing = TRUE, adj = c(0,0.4)) }, bg.border = NA, track.height = 10)
}
dev.off()
