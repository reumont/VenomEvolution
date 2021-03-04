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
setwd("/Users/marianavelasque/Dropbox/Pesquisa/Post-doc/Ivan_research/200722_venom_transcriptome/scripts/")

xyoc_UniPr<- read_delim("../Bee_venom_evo_dataset/Xylocopa/03_Xyvi_Vg.ORP.fasta.pep_UniPr.blastp.outfmt6", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE) %>%
   dplyr::select("transcript_ID" = "X1", "UniPr_protein_ID" = "X2", "UniPr_evalue" = "X11" ) %>%
   mutate(UniPr_protein_ID = str_replace(UniPr_protein_ID, " OX.*", "")) %>%
   mutate(UniPr_protein_ID = str_replace(UniPr_protein_ID, " OS=.*", "")) 

xyoc_ToxPr<- read_delim("../Bee_venom_evo_dataset/Xylocopa/03_Xyvi_Vg.ORP.fasta.pep_ToxPr.blastp.outfmt6", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE)%>%
   dplyr::select("transcript_ID" = "X1", "ToxPr_protein_ID" = "X2", "ToxPr_evalue" = "X11" ) %>%
   mutate(ToxPr_protein_ID = str_replace(ToxPr_protein_ID, " OX.*", "")) %>%
   mutate(ToxPr_protein_ID = str_replace(ToxPr_protein_ID, " OS=.*", "")) 

xyoc_Apoc<- read_delim("../Bee_venom_evo_dataset/Xylocopa/03_Xyvi_Vg.ORP.fasta.pep_Apocrita.blastp.outfmt6", 
                       "\t", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE)%>%
   dplyr::select("transcript_ID" = "X1", "Apoc_protein_ID" = "X2", "Apoc_evalue" = "X11" ) %>%
   mutate(Apoc_protein_ID = str_replace(Apoc_protein_ID, " OX.*", "")) %>%
   mutate(Apoc_protein_ID = str_replace(Apoc_protein_ID, " OS=.*", "")) 

   
## on terminal 

#  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' <03_Xyvi_Vg.ORP.fasta.transdecoder.pep | sed 's/ GENE.*(+)//g' | sed 's/ GENE.*(-)//g' |tr "\n" "\t" |tr ">" "\n"  > 03_Xyvi_Vg.ORP.fasta.transdecoder.pep_singleline


xyoc_trand<- read_delim("../Bee_venom_evo_dataset/Xylocopa/03_Xyvi_Vg.ORP.fasta.transdecoder.pep_singleline", 
                        "\t", escape_double = FALSE, trim_ws = TRUE) %>%
   select("transcript_ID" = "X1", "protein_seq" = "X2")

xyoc_trand <- xyoc_trand %>%
   mutate(transd_ID= str_replace(transcript_ID, ".p.*", "")) 

xyoc_proteins<- xyoc_UniPr %>%
   left_join(xyoc_ToxPr, by ="transcript_ID") %>%
   left_join(xyoc_Apoc, by ="transcript_ID") %>%
   left_join(xyoc_trand, by = "transcript_ID")  %>%
   arrange(rowSums(is.na(.))) %>%
   select("transcript_ID", "transd_ID", "UniPr_protein_ID", "UniPr_evalue", 
          "ToxPr_protein_ID", "ToxPr_evalue", "Apoc_protein_ID",
          "Apoc_evalue", "protein_seq")

write.table( xyoc_proteins , "../Bee_venom_evo_dataset/Xylocopa_blast_combined",  row.names = FALSE, col.names=FALSE, sep="\t", quote = FALSE)


Halictus_ToxPr<- read_delim("../Bee_venom_evo_dataset/Halictus/Ha_Vg6_5.ORP.fasta.pep_ToxPr.blastp.outfmt6", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)%>%
   dplyr::select("transcript_ID" = "X1", "ToxPr_protein_ID" = "X2", "ToxPr_evalue" = "X11" ) %>%
   mutate(ToxPr_protein_ID = str_replace(ToxPr_protein_ID, " OX.*", "")) %>%
   mutate(ToxPr_protein_ID = str_replace(ToxPr_protein_ID, " OS=.*", "")) 

   dplyr::rename("transcript_ID" = "X1", "protein_ID" = "X2", "evalue" = "X11" )
Halictus_Apoc<- read_delim("../Bee_venom_evo_dataset/Halictus/Ha_Vg6_5.ORP.fasta.pep_Apocrita.blastp.outfmt6", 
                           "\t", escape_double = FALSE, col_names = FALSE, 
                           trim_ws = TRUE)%>%
   dplyr::select("transcript_ID" = "X1", "Apoc_protein_ID" = "X2", "Apoc_evalue" = "X11" ) %>%
   mutate(Apoc_protein_ID = str_replace(Apoc_protein_ID, " OX.*", "")) %>%
   mutate(Apoc_protein_ID = str_replace(Apoc_protein_ID, " OS=.*", "")) 

## on terminal 

#   awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' <Ha_Vg6_5.ORP.fasta.transdecoder.pep | sed 's/ GENE.*(+)//g' | sed 's/ GENE.*(-)//g' |tr "\n" "\t" |tr ">" "\n"  > Ha_Vg6_5.ORP.fasta.transdecoder.pep_singleline


Halictus_trand<- read_delim("../Bee_venom_evo_dataset/Halictus/Ha_Vg6_5.ORP.fasta.transdecoder.pep_singleline", 
                            "\t", escape_double = FALSE, trim_ws = TRUE) %>%
   select("transcript_ID" = "X1", "protein_seq" = "X2")

Halictus_trand <- Halictus_trand %>%
   mutate(transd_ID= str_replace(transcript_ID, ".p.*", "")) 

Halictus_proteins<- Halictus_ToxPr %>%
   left_join(Halictus_Apoc, by ="transcript_ID") %>%
   left_join(Halictus_trand, by = "transcript_ID")  %>%
   arrange(rowSums(is.na(.))) %>%
   select("transcript_ID", "transd_ID", "ToxPr_protein_ID", 
          "ToxPr_evalue", "Apoc_protein_ID",
          "Apoc_evalue", "protein_seq")

write.table( Halictus_proteins , "../Bee_venom_evo_dataset/Halictus_blast_combined",  row.names = FALSE, col.names=FALSE, sep="\t", quote = FALSE)


##Apis

Apis_UniPr<- read_delim("../Bee_venom_evo_dataset/Apis/Apis.ORP.fasta.pep_UniPr.blastp.outfmt6", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE) %>%
   dplyr::select("transcript_ID" = "X1", "UniPr_protein_ID" = "X2", "UniPr_evalue" = "X11" ) %>%
   mutate(UniPr_protein_ID = str_replace(UniPr_protein_ID, " OX.*", "")) %>%
   mutate(UniPr_protein_ID = str_replace(UniPr_protein_ID, " OS=.*", "")) 

Apis_ToxPr<- read_delim("../Bee_venom_evo_dataset/Apis/Apis.ORP.fasta.pep_ToxPr.blastp.outfmt6", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE)%>%
   dplyr::select("transcript_ID" = "X1", "ToxPr_protein_ID" = "X2", "ToxPr_evalue" = "X11" ) %>%
   mutate(ToxPr_protein_ID = str_replace(ToxPr_protein_ID, " OX.*", "")) %>%
   mutate(ToxPr_protein_ID = str_replace(ToxPr_protein_ID, " OS=.*", "")) 


Apis_Apoc<- read_delim("../Bee_venom_evo_dataset/Apis/Apis.ORP.fasta.pep_Apocrita.blastp.outfmt6", 
                       "\t", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE)%>%
   dplyr::select("transcript_ID" = "X1", "Apoc_protein_ID" = "X2", "Apoc_evalue" = "X11" ) %>%
   mutate(Apoc_protein_ID = str_replace(Apoc_protein_ID, " OX.*", "")) %>%
   mutate(Apoc_protein_ID = str_replace(Apoc_protein_ID, " OS=.*", "")) 

## on terminal 

      #      awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' <Apis.ORP.fasta.transdecoder.pep | sed 's/ GENE.*(+)//g' | sed 's/ GENE.*(-)//g' |tr "\n" "\t" |tr ">" "\n"  > Apis.ORP.fasta.transdecoder.pep_singleline



Apis_trand<- read_delim("../Bee_venom_evo_dataset/Apis/Apis.ORP.fasta.transdecoder.pep_singleline", 
                        "\t", escape_double = FALSE, trim_ws = TRUE) %>%
   select("transcript_ID" = "X1", "protein_seq" = "X2")

Apis_trand <- Apis_trand %>%
   mutate(transd_ID= str_replace(transcript_ID, ".p.*", "")) 

Apis_proteins<- Apis_UniPr %>%
   left_join(Apis_ToxPr, by ="transcript_ID") %>%
   left_join(Apis_Apoc, by ="transcript_ID") %>%
   left_join(Apis_trand, by = "transcript_ID")  %>%
   arrange(rowSums(is.na(.))) %>%
   select("transcript_ID", "transd_ID", "UniPr_protein_ID", "UniPr_evalue", 
          "ToxPr_protein_ID", "ToxPr_evalue", "Apoc_protein_ID",
          "Apoc_evalue", "protein_seq")


write.table( Apis_proteins , "../Bee_venom_evo_dataset/Apis_blast_combined",  row.names = FALSE, col.names=FALSE, sep="\t", quote = FALSE)
