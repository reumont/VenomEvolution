library(stringr)
library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(tximport)
library(lamisc)
library(tibble)
library(gridExtra)
options(scipen = 999)


#loadLoaded 

filtered_protein_lists <- read_csv("../data/filtered_protein_lists.csv")


Xylocopa_ToxPr<- read_delim("03_Transcriptome_ProtPredictions_Annotations/03_Xyvi_Vg.ORP.fasta.pep_ToxPr.blastp.outfmt6", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE) %>%
   dplyr::select(`Transdecoder protein IDs` = "X1", 
                 "Blast-P toxprot sseqid" = "X2", "Blast-P toxprot e-value" = "X11",
                 "Blast-Ptoxprot bitscore" = X12) %>%
   mutate(`Blast-P toxprot sseqid` = str_replace(`Blast-P toxprot sseqid`, " O.*", "")) %>%
   distinct(`Transdecoder protein IDs`, `Blast-P toxprot sseqid`, .keep_all = TRUE)

Xylocopa_Apoc<- read_delim("03_Transcriptome_ProtPredictions_Annotations/03_Xyvi_Vg.ORP.fasta.pep_Apocrita.blastp.outfmt6", 
                           "\t", escape_double = FALSE, col_names = FALSE, 
                           trim_ws = TRUE)%>%
   dplyr::select(`Transdecoder protein IDs` = "X1", "Blast-P UniProt Apocrita sseqid" = "X2", "Blast-P UniProt Apocrita e-value" = "X11", "Blast-P UniProt Apocrita bitscore" = X12) %>%
   mutate(`Blast-P UniProt Apocrita sseqid` = str_replace(`Blast-P UniProt Apocrita sseqid`, " O.*", "")) %>%
   distinct(`Transdecoder protein IDs`, `Blast-P UniProt Apocrita sseqid`, .keep_all = TRUE)

Xylocopa_UniPr<- read_delim("03_Transcriptome_ProtPredictions_Annotations/03_Xyvi_Vg.ORP.fasta.pep_UniPr.blastp.outfmt6", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE) %>%
   dplyr::select(`Transdecoder protein IDs` = "X1", "Blast-P UniProt All sseqid" = "X2", "Blast-P UniProt All e-value" = "X11", "Blast-P UniProt All bitscore" = X12) %>%
   mutate(`Blast-P UniProt All sseqid` = str_replace(`Blast-P UniProt All sseqid`, " O.*", "")) %>%
   distinct(`Transdecoder protein IDs`, `Blast-P UniProt All sseqid`, .keep_all = TRUE)

Xylocopa_MS_protein <- read_csv("05_Proteomics_MSOutput/Xylocopa_proteins.csv")

Xylocopa_cleaned <- read_csv("06_Proteotranscriptomics_Results/02_IvanResults_checked/Xylocopa_proteome_cleaned_bmvr.csv") 


## Import transcripts 


Xylocopa_path<-file.path("../200722_venom_transcriptome/Bee_venom_evo_dataset/Xylocopa",
                         "Kallisto_output", 'abundance.h5')

#Rename reads 
names(Xylocopa_path) <- "Expression value"

# import gene expression from Kallisto

Xylocopa_list<- tximport(Xylocopa_path, type = "kallisto",
                         txOut = TRUE,  #tx2gene = txt2gene, 
                         countsFromAbundance = "scaledTPM")

# select abundance table (TPM results needs to be summarized to transcript level)
Xylocopa_expre<- Xylocopa_list$abundance %>%
   as.data.frame() %>%
   rownames_to_column(var=`Transcript ID`)

#Import peptide table

Xylocopa_pept<- read_delim("../03_Xyvi_Vg.ORP.fasta.transdecoder.pep_singleline", 
                           "\t", escape_double = FALSE, trim_ws = TRUE) %>%
   select(`Transdecoder protein IDs` = "X1") %>%
   mutate(`Transdecoder protein IDs` = str_replace(`Transdecoder protein IDs`, " .*", "")) %>%
   mutate(`Transcript ID` = str_replace(`Transdecoder protein IDs`, ".p.*", ""))  %>%
   distinct(`Transdecoder protein IDs`, .keep_all = TRUE)


Xylocopa_proteins<- Xylocopa_pept %>%
   full_join(Xylocopa_ToxPr,by ="Transdecoder protein IDs") %>%
   full_join(Xylocopa_UniPr,by ="Transdecoder protein IDs") %>%
   full_join(Xylocopa_Apoc,by ="Transdecoder protein IDs") %>%
   mutate(`Transcript ID`= str_replace(`Transdecoder protein IDs`, ".p.*", "")) %>%
   #left_join(Xylocopa_trand, by = `Transdecoder protein IDs`)  %>%
   left_join(Xylocopa_expre, by  = "Transcript ID")  %>%
   arrange(rowSums(is.na(.))) %>%
   select(`Transdecoder protein IDs`, `Transcript ID`, 'Expression value', 
          `Blast-P toxprot sseqid`, 'Blast-P toxprot e-value', 
          'Blast-Ptoxprot bitscore', `Blast-P UniProt Apocrita sseqid`, 
          'Blast-P UniProt Apocrita e-value', 'Blast-P UniProt Apocrita bitscore',
          `Blast-P UniProt All sseqid`, 'Blast-P UniProt All e-value', 
          'Blast-P UniProt All bitscore') %>%
   drop_na(`Transdecoder protein IDs`) 


Xylocopa_clean <- Xylocopa_proteins %>%
   filter(`Transcript ID` %in% Xylocopa_cleaned$gene_id) %>%
   distinct(`Transdecoder protein IDs`, .keep_all = TRUE) %>%
   droplevels


Xylocopa_protein_ms <- Xylocopa_proteins %>%
   filter(`Transdecoder protein IDs` %in% Xylocopa_MS_protein$Accession)  %>%
   droplevels

write.csv( Xylocopa_proteins , "../results/Xylocopa_protein_expression_complete.csv",  
           row.names = FALSE, quote = FALSE)
write.csv( Xylocopa_protein_ms , "../results/Xylocopa_protein_expression_MS.csv",  
           row.names = FALSE, quote = FALSE)
write.csv( Xylocopa_clean , "../results/Xylocopa_protein_expression_cleaned_bmvr.csv",  
           row.names = FALSE, quote = FALSE)



Halictus_ToxPr<- read_delim("03_Transcriptome_ProtPredictions_Annotations/Ha_Vg6_5.ORP.fasta.pep_ToxPr.blastp.outfmt6", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE) %>%
   dplyr::select(`Transdecoder protein IDs` = "X1", 
                 "Blast-P toxprot sseqid" = "X2", "Blast-P toxprot e-value" = "X11",
                 "Blast-Ptoxprot bitscore" = X12) %>%
   mutate(`Blast-P toxprot sseqid` = str_replace(`Blast-P toxprot sseqid`, " O.*", "")) %>%
   distinct(`Transdecoder protein IDs`, `Blast-P toxprot sseqid`, .keep_all = TRUE)


Halictus_Apoc<- read_delim("03_Transcriptome_ProtPredictions_Annotations/Ha_Vg6_5.ORP.fasta.pep_Apocrita.blastp.outfmt6", 
                           "\t", escape_double = FALSE, col_names = FALSE, 
                           trim_ws = TRUE)%>%
   dplyr::select(`Transdecoder protein IDs` = "X1", "Blast-P UniProt Apocrita sseqid" = "X2", "Blast-P UniProt Apocrita e-value" = "X11", "Blast-P UniProt Apocrita bitscore" = X12) %>%
   mutate(`Blast-P UniProt Apocrita sseqid` = str_replace(`Blast-P UniProt Apocrita sseqid`, " O.*", "")) %>%
   distinct(`Transdecoder protein IDs`, `Blast-P UniProt Apocrita sseqid`, .keep_all = TRUE)


Halictus_UniPr<- read_delim("03_Transcriptome_ProtPredictions_Annotations/Ha_Vg6_5.ORP.fasta.pep_UniPr.blastp.outfmt6", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE) %>%
   dplyr::select(`Transdecoder protein IDs` = "X1", "Blast-P UniProt All sseqid" = "X2", "Blast-P UniProt All e-value" = "X11", "Blast-P UniProt All bitscore" = X12) %>%
   mutate(`Blast-P UniProt All sseqid` = str_replace(`Blast-P UniProt All sseqid`, " O.*", "")) %>%
   distinct(`Transdecoder protein IDs`, `Blast-P UniProt All sseqid`, .keep_all = TRUE)


## Import MS data 
Halictus_cleaned <- read_csv("06_Proteotranscriptomics_Results/02_IvanResults_checked/Halictus_proteome_cleaned_bmvr.csv") 

## Import MS data 
Halictus_MS_protein <- read_csv("05_Proteomics_MSOutput/Halictus_MascotResults_CutOff_MascotValue_40.csv") %>%
   drop_na(Accession) %>%
   mutate('Accession'= str_replace(`Accession`, "  ORF .*", "")) %>% ##removing MS information score
   mutate(Accession=strsplit(Accession, " ")) %>% 
   unnest(Accession) %>%
   mutate(Accession=strsplit(Accession, "~~")) %>% 
   unnest(Accession) %>%
   mutate('Accession'= str_replace(`Accession`, "GENE.", "")) %>% 
   distinct(Accession, .keep_all = TRUE) %>% 
   #filter(`Score Mascot: Mascot`  >= 40 & 
   #          `# Peptides` > 2) %>%
   mutate(`Transcript ID`= str_replace(`Accession`, ".p.*", "")) 


## Import transcripts 


Halictus_path<-file.path("../200722_venom_transcriptome/Bee_venom_evo_dataset/Halictus",
                         "Kallisto_output", 'abundance.h5')

#Rename reads 
names(Halictus_path) <- "Expression value"

# import gene expression with Kallisto

Halictus_list<- tximport(Halictus_path, type = "kallisto",
                         txOut = TRUE,  #tx2gene = txt2gene, 
                         countsFromAbundance = "scaledTPM")

# select abundance table (TPM results needs to be summarized to transcript level)
Halictus_expre<- Halictus_list$abundance %>%
   as.data.frame() %>%
   rownames_to_column(var=`Transcript ID`)

## import peptide data

Halictus_pept<- read_delim("../Ha_Vg6_5.ORP.fasta.transdecoder.pep_singleline", 
                           "\t", escape_double = FALSE, trim_ws = TRUE) %>%
   select(`Transdecoder protein IDs` = "X1") %>%
   mutate(`Transdecoder protein IDs` = str_replace(`Transdecoder protein IDs`, " .*", "")) %>%
   mutate(`Transcript ID` = str_replace(`Transdecoder protein IDs`, ".p.*", ""))  %>%
   distinct(`Transdecoder protein IDs`, .keep_all = TRUE)  


Halictus_proteins<- Halictus_pept %>%
   full_join(Halictus_ToxPr,by ="Transdecoder protein IDs") %>%
   full_join(Halictus_UniPr,by ="Transdecoder protein IDs") %>%
   full_join(Halictus_Apoc,by ="Transdecoder protein IDs") %>%
   mutate(`Transcript ID`= str_replace(`Transdecoder protein IDs`, ".p.*", "")) %>%
   #left_join(Halictus_trand, by = `Transdecoder protein IDs`)  %>%
   left_join(Halictus_expre, by  = "Transcript ID")  %>%
   arrange(rowSums(is.na(.))) %>%
   select(`Transdecoder protein IDs`, `Transcript ID`, 'Expression value', 
          `Blast-P toxprot sseqid`, 'Blast-P toxprot e-value', 
          'Blast-Ptoxprot bitscore', `Blast-P UniProt Apocrita sseqid`, 
          'Blast-P UniProt Apocrita e-value', 'Blast-P UniProt Apocrita bitscore',
          `Blast-P UniProt All sseqid`, 'Blast-P UniProt All e-value', 
          'Blast-P UniProt All bitscore') %>%
   drop_na(`Transdecoder protein IDs`) 


Halictus_clean <- Halictus_proteins %>%
   filter(`Transcript ID` %in% Halictus_cleaned$gene_id) %>%
   distinct(`Transdecoder protein IDs`, .keep_all = TRUE) %>%
   droplevels


Halictus_protein_ms <- Halictus_proteins %>%
   filter(`Transcript ID` %in% Halictus_MS_protein$`Transcript ID`)  %>%
   droplevels

write.csv( Halictus_proteins , "../results/Halictus_protein_expression_complete.csv",  
           row.names = FALSE, quote = FALSE)
write.csv( Halictus_protein_ms , "../results/Halictus_protein_expression_MS.csv",  
           row.names = FALSE, quote = FALSE)
write.csv( Halictus_clean , "../results/Halictus_protein_expression_cleaned_bmvr.csv",  
           row.names = FALSE, quote = FALSE)



##Apis

Apis_ToxPr<- read_delim("03_Transcriptome_ProtPredictions_Annotations/Apis.ORP.fasta.pep_ToxPr.blastp.outfmt6", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE) %>%
   dplyr::select(`Transdecoder protein IDs` = "X1", 
                 "Blast-P toxprot sseqid" = "X2", "Blast-P toxprot e-value" = "X11",
                 "Blast-Ptoxprot bitscore" = X12) %>%
   mutate(`Blast-P toxprot sseqid` = str_replace(`Blast-P toxprot sseqid`, " O.*", "")) %>%
   distinct(`Transdecoder protein IDs`, `Blast-P toxprot sseqid`, .keep_all = TRUE)


Apis_Apoc<- read_delim("03_Transcriptome_ProtPredictions_Annotations/Apis.ORP.fasta.pep_Apocrita.blastp.outfmt6", 
                       "\t", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE)%>%
   dplyr::select(`Transdecoder protein IDs` = "X1", "Blast-P UniProt Apocrita sseqid" = "X2", "Blast-P UniProt Apocrita e-value" = "X11", "Blast-P UniProt Apocrita bitscore" = X12) %>%
   mutate(`Blast-P UniProt Apocrita sseqid` = str_replace(`Blast-P UniProt Apocrita sseqid`, " O.*", "")) %>%
   distinct(`Transdecoder protein IDs`, `Blast-P UniProt Apocrita sseqid`, .keep_all = TRUE)


Apis_UniPr<- read_delim("03_Transcriptome_ProtPredictions_Annotations/Apis.ORP.fasta.pep_UniPr.blastp.outfmt6", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE) %>%
   dplyr::select(`Transdecoder protein IDs` = "X1", "Blast-P UniProt All sseqid" = "X2", "Blast-P UniProt All e-value" = "X11", "Blast-P UniProt All bitscore" = X12) %>%
   mutate(`Blast-P UniProt All sseqid` = str_replace(`Blast-P UniProt All sseqid`, " OX.*", "")) %>%
   mutate(`Blast-P UniProt All sseqid` = str_replace(`Blast-P UniProt All sseqid`, " OS=.*", "")) %>%
   distinct(`Transdecoder protein IDs`, `Blast-P UniProt All sseqid`, .keep_all = TRUE)


## Import transcripts 

Apis_path<-file.path("../200722_venom_transcriptome/Bee_venom_evo_dataset/Apis",
                     "Kallisto_output", 'abundance.h5')

#Rename reads 
names(Apis_path) <- "Expression value"

# import gene expression with Kallisto

Apis_list<- tximport(Apis_path, type = "kallisto",
                     txOut = TRUE,  #tx2gene = txt2gene, 
                     countsFromAbundance = "scaledTPM")

# select abundance table (TPM results needs to be summarized to transcript level)
Apis_expre<- Apis_list$abundance %>%
   as.data.frame() %>%
   rownames_to_column(var=`Transcript ID`)

#import MS data

Apis_MS_protein <- read_csv("05_Proteomics_MSOutput/Apis_proteins.csv")

# import cleaned data 

Apis_cleaned <- read_csv("06_Proteotranscriptomics_Results/02_IvanResults_checked/Apis_proteome_cleaned_bmvr.csv")


## import peptide data

Apis_pept<- read_delim("../Apis.ORP.fasta.transdecoder.pep_singleline", 
                       "\t", escape_double = FALSE, trim_ws = TRUE) %>%
   select(`Transdecoder protein IDs` = "X1") %>%
   mutate(`Transdecoder protein IDs` = str_replace(`Transdecoder protein IDs`, " .*", "")) %>%
   mutate(`Transcript ID` = str_replace(`Transdecoder protein IDs`, ".p.*", "")) %>%
   distinct(`Transdecoder protein IDs`, .keep_all = TRUE) 


Apis_proteins<- Apis_pept %>%
   full_join(Apis_ToxPr,by ="Transdecoder protein IDs") %>%
   full_join(Apis_UniPr,by ="Transdecoder protein IDs") %>%
   full_join(Apis_Apoc,by ="Transdecoder protein IDs") %>%
   mutate(`Transcript ID`= str_replace(`Transdecoder protein IDs`, ".p.*", "")) %>%
   #left_join(Apis_trand, by = `Transdecoder protein IDs`)  %>%
   left_join(Apis_expre, by  = "Transcript ID")  %>%
   arrange(rowSums(is.na(.))) %>%
   select(`Transdecoder protein IDs`, `Transcript ID`, 'Expression value', 
          `Blast-P toxprot sseqid`, 'Blast-P toxprot e-value', 
          'Blast-Ptoxprot bitscore', `Blast-P UniProt Apocrita sseqid`, 
          'Blast-P UniProt Apocrita e-value', 'Blast-P UniProt Apocrita bitscore',
          `Blast-P UniProt All sseqid`, 'Blast-P UniProt All e-value', 
          'Blast-P UniProt All bitscore') %>%
   drop_na(`Transdecoder protein IDs`) 


Apis_clean <- Apis_proteins %>%
   filter(`Transcript ID` %in% Apis_cleaned$gene_id) %>%
   distinct(`Transdecoder protein IDs`, .keep_all = TRUE) %>%
   droplevels


Apis_protein_ms <- Apis_proteins %>%
   filter(`Transdecoder protein IDs` %in% Apis_MS_protein$Accession)  

write.csv( Apis_proteins , "../results/Apis_protein_expression_complete.csv",  
           row.names = FALSE, quote = FALSE)
write.csv( Apis_protein_ms , "../results/Apis_protein_expression_MS.csv",  
           row.names = FALSE, quote = FALSE)

write.csv( Apis_clean , "../results/Apis_protein_expression_cleaned_bmvr.csv",  
           row.names = FALSE, quote = FALSE)

