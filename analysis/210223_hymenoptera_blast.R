setwd("~/Dropbox/Pesquisa/Post-doc/Ivan_research/200722_venom_transcriptome/scripts")
xyoc_UniPr<- read_delim("../Bee_venom_evo_dataset/Xylocopa/03_Xyvi_Vg.ORP.fasta.pep_UniPr.blastp.outfmt6", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE) %>%
   dplyr::rename("gene_id" = "X1", "protein_id" = "X2" )
xyoc_ToxPr<- read_delim("../Bee_venom_evo_dataset/Xylocopa/03_Xyvi_Vg.ORP.fasta.pep_ToxPr.blastp.outfmt6", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE)%>%
   dplyr::rename("gene_id" = "X1", "protein_id" = "X2" )
xyoc_Apoc<- read_delim("../Bee_venom_evo_dataset/Xylocopa/03_Xyvi_Vg.ORP.fasta.pep_Apocrita.blastp.outfmt6", 
                       "\t", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE)%>%
   dplyr::rename("gene_id" = "X1", "protein_id" = "X2" )


clean_blast_results <- function(blast_table){
   blast_table<- blast_table %>%
      separate(protein_id, c("protein_id", "species"), "OS=") 
   blast_table$species <-  sub("OX=.*", " ", blast_table$species)
   blast_table$gene_id <-  sub(".p.*", "", blast_table$gene_id)
   blast_table$protein_id<- trimws(blast_table$protein_id) ##remove space after last string
   blast_table$species<- trimws(blast_table$species)
   blast_table$gene_id<- trimws(blast_table$gene_id)
   blast_table <- as_tibble(blast_table)
}

## negative function
'%ni%' <- Negate('%in%')  

xyoc_UniPr<- clean_blast_results(xyoc_UniPr)
xyoc_ToxPr<- clean_blast_results(xyoc_ToxPr)
xyoc_Apoc<- clean_blast_results(xyoc_Apoc)


# Remove genes present in Apocrita from ToxProt and Uniprot.

xyoc_ToxPr<- xyoc_ToxPr %>%
   filter(gene_id %ni% xyoc_Apoc$gene_id)

xyoc_UniPr<- xyoc_UniPr %>%
   filter(gene_id %ni% xyoc_Apoc$gene_id)

# Remove genes present in ToxProt from Uniprot.
xyoc_UniPr<- xyoc_UniPr %>%
   filter(gene_id %ni% xyoc_ToxPr$gene_id)

## Create the protein data base by merging all data sets

xyoc_proteins<- bind_rows(xyoc_Apoc, xyoc_ToxPr, xyoc_UniPr)

write.table( xyoc_proteins , "../Bee_venom_evo_dataset/Xylocopa_blast_combined",  row.names = FALSE, col.names=FALSE, sep="\t", quote = FALSE)


Halictus_UniPr<- read_delim("../Bee_venom_evo_dataset/Halictus/Ha_Vg6_5.ORP.fasta.transdecoder.pep", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE) %>%
   dplyr::rename("gene_id" = "X1", "protein_id" = "X2" )
Halictus_ToxPr<- read_delim("../Bee_venom_evo_dataset/Halictus/Ha_Vg6_5.ORP.fasta.pep_ToxPr.blastp.outfmt6", 
                            "\t", escape_double = FALSE, col_names = FALSE, 
                            trim_ws = TRUE)%>%
   dplyr::rename("gene_id" = "X1", "protein_id" = "X2" )
Halictus_Apoc<- read_delim("../Bee_venom_evo_dataset/Halictus/Ha_Vg6_5.ORP.fasta.pep_Apocrita.blastp.outfmt6", 
                           "\t", escape_double = FALSE, col_names = FALSE, 
                           trim_ws = TRUE)%>%
   dplyr::rename("gene_id" = "X1", "protein_id" = "X2" )


## negative function
'%ni%' <- Negate('%in%')  

Halictus_UniPr<- clean_blast_results(Halictus_UniPr)
Halictus_ToxPr<- clean_blast_results(Halictus_ToxPr)
Halictus_Apoc<- clean_blast_results(Halictus_Apoc)

Halictus_Apoc_unc<- Halictus_Apoc %>%
   filter(str_detect(protein_id, pattern="Uncharacterized protein")) ## no Uncharacterized protein
Halictus_ToxPr_unc<- Halictus_ToxPr %>%
   filter(str_detect(protein_id, pattern="Uncharacterized protein"))# # no Uncharacterized protein

# Remove genes present in Apocrita from ToxProt and Uniprot.

Halictus_ToxPr<- Halictus_ToxPr %>%
   filter(gene_id %ni% Halictus_Apoc$gene_id)

Halictus_UniPr<- Halictus_UniPr %>%
   filter(gene_id %ni% Halictus_Apoc$gene_id)

# Remove genes present in ToxProt from Uniprot.
Halictus_UniPr<- Halictus_UniPr %>%
   filter(gene_id %ni% Halictus_ToxPr$gene_id)

## Create the protein data base by merging all data sets

Halictus_proteins<- bind_rows(Halictus_Apoc, Halictus_ToxPr, Halictus_UniPr)

write.table( Halictus_proteins , "../Bee_venom_evo_dataset/Halictus_blast_combined",  row.names = FALSE, col.names=FALSE, sep="\t", quote = FALSE)


Apis_UniPr<- read_delim("../Bee_venom_evo_dataset/Apis/Apis.ORP.fasta.pep_UniPr.blastp.outfmt6", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE) %>%
   dplyr::rename("gene_id" = "X1", "protein_id" = "X2" )
Apis_ToxPr<- read_delim("../Bee_venom_evo_dataset/Apis/Apis.ORP.fasta.pep_ToxPr.blastp.outfmt6", 
                        "\t", escape_double = FALSE, col_names = FALSE, 
                        trim_ws = TRUE)%>%
   dplyr::rename("gene_id" = "X1", "protein_id" = "X2" )
Apis_Apoc<- read_delim("../Bee_venom_evo_dataset/Apis/Apis.ORP.fasta.pep_Apocrita.blastp.outfmt6", 
                       "\t", escape_double = FALSE, col_names = FALSE, 
                       trim_ws = TRUE)%>%
   dplyr::rename("gene_id" = "X1", "protein_id" = "X2" )

## negative function
'%ni%' <- Negate('%in%')  

Apis_UniPr<- clean_blast_results(Apis_UniPr)
Apis_ToxPr<- clean_blast_results(Apis_ToxPr)
Apis_Apoc<- clean_blast_results(Apis_Apoc)

Apis_Apoc_unc<- Apis_Apoc %>%
   filter(str_detect(protein_id, pattern="Uncharacterized protein")) ## no Uncharacterized protein
Apis_ToxPr_unc<- Apis_ToxPr %>%
   filter(str_detect(protein_id, pattern="Uncharacterized protein"))# # no Uncharacterized protein

# Remove genes present in Apocrita from ToxProt and Uniprot.

Apis_ToxPr<- Apis_ToxPr %>%
   filter(gene_id %ni% Apis_Apoc$gene_id)

Apis_UniPr<- Apis_UniPr %>%
   filter(gene_id %ni% Apis_Apoc$gene_id)

# Remove genes present in ToxProt from Uniprot.
Apis_UniPr<- Apis_UniPr %>%
   filter(gene_id %ni% Apis_ToxPr$gene_id)

## Create the protein data base by merging all data sets

Apis_proteins<- bind_rows(Apis_Apoc, Apis_ToxPr, Apis_UniPr)

write.table( Apis_proteins , "../Bee_venom_evo_dataset/Apis_blast_combined",  row.names = FALSE, col.names=FALSE, sep="\t", quote = FALSE)
