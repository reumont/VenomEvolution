library(tidyr)
library(dplyr)
library(stringr)
library(readr)



setwd("/Users/marianavelasque/Dropbox/Pesquisa/Post-doc/Ivan_research/200722_venom_transcriptome/201215_ortholog_synteny/201230_MCSscanX")


blast_results<- read_table2("201230_apamin_flank/blast_results.blast", col_names = FALSE)

gff_file<- read_table2("201231_gff_files/merged_gff_files.gff", col_names = FALSE) %>%
   mutate(X2 = str_replace(X2, "exon-", "")) %>%
   mutate(X2 = str_replace(X2, "rna-", "")) %>%
   mutate(X2 = str_replace(X2, "cds-", "")) %>%
   separate( X2, into = c("gene_desc", "gene", "duplicate"), sep = "_") %>%
   select("gene_desc", "gene", "duplicate")

gff_file$gene <- paste(gff_file$gene_desc, "_", gff_file$gene)

gff_file$gene<- sub(" ", "", gff_file$gene)
gff_file$gene<- sub(" ", "", gff_file$gene)
gff_file$gene_desc = NULL

blast_results<- blast_results %>%
   left_join(gff_file, by = c("X1" = "gene"))

blast_results$X1 <- paste(blast_results$X1, "_", blast_results$duplicate)

blast_results$X1<- sub(" _ NA", "", blast_results$X1)
blast_results$X1<- sub(" ", "", blast_results$X1)
blast_results$X1<- sub(" ", "", blast_results$X1)


blast_results$duplicate = NULL


blast_results<- blast_results %>%
   left_join(gff_file, by = c("X2" = "gene"))

blast_results$X2 <- paste(blast_results$X2, "_", blast_results$duplicate)

blast_results$X2<- sub(" _ NA", "", blast_results$X2)
blast_results$X2<- sub(" ", "", blast_results$X2)
blast_results$X2<- sub(" ", "", blast_results$X2)


blast_results$duplicate = NULL

write.table( blast_results, "201230_apamin_flank/blast_results_duplicated.blast",  row.names = FALSE, col.names=FALSE, sep="\t", quote = FALSE)

