# Identifying apamin orthologs using synteny



First Obtain flanking genes from the bed files


    gff2bed < ../GCF_003254395.2_Amel_HAv3.1_genomic.gff >  \
    GCF_003254395.2_Amel_HAv3.1_genomic.bed

    conda activate

    cgat bed2bed --method=merge --merge-by-name -I  GCF_003254395.2_Amel_HAv3.1_genomic.bed >   \
    GCF_003254395.2_Amel_HAv3.1_genomic_merged.bed

    grep -e "cds-" GCF_003254395.2_Amel_HAv3.1_genomic_merged.bed |  \
        sed 's/cds-//g' > GCF_003254395.2_Amel_HAv3.1_genomic_merged_flank.bed

    grep -C18 -E 'NP_001011611.2|NP_001011612.1' GCF_003254395.2_Amel_HAv3.1_genomic_merged_flank.bed  \
    > GCF_003254395.2_Amel_HAv3.1_genomic_merged_flanks.bed            ##grep the flanking genes


    awk '{print $4}' GCF_003254395.2_Amel_HAv3.1_genomic_merged_flanks.bed |  \
    cut -d , -f2 | sort | uniq > apis_flanking_genes

Filter the select genes from the fasta file

    seqkit grep -f apis_flanking_genes ../GCF_003254395.2_Amel_HAv3.1_protein.faa -o GCF_003254395.2_Amel_HAv3.1.flank.fasta

    mkdir 210120_blast_results_orthologs
    cd 210120_blast_results_orthologs

    makeblastdb -in ../GCF_003254395.2_Amel_HAv3.1.flank.fasta \
    -input_type fasta -dbtype prot \
    -title apis_apamin_flanks -out apis_apamin_flanks


    cd ../
    mkdir 201230_apamin_flank

    ## the program recommends using the outdated NCBI BLAST. However, BLAST+ has several improvements over the legacy BLAST applications. 

    blastp -num_threads 8 -query ../GCF_003254395.2_Amel_HAv3.1_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Apis_mellifera \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_000184785.3_Aflo_1.1_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Apis_florea \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_000188075.2_Si_gnH_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Solenopsis_invicta \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_000214255.1_Bter_1.0_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Bombus_terrestris \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_000220905.1_MROT_1.0_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Megachile_rotundata \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_000469605.1_Apis_dorsata_1.3_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Apis_dorsata \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_001263275.1_ASM126327v1_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Habropoda_laboriosa \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_001313835.1_ASM131383v1_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Polistes_canadensis \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_001442555.1_ACSNU-2.0_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Apis_cerana \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_001465965.1_Pdom_r1.2_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Polistes_dominula \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_003227725.1_Cflo_v7.5_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Camponotus_floridanus \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_004153925.1_Obicornis_v3_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Osmia_bicornis \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_009193385.2_Nvit_psr_1.1_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Nasonia_vitripennis \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_011952255.1_Bvos_JDL3184-5_v1.1_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Bombus_vosnesenskii \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

    blastp -num_threads 8 -query ../GCF_014083535.2_V.mandarinia_Nanaimo_p1.0_protein.faa \
    -db 201230_blast_results_orthologs/apis_apamin_flanks \
    -out 201230_apamin_flank/Vespa_mandarinia \
    -outfmt 6 \
    -evalue 0.05  -max_target_seqs 6

Concatenate all blast results

    rm 201230_apamin_flank/blast_results.blast
    cat 201230_apamin_flank/* > 201230_apamin_flank/blast_results.blast

Get orthologs for apamin

      awk '{print $2}' 201230_apamin_flank/blast_results.blast > apamin_apis
      awk '{print $1}' 201230_apamin_flank/blast_results.blast > apamin_orthologs

Edit the gff file. The porogram only allows one entry per gene. So merge all genes that overlap with the same name using cagt

     conda install -c conda-forge -c bioconda cgat-apps

Activate bioconda

     conda activate
     rm 201231_gff_files/* # make sure other files are deleted


Most insect genomes are just scafolds, and their names are not compatible with MCScanX, so thelables were replaced by a fixed number

_Apis mellifera_ preparation


     gff2bed < ../GCF_003254395.2_Amel_HAv3.1_filtered.gff >  \
     GCF_003254395.2_Amel_HAv3.bed

     cgat bed2bed --method=merge --merge-by-name -I  GCF_003254395.2_Amel_HAv3.bed > GCF_003254395.2_Amel_HAv3_merged.bed


Filter bed file to include only the chromosomes where apamin orthologs are present

    grep -Fwf apamin_apis GCF_003254395.2_Amel_HAv3_merged.bed > GCF_003254395.2_Amel_HAv3_merged_filtered.bed

     ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

     awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_003254395.2_Amel_HAv3_merged_filtered.bed| \
     sed '/^[[:blank:]]*#/d;s/#.*//' | \
     awk '$6 = "am1" $(NF+1)' | \
     awk '{print $6,$4,$2,$3}'  OFS='\t' > 201231_gff_files/GCF_003254395.2_Amel_HAv3_merged_filtered_duplicates.bed

_Apis mellifera_

     gff2bed < ../GCF_003254395.2_Amel_HAv3.1_genomic.gff >  \
     GCF_003254395.2_Amel_HAv3._1.bed

     cgat bed2bed --method=merge --merge-by-name -I GCF_003254395.2_Amel_HAv3._1.bed > GCF_003254395.2_Amel_HAv3._1_merged.bed

     #Filter bed file to include only the chromosomes where apamin orthologs are present
     grep -Fwf apamin_orthologs GCF_003254395.2_Amel_HAv3._1_merged.bed > GCF_003254395.2_Amel_HAv3._1_merged_filtered.bed

     ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

     awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_003254395.2_Amel_HAv3._1_merged_filtered.bed | \
     sed '/^[[:blank:]]*#/d;s/#.*//' | \
     awk '$6 = "am1" $(NF+1)' | \
     awk '{print $6,$4,$2,$3}' OFS='\t' > 201231_gff_files/GCF_003254395.2_Amel_HAv3._1_merged_filtered_duplicates.bed


_Apis florea_

     gff2bed < ../GCF_000184785.3_Aflo_1.1_genomic.gff >  \
     GCF_000184785.3_Aflo_1.1.bed

     cgat bed2bed --method=merge --merge-by-name -I GCF_000184785.3_Aflo_1.1.bed > GCF_000184785.3_Aflo_1.1_merged.bed

     grep -Fwf apamin_orthologs GCF_000184785.3_Aflo_1.1_merged.bed > GCF_000184785.3_Aflo_1.1_merged_filtered.bed

     ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

     awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_000184785.3_Aflo_1.1_merged_filtered.bed | \
     sed '/^[[:blank:]]*#/d;s/#.*//' | \
     awk '$6 = "af1" $(NF+1)' | \
     awk '{print $6,$4,$2,$3}' OFS='\t' > 201231_gff_files/GCF_000184785.3_Aflo_1.1_merged_filtered_duplicates.bed


_Solenopsis invicta_

    gff2bed < ../GCF_000188075.2_Si_gnH_genomic.gff >  \
    GCF_000188075.2_Si_gnH.bed

    cgat bed2bed --method=merge --merge-by-name -I GCF_000188075.2_Si_gnH.bed > GCF_000188075.2_Si_gnH_merged.bed

    grep -Fwf apamin_orthologs GCF_000188075.2_Si_gnH_merged.bed > GCF_000188075.2_Si_gnH_merged_filtered.bed

    ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

    awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_000188075.2_Si_gnH_merged_filtered.bed | \
    sed '/^[[:blank:]]*#/d;s/#.*//' | \
    awk '$6 = "si1" $(NF+1)' | \
    awk '{print $6,$4,$2,$3}' OFS='\t' > 201231_gff_files/GCF_000188075.2_Si_gnH_merged_filtered_duplicates.bed

_Bombus terrestris_


    gff2bed < ../GCF_000214255.1_Bter_1.0_genomic.gff >  \
    GCF_000214255.1_Bter.bed

    cgat bed2bed --method=merge --merge-by-name -I GCF_000214255.1_Bter.bed > GCF_000214255.1_Bter_merged.bed
    
    grep -Fwf apamin_orthologs GCF_000214255.1_Bter_merged.bed > GCF_000214255.1_Bter_merged_filtered.bed

    ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

    awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_000214255.1_Bter_merged_filtered.bed | \
    sed '/^[[:blank:]]*#/d;s/#.*//' | \
    awk '$6 = "bt1" $(NF+1)' | \
    awk '{print $6,$4,$2,$3}' OFS='\t' > 201231_gff_files/GCF_000214255.1_Bter_merged_filtered_duplicates.bed

_Megachile rotundata_

    gff2bed < ../GCF_000220905.1_MROT_1.0_genomic.gff >  \
    GCF_000220905.1_MROT.bed

    cgat bed2bed --method=merge --merge-by-name -I GCF_000220905.1_MROT.bed > GCF_000220905.1_MROT_merged.bed

    grep -Fwf apamin_orthologs GCF_000220905.1_MROT_merged.bed > GCF_000220905.1_MROT_merged_filtered.bed

    ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

    awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}'  GCF_000220905.1_MROT_merged_filtered.bed | \
    sed '/^[[:blank:]]*#/d;s/#.*//' | \
    awk '$6 = "mr1" $(NF+1)' | \
    awk '{print $6,$4,$2,$3}' OFS='\t' > 201231_gff_files/GCF_000220905.1_MROT_merged_filtered_duplicates.bed

_Apis dorsata_

    gff2bed < ../GCF_000469605.1_Apis_dorsata_1.3_genomic.gff >  \
    GCF_000469605.1_Apis_dorsata.bed

    cgat bed2bed --method=merge --merge-by-name -I GCF_000469605.1_Apis_dorsata.bed > GCF_000469605.1_Apis_dorsata_merged.bed

    grep -Fwf apamin_orthologs GCF_000469605.1_Apis_dorsata_merged.bed > GCF_000469605.1_Apis_dorsata_merged_filtered.bed

    ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

    awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_000469605.1_Apis_dorsata_merged_filtered.bed | \
    sed '/^[[:blank:]]*#/d;s/#.*//' | \
    awk '$6 = "ad1" $(NF+1)' | \
    awk '{print $6,$4,$2,$3}' OFS='\t' > 201231_gff_files/GCF_000469605.1_Apis_dorsata_merged_filtered_duplicates.bed

_Habropoda laboriosa_

    gff2bed < ../GCF_001263275.1_ASM126327v1_genomic.gff >  \
    GCF_001263275.1_ASM126327v1.bed

    cgat bed2bed --method=merge --merge-by-name -I GCF_001263275.1_ASM126327v1.bed > GCF_001263275.1_ASM126327v1_merged.bed

    grep -Fwf apamin_orthologs GCF_001263275.1_ASM126327v1_merged.bed > GCF_001263275.1_ASM126327v1_merged_filtered.bed

    ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

    awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_001263275.1_ASM126327v1_merged_filtered.bed| \
    sed '/^[[:blank:]]*#/d;s/#.*//' | \
    awk '$6 = "hl1" $(NF+1)' | \
    awk '{print $6,$4,$2,$3}' OFS='\t' > 201231_gff_files/GCF_001263275.1_ASM126327v1_merged_filtered_duplicates.bed

_Polistes canadensis_

    gff2bed < ../GCF_001313835.1_ASM131383v1_genomic.gff >  \
    GCF_001313835.1_ASM131383v1.bed

    cgat bed2bed --method=merge --merge-by-name -I GCF_001313835.1_ASM131383v1.bed > GCF_001313835.1_ASM131383v1_merged.bed

    grep -Fwf apamin_orthologs GCF_001313835.1_ASM131383v1_merged.bed > GCF_001313835.1_ASM131383v1_merged_filtered.bed

    ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

    awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_001313835.1_ASM131383v1_merged_filtered.bed | \
    sed '/^[[:blank:]]*#/d;s/#.*//' | \
    awk '$6 = "pc1" $(NF+1)' | \
    awk '{print $6,$4,$2,$3}' OFS='\t' > 201231_gff_files/GCF_001313835.1_ASM131383v1_merged_filtered_duplicates.bed

_Apis cerana_

    gff2bed < ../GCF_001442555.1_ACSNU-2.0_genomic.gff >  \
    GCF_001442555.1_ACSNU.bed

    cgat bed2bed --method=merge --merge-by-name -I GCF_001442555.1_ACSNU.bed > GCF_001442555.1_ACSNU_merged.bed

    grep -Fwf apamin_orthologs GCF_001442555.1_ACSNU_merged.bed > GCF_001442555.1_ACSNU_merged_filtered.bed

    ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

    awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_001442555.1_ACSNU_merged_filtered.bed | \
    sed '/^[[:blank:]]*#/d;s/#.*//' | \
    awk '$6 = "ac1" $(NF+1)' | \
    awk '{print $6,$4,$2,$3}' OFS='\t' > 201231_gff_files/GCF_001442555.1_ACSNU_merged_filtered_duplicates.bed

_Polistes dominula_

    gff2bed < ../GCF_001465965.1_Pdom_r1.2_genomic.gff >  \
    GCF_001465965.1_Pdom.bed

    cgat bed2bed --method=merge --merge-by-name -I GCF_001465965.1_Pdom.bed > GCF_001465965.1_Pdom_merged.bed

    grep -Fwf apamin_orthologs GCF_001465965.1_Pdom_merged.bed > GCF_001465965.1_Pdom_merged_filtered.bed

    ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

    awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_001465965.1_Pdom_merged_filtered.bed | \
    sed '/^[[:blank:]]*#/d;s/#.*//' | \
    awk '$6 = "pd1" $(NF+1)' | \
    awk '{print $6,$4,$2,$3}' OFS='\t'  > 201231_gff_files/GCF_001465965.1_Pdom_merged_filtered_duplicates.bed

_Camponotus floridanus_

    gff2bed < ../GCF_003227725.1_Cflo_v7.5_genomic.gff >  \
    GCF_003227725.1_Cflo.bed

    cgat bed2bed --method=merge --merge-by-name -I GCF_003227725.1_Cflo.bed > GCF_003227725.1_Cflo_merged.bed

    grep -Fwf apamin_orthologs GCF_003227725.1_Cflo_merged.bed > GCF_003227725.1_Cflo_merged_filtered.bed

    ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

    awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_003227725.1_Cflo_merged_filtered.bed | \
    sed '/^[[:blank:]]*#/d;s/#.*//' | \
    awk '$6 = "cf1" $(NF+1)' | \
    awk '{print $6,$4,$2,$3}' OFS='\t'  > 201231_gff_files/GCF_003227725.1_Cflo_merged_filtered_duplicates.bed

_Osmia bicornis_

    gff2bed < ../GCF_004153925.1_Obicornis_v3_genomic.gff >  \
    GCF_004153925.1_Obicornis.bed

    cgat bed2bed --method=merge --merge-by-name -I GCF_004153925.1_Obicornis.bed > GCF_004153925.1_Obicornis_merged.bed

    grep -Fwf apamin_orthologs GCF_004153925.1_Obicornis_merged.bed > GCF_004153925.1_Obicornis_merged_filtered.bed

    ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

    awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_004153925.1_Obicornis_merged_filtered.bed | \
    sed '/^[[:blank:]]*#/d;s/#.*//' | \
    awk '$6 = "ob1" $(NF+1)' | \
    awk '{print $6,$4,$2,$3}' OFS='\t'  > 201231_gff_files/GCF_004153925.1_Obicornis_merged_filtered_duplicates.bed

_Nasonia vitripennis_

    gff2bed < ../GCF_009193385.2_Nvit_psr_1.1_genomic.gff >  \
    GCF_009193385.2_Nvit_psr.bed

    cgat bed2bed --method=merge --merge-by-name -I GCF_009193385.2_Nvit_psr.bed > GCF_009193385.2_Nvit_psr_merged.bed

    grep -Fwf apamin_orthologs GCF_009193385.2_Nvit_psr_merged.bed > GCF_009193385.2_Nvit_psr_merged_filtered.bed

    ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

    awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_009193385.2_Nvit_psr_merged_filtered.bed | \
    sed '/^[[:blank:]]*#/d;s/#.*//' | \
    awk '$6 = "nv1" $(NF+1)' | \
    awk '{print $6,$4,$2,$3}' OFS='\t'  > 201231_gff_files/GCF_009193385.2_Nvit_psr_merged_filtered_duplicates.bed

_Bombus vosnesenskii_

    gff2bed < ../GCF_011952255.1_Bvos_JDL3184-5_v1.1_genomic.gff >  \
    GCF_011952255.1_Bvos_JDL3184-5.bed

    cgat bed2bed --method=merge --merge-by-name -I GCF_011952255.1_Bvos_JDL3184-5.bed > GCF_011952255.1_Bvos_JDL3184-5_merged.bed

    grep -Fwf apamin_orthologs GCF_011952255.1_Bvos_JDL3184-5_merged.bed > GCF_011952255.1_Bvos_JDL3184-5_merged_filtered.bed

    ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

    awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_011952255.1_Bvos_JDL3184-5_merged_filtered.bed | \
    sed '/^[[:blank:]]*#/d;s/#.*//' | \
    awk '$6 = "bv1" $(NF+1)' | \
    awk '{print $6,$4,$2,$3}' OFS='\t'  > 201231_gff_files/GCF_011952255.1_Bvos_JDL3184-5_merged_filtered_duplicates.bed

_Vespa mandarinia_

    gff2bed < ../GCF_014083535.2_V.mandarinia_Nanaimo_p1.0_genomic.gff >  \
    GCF_014083535.2_V.mandarinia_Nanaimo.bed

    cgat bed2bed --method=merge --merge-by-name -I GCF_014083535.2_V.mandarinia_Nanaimo.bed > GCF_014083535.2_V.mandarinia_Nanaimo_merged.bed

    grep -Fwf apamin_orthologs GCF_014083535.2_V.mandarinia_Nanaimo_merged.bed > GCF_014083535.2_V.mandarinia_Nanaimo_merged_filtered.bed

    ## rename with array the gene duplicates and rename chromosomes to add the species id and clean up

    awk '$4 in a {$4=$4 "_" ++a[$4]}{a[$4];print}' GCF_014083535.2_V.mandarinia_Nanaimo_merged_filtered.bed | \
    sed '/^[[:blank:]]*#/d;s/#.*//' | \
    awk '$6 = "vm1" $(NF+1)' | \
    awk '{print $6,$4,$2,$3}' OFS='\t' > 201231_gff_files/GCF_014083535.2_V.mandarinia_Nanaimo_merged_filtered_duplicates.bed

Concatenate all cleaned bed files

    cat 201231_gff_files/*filtered_duplicates.bed > 201231_gff_files/merged_gff_files.gff

There is a lot of duplicated genes on the gff file that have been renamed, so the BLAST results needs to be ajusted to reflect the change

    conda deactivate

    R < prepare_blast_MCScanX.R --no-save

Put all MCScanX files on the same folder

    sed 's/cds-//g' 201231_gff_files/merged_gff_files.gff | \
    sed 's/rna-//g'  | \
    sed 's/exon-//g' > 201231_gff_files/merged_gff_files_eddited.gff


    rm -rf 201231_mcscan_files
    mkdir 201231_mcscan_files

    cp 201230_apamin_flank/blast_results_duplicated.blast 201231_mcscan_files/apamin.blast
    cp 201231_gff_files/merged_gff_files_eddited.gff 201231_mcscan_files/apamin.gff

Run MCScanX

    MCScanX -s 2 -e 0.05 -b 2 201231_mcscan_files/apamin

Reformulate the results file to use on R

    less 201231_mcscan_files/apamin.collinearity | \
    cut -f 2-  | \
    awk 'NR > 11 { print }' | \
    sed 's/## Alignment .*am1&/###Apis mellifera /g'  | \
    sed 's/am1.*/Apis_mellifera/g'  | \
    sed 's/af1.*/Apis_florea/g'  | \
    sed 's/si1.*/Solenopsis_invicta/g'  | \
    sed 's/mr1.*/Megachile_rotundata/g'  | \
    sed 's/bt1.*/Bombus_terrestris/g'  | \
    sed 's/ad1.*/Apis_dorsata/g'  | \
    sed 's/hl1.*/Habropoda_laboriosa/g'  | \
    sed 's/pc1.*/Polistes_canadensis/g'  | \
    sed 's/ac1.*/Apis_cerana/g'  | \
    sed 's/pd1.*/Polistes_dominula/g'  | \
    sed 's/cf1.*/Camponotus_floridanus/g'  | \
    sed 's/ob1.*/Osmia_bicornis/g'  | \
    sed 's/nv1.*/Nasonia_vitripennis/g'  | \
    sed 's/bv1.*/Bombus_vosnesenskii/g'  | \
    sed 's/vm1.*/Vespa_mandarinia/g'  > \
    201231_mcscan_files/apamin.collinearity_filtered

    grep '###' 201231_mcscan_files/apamin.collinearity_filtered  | \
    sed 's/###//g' > 201231_mcscan_files/apamin.collinearity_filtered_species


Prepare the bed file to import on r


    cgat bed2bed --method=merge --merge-by-name -I  GCF_003254395.2_Amel_HAv3.bed > GCF_003254395.2_Amel_HAv3_cds_filtered_merged.bed

Obtain gene list adjacent to apamin and MCD

    grep -A3 -A3 -E 'NP_001011611.2|NP_001011612.1' GCF_003254395.2_Amel_HAv3.1_genomic_merged_flank.bed> GCF_003254395.2_Amel_HAv3.1_genomic_merged_flanks_filtered_mcscan.bed  ##grep the flanking genes

    awk '{print $4}' GCF_003254395.2_Amel_HAv3.1_genomic_merged_flanks_filtered_mcscan.bed |  \
    cut -d , -f2 | sort | uniq > apis_flanking_genes_filtered_mcscan

     grep -Fwf apis_flanking_genes_filtered_mcscan 201230_apamin_flank/blast_results.blast > 201230_apamin_flank/blast_results_flank_filtered.blast

        #_Apis florea_

        cgat bed2bed --method=merge --merge-by-name -I GCF_000184785.3_Aflo_1.1.bed > GCF_000184785.3_Aflo_1.1_cds_filtered_merged.bed

        #_Solenopsis invicta_

        cgat bed2bed --method=merge --merge-by-name -I GCF_000188075.2_Si_gnH.bed > GCF_000188075.2_Si_gnH_cds_filtered_merged.bed

        #_Bombus terrestris_

        cgat bed2bed --method=merge --merge-by-name -I GCF_000214255.1_Bter.bed > GCF_000214255.1_Bter_cds_filtered_merged.bed


        #_Megachile rotundata_

        cgat bed2bed --method=merge --merge-by-name -I GCF_000220905.1_MROT.bed > GCF_000220905.1_MROT_cds_filtered_merged.bed

        #_Apis dorsata_

        cgat bed2bed --method=merge --merge-by-name -I GCF_000469605.1_Apis_dorsata.bed > GCF_000469605.1_Apis_dorsata_cds_filtered_merged.bed

        #_Habropoda laboriosa_

        cgat bed2bed --method=merge --merge-by-name -I GCF_001263275.1_ASM126327v1.bed > GCF_001263275.1_ASM126327v1_cds_filtered_merged.bed

        #_Polistes canadensis_

        cgat bed2bed --method=merge --merge-by-name -I GCF_001313835.1_ASM131383v1.bed > GCF_001313835.1_ASM131383v1_cds_filtered_merged.bed

        #_Apis cerana_

        cgat bed2bed --method=merge --merge-by-name -I GCF_001442555.1_ACSNU.bed > GCF_001442555.1_ACSNU_cds_filtered_merged.bed

        #_Polistes dominula_

        cgat bed2bed --method=merge --merge-by-name -I GCF_001465965.1_Pdom.bed > GCF_001465965.1_Pdom_cds_filtered_merged.bed

        #_Camponotus floridanus_

        cgat bed2bed --method=merge --merge-by-name -I GCF_003227725.1_Cflo.bed > GCF_003227725.1_Cflo_cds_filtered_merged.bed

        #_Osmia bicornis_

        cgat bed2bed --method=merge --merge-by-name -I GCF_004153925.1_Obicornis.bed > GCF_004153925.1_Obicornis_cds_filtered_merged.bed

        #_Nasonia vitripennis_

        cgat bed2bed --method=merge --merge-by-name -I GCF_009193385.2_Nvit_psr.bed > GCF_009193385.2_Nvit_psr_cds_filtered_merged.bed

        #_Bombus vosnesenskii_

        cgat bed2bed --method=merge --merge-by-name -I GCF_011952255.1_Bvos_JDL3184-5.bed > GCF_011952255.1_Bvos_JDL3184-5_cds_filtered_merged.bed

        #_Vespa mandarinia_

        cgat bed2bed --method=merge --merge-by-name -I GCF_014083535.2_V.mandarinia_Nanaimo.bed > GCF_014083535.2_V.mandarinia_Nanaimo_cds_filtered_merged.bed
