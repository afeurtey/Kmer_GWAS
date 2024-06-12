#!/bin/bash

# Usage : sh Run_blast_kmers_and_plots.sh ./Phenotypes/Phenotype_temperature_merged_lofoutrem.csv /cluster/scratch/afeurtey/Kmer_GWAS/5_GWAS_results/Silvias_GWAS_output_dir_ /cluster/scratch/afeurtey/Kmer_GWAS/all_PacBio.fa all_PacBio

# Set environment 
module load gcc/11.4.0
source /cluster/home/afeurtey/.bashrc
conda activate kmersGWAS
cd /cluster/home/afeurtey/Silvia/KmerGWAS

# Define variables
phenotype_file=$1
GWAS_dir_path=$2 #/cluster/scratch/afeurtey/Kmer_GWAS/5_GWAS_results/Silvias_GWAS_output_dir_
blastdb_path=$3 #/cluster/scratch/afeurtey/Kmer_GWAS/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa
blastdb_suffix=$4 #IPO323


working_dir=/cluster/scratch/afeurtey/Kmer_GWAS/4_All_kmers_GWAS/
#GWAS_dir=/cluster/scratch/afeurtey/Kmer_GWAS/5_GWAS_results/

# Write a list of all single ref
uniq_ref=("IPO323")

echo -e "\n\n ---- \n\n Reading in" $phenotype_file

#Per condition run kmers blast, merge results with association and plot
condition_list=$(head -n 1 $phenotype_file | tr "," "\n" | tail -n +2)
c=1
for i in $condition_list;
do
    echo -e "\n\n ---- \n... Condition is : " $i
    module load gcc/11.4.0
    #Increment the column index
    c=$((c+1))

    # if [ $c -lt 57 ]; then
    #     continue
    # fi
    echo "... Running blast"
    sh Blast_all_kmers.sh $i $GWAS_dir_path $blastdb_path $blastdb_suffix

    nb_significant_kmers=$(cat ${GWAS_dir_path}${i}/kmers/pass_threshold_10per | tail -n +2 | wc -l )
    nb_mapping_kmers=$(wc -l ${GWAS_dir_path}${i}/kmers/pass_threshold_10per_for_fetch.blast.tsv )

    echo "... Found $nb_significant_kmers significant kmers for condition $i"
    echo "... Found $nb_mapping_kmers mappings for significant kmers for condition $i"

    # Check if the blastdb_suffix is in the list of unique references
    found=0
    for name in "${uniq_ref[@]}"; do
    if [ "$name" == "$blastdb_suffix" ]; then
        found=1
        break
    fi
    done

    if [ $found -gt 0 ]; then
        python Manhattan_plots_kmer_GWAS.py $i  $blastdb_suffix $GWAS_dir_path 
    else
        ref_list=$(cut -f 2 ${GWAS_dir_path}${i}/kmers/output/phenotype_value.assoc_for_fetch.blast_all_PacBio.tsv | cut -f 1 -d "." | sort | uniq )
        for ref in $ref_list; do
            echo "Working on reference genome " $ref " for condition " $i
            grep $ref ${GWAS_dir_path}${i}/kmers/output/phenotype_value.assoc_for_fetch.blast_${blastdb_suffix}.tsv > \
                      ${GWAS_dir_path}${i}/kmers/output/phenotype_value.assoc_for_fetch.blast_${ref}.tsv
            python Manhattan_plots_kmer_GWAS.py $i  $ref  $GWAS_dir_path 
        done

    fi

done