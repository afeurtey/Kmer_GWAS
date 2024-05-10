#/bin/bash

# Usage : sh Run_blast_kmers_and_plots.sh ./Phenotypes/Phenotype_temperature_merged_lofoutrem.csv

# Set environment 
module load gcc/11.4.0
source /cluster/home/afeurtey/.bashrc
conda activate kmersGWAS
cd /cluster/home/afeurtey/Silvia/KmerGWAS

# Define variables
phenotype_file=$1
working_dir=/cluster/scratch/afeurtey/Kmer_GWAS/4_All_kmers_GWAS/
GWAS_dir=/cluster/scratch/afeurtey/Kmer_GWAS/5_GWAS_results/


echo -e "\n\n\n\n ---- \n\n Reading in" $phenotype_file

#Per condition run kmers blast, merge results with association and plot
condition_list=$(head -n 1 $phenotype_file | tr "," "\n" | tail -n +2)
c=1
for i in $condition_list;
do
    echo "... Condition is : " $i
    module load gcc/11.4.0
    #Increment the column index
    c=$((c+1))

    sh Blast_all_kmers.sh $i

    nb_significant_kmers=$(cat ${GWAS_dir}GWAS_output_dir_${i}/kmers/pass_threshold_10per | tail -n +2 | wc -l )
    nb_mapping_kmers=$(wc -l ${GWAS_dir}GWAS_output_dir_${i}/kmers/pass_threshold_10per_for_fetch.blast.tsv )

    echo "... Found $nb_significant_kmers significant kmers and $nb_mapping_kmers mapped kmers for condition $i"

    if [ -s ${GWAS_dir}GWAS_output_dir_${i}/kmers/pass_threshold_10per_for_fetch.blast.tsv ]; then
        #
        python Manhattan_plots_kmer_GWAS.py $i
    fi
done