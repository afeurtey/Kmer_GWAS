#!/bin/bash
#SBATCH --time=120:00:00          # Walltime
#SBATCH --mem-per-cpu=3000

# Usage exemple: sbatch  Run_GWAS.sh ./Phenotypes/Phenotype_temperature_merged_lofoutrem.csv 

# Note:
#wget https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-59.tar.gz
#R CMD INSTALL MASS_7.3-59.tar.gz
#Read first argument of the script
phenotype_file=$1
echo -e "\nReading in" $phenotype_file

module load gcc/11.4.0
source /cluster/home/afeurtey/.bashrc
#mamba init
conda activate kmersGWAS
cd /cluster/home/afeurtey/Silvia/KmerGWAS
working_dir=/cluster/scratch/afeurtey/Kmer_GWAS/4_All_kmers_GWAS/
GWAS_dir=/cluster/scratch/afeurtey/Kmer_GWAS/5_GWAS_results/

condition_list=$(head -n 1 $phenotype_file | tr "," "\n" | tail -n +2)

c=1
for i in $condition_list;
do
    echo "... Condition is : " $i
    module load gcc/11.4.0
    #Increment the column index
    c=$((c+1))

    echo -e "accession_id\tphenotype_value" > ${GWAS_dir}phenotypes_${i}.pheno
    #cut -f 1,$c -d "," $phenotype_file | tail -n +2 | grep -v ",NA" | sed 's/,/\t/' >> ${GWAS_dir}phenotypes_${i}.pheno
    cut -f 1,$c -d "," $phenotype_file | tail -n +2 | grep -v ",NA" | grep -v ",$" | sed 's/,/\t/' >> ${GWAS_dir}phenotypes_${i}.pheno
    python2 /cluster/home/afeurtey/Silvia/KmerGWAS/voichek_kmersGWAS/kmers_gwas.py \
    --pheno ${GWAS_dir}phenotypes_${i}.pheno \
    --kmers_table ${working_dir}kmers_table \
    -l 31 -p 8 --mac 10 \
    --outdir ${GWAS_dir}GWAS_output_dir_${i}

  cat ${GWAS_dir}GWAS_output_dir_${i}/kmers/pass_threshold_10per | tail -n +2 | cut -f 2 | sed 's/_/\t/' | awk '{print ">"$2"\n"$1}' > \
    ${GWAS_dir}GWAS_output_dir_${i}/kmers/pass_threshold_10per_for_fetch.fasta ; 

  if [ -s ${GWAS_dir}GWAS_output_dir_${i}/kmers/pass_threshold_10per_for_fetch.fasta ]; then
    # Blast the kmers on the reference genome
    module load gcc/8.2.0 blast-plus/2.12.0
    blastn \
        -db /cluster/scratch/afeurtey/Kmer_GWAS/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa \
        -query ${GWAS_dir}GWAS_output_dir_${i}/kmers/pass_threshold_10per_for_fetch.fasta -outfmt 6 > \
        ${GWAS_dir}GWAS_output_dir_${i}/kmers/pass_threshold_10per_for_fetch.blast.tsv ; 
  fi
done

