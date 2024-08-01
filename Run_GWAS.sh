#!/bin/bash
#SBATCH --time=120:00:00          # Walltime
#SBATCH --mem-per-cpu=3000

# Usage exemple: sbatch  Run_GWAS.sh ./Phenotypes/Phenotype_temperature_merged_lofoutrem.csv 

# Note on installing MASS:
# ------------------------
# Sometimes it is difficult to install this package. Here is how to do it if conda etc does not solve it.
#wget https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-59.tar.gz
#R CMD INSTALL MASS_7.3-59.tar.gz


# Paths to set for your specific environment
# ------------------------------------------
cd /cluster/home/afeurtey/Silvia/KmerGWAS
working_dir=/cluster/scratch/afeurtey/Kmer_GWAS/4_All_kmers_GWAS/
GWAS_dir=/cluster/scratch/afeurtey/Kmer_GWAS/5_GWAS_results/
new_GEMMA_path=/cluster/home/afeurtey/Silvia/KmerGWAS/gemma-0.98.5-linux-static-AMD64
blastDB=/cluster/scratch/afeurtey/Kmer_GWAS/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa # To create blast reference DB: ~/Silvia/KmerGWAS$ makeblastdb  -dbtype nucl -in /cluster/scratch/afeurtey/Kmer_GWAS/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa
kmer_script_path=/cluster/home/afeurtey/Silvia/KmerGWAS/voichek_kmersGWAS/kmers_gwas.py
kmer_table=${working_dir}Silvias_kmers_table
GWAS_output_prefix=Silvias_GWAS_output_dir_

# Loading modules for new Euler running UBUNTU
# --------------------------------------------
#module load gcc/11.4.0
#mamba init
source /cluster/home/afeurtey/.bashrc
conda activate kmersGWAS
source /cluster/project/gdc/shared/stack/GDCstack.sh #new
module load python/2.7.18 #new


#Read first argument of the script
phenotype_file=$1
echo -e "\nReading in" $phenotype_file

condition_list=$(head -n 1 $phenotype_file | tr "," "\n" | tail -n +2)

c=1
for i in $condition_list;
do
    echo "... Condition is : " $i
    #Increment the column index
    c=$((c+1))

    echo -e "accession_id\tphenotype_value" > ${GWAS_dir}phenotypes_${i}.pheno
    #cut -f 1,$c -d "," $phenotype_file | tail -n +2 | grep -v ",NA" | sed 's/,/\t/' >> ${GWAS_dir}phenotypes_${i}.pheno
    cut -f 1,$c -d "," $phenotype_file | tail -n +2 | grep -v ",NA" | grep -v ",$" | sed 's/,/\t/' >> ${GWAS_dir}phenotypes_${i}.pheno


    python ${kmer_script_path} \
    --pheno ${GWAS_dir}phenotypes_${i}.pheno \
    --kmers_table ${kmer_table} \
    -l 31 -p 8 --mac 10 \
    --outdir ${GWAS_dir}${GWAS_output_prefix}${i} \
    --gemma_path ${new_GEMMA_path} #new

  cat ${GWAS_dir}${GWAS_output_prefix}${i}/kmers/pass_threshold_10per | tail -n +2 | cut -f 2 \
    | sed 's/_/\t/' | awk '{print ">"$2"\n"$1}' > \
    ${GWAS_dir}${GWAS_output_prefix}${i}/kmers/pass_threshold_10per_for_fetch.fasta ; 

  # Blast the significant kmers on the reference genome
  if [ -s ${GWAS_dir}${GWAS_output_prefix}${i}/kmers/pass_threshold_10per_for_fetch.fasta ]; then
    
    #module load gcc/8.2.0 blast-plus/2.12.0
    module load stack/2024-06  gcc/12.2.0 blast-plus/2.14.1
    blastn \
        -db ${blastDB} \
        -query ${GWAS_dir}${GWAS_output_prefix}${i}/kmers/pass_threshold_10per_for_fetch.fasta -outfmt 6 > \
        ${GWAS_dir}${GWAS_output_prefix}${i}/kmers/pass_threshold_10per_for_fetch.blast.tsv ; 

  fi
done

