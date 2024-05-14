#!/bin/bash
#SBATCH --mem-per-cpu=4000 
#SBATCH --cpus-per-task=3 

#Usage: sbatch --array=1-123%20 ./KmerGWAS_count_kmers.sh ./Reseq_Silvia_full_read_name.txt
# temp_read_name.txt is a file with the following format:
# sample_name read_prefix
sample_list=$1 #Space separated


script_dir=/cluster/home/afeurtey/Silvia/KmerGWAS/
raw_reads_dir=/cluster/scratch/afeurtey/Kmer_GWAS/0_Raw_data/
filtered_reads_dir=/cluster/scratch/afeurtey/Kmer_GWAS/1_Filtered_data/
kmer_count_dir=/cluster/scratch/afeurtey/Kmer_GWAS/3_Kmer_count/
adapters=/cluster/home/afeurtey/Gene_coexpression_network/Adapters.fasta


# Load modules
# ------------
source /cluster/apps/local/env2lmod.sh
module load gcc/6.3.0 python/3.6.5
module load trimmomatic/0.38


# Define needed variables
# -----------------------
IDX=${SLURM_ARRAY_TASK_ID}
echo ${IDX}
sample=$(sed -n ${IDX}p < ${sample_list} | cut -f 1 -d " ")
read_name=$(sed -n ${IDX}p < ${sample_list} | cut -f 2 -d " ")


echo "Read file is: ${raw_reads_dir}${read_name}_1.fastq.gz"

trimmomatic PE \
  ${raw_reads_dir}${read_name}_1.fq.gz ${raw_reads_dir}${read_name}_2.fq.gz \
  ${filtered_reads_dir}${sample}_1P.fq.gz ${filtered_reads_dir}${sample}_1U.fq.gz \
  ${filtered_reads_dir}${sample}_2P.fq.gz ${filtered_reads_dir}${sample}_2U.fq.gz \
  ILLUMINACLIP:${adapters}:2:30:10 \
  LEADING:15 TRAILING:15 SLIDINGWINDOW:5:15 MINLEN:50



module load gcc/11.4.0


#name=`sed -n ${IDX}p <read.list`
mkdir ${kmer_count_dir}${sample}
cd ${kmer_count_dir}${sample}
printf '%s\n' ${filtered_reads_dir}${sample}_1P.fq.gz  ${filtered_reads_dir}${sample}_2P.fq.gz > FILES_${sample} 
sed -i 's/\s\+/\n/g' FILES_${sample} 

${script_dir}voichek_kmersGWAS/external_programs/kmc_v3 \
    -k31 -t8 -ci2 -cs10000 \
    @FILES_${sample}  \
    output_kmc_cannon ./ 1> kmc_cannon.1 2> kmc_cannon.2 

${script_dir}voichek_kmersGWAS/external_programs/kmc_v3 \
    -k31 -t8 -ci0 -cs10000 -b \
    @FILES_${sample}  \
    output_kmc_all ./ 1> kmc_all.1 2> kmc_all.2 

${script_dir}voichek_kmersGWAS/bin/kmers_add_strand_information \
    -c output_kmc_cannon \
    -n output_kmc_all -k 31 \
    -o kmers_with_strand 1> kmc_strand.1 

#rm ${sample}/*.kmc*
#rm FILES_${name}

#Notes: 1-2 typo for last command, 3-4 no option on top, 5-6 options