#!/bin/bash
#SBATCH --time=120:00:00          # Walltime

module load gcc/11.4.0
cd /cluster/home/afeurtey/Silvia/KmerGWAS
working_dir=/cluster/scratch/afeurtey/Kmer_GWAS/4_All_kmers_GWAS/

# Format of the kmers_list_paths file is two tab-separated columns, the first one with the full path, the other with the sample name
# /cluster/scratch/afeurtey/Kmer_GWAS/3_Kmer_count/Algeria_1992_ALA2a/kmers_with_strand   Algeria_1992_ALA2a
# Examples of commands to produce the file
# Example 1
# with awk split the path and get last directory name
# ls -1 /cluster/scratch/afeurtey/Kmer_GWAS/3_Kmer_count/Only_Silvias/*/kmers_with_strand |  \
#   awk 'BEGIN{OFS="\t"} {  split($0, a, "/"); print $0, a[8] }'  \
#   > /cluster/scratch/afeurtey/Kmer_GWAS/4_All_kmers_GWAS/Silvias_kmers_list_paths.txt #no_kmers.list
# Example 2
#ll /cluster/scratch/afeurtey/Kmer_GWAS/3_Kmer_count/ | tail -n +2 | \
#  awk '{printf "/cluster/scratch/afeurtey/Kmer_GWAS/3_Kmer_count/%s/kmers_with_strand\t%s\n", $NF,$NF}' \
#  > kmers_list_paths.txt

#diff -y no_kmers.list kmers_list_paths.txt

./voichek_kmersGWAS/bin/list_kmers_found_in_multiple_samples \
  -l ${working_dir}Silvias_kmers_list_paths.txt \
  -o ${working_dir}Silvias_kmers_to_use \
  --mac 5 -p 0.2 -k 31

./voichek_kmersGWAS/bin/build_kmers_table \
  -l ${working_dir}Silvias_kmers_list_paths.txt \
  -a ${working_dir}Silvias_kmers_to_use \
  -k 31 \
  -o ${working_dir}Silvias_kmers_table

# ./voichek_kmersGWAS/bin/kmers_table_to_bed \
#   -t ${working_dir}kmers_table \
#   -k 31 \
#   -p ${working_dir}phenotypes.pheno \
#   --maf 0.05 --mac 5 -b 10000000 \
#   -o ${working_dir}output_file

#Calculation of kinship matrix
./voichek_kmersGWAS/bin/emma_kinship_kmers \
  -t ${working_dir}Silvias_kmers_table \
  -k 31 --maf 0.05 \
  > ${working_dir}Silvias_kmers_table.kinship


