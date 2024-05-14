#/bin/bash

#Usage:  bash Blast_all_kmers.sh 25total_pH /cluster/scratch/afeurtey/Kmer_GWAS/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa IPO323
GWAS_dir=/cluster/scratch/afeurtey/Kmer_GWAS/5_GWAS_results/

i=$1 #q25total_pH
blastdb_path=$2 #/cluster/scratch/afeurtey/Kmer_GWAS/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa
blastdb_suffix=$3 #IPO323

zcat ${GWAS_dir}GWAS_output_dir_${i}/kmers/output/phenotype_value.assoc.txt.gz | tail -n +2 | cut -f 2 | \
    sed 's/_/\t/' | awk '{print ">"$1"_"$2"\n"$1}' > \
    ${GWAS_dir}GWAS_output_dir_${i}/kmers/output/phenotype_value.assoc_for_fetch.fasta ; 


if [ -s ${GWAS_dir}GWAS_output_dir_${i}/kmers/output/phenotype_value.assoc_for_fetch.fasta ]; then
  # Blast the kmers on the reference genome
  module load gcc/8.2.0 blast-plus/2.12.0
  blastn \
      -db $blastdb_path \
      -query ${GWAS_dir}GWAS_output_dir_${i}/kmers/output/phenotype_value.assoc_for_fetch.fasta -outfmt 6 > \
      ${GWAS_dir}GWAS_output_dir_${i}/kmers/output/phenotype_value.assoc_for_fetch.blast_${blastdb_suffix}.tsv ; 
fi

