#/bin/bash

#Usage:  bash Blast_all_kmers.sh 25total_pH /cluster/scratch/afeurtey/Kmer_GWAS/5_GWAS_results/Silvias_GWAS_output_dir_ /cluster/scratch/afeurtey/Kmer_GWAS/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa IPO323
#GWAS_dir=/cluster/scratch/afeurtey/Kmer_GWAS/5_GWAS_results/

i=$1 #q25total_pH
GWAS_dir_path=$2 
blastdb_path=$3 #/cluster/scratch/afeurtey/Kmer_GWAS/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa
blastdb_suffix=$4 #IPO323


fasta_file=${GWAS_dir_path}${i}/kmers/output/phenotype_value.assoc_for_fetch.fasta

zcat ${GWAS_dir_path}${i}/kmers/output/phenotype_value.assoc.txt.gz | tail -n +2 | cut -f 2 | \
    sed 's/_/\t/' | awk '{print ">"$1"_"$2"\n"$1}' > \
    $fasta_file ; 


if [ -s $fasta_file ]; then
  # Blast the kmers on the reference genome
  module load gcc/8.2.0 blast-plus/2.12.0
  blastn \
      -db $blastdb_path \
      -query $fasta_file -outfmt 6 > \
      ${GWAS_dir_path}${i}/kmers/output/phenotype_value.assoc_for_fetch.blast_${blastdb_suffix}.tsv ; 
else
    echo "Error: File $fasta_file does not exist."
fi
