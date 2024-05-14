# Kmer_GWAS

Running kmers count per sample
```
sbatch --array=1-435%20 ./KmerGWAS_count_kmers.sh ./Sample_lists/Eschikon_batch2.txt 
```
The step above is very fast. It takes around 10 minutes per sample. 

However, the next step is very very inefficient. For 700 samples, determining which kmers to keep took 24 hours and creating the large table took 5 days. Be prepared to set very long times if needed!

```
sbatch Commands_kmers_to_GWAS.sh
```

Creating reference with all PacBio
```
wget https://github.com/crolllab/datasets/archive/refs/heads/master.zip
unzip master.zip 
rm -r datasets-master/Gene_conservation_Tralamazza_2021
rm -r datasets-master/Structural_variation_scripts_Badet_2020/
cd datasets-master/
ls 19-isolate_pangenome_BMC_Biology_2020/*/*.fa | grep -v "prot" | grep -v "trans" > all_genomes.list
ls European_Z.tritici_genomes_2023/*/*.fa | grep -v "prot" | grep -v "trans" >> all_genomes.list 
{ xargs cat < all_genomes.list ; } > /cluster/scratch/afeurtey/Kmer_GWAS/all_PacBio.fa
sed -i 's/AAAA>IR01/AAAA\n>IR01/' /cluster/scratch/afeurtey/Kmer_GWAS/all_PacBio.fa
sed -i 's/CCTTCG>3D1.chr_1/CCTTCG\n>3D1.chr_1/' /cluster/scratch/afeurtey/Kmer_GWAS/all_PacBio.fa 
sed -i -r -e 's/chr([0-9]{1}$)/chr0\1/' ../all_PacBio.fa
sed -i -r -e 's/chr_([0-9]{1}$)/chr_0\1/' ../all_PacBio.fa

sed -i -r -e 's/>([0-9]{1,2}) />IPO323.chr_0\1 /' /cluster/scratch/afeurtey/Kmer_GWAS/all_PacBio.fa
sed -i -r -e 's/>IPO323.chr_0([0-9]{2})/>IPO323.chr_\1/' /cluster/scratch/afeurtey/Kmer_GWAS/all_PacBio.fa

module load gcc/8.2.0 blast-plus/2.12.0
makeblastdb -in /cluster/scratch/afeurtey/Kmer_GWAS/all_PacBio.fa    -input_type fasta -dbtype nucl

#Only reference 
makeblastdb -in /cluster/scratch/afeurtey/Kmer_GWAS/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa   -input_type fasta -dbtype nucl
```

Running blast on all references
```
sh Run_blast_kmers_and_plots.sh ./Phenotypes/Phenotype_temperature_merged_lofoutrem.csv /cluster/scratch/afeurtey/Kmer_GWAS/all_PacBio.fa all_PacBio
```

for i in eAUC12C eAUC15C eAUC18C eAUC22C eAUC25C ; do   cat ${GWAS_dir}GWAS_output_dir_${i}/kmers/pass_threshold_5per | tail -n +2 | cut -f 2 | sed 's
/_/\t/' | awk '{print ">"$2"\n"$1}' >     ${GWAS_dir}GWAS_output_dir_${i}/kmers/pass_threshold_5per_for_fetch.fasta ;    if [ -s ${GWAS_dir}GWAS_output_dir_$
{i}/kmers/pass_threshold_5per_for_fetch.fasta ]; then     blastn         -db /cluster/scratch/afeurtey/Kmer_GWAS/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa
         -query ${GWAS_dir}GWAS_output_dir_${i}/kmers/pass_threshold_5per_for_fetch.fasta -outfmt 6 >         ${GWAS_dir}GWAS_output_dir_${i}/kmers/pass_thre
shold_5per_for_fetch.blast.tsv ;    fi; done


ls /cluster/work/gdc/shared/p723/Ztritici_Kmer_GWAS/Kmer_counts/ -1 | grep -v "dat" > List_done.list
grep -v -f List_done.list ./KmerGWAS/Sample_lists/allthereads.csv > ./KmerGWAS/Sample_lists/missed_all_reads.csv