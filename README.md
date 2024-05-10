# Kmer_GWAS



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