# Kmer_GWAS in Zymoseptoria tritici

This repo is intended to contain the code needed to do a genome-wide association study using kmers instead of SNPs and Silvia's (and perhaps Cécile's) phenotypic data, collected with *Z. tritici* strains in *in vitro* experiments.

The scripts are very much based on the kmersGWAS pipeline from [Voichek](https://github.com/voichek/kmersGWAS). They were written to run on the ETH cluster Euler in the spring of 2024. Modules are loaded from the corresponding environment and there is no guarantee that they will work in any other systems or even an updated Euler cluster. All other packages used belong to a conda environment and can be recreated from a yaml file included here with the following command: 

```
conda env create -f kmersGWAS_environment.yml
```


## Counting kmers and finding common kmers

The first step is to run the kmer counts per sample. I have used the following script to do this:
```
sbatch --array=1-435%20 ./KmerGWAS_count_kmers.sh ./Sample_lists/Eschikon_batch2.txt 
```
The format of the file passed as argument is "sample_name read_prefix". Warning: the fasta extension is hard-coded in the script! you might need to change it for your own batch.

The step above is very fast. It takes around 10 minutes per sample. However, the next step is very very inefficient. For 700 samples, determining which kmers to keep took 24 hours and creating the large table took 5 days. Be prepared to set very long times if needed! I believe that Jess mentioned that she divided this file into smaller pieces to have them run in parallel? 

```
sbatch Commands_kmers_to_GWAS_Silvias_only.sh
```

The script above requires a file that lists the files created in the previous step. The kmers_list_paths file format is two tab-separated columns, the first one with the full path, the other with the sample name. (see file content for example commands of how to produce this file) Warning: output names are hard coded in the script, feel free to change them.

At the end of these two steps, you will have files per samples and most importantly two files that include all samples: one is the the kmer table and the other is the kinship matrix.  These are needed for the GWAS. 


## Running the GWAS and interpreting the results

```
sbatch  Run_GWAS.sh ./Phenotypes/Phenotype_temperature_merged_lofoutrem.csv 
```
The format of the phenotype file is a csv file with a first column including the sample name and the following ones being the phenotypic values. Examples can be found in the [Phenotypes](./Phenotypes) folder. The outputs from this steps are several phenotype files (starting with "pheno"), some logs files and a kmers folder that includes the files we are interested in. The GWAS results are fiven with two different thresholds (5% and 10%) estimated from permutations. The threshold values can be found in the two corresponding files (e.g. "threshold_10per"). The significant kmers are also recorded in the same folder (e.g. "pass_threshold_10per") as a table including the kmer sequence, the p-value and other info. The output subfolder contains similar information for all tested kmers and not only for the significant ones.

One of the advantages of the kmer GWAS is that it is not reference-based and can thus be used to detect regions or variants which are not in the reference genome. However, in order to visualize or interpret results, we are often forced to use reference genomes. I do not want us to be limited to a single reference genome, though, so I create a "reference" fasta file that includes many PacBio genomes from isolate of different origins. 

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

# Format the fasta file as a database for blast
module load gcc/8.2.0 blast-plus/2.12.0
makeblastdb -in /cluster/scratch/afeurtey/Kmer_GWAS/all_PacBio.fa    -input_type fasta -dbtype nucl

#Only reference 
makeblastdb -in /cluster/scratch/afeurtey/Kmer_GWAS/Zymoseptoria_tritici.MG2.dna.toplevel.mt+.fa   -input_type fasta -dbtype nucl
```

Once the fasta file with all reference genomes was generated, I want to obtain summary tables, blast of the kmers on all reference genomes, and (very imperfect) data visualization of the results. I have written a wrapper [Run_blast_kmers_and_plots.sh](./Run_blast_kmers_and_plots.sh) (see below) that calls two different scripts: 

 * the first one is [Blast_all_kmers.sh](./Blast_all_kmers.sh) that transforms the output from the pipeline into fasta files for blast to read and identifies locations for kmers. It reads the file *kmers/output/phenotype_value.assoc.txt.gz* produces a fasta file from the kmer sequences **kmers/output/phenotype_value.assoc_for_fetch.fasta** and outputs the blast results in **kmers/output/phenotype_value.assoc_for_fetch.blast_${blastdb_suffix}.tsv**, where blastdb_suffix is for example all_PacBio.
 
 * the second one [Manhattan_plots_kmer_GWAS.py](./Manhattan_plots_kmer_GWAS.py) produces a simple visualization of the blast results. This is good for a quick exploration but for each chromosome the width can be misleading as itdepends on the kmer alignment and not on the real start and stop. A proper data viz need to be recreated for publication. If I have a list of references, I produce subfiles of the blast as well as individual graphs.

```
sh Run_blast_kmers_and_plots.sh ./Phenotypes/Phenotype_temperature_merged_lofoutrem.csv /cluster/scratch/afeurtey/Kmer_GWAS/all_PacBio.fa all_PacBio
```

I noticed that often many kmers are found to be significant but only a small fraction are actually found on any of the references. I wanted to see if assembling the kmers into longer fragments would help, but could not find an automated way to look at this. So I simply did it manually! In some folders, you will find a file named **kmers/pass_threshold_10per_for_fetch_assemble.fasta** that contains this manual. For some, most kmers neatly fit into longer sequences (ex. eAUC_SA) but for some conditions, many kmers are still alone or with one or two other kmers.

    Note: This was actually quite interesting for Std_grey18C_08dpi as the fragments revealed complex polymorphism patterns that would have been really difficult to capture with SNP calling. We also capture the forward and reverse sequences which is nice. Perhaps something to highlight in the paper. 


I wanted to find kmers that were found is more than one GWAS, so I ran the following command: 

```
cat /cluster/work/gdc/shared/p723/Ztritici_Kmer_GWAS/3_Results_GWAS/Silvias_GWAS_output_dir_*/kmers/pass_threshold_10per | cut -f 2 | cut -f 1 -d "_" | sort | uniq -c | sort -n 
```
The output was then manually manipulated and analysed to obtain the file [Multiple_hits_kmers.txt](./Multiple_hits_kmers.txt)

Useful command to find the number of significant kmers per GWAS

```
wc -l /cluster/work/gdc/shared/p723/Ztritici_Kmer_GWAS/3_Results_GWAS/Silvias_GWAS_output_dir_*/kmers/pass_threshold_10per_for_fetch.fasta | sort -n 
```