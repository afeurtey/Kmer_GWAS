**Material and methods**

We conducted a k-mer GWAS as it was shown previously that this approach can be more effective at detecting associations (cite Voichek paper). Furthermore, this approach is not based on a reference genome and can detect association with genetic regions which are not found in the reference genome. We followed the method described by Voichek et al (cite paper and GitHub) which briefly starts by counting 31bp-long k-mers in the sequencing reads without any need for alignment to a reference. The k-mers counts per samples are then combined, significance thresholds are estimated by random permutations, and a GWAS is performed for each phenotype using XXX from within the kmer-to-GWAS pipeline. We used the same samples and phenotypic values as for the SNP based GWAS.
For visualisation and identifying potential gene candidates, we used blast to map the k-mers to XXX fully assembled genomes of Z.tritici (cite sup table with accession numbers or papers). 


**Results**

In addition to the initial GWAS approach, which is based on short variants, we also used a k-mer approach. This approach, although less well-used, has several advantages: it has been shown in some cases to be more efficient in detecting associations and is not based on mapping to a reference genome. This second point allows for the detection of regions which are not found in the reference genome and is able to include in the analysis a complete vision of the genetic variants contained within the GWAS population. This is especially relevant in Z. tritici, a species known to have a large pan genome (cite Badet paper), including hyper-variable accessory chromosomes. 

Using this approach, we obtained a total of XXX k-mers associated with any of the phenotypes, including XXX that were found to be associated with at least two phenotypes. Amongst the XXX phenotypes, XXX has no associated k-mers at all. For the other phenotypes, the number of associated k-mers ranged from XXX (phenotype) to XXX (phenotype). 

To visualise the results and identify potential gene candidates, we mapped the significant k-mers to several fully assembled genomes of *Z. tritici*. We used a collection of XXX chromosome-level assemblies, including the reference IPO323, that covers the same geographical range as our samples, including isolates from Europe, the Middle-East, the Americas and Oceania, thus encompassing a larger part of the global genetic diversity. We defined significant loci when at least 2 significant k-mers were aligning less than XXX bp apart. We observed XXX significant loci across all phenotypes including XXX with k-mers in common between at least 2. 

Comparing the results of the k-mer-based GWAS to that of the SNPs-based GWAS, revealed limited overlap. However, the only locus significantly associated when using the strict Bonferoni threshold in the SNPs GWAS is recovered in the kmer GWAS. We generally expect the 



* General stats: 
    * number of kmers associated in general
    * number of kmers associated per phenotype & phenotypes without any kmers
    * number of associated regions and genes
* Comparison between phenotypes
    * number of kmers and regions found in common between phenotypes
    * phenotypes which have appearing/disappearing peaks
* Comparison with SNPs
    * Regions in common
    * Comparative power of detection