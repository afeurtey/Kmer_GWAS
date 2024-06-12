
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse


# Define the condition as an argument
# -----------------------------------
parser = argparse.ArgumentParser(description='Create Manhattan plot for kmer GWAS. Usage: python Manhattan_plots_kmer_GWAS.py eAUC25C')
parser.add_argument('condition', type=str, help='Condition to create the Manhattan plot')
parser.add_argument('genome_suffix', type=str, default = "IPO323")
parser.add_argument('GWAS_dir_path', type=str, default = "Silvias_GWAS_output_dir_")
A = parser.parse_args()

condition = A.condition #"eAUC22C"
print(f"\nCreating Manhattan plot for condition: {condition}")
all_gwas_path =  A.GWAS_dir_path + condition + "/"
dir_path = A.GWAS_dir_path + condition + "/kmers/"

#Define file names
blast_results = dir_path + 'output/phenotype_value.assoc_for_fetch.blast_' + A.genome_suffix + '.tsv'
all_kmers_merged_table = all_gwas_path + 'all_kmers_blast_' + A.genome_suffix + '_' + condition + '.tsv'
signif_kmers_merged_table = all_gwas_path + 'significant_kmers_blast_' + A.genome_suffix + '_' + condition + '.tsv'
manhattan_path = all_gwas_path + 'manhattan_plot_' + A.genome_suffix + '_' + condition + '.png'


# Read the association file for best 10000 kmers
df = pd.read_csv(dir_path + 'output/phenotype_value.assoc.txt.gz', sep='\t')
df['log10_pval'] = -np.log10(df['p_lrt'])

# Significance threshold
single_value_df = pd.read_csv(dir_path + 'threshold_10per', nrows=1, header=None, sep='\t')
single_value = single_value_df.iloc[0, 0]
print(f"Significance threshold: {single_value}")

# Plot histogram
plt.hist(df['log10_pval'], bins=50, edgecolor='black')

# Add vertical line
plt.axvline(single_value, color='r', linestyle='dashed', linewidth=2)
plt.title('Histogram of log10_pval')
plt.xlabel('log10_pval')
plt.ylabel('Frequency')
plt.show()

plt.savefig(all_gwas_path + 'histogram_pval_' + condition + '.png')


# Read the results from the blast search
columns_names = ["rs", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue bitscore", "lala"]
pos_df = pd.read_csv(blast_results, sep='\t', 
                     header=None,  index_col = None, names = columns_names)

pos_df = pos_df[["rs", "sseqid", "sstart", "send"]]
pos_df["position"] = (pos_df["sstart"] + pos_df["send"]) / 2
pos_df["position"] = pos_df["position"].astype(int)
pos_df  = pos_df[["rs", "sseqid", "position"]]
pos_df



# Merging information from the association and the blast results
# -------------------------------------------------------------
df2 = pd.merge(pos_df, df, on='rs', how='inner')
df2 = df2[df2['sseqid'] != 'mt']
# Create a new column that is the rank of the order of the chromosomes
df2['chromosome'] = df2['sseqid'].astype('category')
try :
    df2['chromosome'] = df2['chromosome'].astype(int).astype('category')
except:
    pass
df2['chr_index'] = df2['chromosome'].cat.codes

# Calculate cumulative position along the genome
# Sort the DataFrame by chromosome and position
df2 = df2.sort_values(['chromosome', 'position'])

# Calculate the difference between the current and previous position
df2['position_diff'] = df2.groupby('chromosome')['position'].diff()

# Replace NaN values with the original position
df2['position_diff'] = df2['position_diff'].fillna(df2['position'])

df2['cumulative_position'] = df2['position_diff'].transform(lambda x: x.cumsum())
df2['global_position'] = df2['cumulative_position'] + df2['chr_index'] * 1e7  # add gaps between chromosomes



# Create the Manhattan plot
# -------------------------

colors = ['skyblue', 'yellow']
plt.figure(figsize=(10, 5))
for chrom, data in df2.groupby('chr_index'):
    plt.scatter(data['cumulative_position'], data['log10_pval'], color=colors[chrom % len(colors)])
plt.axhline(single_value, color='grey', linestyle='dashed', linewidth=2)
plt.ylabel('-Log10(p-value)')
plt.xlabel('Position along the genome')
plt.title('Manhattan plot for '+ condition + " on reference " + A.genome_suffix)
plt.savefig(manhattan_path)

# Safe file withall and only significant kmers
df2.to_csv(all_kmers_merged_table, sep='\t', index=False)
print(f" -- Found {df[df['log10_pval'] > single_value]['rs'].nunique()} significant kmers.")
print(f" -- Found {df2['rs'].nunique()} kmers aligning on the genome.")

df2 = df2[df2['log10_pval'] > single_value]
df2.to_csv(signif_kmers_merged_table, sep='\t', index=False)
print(f" -- Found {df2['rs'].nunique()} significant kmers aligning on the genome.")

print(f"All kmers table saved to {all_kmers_merged_table}")
print(f"Significant kmers table saved to {signif_kmers_merged_table}")
print(f"Manhattan plot saved to {manhattan_path}")