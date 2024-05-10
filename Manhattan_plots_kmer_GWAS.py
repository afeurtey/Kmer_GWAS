
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse


# Define the condition as an argument
# -----------------------------------
parser = argparse.ArgumentParser(description='Create Manhattan plot for kmer GWAS. Usage: python Manhattan_plots_kmer_GWAS.py eAUC25C')
parser.add_argument('condition', type=str, help='Condition to create the Manhattan plot')
A = parser.parse_args()

condition = A.condition #"eAUC22C"
print(f"\nCreating Manhattan plot for condition: {condition}")
all_gwas_path =  "/cluster/scratch/afeurtey/Kmer_GWAS/5_GWAS_results/"
dir_path = all_gwas_path + "GWAS_output_dir_" + condition + "/kmers/"


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

columns_names = ["rs", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue bitscore", "lala"]
pos_df = pd.read_csv(dir_path + 'output/phenotype_value.assoc_for_fetch.blast.tsv', sep='\t', 
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
df2['chromosome'] = df2['chromosome'].astype(int).astype('category')
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
plt.title('Manhattan plot')
plt.savefig(all_gwas_path + 'manhattan_plot_IPO323_' + condition + '.png')

# Safe file with only significant kmers
df2 = df2[df2['log10_pval'] > single_value]
df2.to_csv(all_gwas_path + 'significant_kmers_blast_IPO323_' + condition + '.tsv', sep='\t', index=False)


print(f" -- Found {df[df['log10_pval'] > single_value]['rs'].nunique()} significant kmers.")
print(f" -- Found {df2['rs'].nunique()} significant kmers aligning on the genome.")


print(f"Significant kmers saved to {all_gwas_path + 'significant_kmers_blast_IPO323_' + condition + '.tsv'}")
print(f"Manhattan plot saved to {all_gwas_path + 'manhattan_plot_IPO323_' + condition + '.png'}")