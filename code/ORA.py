import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import gc
import os
import networkx as nx
from gseapy import enrichr, prerank

class ORAAnalysis:
    def __init__(self, df_deg, organism='human'):
        self.df_deg= df_deg
        self.organism = organism
        self.ORA_results = {}


    def compute_ORA(self, key='ORA_results', gene_column='Gene', 
        collection='GO_Biological_Process_2023', n_out=50):
        """
        Perform ORA (Over-Representation Analysis)
        """
        if gene_column not in self.df_deg.columns:
            raise ValueError(f"Column '{gene_column}' not found in dataframe")
        
        gene_list = self.df_deg[gene_column].dropna().tolist()
        
        # Ensure we have enough genes for ORA
        if len(gene_list) < 5:
            raise ValueError("ORA requires at least 5 genes for meaningful results.")

        results = enrichr(
            gene_list=gene_list,
            gene_sets=[collection],
            organism=self.organism, 
            outdir=None, 
        ).results

        df = results.loc[:, 
            [ 'Term', 'Overlap', 'Adjusted P-value', 'Genes' ]
        ]

        # Select top `n_out` pathways
        filtered_df = df.nsmallest(n_out, 'Adjusted P-value')
        filtered_df = filtered_df.set_index('Term')

        # Store results
        self.ORA_results[key] = filtered_df

        # Free memory
        gc.collect()

        return filtered_df 

file_path="/Users/ieo7295/Desktop/BC_sh/results/pca_plot"
file=os.path.join(file_path,"Degs_paepvsscr.xlsx")
df_deg=pd.read_excel(file)
df_deg.rename(columns={'Unnamed: 0':'Gene'}, inplace=True)
ora_analysis= ORAAnalysis(df_deg)

ora_results = ora_analysis.compute_ORA(gene_column='Gene')
# ora_results_2=ora_analysis.computeORA(covariate='PC2_loading', collection='GO_Biological_Process_2023')
ora_results.to_excel(os.path.join(file_path, "ORA_results_deg.xlsx"))


def plot_ORA_bar(ora_results, title="Top Enriched Pathways", top_n=50):
    """
    Plots a bar chart of the top enriched pathways based on Adjusted P-value.

    Parameters:
    - ora_results: DataFrame with ORA results (must contain 'Adjusted P-value')
    - title: Title of the plot
    - top_n: Number of top pathways to show
    """
    # Select top pathways
    df_plot = ora_results.nsmallest(top_n, 'Adjusted P-value')

    plt.figure(figsize=(13, 9))
    sns.barplot(
        data=df_plot,
        y=df_plot.index,
        x='Adjusted P-value',
        palette='viridis'
    )
    
    plt.xlabel('Adjusted P-value')
    plt.ylabel('Pathway')
    plt.xlim(left=0.07)
    plt.title(title) 
    plt.tight_layout()
    plt.show()

# Plot the ORA results
plot_ORA_bar(ora_results)



def plot_ORA_dot(ora_results, title="Enrichment Dot Plot", top_n=50):
    """
    Plots a dot plot of enriched pathways showing adjusted P-value and gene overlap.

    Parameters:
    - ora_results: DataFrame with ORA results (must contain 'Adjusted P-value' & 'Overlap')
    - title: Title of the plot
    - top_n: Number of top pathways to show
    """
    df_plot = ora_results.nsmallest(top_n, 'Adjusted P-value')
    
    plt.figure(figsize=(14, 9))
    scatter = plt.scatter(
        x=df_plot['Adjusted P-value'],
        y=df_plot.index,
        s=df_plot['Overlap'].apply(lambda x: int(x.split('/')[0])) * 50,  # Bubble size = overlapping genes
        c=np.log10(df_plot['Adjusted P-value']),
        cmap='coolwarm',
        edgecolors='black'
    )
    
    plt.xscale('log')
    plt.xlabel('Adjusted P-value (log scale)')
    plt.ylabel('Pathway')
    plt.colorbar(scatter, label="Log10(Adjusted P-value)")
    plt.title(title)
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.show()
    plt.savefig(os.path.join(file_path, "ORA_deg_paepvsscr.png"),dpi=300)
# Plot dot plot for ORA results
plot_ORA_dot(ora_results)





def plot_enrichment_map(ora_results, top_n=50):
    """
    Creates a network graph of enriched pathways.
    """
    df_plot = ora_results.nsmallest(top_n, 'Adjusted P-value')

    G = nx.Graph()
    
    for _, row in df_plot.iterrows():
        pathway = row.name  # Pathway name
        genes = row['Genes'].split(";")  # Split gene names
        
        G.add_node(pathway, type='pathway')
        
        for gene in genes:
            G.add_node(gene, type='gene')
            G.add_edge(pathway, gene)
    
    # Draw the network
    plt.figure(figsize=(12, 8))
    pos = nx.spring_layout(G, seed=42)
    node_colors = ['red' if G.nodes[n]['type'] == 'pathway' else 'blue' for n in G.nodes()]
    
    nx.draw(
        G, pos, with_labels=True, node_color=node_colors, edge_color='gray',
        node_size=500, font_size=8, alpha=0.8
    )
    plt.title("Enrichment Map")
    plt.show()

# Generate enrichment map
plot_enrichment_map(ora_results)
