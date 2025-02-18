import pandas as pd
import gseapy as gp
import gc
import pandas as pd 
from gseapy import prerank
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.pyplot as plt
from plotting_utils._plotting import *
matplotlib.use('macOSX')

class GSEAAnalysis:
    def __init__(self, pca_loadings_df, organism='human'):
        """
        Initialize the GSEA Analysis class with PCA loadings.
        """
        self.pca_loadings_df = pca_loadings_df  
        self.organism = organism
        self.GSEA = {}

    def compute_GSEA(self, covariate='PC2_loading', collection='GO_Biological_Process_2023'):
        """
        Perform GSEA using PCA loadings (e.g., PC1 or PC2) as covariate.
        """
        if covariate not in self.pca_loadings_df.columns:
            raise ValueError(f"{covariate} not found in the PCA loadings dataframe")

        self.pca_loadings_df['gene'] = self.pca_loadings_df['gene'].astype(str)

        self.pca_loadings_df = self.pca_loadings_df.dropna(subset=[covariate])

        ranked_gene_list = self.pca_loadings_df[['gene', covariate]].sort_values(by=covariate, ascending=False)
        ranked_gene_list.set_index('gene', inplace=True)

        results = prerank(
            rnk=ranked_gene_list,
            gene_sets=[collection],
            threads=4,
            min_size=15,
            max_size=500,
            permutation_num=200, 
            outdir=None, 
            seed=1234,
            verbose=True,
        )

        df = results.res2d.loc[:, 
            [ 'Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes' ]
        ].rename(columns={'FDR q-val' : 'Adjusted P-value'})
        
        idx = df.sort_values(by='Adjusted P-value').index
        filtered_df = df.loc[idx, :]

        filtered_df['Term'] = filtered_df['Term'].map(lambda x: str(x).split('__')[1] if isinstance(x, str) else str(x))
        filtered_df = filtered_df.set_index('Term')

        self.GSEA['original'] = filtered_df
        gc.collect()
        
        return filtered_df

pca_loadings_df = pd.read_csv("pca_loadings.csv")
gsea_analysis = GSEAAnalysis(pca_loadings_df)

#Run the GSEA
gsea_results_pc1 = gsea_analysis.compute_GSEA(covariate='PC1_loading', collection='GO_Biological_Process_2023')
gsea_results_pc2 = gsea_analysis.compute_GSEA(covariate='PC2_loading', collection='GO_Biological_Process_2023')

gsea_results_pc1.to_excel('GSEA_PC1.xlsx')
gsea_results_pc2.to_excel('GSEA_PC2.xlsx')


"""
Plot stem-plots of top 2 PCs GSEA-enriched pathways
"""
i=1 
pc_column = f'PC{i}_loading'  
loading_data = pca_loadings_df[['gene', pc_column]].sort_values(by=pc_column, ascending=False)
organism='human'
collection='GO_Biological_Process_2023'
gsea = GSEAAnalysis(loading_data, organism=organism)
gsea_results = gsea.compute_GSEA(covariate=pc_column, collection=collection)

   
fig, axs = plt.subplots(1, 2, figsize=(15, 8))
stem_plot(
        gsea_results.iloc[:, [0, 1, 3]].sort_values('NES', ascending=False).head(50),
        'NES', 
        ax=axs[0]
    )
format_ax(axs[0], title='GSEA', xlabel='NES')
  
stem_plot(
        loading_data[f'PC{i}_loading'].to_frame('es').sort_values('es', ascending=False).head(50),
        'es', 
        ax=axs[1]
    )
format_ax(axs[1], title=f'PC{i} loadings', xlabel='Loadings')
fig.suptitle(f'PC{i}, scaled layer, original representation')
fig.tight_layout()
fig.savefig(("GSEA_PC1.png"),dpi=300)
plt.show()



### GSEA for regulatory network genes ###

def network_gene_gsea(network_gene_3 ,collection='GO_Biological_Process_2023'):

     network_gene_3['gene'] = network_gene_3['gene'].astype(str)

     ranked_gene_list = network_gene_3.sort_values(by='logFC', ascending=False)
     ranked_gene_list = ranked_gene_list.set_index('gene')['logFC']

     results = prerank(
            rnk=ranked_gene_list,
            gene_sets=[collection],
            threads=4,
            min_size=1,
            max_size=500,
            permutation_num=200, 
            outdir=None, 
            seed=1234,
            verbose=True,
        )

     df = results.res2d.loc[:, 
            [ 'Term', 'ES', 'NES', 'FDR q-val', 'Lead_genes' ]
        ].rename(columns={'FDR q-val' : 'Adjusted P-value'})
        
     idx = df.sort_values(by='Adjusted P-value').index
     filtered_df = df.loc[idx, :]

     filtered_df['Term'] = filtered_df['Term'].map(lambda x: str(x).split('__')[1] if isinstance(x, str) else str(x))
     filtered_df = filtered_df.set_index('Term')
     gc.collect()
        
     return filtered_df


network_gene=pd.read_csv("PAEP1vsPAEP2_gsea.csv")
network_gene= network_gene.rename(columns={'Unnamed: 0':'gene'})

network_gene_2=pd.read_csv("PAEP1vsSCR_gsea.csv")
network_gene_2= network_gene_2.rename(columns={'Unnamed: 0':'gene'})

network_gene_3=pd.read_csv("PAEP2vsSCR_gsea.csv")
network_gene_3= network_gene_3.rename(columns={'Unnamed: 0':'gene'})

gsea = network_gene_gsea(network_gene_3,collection='GO_Biological_Process_2023')
gsea.to_excel("PAEP1vsPAEP2_gsea.xlsx", index=True)
gsea.to_excel("PAEP1vsSCR_gsea.xlsx", index=True)
gsea.to_excel("PAEP2vsSCR_gsea.xlsx", index=True)

fig, ax = plt.subplots(figsize=(10, 8))
stem_plot(
    gsea.iloc[:, [0, 1, 3]].sort_values('NES', ascending=False).head(50),
    'NES', 
    ax=ax 
)
format_ax(ax, title='GSEA NES', xlabel='NES')

plt.tight_layout()
plt.show()
fig.savefig(("GSEA_PAEP1vsPAEP2_hybrid.png"),dpi=300)


### gene enrichment plot ###
df_pc1=pd.read_excel('GSEA_PC1.xlsx')
df_pc2=pd.read_excel('GSEA_PC2.xlsx')


def plot_gsea_results(df_pc1, title='Gene Set Enrichment Analysis'):
    df = pd.read_excel(df_pc1, index_col=0)
    
    df = df.reindex(df['NES'].abs().sort_values(ascending=False).index).head(20)
    
    plt.figure(figsize=(10, 6))
    sns.barplot(y=df.index, x=df['NES'], palette='viridis', edgecolor='black')
    
    plt.axvline(x=0, color='black', linestyle='--', linewidth=1)
    plt.xlabel('Normalized Enrichment Score (NES)')
    plt.ylabel('Gene Sets')
    plt.title(title)
    
    plt.tight_layout()
    plt.show()

# Example usage
plot_gsea_results('GSEA_PC1.xlsx', title='GSEA Results for PC1 Loading')