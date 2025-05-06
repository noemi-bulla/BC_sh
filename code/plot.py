import matplotlib.pyplot as plt
import pandas as pd
import textwrap
import os

path_results="/Users/ieo7295/Desktop/BC_sh/results/res_no_out"  
df_corr=os.path.join(path_results, "Degs_corr.xlsx")
df_chemor=os.path.join(path_results, "Degs_chemor.xlsx")
df_degs=os.path.join(path_results, "Degs_paepvsscr.xlsx")
df_gsea=os.path.join(path_results,"PAEPvsshSCR_gsea.xlsx")
df_d=pd.read_excel(df_degs).rename(columns={'Unnamed: 0':'Genes'}).sort_values("FDR")
df_c=pd.read_excel(df_corr).rename(columns={'Unnamed: 0':'Genes'}).sort_values("FDR")
df_ch=pd.read_excel(df_chemor).rename(columns={'Unnamed: 0':'Genes'}).sort_values("FDR")
df_gsea=pd.read_excel(df_gsea)
df_d=df_d.drop(columns={"F","logCPM"})
df_c=df_c.drop(columns={"F","logCPM"})
df_ch=df_ch.drop(columns={"F","logCPM"})
path_data="/Users/ieo7295/Desktop/BC_chemo_reproducibility/results/MDA"
df_gsea=pd.read_csv(os.path.join(path_data,"Gsea_promet_AC_vs_promet_NT.csv"))


def wrap_text(text, width):
    return '\n'.join(textwrap.wrap(str(text), width=width))

df_gsea['Lead_genes'] = df_gsea['Lead_genes'].apply(lambda x: wrap_text(x.replace(';', ', '), 40))
df_gsea['Term'] = df_gsea['Term'].apply(lambda x: wrap_text(x, 40))  # Optional: wrap 'Term' column too

# Estimate max lines per row based on wrapped 'Lead_genes'
row_line_counts = df_gsea['Lead_genes'].apply(lambda x: x.count('\n') + 1).tolist()

# Plot
fig_height = sum([0.1 * lines for lines in row_line_counts]) + 1  # total height
fig, ax = plt.subplots(figsize=(20, fig_height))
ax.axis("off")

table = ax.table(
    cellText=df_gsea.values,
    colLabels=df_gsea.columns,
    cellLoc='left',
    loc='center'
)

table.auto_set_font_size(False)
table.set_fontsize(10)
table.auto_set_column_width(col=list(range(len(df_gsea.columns))))

# Set header styling
for j in range(len(df_gsea.columns)):
    cell = table[(0, j)]
    cell.set_text_props(weight='bold', color='black')
    cell.set_facecolor('#ADD8E6')

# Apply uniform height to each row based on wrapped line count
for i in range(1, len(df_gsea) + 1):
    lines = row_line_counts[i - 1]
    row_height = 0.05 * lines  # tweak multiplier if needed
    for j in range(len(df_gsea.columns)):
        cell = table[(i, j)]
        cell.set_height(row_height)
        cell.set_text_props(va='center', ha='left', wrap=True)

plt.tight_layout()
plt.savefig(os.path.join(path_data,"DEG_gsea_prometAC_vs_prometNt.png"), dpi=300)