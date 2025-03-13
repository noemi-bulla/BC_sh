import matplotlib.pyplot as plt
import pandas as pd
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

fig, ax = plt.subplots(figsize=(20, 13)) 
ax.axis("tight")
ax.axis("off")

df_gsea_formatted = df_gsea.iloc[0:50].copy()  
table = ax.table(cellText=df_gsea_formatted.values, colLabels=df_gsea_formatted.columns, cellLoc="center", loc="center")


table.auto_set_font_size(False)
table.set_fontsize(10)
table.auto_set_column_width([0, 1, ])  


for (i, j), cell in table.get_celld().items():
    if i == 0: 
        cell.set_text_props(weight="bold", color="black")
        cell.set_facecolor("#ADD8E6") 
    else:
        cell.set_text_props(weight="bold")
        cell.set_facecolor("white") 


plt.savefig(os.path.join(path_results, "DEG_PAEP_gsea_table_.png"), dpi=300)