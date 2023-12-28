import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

def import_dataframe(p_xlsx_1: str, p_xlsx_2: str, targets: list, filter=True):
    if p_xlsx_1.endswith(".xlsx") != True or p_xlsx_2.endswith(".xlsx") != True:
        print("One of the files is not an excel file LOL")
        return RuntimeError()
    
    xlsx_1 = pd.read_excel(p_xlsx_1)
    xlsx_2 = pd.read_excel(p_xlsx_2)

    if p_xlsx_1.endswith("Destructivum_Expression.xlsx"):
        xlsx_1["Species"] = ["Destructivum"] *len(xlsx_1)
        xlsx_2["Species"] = ["Higginsianum"] * len(xlsx_2)
        xlsx_1["GeneID.1"] = xlsx_1["Gene ID"]
        xlsx_2["GeneID.1"] = xlsx_2["Gene ID"].astype(str) + ".1"

    elif p_xlsx_1.endswith("Higginsianum_Expression.xlsx"):
        xlsx_1["GeneID.1"] = xlsx_1["Gene ID"].astype(str) + ".1"
        xlsx_2["GeneID.1"] = xlsx_2["Gene ID"]
        xlsx_1["Species"] = ["Higginsianum"] *len(xlsx_1)
        xlsx_2["Species"] = ["Destructivum"] * len(xlsx_2)


    target_rows_1 = []
    target_rows_2 = []
    xlsx_1 = xlsx_1.dropna()
    xlsx_2 = xlsx_2.dropna()

    #print(f"Dataframe 1 Spalten:{xlsx_1.columns}")
    #print(f"Dataframe 2 Spalten:{xlsx_2.columns}")

    if xlsx_1["Species"].iloc[0] == "Destructivum":
        for i in targets:
            target_rows_1.append(xlsx_1[xlsx_1["GeneID.1"] == i])
            target_rows_2.append(xlsx_2[xlsx_2["GeneID.1"] == i])

    elif xlsx_1["Species"].iloc[0] == "Higginsianum":
        for i in targets:
            target_rows_1.append(xlsx_1[xlsx_1["GeneID.1"] == i])
            target_rows_2.append(xlsx_2[xlsx_2["GeneID.1"] == i])

    df_1 = pd.concat(target_rows_1)
    df_2 = pd.concat(target_rows_2)

    return df_1, df_2

def effectors_from_excel(e_path):
    df = pd.read_excel(e_path, skiprows=1, sheet_name="Predicted effectors")
    df = df.dropna()
    return df["Chig"].to_list()

def prefilter(df1: pd.DataFrame,df2: pd.DataFrame, min_TPM):
    if df1["Species"].iloc[0] == "Destructivum":
        df2 = df2[(df2["BP"] > min_TPM) | (df2["NP"] > min_TPM)]
        df1 = df1[(df1["48h"] > min_TPM) | (df1["72h"] > min_TPM)]

        # First, find the intersection of 'GeneID.1' values in both DataFrames
        common_genes = df1['GeneID.1'][df1['GeneID.1'].isin(df2['GeneID.1'])].unique()

        # Filter both DataFrames to only include rows with common 'GeneID.1' values
        df1_common = df1[df1['GeneID.1'].isin(common_genes)].copy()
        df2_common = df2[df2['GeneID.1'].isin(common_genes)].copy()

        # Sort both DataFrames by 'GeneID.1' to ensure they are in the same order
        df1_common.sort_values('GeneID.1', inplace=True)
        df2_common.sort_values('GeneID.1', inplace=True)
        
        dest_ar = df1_common[['48h', '72h']].to_numpy()
        higg_ar = df2_common[["BP", "NP"]].to_numpy()

    else:
        df1 = df1[(df1["BP"] > min_TPM) | (df1["NP"] > min_TPM)]
        df2 = df2[(df2["48h"] > min_TPM) | (df2["72h"] > min_TPM)]

        # First, find the intersection of 'GeneID.1' values in both DataFrames
        common_genes = df1['GeneID.1'][df1['GeneID.1'].isin(df2['GeneID.1'])].unique()

        # Filter both DataFrames to only include rows with common 'GeneID.1' values
        df1_common = df1[df1['GeneID.1'].isin(common_genes)].copy()
        df2_common = df2[df2['GeneID.1'].isin(common_genes)].copy()

        # Sort both DataFrames by 'GeneID.1' to ensure they are in the same order
        df1_common.sort_values('GeneID.1', inplace=True)
        df2_common.sort_values('GeneID.1', inplace=True)

        dest_ar = df2_common[["48h", "72h"]].to_numpy()
        higg_ar = df1_common[["BP", "NP"]].to_numpy()
    return dest_ar, higg_ar, df1_common['GeneID.1'].tolist()

def make_covariance(dest_ar: np.array, higg_ar: np.array, top_n_print = None, export_xlsx_cols=None, export_path=None):
    corr = np.corrcoef(dest_ar, higg_ar, rowvar=True)[0:dest_ar.shape[0], dest_ar.shape[0]:]  
    
    if top_n_print:
        top_k = 10  # for example, top 10 correlations

# Flatten the matrix and sort by absolute value
        flat_corr = corr.flatten()
        sorted_indices = np.argsort(np.abs(flat_corr))[::-1]  # Descending sort

        # Get the top k indices
        top_k_indices = sorted_indices[:top_k]

        # Map indices to matrix coordinates and gene names
        for idx in top_k_indices:
            row_idx, col_idx = np.unravel_index(idx, corr.shape)
            gene1 = export_xlsx_cols[row_idx]
            gene2 = export_xlsx_cols[col_idx]
            corr_value = corr[row_idx, col_idx]
            print(f"{gene1} - {gene2}: {corr_value}")

    if export_xlsx_cols:
        if not export_path:
            print("Specify Output Path por favor")
        corr_df = pd.DataFrame(corr, index=export_xlsx_cols, columns=export_xlsx_cols)
        
        corr_df_diag = pd.DataFrame([np.diagonal(corr)], columns=export_xlsx_cols)
        corr_df.to_excel(export_path, index_label="Columns are C_dest, Rows are C_higginsianum", sheet_name="Heatmap")
        corr_df_diag.to_excel(export_path, sheet_name="Correlations between each identical Gene (Diagonal of Heatmap)")
    
    return corr


    
    return 

def main():
    e_list = effectors_from_excel("/Users/floriangrun/Desktop/Bachelorarbeit_Koch/Expressionsmuster/Supplementary_Table_6_AllGenes_vf.xlsx")
    df1, df2 = import_dataframe("/Users/floriangrun/Desktop/Bachelorarbeit_Koch/Expressionsmuster/Destructivum_Expression.xlsx", 
                                "/Users/floriangrun/Desktop/Bachelorarbeit_Koch/Expressionsmuster/Higginsianum_Expression.xlsx",
                                e_list)
                                #["CH63R_14384", "CH63R_14389", "CH63R_14516"])
    d_ar, h_ar, geneIDs = prefilter(df1, df2, 20)
    corr = make_covariance(d_ar, h_ar, 10, geneIDs,"/Users/floriangrun/Desktop/Bachelorarbeit_Koch/Expressionsmuster/output.xlsx")
    plt.imshow(corr)
    plt.colorbar()
    plt.show()

    return

if __name__ == "__main__":
    main()
