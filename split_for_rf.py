import pandas as pd
import numpy as np
import pickle

chrom = "22"

def get_gene_annotation (gene_anot_filepath, chrom, gene_types=["protein_coding"]):
     gene_anot = pd.read_csv(gene_anot_filepath, sep="\t")
     gene_anot = gene_anot[(gene_anot["chr"]==str(chrom)) &
                           (gene_anot["gene_type"].isin(gene_types))]
     return gene_anot

def get_gene_expression(gene_expression_file_name, gene_annot):
	expr_df = pd.read_csv(gene_expression_file_name, header = 0, index_col = 0, delimiter='\t')
	expr_df = expr_df.T
	inter = list(set(gene_annot['gene_id']).intersection(set(expr_df.columns)))
	#print(len(inter))
	expr_df = expr_df.loc[:, inter ]
	return expr_df

gene_expression_file = "Z:/data/mesa_models/meqtl_sorted_AFA_MESA_Epi_GEX_data_sidno_Nk-10.txt"
gene_annotation_file = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt"
geneannot = get_gene_annotation(gene_annotation_file, chrom)
expr_df = get_gene_expression(gene_expression_file, geneannot)
genes = list(expr_df.columns)

len(genes)

a = 0
for i in range(1,23,1):
	if i == 22:
		ls1 = genes[231:244]
		with open("Z:/data/paper_hyperopt/RF/chr22_"+str(i), "wb") as chunk:
			pickle.dump(ls1, chunk)
	else:
		ls1 = genes[11+a:22+a]
		with open ("Z:/data/paper_hyperopt/RF/chr22_"+str(i), "wb") as chunk:
			pickle.dump(ls1, chunk)
		a = a + 11
