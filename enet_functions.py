import pandas as pd
import numpy as np
import gzip 

def get_maf_filtered_genotype(genotype_file_name,  maf, samples):
	gt_df = pd.read_csv(genotype_file_name, 'r', compression='gzip', header = 0, index_col = 0,delimiter='\t')
	effect_allele_freqs = gt_df.mean(axis=1)
	effect_allele_freqs = [ x / 2 for x in effect_allele_freqs ]
	effect_allele_boolean = pd.Series([ ((x >= maf) & (x <= (1 - maf))) for x in effect_allele_freqs ]).values
	gt_df = gt_df.loc[ effect_allele_boolean ]
	gt_df = gt_df.T
	return gt_df

def get_gene_expression(gene_expression_file_name, gene_annot):
	expr_df = pd.read_csv(gene_expression_file_name, compression='gzip', header = 0, index_col = 0, delimiter='\t')
	expr_df = expr_df.T
	inter = list(set(gene_annot['gene_id']).intersection(set(expr_df.columns)))
	print(len(inter))
	expr_df = expr_df.loc[:, inter ]
	return expr_df

##still working on this one
def get_cis_genotype(gt_df, snp_annot, coords, cis_window)
	snp_boolean = snp_annot[(snp_annot['pos'] >= (coords[0] - cis_window)) & (snp_annot['rsid'].notna()) & (snp_annot['pos'] <= (coords[1] + cis_window))]
	snp_info = snp_anno[snp_boolean
	if (len(snp_info.index) == 0)
		return(NaN)
	cis_gt <- gt_df %>% select(one_of(intersect(snp_info$varID, colnames(gt_df))))
	column_labels <- colnames(cis_gt)
	row_labels <- rownames(cis_gt)
	cis_gt <- matrix(as.matrix(cis_gt), ncol=ncol(cis_gt)) # R is such a bad language.
	colnames(cis_gt) <- column_labels
	rownames(cis_gt) <- row_labels
	cis_gt

