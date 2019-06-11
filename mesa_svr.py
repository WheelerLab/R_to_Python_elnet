import numpy as np
from sklearn.svm import SVR
#from sklearn.model_selection import train_test_split
import pandas as pd
from sklearn.model_selection import cross_val_score
from statistics import mean

#MESA
cis_gt = pd.read_csv("/Users/okoro/OneDrive/Desktop/svr/svr_cis_gt_chr6_HLA-DRB1.csv")
adj_exp = pd.read_csv("/Users/okoro/OneDrive/Desktop/svr/svr_adj_expression_chr6_HLA-DRB1.csv")


#convert the dataframe to numpy array
cis_gt = cis_gt.values
adj_exp = adj_exp.values

#split data into train and test

#X_train, X_test, y_train, y_test = train_test_split(cis_gt, adj_exp, test_size=0.1, random_state=42)

#print()
print("Random Forest")
#run a random forest regression
from sklearn.ensemble import RandomForestRegressor
#CV
rf = RandomForestRegressor(max_depth=None, random_state=1234, n_estimators=100)

t0 = time.time() # or time.process_time()
print(mean(cross_val_score(rf, cis_gt, adj_exp.ravel(), cv=5)))
t1 = time.time() #or time.process_time()
total_time = t1 - t0

#fit the rf with all data
rf.fit(cis_gt, adj_exp.ravel())

std = np.std(rf.feature_importances_)
importances = rf.feature_importances_
indices = np.argsort(importances)[::-1]
print("Feature ranking:")
for f in range(cis_gt.shape[1]):
    print("%d. feature %d (%f)" % (f + 1, indices[f], importances[indices[f]]))

#select only the most important features (SNPs)
imp_indices = indices[0:500] #first 500 most imprtant indices
imp_cis_gt = cis_gt[:,imp_indices]

#Try simple linear regression with the feature importance
from sklearn import linear_model
reg = linear_model.LinearRegression() #it performed very poor
mean(cross_val_score(reg, cis_gt, adj_exp, cv=5))
mean(cross_val_score(reg, imp_cis_gt, adj_exp, cv=5)

from sklearn.svm import SVR
svrl = SVR(kernel="linear", gamma="auto")
mean(cross_val_score(svrl, cis_gt, adj_exp.ravel(), cv=5))
mean(cross_val_score(svrl, imp_cis_gt, adj_exp.ravel(), cv=5))
svr = SVR(kernel="rbf", gamma="auto")
mean(cross_val_score(svr, cis_gt, adj_exp.ravel(), cv=5))
mean(cross_val_score(svr, imp_cis_gt, adj_exp.ravel(), cv=5))

from sklearn.neighbors import KNeighborsRegressor
knn = KNeighborsRegressor(n_neighbors=10, weights = "distance")
mean(cross_val_score(knn, cis_gt, adj_exp, cv=5))
mean(cross_val_score(knn, imp_cis_gt, adj_exp, cv=5))

from sklearn.linear_model import ElasticNet
elnet = ElasticNet(alpha=0.1, random_state=1234)
mean(cross_val_score(elnet, cis_gt, adj_exp, cv=5))
mean(cross_val_score(elnet, imp_cis_gt, adj_exp, cv=5))
#KNN performed the best with feature importance R2 0.89
#followed by elasticNet R2 0.78
#svm performd at R2 0.53

imp_snps = cis_gt.columns[imp_indices]
pruned_cis_gt = cis_gt[imp_snps]
pruned_cis_gt = pruned_cis_gt.values

#this is how to index the pandas dataframe to get the snp names of the
# important features after running random forest feature importance
# this will be done so as to get the snp names, and then look for those snp
# names in the test dataframe. then take the overlaps and create a new test
# dataframe. Then fit the ML algorithm with the old important features and then
# predict expression in the new pruned test set

#to get column names of pandas dataframe
# for col in cis_gt.columns:
#		print(col)

# to access specific column name by col index e.g cis_gt.columns[1]
# to access list of col names by col indices e.g cis_gt.columns[indices]
# to the get the important colnames as selected by RF e.g cis_gt.columns[imp_indices]
# to finally create or select out only important features columns
# is df[colname] or df[list of colnames] e.g cis_gt[cis_gt.columns[imp_indices]]
# store the pruned cis_dataframe in a new variable
# e.g   pruned_cis = cis_gt[cis_gt.columns[imp_indices]]
#then fit your ML with the pruned cis
#note however, that your measured expression should be kept intact
#use these same steps to select out the important snps from the test data
#then use the fitted ML model to predict expression on the overlap snps of test
#data. then do the correlation of the

#make sure that the imp_snps are also present in the test snps
#therefore take the intersect of the two snps lists
# example function to take intersection of two lists

def snps_intersect(list1, list2):
     return list(set(list1) & set(list2))

def snps_intersect(list1, list2):
     return set(list1).intersection(list2)

#take snps intersect
snp_intersect = snps_intersect(imp_snps, yri_snps)

pruned_test = yri_cis[snp_intersect]
pruned_test = pruned_test.values
pruned_train = afa_cis[snp_intersect]
pruned_train = pruned_train.values

mean(cross_val_score(rf, pruned_train, afaadj_np.ravel(), cv=5))
rf.fit(pruned_train, afaadj_np.ravel())
rf.score(pruned_test, yriadj_np.ravel())

ypred = rf.predict(pruned_test)
calc_corr(yriadj_np, ypred)

#turn all these code and data file into a wrapper, and or a docker file that can be installed and run by users... with database capability
#even better, without having to save the data, save the fitted model objects for each algorithm into a file, and reload from there.
#this object savings can be done with python pickle module. import pickle


"""
import math

def calc_R2 (y, y_pred):
    tss = 0
    rss = 0
    for i in range(len(y)):
        tss = tss + (y[i])**2
        rss = rss + (((y[i]) - (y_pred[i]))**2)
    tss = float(tss)
    rss = float(rss)
    r2 = 1 - (rss/tss)
    return r2


def calc_corr (y, y_pred):
    num = 0
    dem1 = 0
    dem2 = 0
    for i in range(len(y)):
        num = num + ((y[i]) * (y_pred[i]))
        dem1 = dem1 + (y[i])**2
        dem2 = dem2 + (y_pred[i])**2
    num = float(num)
    dem1 = math.sqrt(float(dem1))
    dem2 = math.sqrt(float(dem2))
    rho = num/(dem1*dem2)
    return rho
"""


get_filtered_snp_annot <- function(snp_annot_file_name) {
  snp_annot <- read.table(snp_annot_file_name, header = T, stringsAsFactors = F) %>%
    filter(!((refAllele == 'A' & effectAllele == 'T') |
               (refAllele == 'T' & effectAllele == 'A') |
               (refAllele == 'C' & effectAllele == 'G') |
               (refAllele == 'G' & effectAllele == 'C')) &
             !(is.na(rsid))) %>%
    distinct(varID, .keep_all = TRUE)
  snp_annot
}

snpfilepath = "/Users/okoro/OneDrive/Desktop/svr/AFA_1_annot.txt"
snpfile = "/Users/okoro/OneDrive/Desktop/svr/YRI_annot.chr1.txt"
snpfile = "/Users/okoro/OneDrive/Desktop/svr/yri_dummy_anot.txt"
     
def get_filtered_snp_annot (snpfilepath):
     snpanot = pd.read_csv(snpfilepath, sep="\t")
     #snpanot = snpanot[((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="C")) | ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="A")) | ((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="G")) | ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="A")) | ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="G")) | ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="T")) | ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="C")) | ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="T"))]
     snpanot = snpanot[(((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="C")) | ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="A")) | ((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="G")) | ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="A")) | ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="G")) | ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="T")) | ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="C")) | ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="T"))) & (snpanot["rsid"].notna())]
     snpanot = snpanot.drop_duplicates(["varID"])
     return snpanot
     

get_gene_annotation <- function(gene_annot_file_name, chrom, gene_types=c('protein_coding')){
  gene_df <- read.table(gene_annot_file_name, header = TRUE, stringsAsFactors = FALSE) %>%
    filter((chr == chrom) & gene_type %in% gene_types)
  gene_df
}

geneanotfile = "/Users/okoro/OneDrive/Desktop/svr/gencode.v18.annotation.parsed.txt"

def get_gene_annotation (gene_anot_filepath, chrom, gene_types=["protein_coding"]):
     gene_anot = pd.read_csv(gene_anot_filepath, sep="\t")
     gene_anot = gene_anot[(gene_anot["chr"]==str(chrom)) & (gene_anot["gene_type"].isin(gene_types))]
     return gene_anot


get_gene_type <- function(gene_annot, gene) {
  filter(gene_annot, gene_id == gene)$gene_type
}


def get_gene_type (gene_anot, gene):
     gene_type = gene_anot[gene_anot["gene_id"]==gene]
     gene_type = gene_type.iloc[0,5]
     return gene_type


get_gene_coords <- function(gene_annot, gene) {
  row <- gene_annot[which(gene_annot$gene_id == gene),]
  c(row$start, row$end)

def get_gene_coords (gene_anot, gene):
     gene_type = gene_anot[gene_anot["gene_id"]==gene]
     gene_coord = [gene_type.iloc[0,3], gene_type.iloc[0,4]]
     return gene_coord

get_covariates <- function(covariate_file_name, samples) {
  cov_df <- read.table(covariate_file_name, header = TRUE, stringsAsFactors = FALSE, row.names = 2)
  cov_df <- cov_df[,3:5] %>% as.data.frame()
  cov_df
}#cov_df <- cov_df[,-1:-2] %>% as.data.frame()


cov_file = "/Users/okoro/OneDrive/desktop/svr/AFA_3_PCs.txt"
yri_covfile = "/Users/okoro/OneDrive/Desktop/svr/YRI_final_pcs.txt"
#cov file is space separated " ".
# the column names contain IID which is sample ID
def get_covariates (cov_filepath):
      cov = pd.read_csv(cov_filepath, sep=" ")
      cov = cov.set_index("IID") #make IID to be the row names
      cov.index.names = [None] # remove the iid name from the row
      pc = ["PC1", "PC2", "PC3"] #a list of the PCs to retain
      cov = cov[pc]
      return cov
  
gex = "/Users/okoro/OneDrive/Desktop/svr/meqtl_sorted_AFA_MESA_Epi_GEX_data_sidno_Nk-10.txt"
gex_yri = "/Users/okoro/OneDrive/Desktop/svr/YRI_expression_ens.txt"
  
def get_gene_expression(gene_expression_file_name, gene_annot):
	expr_df = pd.read_csv(gene_expression_file_name, compression='gzip', header = 0, index_col = 0, delimiter='\t')
	expr_df = expr_df.T
	inter = list(set(gene_annot['gene_id']).intersection(set(expr_df.columns)))
	print(len(inter))
	expr_df = expr_df.loc[:, inter ]
	return expr_df


adjust_for_covariates <- function(expression_vec, cov_df) {
  combined_df <- cbind(expression_vec, cov_df)
  expr_resid <- summary(lm(expression_vec ~ ., data=combined_df))$residuals
  expr_resid <- scale(expr_resid, center = TRUE, scale = TRUE)
  expr_resid
}

import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import scale
from pandas import DataFrame

#newdf = DataFrame(scale(df), index=df.index, columns=df.columns) # if we want
#to retain the scale array to dataframe with rownames and colnames  
# yri_exp["ENSG00000116791.9"] # test gene expression vector

def adjust_for_covariates (expr_vec, cov_df):
      reg = LinearRegression().fit(cov_df, expr_vec)
      ypred = reg.predict(cov_df)
      residuals = expr_vec - ypred
      residuals = scale(residuals)
      return residuals

  
  
  
  
