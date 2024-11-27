# installing basic libraries
import datatable as dt
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold, StratifiedKFold, cross_val_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import random
import time
import sys

import pickle

###load in needed data

# Get the trait as a command line argument

trait = sys.argv[1]

print(trait)

# Load phenotype data

file_all = f"/usr/users/osatohanmwen/Simulation/Sim_hybrid_oligogenic/Geno_pheno_hybrid.csv"

Full_data = dt.fread(file_all)
Full_data = Full_data.to_pandas()

### for Normal

print("Start of Data Processing for Normal")

features = Full_data.iloc[:,1:30001].values
outcome = Full_data[trait].values


##split into training and testing test and reapet 20 times

# Set the number of times you want to run k-fold cross-validation
n_repeats = 20

# Initialize a dictionary to store the data splits
Data_OP = {}

for repeat in range (n_repeats):

    random.seed(repeat)  # Set a different random seed for each repeat
    shuffled_indices = random.sample(range(len(features)), len(features))
    shuffled_features = features[shuffled_indices]
    shuffled_outcome = outcome[shuffled_indices]

    ##standardise the data
    scaler=StandardScaler()

    #apply to both training and testing set
    features_scaled = scaler.fit_transform(shuffled_features)

    # ### using the eingen vector of SNPs instead of the markers
    #
    # pca = PCA(n_components=0.99)
    #
    # features_pca= pca.fit_transform(features_scaled)

    kf =KFold(n_splits=5, shuffle=True, random_state=repeat)

    # split data into training and test set, then run Bayessearch to determine the best parameters and run the Random forest for prediction.

    cnt = 1

    for train_index, test_index in kf.split(features_scaled, shuffled_outcome):

        xtrain = features_scaled[train_index]
        xtest = features_scaled[test_index]
        ytrain = shuffled_outcome[train_index]
        ytest = shuffled_outcome[test_index]

        ## Data for analysis
        key = f'Repeat_{repeat + 1}_Fold_{cnt}'
        print(key)
        Data_OP[key] = [xtrain,xtest,ytrain,ytest]

        cnt += 1


with open(f'{trait}_ML_data.pkl', "wb") as fp:
    pickle.dump(Data_OP, fp)

### for Combined Dominance

# print("Start of Data Processing for combined")
# 
# features = Full_data.iloc[:,1:30001].values
# outcome = Full_data[trait].values
# 
# ##split into training and testing test and reapet 20 times
# 
# # Set the number of times you want to run k-fold cross-validation
# n_repeats = 20
# 
# # Initialize a dictionary to store the data splits
# Data_OP = {}
# 
# for repeat in range (n_repeats):
# 
#     random.seed(repeat)  # Set a different random seed for each repeat
#     shuffled_indices = random.sample(range(len(features)), len(features))
#     shuffled_features = features[shuffled_indices]
#     shuffled_outcome = outcome[shuffled_indices]
# 
#     kf =KFold(n_splits=5, shuffle=True, random_state=repeat)
# 
#     # split data into training and test set, then run Bayessearch to determine the best parameters and run the Random forest for prediction.
# 
#     cnt = 1
# 
#     for train_index, test_index in kf.split(shuffled_features,shuffled_outcome):
# 
#         xtrain = shuffled_features[train_index]
#         xtest  = shuffled_features[test_index]
#         ytrain = shuffled_outcome[train_index]
#         ytest  = shuffled_outcome[test_index]
# 
#         snp_mean = np.zeros((xtrain.shape[1], 3))
# 
#         # Loop over each SNP
#         for snp in range(xtrain.shape[1]):
#             for genotype in range(3):
#                 # Calculate the mean of y_train for each genotype category
#                 snp_mean[snp, genotype] = np.nanmean(ytrain[xtrain[:, snp] == genotype])
#                 # Replace NaN values with 0
#                 snp_mean[np.isnan(snp_mean)] = 0.0
# 
#         # Create a condition to select rows where snpMean[:, 0] > snpMean[:, 2]
#         my_choose = snp_mean[:, 0] > snp_mean[:, 2]
# 
#         # Replace values in snpMean based on my_choose condition
#         snp_mean[my_choose, 0] = snp_mean[my_choose, 2]
#         snp_mean[my_choose, 2] = snp_mean[my_choose, 0]
# 
#         xtrain[:, my_choose] = np.abs(xtrain[:, my_choose] - 2)
# 
#         # Calculate d, which is the degree of dominance
#         d = (snp_mean[:, 1] - snp_mean[:, 0]) / (snp_mean[:, 2] - snp_mean[:, 0]) * 2
# 
#         # Replace missing values with 1
#         my_choose_na = np.isnan(d)
#         d[my_choose_na] = 1
# 
# 
#         # Set boundaries for d
#         my_choose1 = d >= 2
#         d[my_choose1] = 2
# 
#         my_choose2 = d <= 0
#         d[my_choose2] = 0
# 
#         # Weighted heterozygous genotypes
#         for snp in range(shuffled_features.shape[1]):
#             my_choose = shuffled_features[:, snp] == 1
#             shuffled_features[my_choose, snp] = d[snp]
# 
# 
#         # d now contains the degree of dominance, and SNP_data_G has been modified as described
# 
#         ##standardise the data
# 
#         scaler=StandardScaler()
# 
#         #apply to both training and testing set
#         features_scaled = scaler.fit_transform(shuffled_features)
# 
#         ### using the eingen vector of SNPs instead of the markers
#         # from sklearn.decomposition import PCA
#         #
#         # pca = PCA(n_components=0.99)
#         # features_pca=pca.fit_transform(features_scaled)
#         #
#         # print(features_pca.shape)
# 
#         xtrain_pca = features_scaled[train_index]
#         xtest_pca  = features_scaled[test_index]
# 
#         ## Data for analysis
#         key = f'Repeat_{repeat + 1}_Fold_{cnt}'
#         print(key)
#         Data_OP[key] = [ xtrain_pca,xtest_pca,ytrain,ytest]
# 
#         cnt += 1
# 
# 
# with open(f'{trait}_ML_domdata.pkl', "wb") as fp:
#     pickle.dump(Data_OP, fp)


