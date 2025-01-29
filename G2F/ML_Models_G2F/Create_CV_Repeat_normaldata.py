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

file_pheno = f"/usr/users/osatohanmwen/GP_Dom_G2Fdata/emmreml/{trait}_BLUES_All_final.csv"
Pheno_Blues_2018_2021 = dt.fread(file_pheno)
Pheno_Blues_2018_2021 = Pheno_Blues_2018_2021.to_pandas()


# Remove "Hybrid" from the Hybrid column
Pheno_Blues_2018_2021['Hybrid'] = Pheno_Blues_2018_2021['Hybrid'].str.replace("Hybrid", "", regex=False)

# Load filtered and recoded genotype data
file_genotype = "/usr/users/osatohanmwen/GP_Dom_G2Fdata/emmreml/Genotype_True1.csv"

SNP_data = dt.fread(file_genotype)
SNP_data = SNP_data.to_pandas()

#SNP_data.to_csv("/Users/bright/Downloads/Genotype_True1.csv", index=False)

# Filter genotype data based on Hybrids present in phenotype data
SNP_data = SNP_data[SNP_data['Hybrid'].isin(Pheno_Blues_2018_2021['Hybrid'])]

# Filter phenotype data based on Hybrids present in genotype data
Pheno_Blues_2018_2021 = Pheno_Blues_2018_2021[Pheno_Blues_2018_2021['Hybrid'].isin(SNP_data['Hybrid'])]

Full_data = SNP_data.merge(Pheno_Blues_2018_2021,how='inner', on="Hybrid")

### for Normal

print("Start of Data Processing for Normal")

features = Full_data.iloc[:,1:-1].values
outcome = Full_data.iloc[:,-1].values


##split into training and testing test and repeat 10 times

# Set the number of times you want to run k-fold cross-validation
n_repeats = 10

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
  
