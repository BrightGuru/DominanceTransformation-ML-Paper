import sys
import pandas as pd
from skopt import BayesSearchCV
import sklearn.metrics as metrics
from scipy.stats import pearsonr
import numpy as np
import time
import xgboost as xgb

i = sys.argv[2]

trait = sys.argv[1]

print(i)
print(trait)

t = time.time()

file_pheno = f"/usr/users/osatohanmwen/Simulation/ML_data/{trait}_ML_domdata.pkl"

XGB_data = pd.read_pickle(file_pheno)

data = XGB_data[str(i)]

xtrain = data[0]
xtest  = data[1]
ytrain = data[2]
ytest  = data[3]

print(xtrain.shape)
print(xtest .shape)
print(ytrain.shape)
print(ytest.shape)

###initialise XBM
XGBr_clf= xgb.XGBRegressor(random_state = 99)

# Bayesian optimization using an iterative Gaussian process

params_lg = {
  'n_estimators': (3000,8000, "log-uniform"), # No of trees# 220,520,620,
  'learning_rate' : (0.005, 0.3,"log-uniform"),

  'max_depth':(2,15,"log-uniform") , # maximum depth to explore (meaning the longest path between the root node and the leaf node.'max_depth': [10,15],
  'min_child_weight' : (1,20, "log-uniform"),
  'subsample': (0.1, 1),
  'colsample_bytree' : (0.1, 1),

  'gamma' : (0.01, 0.5, "log-uniform"),
  'reg_alpha' : (10, 200, "log-uniform"),
  'reg_lambda' : (1, 10, "log-uniform"),

  'n_jobs' : [-1] ,
  'seed' : [99]
  }

Bayes_search = BayesSearchCV(
    XGBr_clf,params_lg,n_iter=25,
    n_jobs = -1,cv=3,scoring='neg_mean_squared_error',
    random_state=99,
    verbose=0)

print(Bayes_search.total_iterations)

Bayes_search.fit(xtrain,ytrain)

print(Bayes_search.best_estimator_)

#update XGBM parameters usingb the best estimator which was obtained by using BayesSearchCV.

XGBr_clf_best = Bayes_search.best_estimator_

XGBr_clf_best.fit(xtrain, ytrain)

y_pred_grid = XGBr_clf_best.predict(xtest)

##observed and predicted dataframe
data_pre_observed = pd.DataFrame(
    {'Observed' : ytest,
        'Predicted': y_pred_grid
        }
)

data_accuracy = pd.DataFrame(
        {'Fold' : i,
        'r2_score' : metrics.r2_score(ytest, y_pred_grid),
         'pearsonr' : pearsonr(ytest, y_pred_grid)[0],
         'rmse'     : np.sqrt(metrics.mean_squared_error(ytest, y_pred_grid))
         }, index=[0]
    )

data_pre_observed.to_csv(str(trait) +'XBM_dom_obsv_pred.csv' + str(i) + '.csv',index=False)

data_accuracy.to_csv(str(trait) +'XBM_dom_acc.csv' + str(i) + '.csv' ,index=False)

print(i)
print(time.time() - t)
