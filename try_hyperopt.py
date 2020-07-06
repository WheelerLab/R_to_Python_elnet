import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import ElasticNet
from sklearn.svm import SVR
from sklearn.metrics import make_scorer
from sklearn.metrics import r2_score
 
r2 = make_scorer(r2_score, greater_is_better=True)
from hyperopt import fmin, tpe, hp, SparkTrials, Trials, STATUS_OK

trials  = Trials()

x = pd.read_csv("Z:/svr_cis_gt_chr1.csv", sep=",")
y = pd.read_csv("Z:/svr_adj_expression_chr1.csv", sep=",")

x_train = x.values

y_train = y.values
en = ElasticNet()
cv = cross_val_score(en, x_train, y_train, scoring=r2, cv=5)
cv.mean()

svr = SVR(kernel="linear", gamma="scale")
cross_val_score(svr, x_train, y_train.ravel(), scoring=r2, cv=5).mean()

rf = RandomForestRegressor(random_state=1234)
knn = KNeighborsRegressor()

#try hyperopt

#1 Define an objective function

def objective(params):
    regressor_type = params["type"]
    del params["type"]
    if regressor_type == "elastic_net":
        regressor = ElasticNet(**params)
    elif regressor_type == "svm":
        regressor = SVR(gamma="scale", **params)
    else:
        return 0
    r2_mean = cross_val_score(regressor, x_train, y_train.ravel(), scoring=r2, cv=5).mean()

    return {"loss": -r2_mean, "status": STATUS_OK}


def objective(params):
    regressor_type = params["type"]
    del params["type"]
    if regressor_type == "elastic_net":
        regressor = ElasticNet(**params)
    elif regressor_type == "rf":
        regressor = RandomForestRegressor(random_state=1234, **params)
    elif regressor_type == "svm":
        regressor = SVR(gamma="scale", **params)
    elif regressor_type == "knn":
        regressor = KNeighborsRegressor(**params)
    else:
        return 0
    r2_mean = cross_val_score(regressor, x_train, y_train.ravel(), scoring=r2, cv=5).mean()

    return {"loss": -r2_mean, "status": STATUS_OK}


#2 Define search space

search_space = hp.choice("regressor_type", [
    {
        "type": "elastic_net",
    },
    {
        "type": "svm",
        "C": hp.lognormal("C", 0, 1.0),
        "kernel": hp.choice("kernel", ["linear", "rbf"])
    },
])


#3 Choose search algorithm
algo = tpe.suggest

"""
#4 apply hyperopt fmin()
#set max_evals, parallelism, and timeout

spark_trials = SparkTrials(parallelism=2, timeout=100)

#run fmin()
best_result = fmin(
    fn = objective,
    space = search_space,
    algo = algo,
    max_evals = 16,
    trials = spark_trials)

"""






#############
#Do hyperopt per algorithm (That is create search space per algorithm)

#Elastic Net

en_space = hp.choice("regressor", [
    {
        "type": "elastic_net",
        "alpha": hp.lognormal("alpha", 0, 10.0)
    }
])

svm_space = hp.choice("regressor", [
    {
        "type": "svm",
        "C": hp.lognormal("C", 0, 1.0),
        "kernel": hp.choice("kernel", ["linear", "rbf", "sigmoid", "poly"])
    }
])

svm_space = hp.choice("regressor", [
    {
        "type": "svm",
        "C": hp.lognormal("C", 0, 1.0),
        "kernel": hp.choice("kernel",
                            ["linear", "rbf", "sigmoid", {
                                "kernel": "poly",
                                "degree": hp.lognormal("degree",
                                                       2, 7)}])
    }
])

svm_space = hp.choice("regressor", [
    {
        "type": "svm",
        "C": hp.lognormal("C", 0, 1.0),
        "kernel": hp.choice("kernel", ["linear", "rbf", "sigmoid"]),
        "kernel": "poly", "poly": hp.lognormal("degree", 2, 7)
    }
])
#import hyperopt
#hyperopt.space_eval(svm_space,best_result)
#best_result = fmin(fn=objective, space=search_space, algo=algo, max_evals=25,trials=trials)


#Search spaces tryout

svm_space = {
    "type": "svm",
    "C": hp.lognormal("C", 0, 1.0),
    "kernel": hp.choice("kernel", ["linear", "rbf", "sigmoid", "poly"]),
    "degree": hp.choice("degree", range(2,8,1))
}

rf_space = {
    "type": "rf",
    "n_estimators": hp.choice("trees", range(50, 550, 50))
}

knn_space = {
    "type": "knn",
    "n_neighbors": hp.choice("neighbors", range(3, 33, 2)),
    "weights": hp.choice("weights", ["uniform", "distance"]),
    "p": hp.choice("p", range(1, 4, 1))
}

#concern
#why does hyperopt fmin keep changing results
#the optimum hyperparameters keep changing at each run
#Answer = Because max_evals less than 20 is done randomly
#solution = Increase max_evals to about 50


#trials object is where the loss (ie the minimized negative r2) values are.
#result class is where the loss are in the trials object
#thus
#table = pd.DataFrame(trials.results)
#table

#I can also just go straight to the best loss (which is just negative r2) with
#trials.best_trial["result"]["loss"]
#thus 5 fold CV R2
#-1 * trials.best_trial["result"]["loss"]
