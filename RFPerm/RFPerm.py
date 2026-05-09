#RFPerm:
import os
from dataclasses import dataclass
from typing import Optional, Dict, Any
import numpy as np
import pandas as pd
import torch
import random
from datasets import Dataset
from transformers import (
    AutoTokenizer,
    AutoModel,
    GPT2Tokenizer,
    GPT2Model,
    AutoModelForCausalLM,
    AutoModelForMaskedLM,
    DataCollatorForLanguageModeling,
    Trainer,
    TrainingArguments,
    set_seed,
)
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.preprocessing import OrdinalEncoder

def perm_oobtest(df_exist, df_new, n_B = 500, B = 5000, tuning = False,
    random_state = 42, regressor = True):
    X_exist = df_exist.iloc[:, :-1]
    Y_exist = df_exist.iloc[:, -1]
    X_new = df_new.iloc[:, :-1]
    Y_new = df_new.iloc[:, -1]
    pval_list = []
    n_minnode = round(np.sqrt(X_exist.shape[0]))
    for i in range(n_B):
        X_train, X_valid, Y_train, Y_valid = train_test_split(X_exist, Y_exist, test_size = 0.25, random_state = i)
        if regressor:
            model = RandomForestRegressor(
                    min_samples_leaf = n_minnode,
                    max_features = 0.33,
                    n_estimators = 150,
                    oob_score = True,
                    bootstrap = True,
                    random_state = random_state)
        else:
            model = RandomForestClassifier(
                    min_samples_leaf = n_minnode,
                    max_features = 0.33,
                    n_estimators = 150,
                    oob_score = True,
                    bootstrap = True,
                    random_state = random_state
                )
        model.fit(X_train, Y_train)
        Y_hat_oob = model.oob_decision_function_[:,1]
        Y_hat_test = model.predict_proba(X_new)[:,1]
        pval_list.append(permTest_(l2_loss(Y_hat_oob, Y_train), l2_loss(Y_hat_test, Y_new)))
        print(i)
    return statistics.mean([1 if t < 0.05 else 0 for t in pval_list])




def permOOB_randomSearch(df_exist, df_new, n_B = 500, B = 5000, randomSearch = True):
    X_exist = df_exist.iloc[:, :-1]
    Y_exist = df_exist.iloc[:, -1]
    X_new = df_new.iloc[:, :-1]
    Y_new = df_new.iloc[:, -1]
    for i in range(n_B):
        pipeline = Pipeline([
        ('rf', RandomForestRegressor())
        ])
        param_grid = {
        'rf__min_samples_leaf': [5, 10, 20],
        'rf__max_depth': [3, 5, 8],
        'rf__min_child_weight': [2, 4, 7],
        'rf__colsample_bytree': [0.5, 0.8, 0.99]
        }
        grid = RandomizedSearchCV(pipeline, param_grid, cv = 5, 
            n_jobs = -1, scoring = 'neg_mean_squared_error')
        

#Adding another paragraph here:



























