#!/usr/bin/env python
# coding: utf-8

# # Random Forest Model

# In[ ]:


import pandas as pd
import numpy as np
import argparse

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score
from sklearn import metrics


# In[ ]:


argparser = argparse.ArgumentParser(description = 'Ancestry inference via random forest model')
argparser.add_argument('-r', '--reference_PCs', required = True, help = 'Reference PCs in .coord file')
argparser.add_argument('-s', '--study_PCs', required = True, help = 'Study projected PCs in ProPC.coord file')
argparser.add_argument('-k', '--nPCs', type = int, required = False, help = 'Number of PCs to use', default=20)
argparser.add_argument('-p', '--reference_labels', required = True, help = 'Input file containing reference population labels')
argparser.add_argument('-m','--min_prob', required = False, default = 0, type = float, help = 'Minimum probability threshold for assigning population label')
argparser.add_argument('--seed', required = False, default = 11, type = int, help = 'Random seed')
#argparser.add_argument('-i', '--ID-col', help = 'Name of ID column in input files', default='ID')
#argparser.add_argument('-c', '--population_col', help = 'Name of population label column in input reference file', default='genetic_region')


# In[ ]:


if __name__ == '__main__':
    args = argparser.parse_args()
    
    # import files
    refPC = pd.read_csv(args.reference_PCs, sep='\t')
    studyProPC = pd.read_csv(args.study_PCs, sep='\t')
    ref_ancestry = pd.read_csv(args.reference_labels)
    
    # merge reference labels and PCs
    ref = pd.merge(left=refPC, right=ref_ancestry, left_on='indivID', right_on='ID', how='right')
    
    # extract PCs as features
    k = args.nPCs
    feature_list = ["PC"+str(i) for i in range(1,k+1)]
    labels = np.array(ref['genetic_region'])
    features = ref[feature_list]
    features = np.array(features)
    
    # split data into test and training sets
    seed = args.seed
    train_features, test_features, train_labels, test_labels = train_test_split(features, labels, test_size = 0.25, random_state = seed)
    # build RF model
    clf = RandomForestClassifier(n_estimators=100, random_state=seed)
    clf.fit(train_features, train_labels)
    
    # check model accuracy w/ cross-validation
    scores = cross_val_score(clf, train_features, train_labels, cv=5)
    # check model accuracy w/ separated test data
    test_pred=clf.predict(test_features)
    
    print('Average cross-validated accuracy:', np.mean(scores), ', with SE=', np.std(scores) / np.sqrt(len(scores)))
    print("Test data prediction accuracy:",metrics.accuracy_score(test_labels, test_pred))
    
    
    # predict ancestry for study samples
    studyProPC['predicted_ancestry'] = clf.predict(studyProPC[feature_list].values)
    # calculate and append probability for each prediction
    min_prob = args.min_prob
    probs = clf.predict_proba(studyProPC[feature_list].values)
    probs = pd.DataFrame(probs, columns=[f"prob_{p}" for p in clf.classes_])
    probs["max_prob"] = probs.max(axis=1)
    studyProPC.loc[probs["max_prob"] < min_prob, 'predicted_ancestry'] = "undefined"

    output = pd.concat([studyProPC[['indivID','predicted_ancestry']], probs],axis=1)
    output.to_csv('./predicted_ancestry.txt', index=False, sep='\t')
    
    # print out results summary
    result = pd.DataFrame(output.predicted_ancestry.value_counts())
    print("\n")
    for i in range(len(result.index)):
        print(result.predicted_ancestry[i]," samples classified as ",result.index[i])
    

