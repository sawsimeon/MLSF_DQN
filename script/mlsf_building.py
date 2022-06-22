#### selecting optimal features

from random import Random
import sys
sys.path.append("script/")
from mlsf_optimization import featurizers
from sklearn.tree import DecisionTreeRegressor
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor
from xgboost.sklearn import XGBRegressor
from sklearn.svm import SVR
from sklearn.neural_network import MLPRegressor
from tensorflow import keras  
import tensorflow as tf
from tensorflow.keras import regularizers, layers
from tensorflow.keras.layers import Dense, BatchNormalization, Dropout, Activation
from tensorflow.keras.optimizers import Adadelta, Adam, RMSprop
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr
from sklearn import preprocessing
protein = '/home/simeons/Desktop/MLSF_DQN/script/data/protein/1M17_insert_addH_modelledHis_Pro.mol2'
sdf = '/home/simeons/Desktop/MLSF_DQN/script/data/compound/egfr_ligand_library.sdf'
bioactivity_path = '/home/simeons/Desktop/MLSF_DQN/script/data/bioactivity.csv'


class DNNREgressor():
  def __init__(self):
    None
  def build(self):
    dnn = keras.Sequential([
      layers.Dense(units = 8192,kernel_regularizer = regularizers.l2(0), activation = "relu"),
      layers.BatchNormalization(),
      layers.Dropout(0),
      layers.Dense(units = 4396, kernel_regularizer = regularizers.l2(0), activation = "relu"),
      layers.BatchNormalization(),
      layers.Dense(units = 2048, kernel_regularizer = regularizers.l2(0), activation = "relu"),
      layers.Dropout(0),
      layers.Dense(units = 1, activation = "linear")
    ])
    dnn.compile(optimizer = "RMSprop", loss = "mean_squared_error", metrics = ['mean_squared_error'])
    return dnn 

def predict_function(model, train_features, train_pIC50, test_features, test_pIC50, name):
  if name == 'sklearn': 
    model.fit(np.array(train_features), train_pIC50)
    predict_pIC50 = model.predict(np.array(test_features))
  else:
    train_pIC50 = np.array(train_pIC50).astype("float32")
    model.fit(np.array(train_features) / train_features.shape[1], train_pIC50, epochs = 100, batch_size = 500, verbose = False)
    predict_pIC50 = model.predict(np.array(test_features))[:, 0]

  prediction_df = pd.DataFrame({"Observed_pIC50": test_pIC50,
                                "Predicted_pIC50": predict_pIC50})
  rmse = mean_squared_error(prediction_df['Predicted_pIC50'], prediction_df['Observed_pIC50'], squared= False)
  #r2 = r2_score(prediction_df['Predicted_pIC50'], prediction_df['Observed_pIC50'])
  r = pearsonr(x = np.array(prediction_df['Predicted_pIC50']), y = np.array(prediction_df['Observed_pIC50']))[0]
  return rmse, r

min_max_scaler = preprocessing.MinMaxScaler()

def load_data(feature_name, bioactivity_path):
  global features
  df_bioactivity = pd.read_csv(bioactivity_path)[['ChEMBLID', 'pIC50', 'Label']]
  if feature_name == "IF":
    features = featurizers(sdf = sdf, protein = protein).IF()
  if feature_name == 'SFP': 
    features = featurizers(sdf = sdf, protein = protein).SFP()
  if feature_name == 'BINANA':
    features = featurizers(sdf = sdf, protein = protein).BINANA()
  if feature_name == "RFSCORE V1":
    features = featurizers(sdf = sdf, protein = protein).RFSCORE(version = 1)
  if feature_name == "RFSCORE V2": 
    features = featurizers(sdf = sdf, protein = protein).RFSCORE(version = 2)
  if feature_name == 'RFSCORE V3': 
    features = featurizers(sdf = sdf, protein = protein).RFSCORE(version = 3)

  
  data = features.merge(df_bioactivity.drop_duplicates(subset = ['ChEMBLID']), how = 'left')
  train_descriptors = data[data['Label'].str.contains("Training")]
  test_descriptors = data[data['Label'].str.contains("Test")]
  train_pIC50 = list(train_descriptors['pIC50'])
  test_pIC50 = list(test_descriptors['pIC50'])
  train_descriptors.drop(['ChEMBLID', 'pIC50', 'Label'], axis = 1, inplace = True)
  test_descriptors.drop(['ChEMBLID', 'pIC50', 'Label'], axis = 1, inplace = True)
  train_descriptors = min_max_scaler.fit_transform(train_descriptors)
  test_descriptors = min_max_scaler.fit_transform(test_descriptors)

  return train_descriptors, test_descriptors, train_pIC50, test_pIC50

  
    
regressors = (Ridge, DecisionTreeRegressor, RandomForestRegressor, XGBRegressor, SVR, MLPRegressor, DNNREgressor)

methodnames = {Ridge: "Ridge", DecisionTreeRegressor: "DT", RandomForestRegressor: "RF", XGBRegressor: "XGB", 
               SVR: "SVM", MLPRegressor: "ANN", DNNREgressor: "DNN"}

features = ['IF', 'SFP', 'BINANA', 'RFSCORE V1', 'RFSCORE V2', 'RFSCORE V3', 'PLEC']

rmse_list,feature_list,method_list, r_list = [],[],[],[]
for feat in features:
  train_features, test_features, train_pIC50, test_pIC50 = load_data(feature_name = feat, bioactivity_path = bioactivity_path)
  for model in regressors:
    methodname = methodnames[model]
    if methodname == "Ridge":
      estimator = model(random_state = 0)
      rmse, r = predict_function(estimator, train_features, train_pIC50, test_features, test_pIC50, name = "sklearn")
    if methodname == 'DT':
      estimator = model(random_state = 0)
      rmse, r = predict_function(estimator, train_features, train_pIC50, test_features, test_pIC50, name = "sklearn")
    if methodname == 'RF':
      estimator = model(random_state = 0, n_estimators = 500)
      rmse, r = predict_function(estimator, train_features, train_pIC50, test_features, test_pIC50, name = "sklearn")
    if methodname == 'XGB':
      estimator = model(random_state = 0)
      rmse, r = predict_function(estimator, train_features, train_pIC50, test_features, test_pIC50, name = "sklearn")
    if methodname == 'SVM':
      estimator = model()
      rmse, r = predict_function(estimator, train_features, train_pIC50, test_features, test_pIC50, name = "sklearn")
    if methodname == 'ANN':
      estimator = model(random_state = 0)
      rmse, r = predict_function(estimator, train_features, train_pIC50, test_features, test_pIC50, name = "sklearn")
    if methodname == 'DNN':
      dnn = model()
      estimator = dnn.build()
      rmse, r  = predict_function(estimator, train_features, train_pIC50, test_features, test_pIC50, name = "tensorflow")
    print(rmse)
    
    rmse_list.append(rmse)
    feature_list.append(feat)
    method_list.append(methodname)
    r_list.append(r)

result = pd.DataFrame({"Feature": feature_list,
                      "Method": method_list,
                      "RMSE": rmse_list,
                      "R": r_list})

result.to_csv("/home/simeons/Desktop/MLSF_DQN/script/data/result_default.csv", index = False)