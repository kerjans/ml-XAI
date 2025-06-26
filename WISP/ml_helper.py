import numpy as np
import pickle

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator

from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import KFold
from sklearn.model_selection import HalvingRandomSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import statistics
from sklearn.metrics import mean_absolute_error, mean_squared_error, max_error, accuracy_score, precision_score, recall_score, f1_score, confusion_matrix
from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier
from sklearn.svm import SVC
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.neural_network import MLPRegressor
from sklearn.linear_model import BayesianRidge
from sklearn.gaussian_process.kernels import Matern
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import Lasso
from sklearn.svm import SVR

from WISP.plotting_helper import *

def get_features(data, CLOUMS):
    """
    Preprocesses the features according to theyr type.

    Keyword arguments:
    -- data: Table where the features are stored.
    -- CLOUMS: Colum where the feature is stored.

    Returns:
    -- X: Preprocessed feature.
    """
    features = []

    for i in CLOUMS:
        if type(data[i].values[0]) is np.ndarray:
            x = np.stack(data[i].values)
        else:
            x = data[i].values.reshape(-1, 1)
        features.append(x)

    X = np.concatenate(features, axis=1)
    return X

def get_morgan_fingerprint(smiles):
    """
    Calculates the Morgan Fingerprint (2028 bits, radius of 2) for the input smiles.

    Keyword arguments:
    -- smiles: Smiles for the which the fingerprint should be calculated.

    Returns:
    -- fingerprint: Morgan Fingerprint as array to the respective smiles.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol)#to keep the explicit hydrogens
    
    if mol is not None:
        generator = GetMorganGenerator(radius=2, fpSize=2048)
        fingerprint = generator.GetFingerprint(mol)
        #fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        fingerprint = fingerprint.ToBitString()
        fingerprint = np.array(list(fingerprint))
        return fingerprint
    else:
        return None

def get_MACCS_fingerprint(smiles):
    """
    Calculates the MCCS Keys Fingerprint for the input smiles.

    Keyword arguments:
    -- smiles: Smiles for the which the fingerprint should be calculated.

    Returns:
    -- maccs_fp: MCCS Keys Fingerprint as array to the respective smiles.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol)#to keep the explicit hydrogens
    
    if mol is not None:
        maccs_fp = AllChem.GetMACCSKeysFingerprint(mol)
        maccs_fp = maccs_fp.ToBitString()
        maccs_fp = np.array(list(maccs_fp))
        return maccs_fp
    else:
        return None
        

def get_RDK_fingerprint(smiles):
    """
    Calculates the RDK Fingerprint (Maximal pathlenght of 7 and 2048 bits) for the input smiles.

    Keyword arguments:
    -- smiles: Smiles for the which the fingerprint should be calculated.

    Returns:
    -- rdkit_fp: RDK Fingerprint as array to the respective smiles.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol)#to keep the explicit hydrogens
    
    if mol is not None:
        rdkit_fp = AllChem.RDKFingerprint(mol, maxPath=7)
        rdkit_fp = rdkit_fp.ToBitString()
        rdkit_fp = np.array(list(rdkit_fp))
        return rdkit_fp
    else:
        return None
    
def hp_search_helper(model,df_train,target,feature):
    
    PARAM_GRID = {
    'SVC': {'model__C': [0.1, 1, 10, 100], 'model__kernel': ['rbf'], 'model__class_weight': ['balanced'], 'model__gamma': ['scale', 'auto', 1, 0.001, 0.01, 0.1]},
    'RandomForestClassifier': {'n_estimators': [400,700,1000], 'class_weight': ['balanced'], 'min_samples_leaf': [2,3]},
    
    'MLPClassifier': {'model__hidden_layer_sizes': [(50,), (100,), (50, 50)], 'model__activation': ['relu', 'tanh'], 'model__solver': ['adam'], 'model__alpha': [0.0001, 0.001, 0.01], 'model__learning_rate': ['constant', 'adaptive'], 'model__max_iter': [200, 500]},
    'GradientBoostingClassifier': {'model__n_estimators': [100, 300, 500], 'model__learning_rate': [0.01, 0.05, 0.1], 'model__max_depth': [3, 5, 7], 'model__subsample': [0.6, 0.8, 1.0], 'model__min_samples_split': [2, 5, 10]},
    'GaussianProcessClassifier': {'model__n_restarts_optimizer': [0, 2, 5], 'model__max_iter_predict': [100, 200], 'model__multi_class': ['one_vs_rest'], 'model__warm_start': [True, False]},

    'SVR': {'model__C': [0.1, 1, 10, 100], 'model__kernel': ['rbf'],  'model__gamma': ['scale', 'auto', 1, 0.001, 0.01, 0.1]},
    'RandomForestRegressor': {'model__n_estimators': [400,700,1000], 'model__max_depth':[30, 50], 'model__n_jobs':[32], 'model__random_state': [42]},
    
    'MLPRegressor': {'model__hidden_layer_sizes': [(50,), (100,), (50, 50)],'model__activation': ['relu', 'tanh'],'model__solver': ['adam'],'model__alpha': [0.0001, 0.001, 0.01],'model__learning_rate': ['constant', 'adaptive'],'model__max_iter': [200, 500], 'model__random_state': [42]},
    'BayesianRidge': {'model__alpha_1': [1e-7, 1e-6, 1e-5],'model__alpha_2': [1e-7, 1e-6, 1e-5],'model__lambda_1': [1e-7, 1e-6, 1e-5],'model__lambda_2': [1e-7, 1e-6, 1e-5]},
    'Lasso': {'model__alpha': [0.001, 0.01, 0.1, 1.0, 10.0],'model__max_iter': [1000, 5000],'model__tol': [1e-4, 1e-3, 1e-2],'model__selection': ['cyclic', 'random'], 'model__random_state': [42]},
    'GradientBoostingRegressor': {'model__n_estimators': [100, 300, 500],'model__learning_rate': [0.01, 0.05, 0.1],'model__max_depth': [3, 5, 7],'model__subsample': [0.6, 0.8, 1.0],'model__min_samples_split': [2, 5, 10], 'model__random_state': [42]},
    'LinearRegression': {},
    'GaussianProcessRegressor': {'model__alpha': [1e-10, 1e-5, 1e-2],'model__n_restarts_optimizer': [0, 2, 5, 10],'model__normalize_y': [True, False]}
    }

    SCORING = {
        'SVC': 'f1',
        'RandomForestClassifier': 'f1',
        'MLPClassifier': 'f1',
        'GradientBoostingClassifier': 'f1',
        'GaussianProcessClassifier': 'f1',
        'RandomForestRegressor': "r2",
        'SVR': "r2",
        'MLPRegressor': 'r2',
        'BayesianRidge': 'r2',
        'Lasso': 'r2',
        'GradientBoostingRegressor': 'r2',
        'LinearRegression': 'r2',
        'GaussianProcessRegressor': 'r2'
    }

    param_grid = PARAM_GRID[model.__class__.__name__]
    kf = KFold(n_splits=5, shuffle=True, random_state=42) 
    pipe = HalvingRandomSearchCV(Pipeline([('scaler', StandardScaler()), ('model', model.__class__())]), param_distributions=param_grid, random_state=42, refit=True, cv=kf, scoring=SCORING[model.__class__.__name__])
    
    prep_train = get_features(df_train, feature)
    target_train = df_train[target].values

    pipe.fit(prep_train, target_train)

    splitted_data = kf.split(df_train)

    correlation_coff = []
    coef_of_determ = []
    mean_AE = []
    RMSE =[]

    for train_fold, test_fold in splitted_data:
        train_data = df_train.iloc[train_fold]
        test_data = df_train.iloc[test_fold]

        prep_train = get_features(train_data, feature)
        prep_test = get_features(test_data, feature)
        target_train = train_data[target].values
        target_test = test_data[target].values
            
        best_estimator = pipe.best_estimator_
        model = best_estimator.fit(prep_train, target_train)
            
        predictions = model.predict(prep_test)

        #stats
        r2 = np.corrcoef(target_test.flatten(), predictions.flatten())[0,1]**2
        correlation_coff.append(r2)
        
        RMSE_z = np.sqrt(np.mean((target_test.flatten() - predictions.flatten())**2))
        RMSE_n = np.sqrt(np.mean((target_test.flatten() - np.mean(target_test.flatten()))**2))
        R2 = 1 - RMSE_z**2/RMSE_n**2
        coef_of_determ.append(R2)
        
        MAE = mean_absolute_error(target_test.flatten(), predictions.flatten())
        mean_AE.append(MAE)
        
        mse = mean_squared_error(target_test.flatten(), predictions.flatten())
        rmse = np.sqrt(mse)
        RMSE.append(rmse)

    mean_r2 = statistics.mean(correlation_coff)
    mean_R2 = statistics.mean(coef_of_determ)
    mean_MAE = statistics.mean(mean_AE)
    mean_RMSE = statistics.mean(RMSE)
    
    return pipe.best_estimator_, mean_r2, mean_R2, mean_MAE, mean_RMSE

def split_data(data, Target_Column_Name, working_dir):
    nr_test_samples = round(len(data) / 5) # 80/20 split
    test = data.sample(n=nr_test_samples, random_state=6)
    target_test = test[Target_Column_Name].values
    train = data.drop(test.index)
    #plot_histogram(test, train,  Target_Column_Name, 'Target', 'Stucture Count', 'Test Set', 'Train Set',  working_dir + 'Count-Bins-test-train.png')
    return test, target_test, train

def features_and_reg_model_types(data):
    data['Morgan_Fingerprint 2048Bit 2rad'] = data['smiles_std'].apply(get_morgan_fingerprint)
    data['MACCS_Fingerprint'] = data['smiles_std'].apply(get_MACCS_fingerprint)
    data['RDK_Fingerprint'] = data['smiles_std'].apply(get_RDK_fingerprint)

    ALLfeatureCOLUMS = ['Morgan_Fingerprint 2048Bit 2rad',
            'RDK_Fingerprint',
            'MACCS_Fingerprint']
        
    model_types = [MLPRegressor(), 
            BayesianRidge(), 
            Lasso(), 
            GradientBoostingRegressor(), 
            LinearRegression(), 
            RandomForestRegressor(), 
            SVR(),
            GaussianProcessRegressor(kernel=Matern())]
    return data, ALLfeatureCOLUMS, model_types

def features_and_class_model_types(data):
    data['Morgan_Fingerprint 2048Bit 2rad'] = data['smiles_std'].apply(get_morgan_fingerprint)
    data['MACCS_Fingerprint'] = data['smiles_std'].apply(get_MACCS_fingerprint)
    data['RDK_Fingerprint'] = data['smiles_std'].apply(get_RDK_fingerprint)

    ALLfeatureCOLUMS = ['Morgan_Fingerprint 2048Bit 2rad',
            'RDK_Fingerprint',
            'MACCS_Fingerprint']
        
    model_types = [MLPClassifier(), 
            GradientBoostingClassifier(), 
            RandomForestClassifier(), 
            SVC(),
            GaussianProcessClassifier()]
    return data, ALLfeatureCOLUMS, model_types

def get_best_reg_model(model_types, ALLfeatureCOLUMS, train, Target_Column_Name, working_dir):
    results = []

    for model_arc in model_types:          
        for feature in ALLfeatureCOLUMS:
            model, r2, R2, MAE, RMSE = hp_search_helper(model_arc,train,Target_Column_Name,[str(feature)])
            results.append({'Feature': feature,'Model_Type': model_arc,'Model': model,'r2': r2,'R2': R2,'MAE': MAE,'RMSE': RMSE})
    

    results_df = pd.DataFrame(results)
    best_model_row = results_df.loc[results_df['MAE'].idxmin()]
    model = best_model_row['Model']
    print('Best Model: ', best_model_row['Model'])
    print('With a MAE of: ', best_model_row['MAE'])
    print('Feature: ', best_model_row['Feature'])
    results_df.to_csv(working_dir + "Grid-Search.csv", index=False)

    #pick feature function
    if best_model_row['Feature'] == 'Morgan_Fingerprint 2048Bit 2rad':
        feature_function = get_morgan_fingerprint
    if best_model_row['Feature'] == 'RDK_Fingerprint':
        feature_function = get_RDK_fingerprint
    if best_model_row['Feature'] == 'MACCS_Fingerprint':
        feature_function = get_MACCS_fingerprint

    #train model on trining set
    featureCOLUMS = [best_model_row['Feature']]

    return model, feature_function, featureCOLUMS

class ProbabilisticRandomForest:
    def __init__(self, **kwargs):
        self.model = RandomForestClassifier(**kwargs)

    def fit(self, X, y):
        return self.model.fit(X, y)

    def predict(self, X):
        #return self.model.predict_proba(X)
        proba = self.model.predict_proba(X)
        return proba[:, 1]

    def predict_classes(self, X):
        return self.model.predict(X)

    def __getattr__(self, attr):
        # Delegate attribute access to the underlying model
        return getattr(self.model, attr)
    
    def __getstate__(self):
        return self.model

    def __setstate__(self, state):
        self.model = state
    
    @property
    def estimators_(self):
        return self.model.estimators_

    @property
    def n_classes_(self):
        return self.model.n_classes_

    @property
    def classes_(self):
        return self.model.classes_

def get_and_train_class_model(train, test, Target_Column_Name, target_test, working_dir):
    
    #fixed settings
    feature_function = get_morgan_fingerprint
    featureCOLUMS = ['Morgan_Fingerprint 2048Bit 2rad']
    #model = Pipeline(steps=[('scaler', StandardScaler()),
                #('model',
                 #RandomForestClassifier(random_state=42))])
    model = Pipeline(steps=[('scaler', StandardScaler()),
                ('model',
                 ProbabilisticRandomForest(random_state=42))])

    #train on whole training set
    prep_train = get_features(train, featureCOLUMS)
    target_train = train[Target_Column_Name].values
    model.fit(prep_train, target_train)

    #performance on test set
    prep_test = get_features(test, featureCOLUMS)
    predictions = model.predict(prep_test)

    #statistic on testset
    predicted_labels = (predictions >= 0.5).astype(int)
    accuracy = accuracy_score(target_test, predicted_labels)
    precision = precision_score(target_test, predicted_labels, average='weighted')
    recall = recall_score(target_test, predicted_labels, average='weighted')
    f1 = f1_score(target_test, predicted_labels, average='weighted')
    conf_matrix = confusion_matrix(target_test, predicted_labels)
    '''
    accuracy = accuracy_score(target_test, predictions)
    precision = precision_score(target_test, predictions, average='weighted')  # or 'macro', 'micro'
    recall = recall_score(target_test, predictions, average='weighted')
    f1 = f1_score(target_test, predictions, average='weighted')
    conf_matrix = confusion_matrix(target_test, predictions)
    '''

    # Print results
    print("Performance on testset:")
    print(f"Accuracy: {accuracy}")
    print(f"Precision: {precision}")
    print(f"Recall: {recall}")
    print(f"F1-score: {f1}")
    print("Confusion Matrix:\n", conf_matrix)

    pickle.dump(model, open(working_dir + "model.pkl", "wb"))

    return model, feature_function, featureCOLUMS


def train_and_evaluate_reg_model(model, train, test, featureCOLUMS, Target_Column_Name, target_test, working_dir):

    #train on whole training set
    prep_train = get_features(train, featureCOLUMS)
    target_train = train[Target_Column_Name].values
    model.fit(prep_train, target_train)

    #performance on test set
    prep_test = get_features(test, featureCOLUMS)
    predictions = model.predict(prep_test)

    #statistic on testset
    r2 = np.corrcoef(target_test.flatten(), predictions.flatten())[0,1]**2
    RMSE_z = np.sqrt(np.mean((target_test.flatten() - predictions.flatten())**2))
    RMSE_n = np.sqrt(np.mean((target_test.flatten() - np.mean(target_test.flatten()))**2))
    R2 = 1 - RMSE_z**2/RMSE_n**2
    MAE = mean_absolute_error(target_test.flatten(), predictions.flatten())
    max_err = max_error(target_test.flatten(), predictions.flatten())
    mse = mean_squared_error(target_test.flatten(), predictions.flatten())
    rmse = np.sqrt(mse)

    #print/plot results
    print('Performance on testset(r2, R2, MAE, RMSE, Maximal Error, MSE):',r2,';',R2,';',MAE,';',rmse,';',max_err,';', mse)

    r2 = np.corrcoef(predictions, target_test)[0,1]**2
    plot_2D(['r$^2$ = ' + str(f"{r2:.2f}")], 'upper left', predictions , target_test,
            'predicted', 'experimental', working_dir + '20-80-split-true-pred.png', '#A43341', 
            include_line=False, line_style='None')
    
    #save model training results
    pickle.dump(model, open(working_dir + "model.pkl", "wb"))

def add_predictions(data_MMPs, feature_function, model):
    #for 1
    data_MMPs['Feature_1'] = data_MMPs['smiles_1'].apply(feature_function)
    X_data_attributions_1 = get_features(data_MMPs, ['Feature_1'])
    predictions_1 = model.predict(X_data_attributions_1)
    data_MMPs['predictions_1'] = predictions_1

    #for 2
    data_MMPs['Feature_2'] = data_MMPs['smiles_2'].apply(feature_function)
    X_data_attributions_1 = get_features(data_MMPs, ['Feature_2'])
    predictions_1 = model.predict(X_data_attributions_1)
    data_MMPs['predictions_2'] = predictions_1
    
    return data_MMPs

def train_test_dependency(data_MMPs, test, Attribution_Colums):
    data_MMPs['set_1'] = data_MMPs['smiles_1'].isin(test['smiles_std']).map({True: 'test', False: 'train'})
    data_MMPs['set_2'] = data_MMPs['smiles_2'].isin(test['smiles_std']).map({True: 'test', False: 'train'})

    def compute_r2(data, label):
        mask = data['set_1'].str.contains(label) & data['set_2'].str.contains(label)
        filtered = data[mask]
        for attr in Attribution_Colums:
            r2 = np.corrcoef(
                filtered['delta_prediction'],
                filtered[f'delta_sum_{attr}']
            )[0, 1] ** 2
            print(f'{attr}_{label}:', r2)
    
    for split in ['train', 'test']:
        compute_r2(data_MMPs, split)
    
    return data_MMPs

def split_MMPs_by_set(data_MMPs, test):
    data_MMPs['set_1'] = data_MMPs['smiles_1'].isin(test['smiles_std']).map({True: 'test', False: 'train'})
    data_MMPs['set_2'] = data_MMPs['smiles_2'].isin(test['smiles_std']).map({True: 'test', False: 'train'})

    train_set = data_MMPs[(data_MMPs['set_1'] == 'train') & (data_MMPs['set_2'] == 'train')]
    test_set = data_MMPs[(data_MMPs['set_1'] == 'test') & (data_MMPs['set_2'] == 'test')]

    return train_set, test_set