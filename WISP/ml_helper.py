import numpy as np
from standardizer.io import else_none
np.seterr(divide='ignore', invalid='ignore')
import pickle
import shutil

import chemprop

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
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Lasso
from sklearn.svm import SVR

from WISP.plotting_helper import *
from WISP.chemprop import *

def get_features(data, COLUMNS):
    """
    Preprocesses the features according to theyr type.

    Keyword arguments:
    -- data: Table where the features are stored.
    -- COLUMNS: Column where the feature is stored.

    Returns:
    -- X: Preprocessed feature.
    """
    features = []

    for i in COLUMNS:
        if type(data[i].values[0]) is np.ndarray:
            x = np.stack(data[i].values)
        else:
            x = data[i].values.reshape(-1, 1)
        features.append(x)

    X = np.concatenate(features, axis=1)
    return X

from mordred import Calculator, descriptors
def get_mordred_descriptors(smiles):
    calc = Calculator(descriptors, ignore_3D=True)
    mols = [else_none(Chem.MolFromSmiles)(smi) for smi in smiles]
    df_descs = calc.pandas(mols,nproc=1,quiet=True,)
    
    mat_descs = []
    for col in sorted(df_descs.columns):
        mat_descs.append(np.vstack(df_descs[col].tolist()))

    return np.hstack(mat_descs)

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
    
def hp_search_helper(model, df_train, target,feature):
    """
    Perform a halving random‐search for hyperparameters on a scaler+model pipeline,
    then evaluate cross‐validated performance on df_train.

    Parameters:
        model (sklearn estimator):
            An unfitted classifier or regressor instance.
        df_train (pd.DataFrame):
            Training data containing feature columns and a target column.
        target (str):
            Name of the target column in df_train.
        feature (any):
            Argument forwarded to `get_features(df, feature)` to extract model inputs.

    Returns:
        best_estimator (Pipeline):
            The pipeline (StandardScaler + model) refit on the full training set.
        mean_r2 (float):
            Average squared Pearson correlation (r²) across CV folds.
        mean_R2 (float):
            Average coefficient of determination across folds.
        mean_MAE (float):
            Average mean absolute error across folds.
        mean_RMSE (float):
            Average root mean squared error across folds.
    """
    
    PARAM_GRID = {
    'SVC': {'model__C': [0.1, 1, 10, 100], 'model__kernel': ['rbf'], 'model__class_weight': ['balanced'], 'model__gamma': ['scale', 'auto', 1, 0.001, 0.01, 0.1]},
    'RandomForestClassifier': {'model__n_jobs':[6], 'n_estimators': [400,700,1000], 'class_weight': ['balanced'], 'min_samples_leaf': [2,3]},
    
    'MLPClassifier': {'model__hidden_layer_sizes': [(50,), (100,), (50, 50)], 'model__activation': ['relu', 'tanh'], 'model__solver': ['adam'], 'model__alpha': [0.0001, 0.001, 0.01], 'model__learning_rate': ['constant', 'adaptive'], 'model__max_iter': [200, 500]},
    'GradientBoostingClassifier': {'model__n_estimators': [100, 300, 500], 'model__learning_rate': [0.01, 0.05, 0.1], 'model__max_depth': [3, 5, 7], 'model__subsample': [0.6, 0.8, 1.0], 'model__min_samples_split': [2, 5, 10]},
    'HistGradientBoostingClassifier': {'model__learning_rate': [0.01, 0.05, 0.1],'model__max_depth': [3, 5, 7, None], 'model__random_state': [42]},
    'GaussianProcessClassifier': {'model__n_restarts_optimizer': [0, 2, 5], 'model__max_iter_predict': [100, 200], 'model__multi_class': ['one_vs_rest'], 'model__warm_start': [True, False]},

    'SVR': {'model__C': [0.1, 1, 10, 100], 'model__kernel': ['rbf'],  'model__gamma': ['scale', 'auto', 1, 0.001, 0.01, 0.1]},
    'RandomForestRegressor': {'model__n_estimators': [400,700,1000], 'model__max_depth':[30, 50], 'model__n_jobs':[6], 'model__random_state': [42]},
    
    'MLPRegressor': {'model__hidden_layer_sizes': [(50,), (100,), (50, 50)],'model__activation': ['relu', 'tanh'],'model__solver': ['adam'],'model__alpha': [0.0001, 0.001, 0.01],'model__learning_rate': ['constant', 'adaptive'],'model__max_iter': [200, 500], 'model__random_state': [42]},
    'BayesianRidge': {'model__alpha_1': [1e-7, 1e-6, 1e-5],'model__alpha_2': [1e-7, 1e-6, 1e-5],'model__lambda_1': [1e-7, 1e-6, 1e-5],'model__lambda_2': [1e-7, 1e-6, 1e-5]},
    'Lasso': {'model__alpha': [0.001, 0.01, 0.1, 1.0, 10.0],'model__max_iter': [1000, 5000],'model__tol': [1e-4, 1e-3, 1e-2],'model__selection': ['cyclic', 'random'], 'model__random_state': [42]},
    'GradientBoostingRegressor': {'model__n_estimators': [100, 300, 500],'model__learning_rate': [0.01, 0.05, 0.1],'model__max_depth': [3, 5, 7],'model__subsample': [0.6, 0.8, 1.0],'model__min_samples_split': [2, 5, 10], 'model__random_state': [42]},
    'HistGradientBoostingRegressor': {'model__learning_rate': [0.01, 0.05, 0.1],'model__max_depth': [3, 5, 7, None], 'model__random_state': [42]},
    'LinearRegression': {},
    'GaussianProcessRegressor': {'model__alpha': [1e-10, 1e-5, 1e-2],'model__n_restarts_optimizer': [0, 2, 5, 10],'model__normalize_y': [True, False]}
    }

    SCORING = {
        'SVC': 'f1',
        'RandomForestClassifier': 'f1',
        'MLPClassifier': 'f1',
        'GradientBoostingClassifier': 'f1',
        'HistGradientBoostingClassifier': 'f1',
        'GaussianProcessClassifier': 'f1',
        'RandomForestRegressor': "r2",
        'SVR': "r2",
        'MLPRegressor': 'r2',
        'BayesianRidge': 'r2',
        'Lasso': 'r2',
        'GradientBoostingRegressor': 'r2',
        'HistGradientBoostingRegressor': 'r2',
        'LinearRegression': 'r2',
        'GaussianProcessRegressor': 'r2'
    }

    param_grid = PARAM_GRID[model.__class__.__name__]
    kf = KFold(n_splits=5, shuffle=True, random_state=42) 
    pipe = HalvingRandomSearchCV(Pipeline([('scaler', StandardScaler()), ('model', model.__class__())]), param_distributions=param_grid, random_state=42, refit=True, cv=kf, scoring=SCORING[model.__class__.__name__])
    
    prep_train = get_features(df_train, feature)
    target_train = df_train[target].values

    pipe.fit(prep_train, target_train)

    splits_data = kf.split(df_train)

    correlation_coeff = []
    coef_of_determ = []
    mean_AE = []
    RMSE =[]

    for train_fold, test_fold in splits_data:
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
        correlation_coeff.append(r2)
        
        RMSE_z = np.sqrt(np.mean((target_test.flatten() - predictions.flatten())**2))
        RMSE_n = np.sqrt(np.mean((target_test.flatten() - np.mean(target_test.flatten()))**2))
        R2 = 1 - RMSE_z**2/RMSE_n**2
        coef_of_determ.append(R2)
        
        MAE = mean_absolute_error(target_test.flatten(), predictions.flatten())
        mean_AE.append(MAE)
        
        mse = mean_squared_error(target_test.flatten(), predictions.flatten())
        rmse = np.sqrt(mse)
        RMSE.append(rmse)

    mean_r2 = statistics.mean(correlation_coeff)
    mean_R2 = statistics.mean(coef_of_determ)
    mean_MAE = statistics.mean(mean_AE)
    mean_RMSE = statistics.mean(RMSE)
    
    return pipe.best_estimator_, mean_r2, mean_R2, mean_MAE, mean_RMSE

def split_data(data, Target_Column_Name):
    """
    Split a DataFrame into an 80/20 train/test split and return the test set,
    its target values, and the remaining training set.

    Parameters:
        data (pd.DataFrame): Full dataset containing features and target.
        Target_Column_Name (str): Name of the target column in `data`.

    Returns:
        test (pd.DataFrame): Randomly sampled ~20% of the rows.
        target_test (np.ndarray): 1D array of target values for `test`.
        train (pd.DataFrame): The remaining ~80% of the rows.
    """
    nr_test_samples = round(len(data) / 5)
    test = data.sample(n=nr_test_samples, random_state=6)
    target_test = test[Target_Column_Name].values
    train = data.drop(test.index)
    return test, target_test, train

def features_and_reg_model_types(data,fast_run=False):
    """
    Generate molecular fingerprint features and return regression model candidates.

    Parameters:
        data (pd.DataFrame): Must contain a 'smiles_std' column of standardized SMILES strings.

    Returns:
        data (pd.DataFrame):
            The input DataFrame augmented with:
              - 'Morgan_Fingerprint 2048Bit 2rad'
              - 'MACCS_Fingerprint'
              - 'RDK_Fingerprint'
        feature_columns (list of str):
            Names of the fingerprint columns added above.
        model_types (list of sklearn estimators):
            A list of untrained regression model instances to try.
    """
    data['Morgan_Fingerprint 2048Bit 2rad'] = data['smiles_std'].apply(get_morgan_fingerprint)
    data['MACCS_Fingerprint'] = data['smiles_std'].apply(get_MACCS_fingerprint)
    data['RDK_Fingerprint'] = data['smiles_std'].apply(get_RDK_fingerprint)

    descs = get_mordred_descriptors(data['smiles_std'].tolist())
    data['mordred'] = list(descs)

    if fast_run:
        ALLfeatureCOLUMNS = [
            'Morgan_Fingerprint 2048Bit 2rad',
            #"mordred",
            ]
    else:
        ALLfeatureCOLUMNS = ['Morgan_Fingerprint 2048Bit 2rad',
            'RDK_Fingerprint',
            'MACCS_Fingerprint']
    
        
    if fast_run:
        model_types = [
            #RandomForestRegressor(),
             HistGradientBoostingRegressor(),
             ]
    else:
        model_types = [MLPRegressor(), 
            BayesianRidge(), 
            Lasso(), 
            GradientBoostingRegressor(), 
            LinearRegression(), 
            RandomForestRegressor(), 
            SVR(),
            GaussianProcessRegressor(kernel=Matern())]
    return data, ALLfeatureCOLUMNS, model_types

def features_and_class_model_types(data,fast_run=False,):
    """
    Create fingerprint features and return classification model candidates.

    Parameters:
        data (pd.DataFrame): Must have a 'smiles_std' column of standardized SMILES.

    Returns:
        data (pd.DataFrame): The input augmented with:
            - 'Morgan_Fingerprint 2048Bit 2rad'
            - 'MACCS_Fingerprint'
            - 'RDK_Fingerprint'
        feature_columns (list of str): Names of the added fingerprint columns.
        model_types (list): Unfitted classification model instances to try.
    """
    data['Morgan_Fingerprint 2048Bit 2rad'] = data['smiles_std'].apply(get_morgan_fingerprint)
    data['MACCS_Fingerprint'] = data['smiles_std'].apply(get_MACCS_fingerprint)
    data['RDK_Fingerprint'] = data['smiles_std'].apply(get_RDK_fingerprint)
    data['mordred'] = get_mordred_descriptors(data['smiles_std'].tolist())

    if fast_run:
        ALLfeatureCOLUMNS = [
            #'Morgan_Fingerprint 2048Bit 2rad',
            "mordred"]
    else:
        ALLfeatureCOLUMNS = ['Morgan_Fingerprint 2048Bit 2rad',
            'RDK_Fingerprint',
            'MACCS_Fingerprint']
        
    if fast_run:
        model_types = [
            RandomForestClassifier(),HistGradientBoostingClassifier(),
        ]
    else:
        model_types = [MLPClassifier(), 
            GradientBoostingClassifier(), 
            RandomForestClassifier(), 
            SVC(),
            GaussianProcessClassifier()]

    return data, ALLfeatureCOLUMNS, model_types

def get_best_reg_model(model_types, ALLfeatureCOLUMNS, train, Target_Column_Name, working_dir, use_GNN):
    """
    Try GNN (optional) and grid‐searched regressors on fingerprint features, 
    then pick and return the best model and its feature function based on the MAE.

    Parameters:
        model_types (list): Unfitted regressor instances to tune.
        ALLfeatureCOLUMNS (list of str): Names of fingerprint columns in `train`.
        train (pd.DataFrame): Training set with features and target.
        Target_Column_Name (str): Name of the regression target column.
        working_dir (str): Directory for checkpoints and CSV output.
        use_GNN (bool): If True, train/evaluate a Chemprop GNN first.

    Returns:
        best_model:    The fitted model with lowest MAE.
        feature_fn:    Function used to compute features for that model.
        feature_cols:  List containing the single feature column name.
    """
    results = []

    if use_GNN is True:
        checkpoint_dir = working_dir + 'checkpoints/' 
        if os.path.exists(checkpoint_dir) and os.path.isdir(checkpoint_dir):
            shutil.rmtree(checkpoint_dir) #deletes old checkpoint dir
        r2, R2, MAE, RMSE = train_GNN(train, 'smiles_std', Target_Column_Name, working_dir)
        model_GNN = SklChemprop(problem_type="regression", max_epochs=50, Smiles_Column_Name='smiled_std', Target_Column_Name=Target_Column_Name, working_dir=working_dir)
        results.append({'Feature': 'smiles_std','Model_Type': 'SklChemprop','Model': model_GNN,'r2': r2,'R2': R2,'MAE': MAE,'RMSE': RMSE})
    for model_arc in model_types:          
        for feature in ALLfeatureCOLUMNS:
            model, r2, R2, MAE, RMSE = hp_search_helper(model_arc,train,Target_Column_Name,[str(feature)])
            results.append({'Feature': feature,'Model_Type': model_arc,'Model': model,'r2': r2,'R2': R2,'MAE': MAE,'RMSE': RMSE})

    results_df = pd.DataFrame(results)
    best_model_row = results_df.loc[results_df['MAE'].idxmin()]
    model = best_model_row['Model']
    print('Best Model: ', best_model_row['Model'])
    print('With a MAE of: ', best_model_row['MAE'])
    print('Feature: ', best_model_row['Feature'])
    results_df.to_csv(working_dir + "Grid-Search.csv", index=False)

    #delete checkpoints folder
    if model.__class__.__name__ != "SklChemprop":
        checkpoints_path = os.path.join(working_dir, 'checkpoints')
        if os.path.isdir(checkpoints_path):
            shutil.rmtree(checkpoints_path)

    #pick feature function
    if best_model_row['Feature'] == 'Morgan_Fingerprint 2048Bit 2rad':
        feature_function = get_morgan_fingerprint
    if best_model_row['Feature'] == 'RDK_Fingerprint':
        feature_function = get_RDK_fingerprint
    if best_model_row['Feature'] == 'MACCS_Fingerprint':
        feature_function = get_MACCS_fingerprint
    if best_model_row['Feature'] == 'mordred':
        feature_function = get_mordred_descriptors
    if best_model_row['Feature'] == 'smiles_std':
        feature_function = identity

    featureCOLUMNS = [best_model_row['Feature']]

    return model, feature_function, featureCOLUMNS

class ProbabilisticRandomForest:
    """
    A thin wrapper around sklearn’s RandomForestClassifier that:
      - uses predict(X) to return the probability of the positive class (column 1)
      - provides a separate predict_classes(X) for hard labels
      - delegates all other attributes/methods to the underlying classifier
      - supports pickle via __getstate__/__setstate__
    """
    def __init__(self, **kwargs):
        self.model = RandomForestClassifier(**kwargs)

    def fit(self, X, y):
        return self.model.fit(X, y)

    def predict(self, X):
        proba = self.model.predict_proba(X)
        return proba[:, 1]

    def predict_classes(self, X):
        return self.model.predict(X)

    def __getattr__(self, attr):
        # delegate attribute access to the underlying model
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
    """
    Train a probabilistic random forest on Morgan fingerprints, evaluate on the test set,
    print classification metrics, save the fitted model, and return the pipeline and feature info.

    Parameters:
        train (pd.DataFrame): Training data with a 'smiles_std' column and target column.
        test (pd.DataFrame): Test data in the same format.
        Target_Column_Name (str): Name of the target column in train/test.
        target_test (np.ndarray): 1D array of true labels for test.
        working_dir (str): Directory where 'model.pkl' will be saved.

    Returns:
        model (Pipeline): Fitted sklearn Pipeline (StandardScaler + ProbabilisticRandomForest).
        feature_function (callable): get_morgan_fingerprint.
        feature_columns (list of str): ['Morgan_Fingerprint 2048Bit 2rad'].
    """    
    #fixed settings
    feature_function = get_morgan_fingerprint
    featureCOLUMNS = ['Morgan_Fingerprint 2048Bit 2rad']
    model = Pipeline(steps=[('scaler', StandardScaler()),
                ('model',
                 ProbabilisticRandomForest(random_state=42))])

    #train on whole training set
    prep_train = get_features(train, featureCOLUMNS)
    target_train = train[Target_Column_Name].values
    model.fit(prep_train, target_train)

    #performance on test set
    prep_test = get_features(test, featureCOLUMNS)
    predictions = model.predict(prep_test)

    #statistic on test set
    predicted_labels = (predictions >= 0.5).astype(int)
    accuracy = accuracy_score(target_test, predicted_labels)
    precision = precision_score(target_test, predicted_labels, average='weighted')
    recall = recall_score(target_test, predicted_labels, average='weighted')
    f1 = f1_score(target_test, predicted_labels, average='weighted')
    conf_matrix = confusion_matrix(target_test, predicted_labels)

    #print results
    print("Performance on test set:")
    print(f"Accuracy: {accuracy}")
    print(f"Precision: {precision}")
    print(f"Recall: {recall}")
    print(f"F1-score: {f1}")
    print("Confusion Matrix:\n", conf_matrix)

    pickle.dump(model, open(working_dir + "model.pkl", "wb"))

    return model, feature_function, featureCOLUMNS


def train_and_evaluate_reg_model(model, train, test, featureCOLUMNS, Target_Column_Name, target_test, working_dir):
    """
    Train (if not a Chemprop model) and evaluate a regressor on an 80/20 split,
    printing key metrics, plotting true vs. predicted values, and saving the model.

    Parameters:
        model (estimator or Pipeline):
            The regression model or Pipeline to train/evaluate.
        train (pd.DataFrame):
            Training set containing feature columns and Target_Column_Name.
        test (pd.DataFrame):
            Test set containing feature columns and Target_Column_Name.
        featureCOLUMNS (list of str):
            Column names to pass to get_features() for X‐matrix construction.
        Target_Column_Name (str):
            Name of the target column in train/test DataFrames.
        target_test (array‐like):
            True target values for the test set.
        working_dir (str):
            Directory path where the plot and model.pkl are saved.

    Returns:
        None
    """
    #train on whole training set
    if model.__class__.__name__ != "SklChemprop":
        prep_train = get_features(train, featureCOLUMNS)
        target_train = train[Target_Column_Name].values
        model.fit(prep_train, target_train)

    #performance on test set
    prep_test = get_features(test, featureCOLUMNS)
    predictions = model.predict(prep_test)

    #statistic on test set
    r2 = np.corrcoef(target_test.flatten(), predictions.flatten())[0,1]**2
    RMSE_z = np.sqrt(np.mean((target_test.flatten() - predictions.flatten())**2))
    RMSE_n = np.sqrt(np.mean((target_test.flatten() - np.mean(target_test.flatten()))**2))
    R2 = 1 - RMSE_z**2/RMSE_n**2
    MAE = mean_absolute_error(target_test.flatten(), predictions.flatten())
    max_err = max_error(target_test.flatten(), predictions.flatten())
    mse = mean_squared_error(target_test.flatten(), predictions.flatten())
    rmse = np.sqrt(mse)

    #print/plot results
    print('Performance on test set(r2, R2, MAE, RMSE, Maximal Error, MSE):',r2,';',R2,';',MAE,';',rmse,';',max_err,';', mse)

    r2 = np.corrcoef(predictions.flatten(), target_test.flatten())[0,1]**2
    plot_2D(['$r^2$ = ' + str(f"{r2:.2f}")], 'upper left', predictions , target_test,
            'predicted', 'experimental', working_dir + '20-80-split-true-pred.png', '#A43341', 
            include_line=False, line_style='None')
    
    #save model training results
    pickle.dump(model, open(working_dir + "model.pkl", "wb"))

def add_predictions(data_MMPs, feature_function, model):
    """
    Compute model predictions for the two SMILES columns of an MMP DataFrame.

    Parameters:
        data_MMPs (pd.DataFrame): Must contain 'smiles_1' and 'smiles_2' columns.
        feature_function (callable): Converts a SMILES string to a feature vector.
        model: Fitted estimator with a .predict(X) method.

    Returns:
        pd.DataFrame: The input DataFrame augmented with:
            - 'Feature_1' and 'predictions_1' for smiles_1
            - 'Feature_2' and 'predictions_2' for smiles_2
    """
    data_MMPs['Feature_1'] = data_MMPs['smiles_1'].apply(else_none(feature_function))
    data_MMPs['Feature_2'] = data_MMPs['smiles_2'].apply(else_none(feature_function))

    data_MMPs = data_MMPs[(~data_MMPs["Feature_1"].isna()) & (~data_MMPs["Feature_2"].isna())]
    X_data_attributions_1 = get_features(data_MMPs, ['Feature_1'])
    predictions_1 = model.predict(X_data_attributions_1)
    data_MMPs['predictions_1'] = predictions_1

    X_data_attributions_2 = get_features(data_MMPs, ['Feature_2'])
    predictions_2 = model.predict(X_data_attributions_2)
    data_MMPs['predictions_2'] = predictions_2

    return data_MMPs

def split_MMPs_by_set(data_MMPs, test):
    """
    Split an MMP DataFrame into pairs where both molecules are in the train or test split.

    Parameters:
        data_MMPs (pd.DataFrame): Must have 'smiles_1' and 'smiles_2' columns.
        test (pd.DataFrame): Test‐set DataFrame with a 'smiles_std' column.

    Returns:
        train_set (pd.DataFrame): Rows where both smiles_1 and smiles_2 are outside the test set.
        test_set  (pd.DataFrame): Rows where both smiles_1 and smiles_2 are in the test set.
    """
    data_MMPs['set_1'] = data_MMPs['smiles_1'].isin(test['smiles_std']).map({True: 'test', False: 'train'})
    data_MMPs['set_2'] = data_MMPs['smiles_2'].isin(test['smiles_std']).map({True: 'test', False: 'train'})

    train_set = data_MMPs[(data_MMPs['set_1'] == 'train') & (data_MMPs['set_2'] == 'train')]
    test_set = data_MMPs[(data_MMPs['set_1'] == 'test') & (data_MMPs['set_2'] == 'test')]

    return train_set, test_set
