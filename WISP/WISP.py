import inspect
import multiprocessing
import os
import sys
import os
import contextlib

import warnings

warnings.filterwarnings("ignore", module="lightning.pytorch")
warnings.filterwarnings("ignore", module="sklearn")


from WISP.log_step import LogStep, suppress_output
from standardizer.mol_standardizer import *
from WISP.ml_helper import *
from WISP.SHAP_MORGAN_attributor import *
from WISP.plotting_helper import *
from WISP.RDKit_attributor import *
from WISP.create_MMPs import *
from WISP.get_indices import *

def detect_binary_classification(data, Target_Column_Name):
    """
    Determine if the specified target/property of interest column defines a binary classification task or regression.

    Parameters:
        data (pd.DataFrame): The input DataFrame.
        Target_Column_Name (str): Name of the target column to inspect.

    Returns:
        str: 'classification' if the column has exactly two unique values (NaNs counted),
             otherwise 'regression'.
    """
    unique_values = data[Target_Column_Name].nunique(dropna=False)
    return 'classification' if unique_values == 2 else 'regression'

def MMP_accuracy(data, attr_method):
    """
    Compute and print the binary accuracy of an attribution method versus prediction deltas.

    Parameters:
        data (pd.DataFrame): Must contain 'delta_prediction' and the column named by attr_method.
        attr_method (str): Column name of attributions to threshold at zero.

    Returns:
        None: Prints "<attr_method> Accuracy: XX.XX".
    """
    y_true = (data['delta_prediction'] > 0).astype(int)
    y_pred = (data[attr_method] > 0).astype(int)

    accuracy = accuracy_score(y_true, y_pred)
    print(attr_method + f" Accuracy: {accuracy:.2f}")

def normalize_atom_attributions(data, input_col):
    """
    Normalize per-atom attribution arrays by the global standard deviation of the respective dataset.

    Parameters:
        data (pd.DataFrame): Contains column `input_col` of array-like attributions.
        input_col (str): Name of the column to normalize.

    Returns:
        pd.DataFrame: The original DataFrame with a new column
                      `input_col + '_std'` holding the normalized arrays.
    """
    #output_col = input_col + '_std'

    attr = np.concatenate(data[input_col].values)
    std = np.nanstd(attr)

    data[input_col] = data[input_col].apply(lambda x: np.array(x) / std)
    
    return data






def WISP(working_dir, input_dir, ID_Column_Name, Smiles_Column_Name, Target_Column_Name, model_available=None, use_GNN=True, fast_run=False,):
    """
    Executes the main WISP workflow: Load data, standardize SMILES, train or load a model, compute atom‐level attributions,
    generate heatmaps and matched‐molecular‐pair (MMP) analyses, then return the MMP DataFrame as well as the analyses plots.

    Parameters:
        working_dir (str): Directory for outputs (models, plots, CSVs).
        input_dir (str): Path to the input CSV file.
        ID_Column_Name (str): Name of the ID column in the CSV.
        Smiles_Column_Name (str): Name of the SMILES column in the CSV.
        Target_Column_Name (str): Name of the target column (regression or binary class).
        model_available (str or None): If provided, loads “model.pkl” from working_dir
            and skips training; else trains a new model. In this case you need to provide the feature 
            function name, which is the function that converts the input smiles to the 
            respective machine learning feature.
        use_GNN (bool): If True, include graph‐neural‐network models in the search.
        fast_run (bool): If True, only try MorganFingerprints+RandomForestRegressor/Classifier

    Returns:
        pd.DataFrame: The final MMP DataFrame with predictions, attributions, plots, and metrics.
    """
    #interactive questions
    if model_available is not None:
        print('Please provide the name of the function to create the features based on smiles as input:')
        function_name = input().strip()
        feature_function = globals()[function_name]
        print(function_name)

    #load Data
    data = pd.read_csv(input_dir)
    data.rename(columns={ID_Column_Name: 'ID'}, inplace=True)

    #set type
    task_type = detect_binary_classification(data, Target_Column_Name)   


    with LogStep("Standardization"):
        with suppress_output():
            std = Standardizer(max_num_atoms=1000,#tipp jan: 100
                        max_num_tautomers=10,
                        include_stereoinfo=False,
                        keep_largest_fragment=True, 
                        canonicalize_tautomers=True, 
                        normalize=True, 
                        sanitize_mol=True)
        data["smiles_std"] = data[Smiles_Column_Name].apply(lambda smi: std(smi)[0]) 

    data = filter_duplicates(data, 'smiles_std', Target_Column_Name)

    #save data in smi format
    data[['smiles_std', 'ID', Target_Column_Name]].to_csv(working_dir + "data.smi", sep='\t', index=False, header=False)

    if model_available is None:

        with LogStep("Train Model"):
            if task_type == 'regression':

                #calculate fingerprints as descriptors and set settings
                data, ALLfeatureCOLUMNS, model_types = features_and_reg_model_types(data,fast_run=fast_run)

                #test/train split
                test, target_test, train = split_data(data, Target_Column_Name)
                
                #find and train best model
                model, feature_function, featureCOLUMNS  = get_best_reg_model(model_types, ALLfeatureCOLUMNS, train, Target_Column_Name, working_dir, use_GNN)

                train_and_evaluate_reg_model(model, train, test, featureCOLUMNS, Target_Column_Name, target_test, working_dir)

            if task_type == 'classification':
                
                #calculate fingerprints as descriptors and set settings
                data, ALLfeatureCOLUMNS, model_types = features_and_class_model_types(data,fast_run=fast_run)

                #test/train split
                test, target_test, train = split_data(data, Target_Column_Name)

                #find and train standard model
                model, feature_function, featureCOLUMNS  = get_and_train_class_model(train, test, Target_Column_Name, target_test, working_dir)


    #load the provided or just trained model
    model = pickle.load(open(working_dir + "model.pkl", 'rb'))

    #attribute atoms
    Attribution_Columns = ['Atom Attributions']
    color_coding =['#10384f']


    with LogStep("Calc Atom Attributions"):
        VECTORIZE_AA = True
        if VECTORIZE_AA:
            from WISP.atom_attributor_vec import attribute_atoms
            data["Atom Attributions"] = attribute_atoms(data["smiles_std"].tolist(), model, feature_function)
            data = normalize_atom_attributions(data, 'Atom Attributions')
        else:
            from WISP.atom_attributor import attribute_atoms
            if "identity" in feature_function.__name__:
                unvec_feature_fun = feature_function
            else:
                unvec_feature_fun = lambda smi: feature_function([smi])
            data["Atom Attributions"] = data["smiles_std"].apply(lambda smiles: attribute_atoms(smiles, model, unvec_feature_fun))
            data = normalize_atom_attributions(data, 'Atom Attributions')


    if not fast_run:
        #SHAP explainer
        if "Morgan" in inspect.getsource(feature_function):
            if task_type == 'regression':
                data['Morgan_Fingerprint 2048Bit 2rad'] = data['smiles_std'].apply(feature_function)
                
                explainer = pick_shap_explainer(model, data)
                
                data = get_SHAP_Morgan_attributions(data, 'Morgan_Fingerprint 2048Bit 2rad', 'smiles_std', model, explainer)
                data = normalize_atom_attributions(data, 'SHAP Attributions')
                Attribution_Columns.append('SHAP Attributions')
                color_coding.append('#9C0D38')

                print("SHAP Attribution done")
    
        #RDKit
        if "Morgan" in inspect.getsource(feature_function):
            data['RDKit Attributions'] = data['smiles_std'].apply(lambda s: RDKit_attributor(s, SimilarityMaps.GetMorganFingerprint, model))
            data = normalize_atom_attributions(data, 'RDKit Attributions')
            Attribution_Columns.append('RDKit Attributions')
            color_coding.append('#758ECD')
            print("RDKit Attribution done")

        if "RDK" in inspect.getsource(feature_function):
            def fp_func(m, a):
                return SimilarityMaps.GetRDKFingerprint(m, atomId=a, maxPath=7)
            data['RDKit Attributions'] = data['smiles_std'].apply(lambda s: RDKit_attributor(s, fp_func, model))
            data = normalize_atom_attributions(data, 'RDKit Attributions')
            Attribution_Columns.append('RDKit Attributions')
            color_coding.append('#758ECD')
            print("RDKit Attribution done")

    data.to_csv(working_dir + "Attribution_Data.csv", index=False)

    directory = working_dir + "HeatMaps"
    if not os.path.exists(directory):
        os.makedirs(directory)

    with LogStep("Generating Heatmaps"):
        for index, row in data.iterrows():
            for attr_method in Attribution_Columns:
                output_dir = directory + '/'
                generate_heatmap(data, index, output_dir, 'smiles_std', attr_method, 'ID', task_type)

    with LogStep("Creating MMPs"):
        columns_to_keep = Attribution_Columns + [Target_Column_Name]
        data_MMPs = create_MMP_database(working_dir + "data.smi", working_dir ,data, columns_to_keep)
        data_MMPs.to_csv(working_dir + "MMPs_with_attributions.csv", index=False)

    with LogStep("Predict on MMPs"):
        data_MMPs = add_predictions(data_MMPs, feature_function, model)

    data_MMPs[["unmatched_atom_index_1", "unmatched_atom_index_2"]] = data_MMPs.apply(
    lambda row: pd.Series(get_unmatched_atom_indices_fragments(row["smiles_1"], row["smiles_2"], row["constant"])), axis=1)

    #add Plots   
    if task_type == 'regression':
        if model_available is not None:
            for attr_method, color in zip(Attribution_Columns, color_coding): 
                data_MMPs = plot_MMP_correlations(data_MMPs, attr_method, color, working_dir, Target_Column_Name)
                data_MMPs = plot_const_histogram(data_MMPs, attr_method, color, working_dir)

    #test/train dependency and plots
    if model_available is None:

        train_set, test_set = split_MMPs_by_set(data_MMPs, test)

        for attr_method, color in zip(Attribution_Columns, color_coding): 

            print('For the training set:')
            train_set = plot_MMP_correlations(train_set, attr_method, color, working_dir, Target_Column_Name, header='Training Set')
            train_set = plot_const_histogram(train_set, attr_method, color, working_dir, header='Training Set')

            print('For the test set:')
            test_set = plot_MMP_correlations(test_set, attr_method, color, working_dir, Target_Column_Name, header='Test Set')
            test_set = plot_const_histogram(test_set, attr_method, color, working_dir, header='Test Set')
        
        for attr_method, color in zip(Attribution_Columns, color_coding): 
        
            columns = ['delta_sum_' + attr_method, 'delta_sum_fragment_contributions_' + attr_method]
            for col in columns:
                print('For the training set:')
                MMP_accuracy(train_set, col)
                print('For the test set:')
                MMP_accuracy(test_set, col)
        
        #save data
        train_set.to_csv(working_dir + "Complete_Data_Training.csv", index=False)
        test_set.to_csv(working_dir + "Complete_Data_Test.csv", index=False)


    #get MMP accuracy
    if model_available is not None:
        for attr_method in Attribution_Columns:
            
            if task_type == 'classification':
                #to get the sums also for the regression part
                _, _, data_MMPs = get_r2_and_summed_data_attributions(data_MMPs, 'predictions_1', 'predictions_2', attr_method + '_1', attr_method + '_2', attr_method, 'unmatched_atom_index_1' , 'unmatched_atom_index_2')
            
            columns = ['delta_sum_' + attr_method, 'delta_sum_fragment_contributions_' + attr_method]
            for col in columns:
                MMP_accuracy(data_MMPs, col)
        
    #save data
    data_MMPs.to_csv(working_dir + "Complete_Data.csv", index=False)

    print("WISP done")

    return data_MMPs