import inspect
import os
import sys
import os
import contextlib

from standardizer.mol_standardizer import *
from WISP.ml_helper import *
from WISP.SHAP_MORGAN_attributor import *
from WISP.plotting_helper import *
from WISP.atom_attributor import *
from WISP.dropout_attributor import *
from WISP.path_attributor import *
from WISP.RDKit_attributor import *
from WISP.create_MMPs import *
from WISP.get_indices import *


def analyze_locality(data_MMPs, Attribution_Colums, working_dir):
    data_MMPs['unmatched_atom_index_1_with_neighbors'] = data_MMPs.apply(lambda row: get_neighbors(row['smiles_1'], row['unmatched_atom_index_1']), axis=1)
    data_MMPs['unmatched_atom_index_2_with_neighbors'] = data_MMPs.apply(lambda row: get_neighbors(row['smiles_2'], row['unmatched_atom_index_2']), axis=1)

    methods = {
    'Atom': ('Atom Attributions_1_fix', 'Atom Attributions_2_fix', '#10384f', 'ATOM')}

    if 'SHAP Attributions' in Attribution_Colums:
        methods['SHAP'] = ('SHAP Attributions_1_fix', 'SHAP Attributions_2_fix', '#9C0D38', 'SHAP')
    if 'RDKit Attributions' in Attribution_Colums:
        methods['RDKit'] = ('RDKit Attributions_1_fix', 'RDKit Attributions_2_fix', '#758ECD', 'RDKIT')

    for method, (attr1, attr2, color, label) in methods.items():
        key1 = f'summ_unmatched_{method}_contributions_1_sphere'
        key2 = f'summ_unmatched_{method}_contributions_2_sphere'
        delta_key = f'delta_sum_fragment_{method}_contributions_sphere'

        data_MMPs[key1] = get_unmatched_attributions(data_MMPs, attr1, 'unmatched_atom_index_1_with_neighbors')
        data_MMPs[key2] = get_unmatched_attributions(data_MMPs, attr2, 'unmatched_atom_index_2_with_neighbors')
        data_MMPs[delta_key] = data_MMPs[key1] - data_MMPs[key2]

        r2 = np.corrcoef(data_MMPs['delta_prediction'], data_MMPs[delta_key])[0,1]**2

        plot_2D(
            [f'$r^2$({method}) = {r2:.2f}'], 'upper left',
            data_MMPs['delta_prediction'],
            data_MMPs[delta_key],
            '$\Delta$Predictions MMP',
            '$\Delta$Attributions MMP (Fragment+Neig)',
            working_dir + label + 'attributor_NEIG.png',
            color,
            include_line=True,
            line_style='None'
        )
    return data_MMPs

def detect_binary_classification(data, Target_Column_Name):
    unique_values = data[Target_Column_Name].nunique(dropna=False)
    return 'classification' if unique_values == 2 else 'regression'

def MMP_accuracy(data, attr_method):
    y_true = (data['delta_prediction'] > 0).astype(int)
    y_pred = (data[attr_method] > 0).astype(int)

    accuracy = accuracy_score(y_true, y_pred)
    print(attr_method + f" Accuracy: {accuracy:.2f}")

@contextlib.contextmanager
def suppress_output():
    with open(os.devnull, 'w') as devnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        sys.stdout = devnull
        sys.stderr = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr

def WISP(working_dir, input_dir, ID_Column_Name, Smiles_Column_Name, Target_Column_Name, model_available=False):
    '''
    If you have your model available please place it in the working_dir as model.pkl
    Input as comma seperated file
    '''
    #Interactive questions
    if model_available is not None:
        print('Please provide the name of the function to create the features based on smiles as input:')
        function_name = input().strip()
        feature_function = globals()[function_name]
        print(function_name)

    #Load Data
    data = pd.read_csv(input_dir)
    data.rename(columns={ID_Column_Name: 'ID'}, inplace=True)

    #set type
    task_type = detect_binary_classification(data, Target_Column_Name)   

    #Preprocessing/Standadizing
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

    print("Standardization done")

    if model_available is None:

        if task_type == 'regression':

            #Calculate Fingerprints as Descriptors and set settings
            data, ALLfeatureCOLUMS, model_types = features_and_reg_model_types(data)

            #test/train split
            test, target_test, train = split_data(data, Target_Column_Name, working_dir)
            
            #find and train best model
            model, feature_function, featureCOLUMS  = get_best_reg_model(model_types, ALLfeatureCOLUMS, train, Target_Column_Name, working_dir)

            train_and_evaluate_reg_model(model, train, test, featureCOLUMS, Target_Column_Name, target_test, working_dir)

        if task_type == 'classification':
            
            #Calculate Fingerprints as Descriptors and set settings
            data, ALLfeatureCOLUMS, model_types = features_and_class_model_types(data)

            #test/train split
            test, target_test, train = split_data(data, Target_Column_Name, working_dir)

            #find and train standard model ---------------------> could be extended
            model, feature_function, featureCOLUMS  = get_and_train_class_model(train, test, Target_Column_Name, target_test, working_dir)


    #load the provided or just trained model
    model = pickle.load(open(working_dir + "model.pkl", 'rb'))

    #Attribute Atoms
    Attribution_Colums = ['Atom Attributions']#, 'Path Attributions', 'Dropout Attributions'
    color_coding =['#10384f'] #,'#89d329','#00bcff'
    
    #model/descriptor agnostic
    data['Atom Attributions'] = data['smiles_std'].apply(lambda s: attribute_atoms(s, model, feature_function))
    #data['Dropout Attributions'] = data['smiles_std'].apply(lambda s: attribute_atoms_dropout(s, model, feature_function))
    #data['Path Attributions'] = data['smiles_std'].apply(lambda s: attribute_atoms_paths(s, model, feature_function))

    print("XSMILES Attribution done")
    
    #SHAP explainer
    if "Morgan" in inspect.getsource(feature_function):#to only use on Morgan
        if task_type == 'regression': #need to get classification to work 
            data['Morgan_Fingerprint 2048Bit 2rad'] = data['smiles_std'].apply(feature_function)
            
            explainer = pick_shap_explainer(model)
            
            data = get_SHAP_Morgan_attributions(data, 'Morgan_Fingerprint 2048Bit 2rad', 'smiles_std', model, explainer)
            Attribution_Colums.append('SHAP Attributions')
            color_coding.append('#9C0D38')

            print("SHAP Attribution done")
    
    #RDKit
    if "Morgan" in inspect.getsource(feature_function):
        data['RDKit Attributions'] = data['smiles_std'].apply(lambda s: RDKit_attributor(s, SimilarityMaps.GetMorganFingerprint, model))
        Attribution_Colums.append('RDKit Attributions')
        color_coding.append('#758ECD')
        print("RDKit Attribution done")
    if "RDK" in inspect.getsource(feature_function):
        def fp_func(m, a):
            return SimilarityMaps.GetRDKFingerprint(m, atomId=a, maxPath=7)
        data['RDKit Attributions'] = data['smiles_std'].apply(lambda s: RDKit_attributor(s, fp_func, model))
        Attribution_Colums.append('RDKit Attributions')
        color_coding.append('#758ECD')
        print("RDKit Attribution done")

    #create heatmaps
    directory = working_dir + "HeatMaps"
    if not os.path.exists(directory):
        os.makedirs(directory)

    for index, row in data.iterrows():
        for attr_method in Attribution_Colums:
            output_dir = directory + '/'
            generate_heatmap(data, index, output_dir, 'smiles_std', attr_method, 'ID', task_type)
    print("Heatmaps have been created")

    #save attribution data
    data.to_csv(working_dir + "Attribution_Data.csv", index=False)

    #Creat MMP database
    colums_to_keep = Attribution_Colums + [Target_Column_Name]
    data_MMPs = create_MMP_database(working_dir + "data.smi", working_dir ,data, colums_to_keep)
    data_MMPs.to_csv(working_dir + "MMPs_with_attributions.csv", index=False)
    print("MMPs created")

    #add predictions
    data_MMPs = add_predictions(data_MMPs, feature_function, model)

    #Add indices
    data_MMPs[["unmatched_atom_index_1", "unmatched_atom_index_2"]] = data_MMPs.apply(
    lambda row: pd.Series(get_unmatched_atom_indices_fragments(row["smiles_1"], row["smiles_2"], row["constant"])), axis=1)

    #Add Plots   
    if task_type == 'regression':
        if model_available is not None:
            for attr_method, color in zip(Attribution_Colums, color_coding): 
                data_MMPs = plot_MMP_correlations(data_MMPs, attr_method, color, working_dir, Target_Column_Name)
                data_MMPs = plot_const_histogram(data_MMPs, attr_method, color, working_dir)
            
            #Analyse Locality
            data_MMPs = analyze_locality(data_MMPs, Attribution_Colums, working_dir)

    #test/train dependency and plots
    if model_available is None:

        train_set, test_set = split_MMPs_by_set(data_MMPs, test)

        for attr_method, color in zip(Attribution_Colums, color_coding): 

            print('For the training set:')
            train_set = plot_MMP_correlations(train_set, attr_method, color, working_dir, Target_Column_Name, header='Training Set')
            train_set = plot_const_histogram(train_set, attr_method, color, working_dir, header='Training Set')

            print('For the test set:')
            test_set = plot_MMP_correlations(test_set, attr_method, color, working_dir, Target_Column_Name, header='Test Set')
            test_set = plot_const_histogram(test_set, attr_method, color, working_dir, header='Test Set')
        
        #Analyse Locality
        print('For the training set:')
        train_set = analyze_locality(train_set, Attribution_Colums, working_dir)
        print('For the test set:')
        test_set = analyze_locality(test_set, Attribution_Colums, working_dir)
        
        for attr_method, color in zip(Attribution_Colums, color_coding): 
        
            columns = ['delta_sum_' + attr_method, 'delta_sum_fragment_contributions_' + attr_method]
            for col in columns:
                print('For the training set:')
                MMP_accuracy(train_set, col)
                print('For the test set:')
                MMP_accuracy(test_set, col)


    #get MMP accuracy
    if model_available is not None:
        for attr_method in Attribution_Colums:
            
            if task_type == 'classification':
                #to get the sums also for the regression part
                r2_whole, r2_fragment, data_MMPs = get_r2_and_summed_data_attributions(data_MMPs, 'predictions_1', 'predictions_2', attr_method + '_1', attr_method + '_2', attr_method, 'unmatched_atom_index_1' , 'unmatched_atom_index_2')
            
            columns = ['delta_sum_' + attr_method, 'delta_sum_fragment_contributions_' + attr_method]
            for col in columns:
                MMP_accuracy(data_MMPs, col)
        
    
    #save data
    data_MMPs.to_csv(working_dir + "Complete_Data.csv", index=False)

    print("WISP done")

    return data_MMPs