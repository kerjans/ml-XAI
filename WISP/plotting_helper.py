import matplotlib.pyplot as plt
import pandas as pd
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from IPython.display import Image
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import os
import ast
import matplotlib.colors as mcolors
import ast
from matplotlib.colorbar import ColorbarBase
from WISP.get_indices import *



def plot_2D(LEGEND, placeLegend, Xaxes, Yaxes, Xlable, Ylable, SAVEDIR, color, header=None, include_line=True, line_style='-'):
    plt.rcParams['font.family'] = 'arial'
    plt.plot(Xaxes, Yaxes, '.', color=color, markersize=10, linestyle=line_style, alpha=0.5)
    plt.xlabel(Xlable, fontsize=18)
    plt.ylabel(Ylable, fontsize=18)
    if include_line:
        plt.plot((np.min(Xaxes),np.max(Xaxes)), (np.min(Yaxes), np.max(Yaxes)), color='gray')
    if header is not None:
        plt.title(header, fontsize=20)
        #for correct saving
        safe_header = header.replace(" ", "_").replace("/", "-")
        directory, filename = os.path.split(SAVEDIR)
        name, ext = os.path.splitext(filename)
        new_filename = f"{name}_{safe_header}{ext}"
        SAVEDIR = os.path.join(directory, new_filename)
    plt.tick_params(labelsize=16)
    plt.legend(LEGEND, fontsize=16, loc=placeLegend)
    plt.savefig(SAVEDIR, dpi=300, bbox_inches='tight')
    plt.show()

def plot_histogram(dataset1, dataset2, colum_of_interest, x_Label, y_Label, dataset1_name,dataset2_name,  SAVEDIR):
    '''
    Plots a histogram on a property on interest in bins on 10.

    Keyword arguments:
    -- dataset1: first dataset
    -- dataset2: second dataset
    -- colum_of_interest: Colum with the data on which the histogram should be plotted (needs to be present in both datasets).
    -- x_Label: Label which should be displayed on the y-axis.
    -- y_Label: Label which should be displayed on the y-axis.
    -- dataset1_name: name of the first dataset
    -- dataset2_name: name of the second dataset
    -- SAVEDIR: Directory where the plot should be stored.
    '''
    dataset1['data'] = 'dataset1'
    dataset2['data'] = 'dataset2'
    data = pd.concat([dataset1, dataset2], ignore_index=True)

    bins = np.arange(math.floor(min(data[colum_of_interest])), math.ceil(max(data[colum_of_interest])) + 1, 1)

    dataset1_counts, _ = np.histogram(dataset1[colum_of_interest], bins=bins)
    dataset2_counts, _ = np.histogram(dataset2[colum_of_interest], bins=bins)

    plt.rcParams['font.family'] = 'arial'
    bar_width = 0.4  # Width of the bars
    plt.bar(bins[:-1] - bar_width / 2, dataset1_counts, width=bar_width, label=dataset1_name, align='center', color='#A43341')
    plt.bar(bins[:-1] + bar_width / 2, dataset2_counts, width=bar_width, label=dataset2_name, align='center', color='#b2b2b2')

    plt.xlabel(x_Label, fontsize=18)
    plt.ylabel(y_Label, fontsize=18)
    plt.legend(fontsize=16)

    bin_labels = [f'{int(b)}-{int(b + 1)}' for b in bins[:-1]]
    plt.xticks(bins[:-1], labels=bin_labels, rotation=90)  #
    plt.tick_params(labelsize=16)
    plt.xticks(bins[:-1], labels=bin_labels)
    plt.savefig(SAVEDIR, dpi=300, bbox_inches='tight')
    plt.show()

def plot_histogram_one_dataset(data, colum_of_interest, label, color, attr_method, save_dir, header):

    #get standard deviation
    std_dev = np.std(data[colum_of_interest])

    # color
    base_color = mcolors.to_rgba(color) 
    color_with_alpha = base_color[:3] + (0.5,)

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.hist(data[colum_of_interest], bins=30, color=color_with_alpha, edgecolor=color, label=f"{attr_method} \n(Standard Deviation: {std_dev:.2f})")  
    
    if header is not None:
        plt.title(header, fontsize=20)
        #to save correctly
    else:
        header = "plot"

    plt.rcParams['font.family'] = 'arial'
    plt.xlabel(label, fontsize=18)
    plt.ylabel('MMP Count', fontsize=18)
    plt.legend(fontsize=16, loc='upper left')

    # Save plot
    filename = f"{attr_method.replace(' ', '_').lower()}_histogram_{header.replace(' ', '_').lower()}.png"
    save_path = os.path.join(save_dir, filename)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

def replace_nan_with_zero(vector_str):
    #vector_str = vector_str.replace('nan', 'np.nan')
    #vector = eval(vector_str, {"__builtins__": None}, {"np": np})
    return [0 if pd.isna(x) else x for x in vector_str]

def get_unmatched_attributions(data, attributions, indices):
    sum_weights = []

    # Loop through both columns and accumulate the sum of the weights for each index
    for idx, weight in zip(data[indices], data[attributions]):
        total_weight = sum(weight[i] for i in idx)  
        sum_weights.append(total_weight)

    return sum_weights

def get_r2_and_summed_data_attributions(data, pred_column_1, pred_column_2, attribution_column_1, attribution_column_2, attribution_type, unmatched_atom_index_column_1, unmatched_atom_index_column_2):
    data['delta_prediction'] = data[pred_column_1] - data[pred_column_2]

    #convert strings to real lists
    data[attribution_column_1 + "_fix"] = data[attribution_column_1].apply(replace_nan_with_zero)
    data[attribution_column_2 + "_fix"] = data[attribution_column_2].apply(replace_nan_with_zero)

    data["sum_" + attribution_column_1] = data[attribution_column_1 + "_fix"].apply(lambda x: sum(x))
    data["sum_" + attribution_column_2] = data[attribution_column_2 + "_fix"].apply(lambda x: sum(x))
    data['delta_sum_' + attribution_type] = data["sum_" + attribution_column_1] - data["sum_" + attribution_column_2]

    r2_whole_mol = np.corrcoef(data['delta_prediction'], data['delta_sum_' + attribution_type])[0,1]**2

    data['sum_unmatched_contributions_1_' + attribution_type] = get_unmatched_attributions(data, attribution_column_1 + "_fix", unmatched_atom_index_column_1)
    data['sum_unmatched_contributions_2_' + attribution_type] = get_unmatched_attributions(data, attribution_column_2 + "_fix", unmatched_atom_index_column_2)

    data['delta_sum_fragment_contributions_' + attribution_type] = data['sum_unmatched_contributions_1_' + attribution_type] - data['sum_unmatched_contributions_2_' + attribution_type]

    r2_fragment = np.corrcoef(data['delta_prediction'], data['delta_sum_fragment_contributions_' + attribution_type])[0,1]**2

    return r2_whole_mol, r2_fragment, data

def get_r2_and_summed_data_attributions_const(data, pred_column_1, pred_column_2, attribution_column_1, attribution_column_2, attribution_type, unmatched_atom_index_column_1, unmatched_atom_index_column_2):
    data['delta_prediction'] = data[pred_column_1] - data[pred_column_2]

    data['sum_constant_contributions_1_' + attribution_type] = get_unmatched_attributions(data, attribution_column_1 + "_fix", unmatched_atom_index_column_1)
    data['sum_constant_contributions_2_' + attribution_type] = get_unmatched_attributions(data, attribution_column_2 + "_fix", unmatched_atom_index_column_2)

    data['delta_sum_const_contributions_' + attribution_type] = data['sum_constant_contributions_1_' + attribution_type] - data['sum_constant_contributions_2_' + attribution_type]

    r2_const = np.corrcoef(data['delta_prediction'], data['delta_sum_const_contributions_' + attribution_type])[0,1]**2

    return r2_const, data

def plot_MMP_correlations(data_MMPs, attr_method, color, working_dir, Target_Column_Name, header=None):

    r2_whole, r2_fragment, data_MMPs = get_r2_and_summed_data_attributions(data_MMPs, 'predictions_1', 'predictions_2', attr_method + '_1', attr_method + '_2', attr_method, 'unmatched_atom_index_1' , 'unmatched_atom_index_2')

    plot_2D(['$r^2$(' + attr_method + ') = ' + str(f"{r2_whole:.2f}")], 'upper left', data_MMPs['delta_prediction'], data_MMPs['delta_sum_' + attr_method],
            '$\Delta$Predictions MMP', '$\Delta$Attributions MMP (whole Mol)', working_dir + 'PREDvsCONTRIBUTIONSwhole' + attr_method + '.png', 
            color, header,
            include_line=False, line_style='None')

    plot_2D(['$r^2$(' + attr_method + ') = ' + str(f"{r2_fragment:.2f}")], 'upper left', data_MMPs['delta_prediction'], data_MMPs['delta_sum_fragment_contributions_' + attr_method],
            '$\Delta$Predictions MMP', '$\Delta$Attributions MMP (Fragment)', working_dir + 'PREDvsCONTRIBUTIONSfragment' + attr_method + '.png', 
            color, header,
            include_line=False, line_style='None')


    data_MMPs['delta_target'] = data_MMPs[Target_Column_Name + '_1'] - data_MMPs[Target_Column_Name + '_2']

    r2 = np.corrcoef(data_MMPs['delta_target'], data_MMPs['delta_sum_' + attr_method])[0,1]**2

    plot_2D(['$r^2$(' + attr_method + ') = ' + str(f"{r2:.2f}")], 'upper left', data_MMPs['delta_target'], data_MMPs['delta_sum_' + attr_method],
            '$\Delta$Target MMP', '$\Delta$Attributions MMP (whole Mol)', working_dir + 'EXPvsCONTRIBUTIONSwhole' + attr_method + '.png', 
            color, header,
            include_line=False, line_style='None')
    
    return data_MMPs

def plot_const_histogram(data_MMPs, attr_method, color, working_dir, header=None):
    data_MMPs['const_indices_1'] = data_MMPs.apply(lambda row: get_unselected_atom_indices(row['smiles_1'], row['unmatched_atom_index_1']), axis=1)
    data_MMPs['const_indices_2'] = data_MMPs.apply(lambda row: get_unselected_atom_indices(row['smiles_2'], row['unmatched_atom_index_2']), axis=1)

    r2_const, data_MMPs = get_r2_and_summed_data_attributions_const(data_MMPs, 'predictions_1', 'predictions_2', attr_method + '_1', attr_method + '_2', attr_method, 'const_indices_1' , 'const_indices_2')

    plot_histogram_one_dataset(data_MMPs, 'delta_sum_const_contributions_' + attr_method, '$\Delta$Attributions MMP (Constant)', color, attr_method, working_dir, header)
    
    return data_MMPs

def generate_heatmap(data, index, output_dir, smiles_column, attribution_column, ID_column, task_type):
    '''
    Red is positive while blue is negative (to be consistant with PIXIE)
    '''
    smiles = data[smiles_column][index]
    attributions = data[attribution_column][index]
    attributions = [0 if isinstance(x, float) and math.isnan(x) else x for x in attributions]
    attributions = [float(x) for x in attributions]
    mol_id = data[ID_column][index]

    if task_type == 'regression':
        vmax = data[attribution_column].apply(lambda lst: max(map(abs, lst))).max() * 0.7
    if task_type == 'classification':
        vmax = 0.7
    vmin = -vmax

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol)#to keep the explicit hydrogens
    
    if mol.GetNumAtoms() < 2: # to avoid errors later
        return

    # Draw similarity map
    draw2d = Draw.MolDraw2DCairo(900, 900)
    d = GetSimilarityMapFromWeightsWithScale(mol, attributions, draw2d, '#10384f', '#9C0D38', vmin, vmax)
    d.FinishDrawing()

    # Create a color bar
    fig, ax = plt.subplots(figsize=(6, 1))
    fig.subplots_adjust(bottom=0.5)

    colors = ['#10384f', '#ffffff', '#9C0D38'] 
    cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    cbar = ColorbarBase(ax, norm=norm, cmap=cmap, orientation='horizontal')

    tick_positions = [vmin, vmax]
    tick_labels = [vmin, vmax]
    cbar.set_ticks(tick_positions)
    cbar.set_ticklabels(tick_labels)

    plt.savefig(os.path.join(output_dir, attribution_column + '_Legend.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Define and write output path
    filename = f"{mol_id}_{attribution_column}.png"
    output_path = os.path.join(output_dir, filename)
    with open(output_path, "wb") as f:
        f.write(d.GetDrawingText())

    return Image(filename=output_path)


def GetSimilarityMapFromWeightsWithScale(mol, weights, draw2d, hex_color1, hex_color2, overall_min:float, overall_max:float, colorMap=None, *args,**kwargs, ):
    assert colorMap is None, "custom colorMap may not be specified"

    weights = list(weights)
    weights = [min(max(w,overall_min),overall_max) for w in weights]

    this_max = max(weights)
    this_min = min(weights)

    color1 = np.array(mcolors.to_rgb(hex_color1))
    color2 = np.array(mcolors.to_rgb(hex_color2))
    white = np.array([1,1,1])
    
    scale_min = this_min / np.array(overall_min)
    scale_max = this_max / np.array(overall_max)
    
    if scale_min < 0:
        scale_min = 0
    if scale_max < 0:
        scale_max = 0
    if scale_min > 1:
        scale_min = 1
    if scale_max > 1:
        scale_max = 1
    
    color1 = color1 * scale_min + white * (1 - scale_min)
    color2 = color2 * scale_max + white * (1 - scale_max)
       
    custom_colors = [color1, '#ffffff', color2]
    custom_cmap = LinearSegmentedColormap.from_list("custom", custom_colors, N=256)
    
    d = SimilarityMaps.GetSimilarityMapFromWeights(mol,list(weights),draw2d,colorMap=custom_cmap,
                                    *args,**kwargs,
                                   )
    return d

def get_color(weight, neg_color, neut_color, pos_color, vmin, vmax):
    """
    Creates a custom colormap from black to red and colors the bonds accordingly.

    Keyword arguments:
    -- weight: Weight of the bond to be colored.
    -- neg_color: Color with which the decrease of the target variable should be described.
    -- neut_color: Color with which no influence of the target variable should be described.
    -- pos_color: Color with which the increase of the target variable should be described.
    -- vmin: Lower bound to normalize the colourcoding on the bond weights.
    -- vmax: Upper bound to normalize the colourcoding on the bond weights.

    Returns:
    -- tuple(rgba[:3]): Colour of the bond.
    -- cmap: Custom color map.
    """
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    colors = [neg_color, neut_color, pos_color] # 'dimgray'
    cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
    rgba = cmap(norm(weight))
    return tuple(rgba[:3]), cmap