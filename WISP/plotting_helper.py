import pandas as pd
import math
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps

from IPython.display import Image

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.colorbar import ColorbarBase
import matplotlib.colors as mcolors

import os
import ast

from WISP.get_indices import *


DISPLAY_PLOTS = True

def plot_2D(LEGEND, placeLegend, Xaxes, Yaxes, Xlabel, Ylabel, SAVEDIR, color, header=None, include_line=True, line_style='-'):
    """
    Create and save a 2D scatter plot (with optional unity line), legend, and title.

    Parameters:
        LEGEND (list of str):
            Labels to display in the legend.
        placeLegend (str):
            Location string for plt.legend (e.g. 'upper left').
        Xaxes (array‐like):
            X coordinates for the scatter points.
        Yaxes (array‐like):
            Y coordinates for the scatter points.
        Xlabel (str):
            Label for the X axis.
        Ylabel (str):
            Label for the Y axis.
        SAVEDIR (str):
            Full file path where the figure will be saved.
        color (str):
            Matplotlib color code for the points.
        header (str, optional):
            Title text; if provided, also appended to the saved filename.
        include_line (bool, default=True):
            If True, draw the y=x unity line.
        line_style (str, default='-'):
            Line style for points and/or unity line.

    Returns:
        None  (saves the figure to `SAVEDIR` and optionally displays it)
    """
    plt.rcParams['font.family'] = 'arial'
    plt.plot(Xaxes, Yaxes, '.', color=color, markersize=10, linestyle=line_style, alpha=0.5)
    plt.xlabel(Xlabel, fontsize=18)
    plt.ylabel(Ylabel, fontsize=18)
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
    if DISPLAY_PLOTS:
        plt.show()
    else:
        plt.clf()

def plot_histogram(dataset1, dataset2, column_of_interest, x_Label, y_Label, dataset1_name,dataset2_name,  SAVEDIR):
    '''
    Plots a histogram on a property on interest in bins on 10.

    Keyword arguments:
    -- dataset1: first dataset
    -- dataset2: second dataset
    -- column_of_interest: Column with the data on which the histogram should be plotted (needs to be present in both datasets).
    -- x_Label: Label which should be displayed on the y-axis.
    -- y_Label: Label which should be displayed on the y-axis.
    -- dataset1_name: name of the first dataset
    -- dataset2_name: name of the second dataset
    -- SAVEDIR: Directory where the plot should be stored.
    '''
    dataset1['data'] = 'dataset1'
    dataset2['data'] = 'dataset2'
    data = pd.concat([dataset1, dataset2], ignore_index=True)

    bins = np.arange(math.floor(min(data[column_of_interest])), math.ceil(max(data[column_of_interest])) + 1, 1)

    dataset1_counts, _ = np.histogram(dataset1[column_of_interest], bins=bins)
    dataset2_counts, _ = np.histogram(dataset2[column_of_interest], bins=bins)

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
    if DISPLAY_PLOTS:
        plt.show()
    else:
        plt.clf()

def plot_histogram_one_dataset(data, column_of_interest, label, color, attr_method, save_dir, header):
    """
    Plot and save a histogram of a single data column with annotated standard deviation.

    Parameters:
        data (pd.DataFrame):
            DataFrame containing the column to plot.
        column_of_interest (str):
            Name of the numeric column to histogram.
        label (str):
            X-axis label for the plot.
        color (str):
            Matplotlib color code for the bars.
        attr_method (str):
            Name of the attribution method, used in the legend and filename.
        save_dir (str):
            Directory path where the PNG file will be saved.
        header (str, optional):
            Title for the plot. If None, defaults to 'plot'.

    Returns:
        None  (the figure is saved to disk and optionally displayed)
    """
    std_dev = np.std(data[column_of_interest])

    base_color = mcolors.to_rgba(color) 
    color_with_alpha = base_color[:3] + (0.5,)

    plt.hist(data[column_of_interest], bins=30, color=color_with_alpha, edgecolor=color, label=f"{attr_method} \n(Standard Deviation: {std_dev:.2f})")  
    
    if header is not None:
        plt.title(header, fontsize=20)
    else:
        header = "plot"

    plt.rcParams['font.family'] = 'arial'
    plt.xlabel(label, fontsize=18)
    plt.ylabel('MMP Count', fontsize=18)
    plt.legend(fontsize=16, loc='upper left')

    filename = f"{attr_method.replace(' ', '_').lower()}_histogram_{header.replace(' ', '_').lower()}.png"
    save_path = os.path.join(save_dir, filename)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')

    if DISPLAY_PLOTS:
        plt.show()
    else:
        plt.clf()

def replace_nan_with_zero(vector_str):
    return [0 if pd.isna(x) else x for x in vector_str]

def get_unmatched_attributions(data, attributions, indices):
    """
    Sum the atom‐level attribution values for the given index lists.

    Parameters:
        data (pd.DataFrame): Must contain two columns:
            - attributions: each row is an array‐like of atom scores.
            - indices: each row is a list of atom indices to sum.
        attributions (str): Column name of the attribution arrays.
        indices (str): Column name of the index lists.

    Returns:
        list of float: For each row, the sum of attributions at the specified indices.
    """
    sum_weights = []

    for idx, weight in zip(data[indices], data[attributions]):
        total_weight = sum(weight[i] for i in idx)  
        sum_weights.append(total_weight)

    return sum_weights

def get_r2_and_summed_data_attributions(data, pred_column_1, pred_column_2, attribution_column_1, attribution_column_2, attribution_type, unmatched_atom_index_column_1, unmatched_atom_index_column_2):
    """
    Compute delta‐predictions, sum total and fragment attributions, and return r² metrics.

    Parameters:
        data (pd.DataFrame):
            Must contain prediction columns and raw attribution lists.
        pred_column_1 (str):
            Column name for predictions of molecule 1.
        pred_column_2 (str):
            Column name for predictions of molecule 2.
        attribution_column_1 (str):
            Column name of raw attributions for molecule 1 (list per row).
        attribution_column_2 (str):
            Column name of raw attributions for molecule 2.
        attribution_type (str):
            Suffix used in new delta and fragment‐sum column names.
        unmatched_atom_index_column_1 (str):
            Column name containing atom‐index lists for molecule 1.
        unmatched_atom_index_column_2 (str):
            Column name for molecule 2 index lists.

    Returns:
        r2_whole_mol (float):
            r² between Δpredictions and Δtotal attributions.
        r2_fragment (float):
            r² between Δpredictions and Δfragment-only attributions.
        data (pd.DataFrame):
            Input DataFrame augmented with columns:
              - 'delta_prediction'
              - '{attribution_column}_fix' (lists with NaNs replaced)
              - 'sum_{attribution_column}'
              - 'delta_sum_{attribution_type}'
              - 'sum_unmatched_contributions_1_{attribution_type}'
              - 'sum_unmatched_contributions_2_{attribution_type}'
              - 'delta_sum_fragment_contributions_{attribution_type}'
    """
    data['delta_prediction'] = data[pred_column_1] - data[pred_column_2]

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
    """
    Compute r² between Δpredictions and Δ“constant” fragment attributions.

    Parameters:
        data (pd.DataFrame): MMP DataFrame with prediction and attribution‐fix columns.
        pred_column_1 (str): Column name for predictions of molecule 1.
        pred_column_2 (str): Column name for predictions of molecule 2.
        attribution_column_1 (str): Base name of the “fix” attribution column for mol 1.
        attribution_column_2 (str): Base name of the “fix” attribution column for mol 2.
        attribution_type (str): Suffix used in output column names.
        unmatched_atom_index_column_1 (str): Col name of atom‐index lists for mol 1.
        unmatched_atom_index_column_2 (str): Col name of atom‐index lists for mol 2.

    Returns:
        r2_const (float): Squared Pearson correlation of Δpred vs. Δconstant attributions.
        data (pd.DataFrame): Same DataFrame with added columns:
            - sum_constant_contributions_1_{type}
            - sum_constant_contributions_2_{type}
            - delta_sum_const_contributions_{type}
    """
    data['delta_prediction'] = data[pred_column_1] - data[pred_column_2]

    data['sum_constant_contributions_1_' + attribution_type] = get_unmatched_attributions(data, attribution_column_1 + "_fix", unmatched_atom_index_column_1)
    data['sum_constant_contributions_2_' + attribution_type] = get_unmatched_attributions(data, attribution_column_2 + "_fix", unmatched_atom_index_column_2)

    data['delta_sum_const_contributions_' + attribution_type] = data['sum_constant_contributions_1_' + attribution_type] - data['sum_constant_contributions_2_' + attribution_type]

    r2_const = np.corrcoef(data['delta_prediction'], data['delta_sum_const_contributions_' + attribution_type])[0,1]**2

    return r2_const, data

def plot_MMP_correlations(data_MMPs, attr_method, color, working_dir, Target_Column_Name, header=None):
    """
    Plot and save 2D correlations for MMP prediction vs. attributions and experimental Δ.

    Parameters:
        data_MMPs (pd.DataFrame):
            Must contain columns
              - 'predictions_1', 'predictions_2'
              - '{attr_method}_1' and '{attr_method}_2'
              - 'unmatched_atom_index_1', 'unmatched_atom_index_2'
              - '{target_col}_1', '{target_col}_2'
        attr_method (str): Base name of the attribution column (no '_1'/'_2').
        color (str): Color code for the scatter points.
        working_dir (str): Directory path where plots will be saved.
        target_col (str): Base name of the reference/experimental target column.
        header (str, optional): Title to add to each plot and filename.

    Returns:
        pd.DataFrame: The same `data_MMPs` with 'delta_target' added.
    """
    r2_whole, r2_fragment, data_MMPs = get_r2_and_summed_data_attributions(data_MMPs, 'predictions_1', 'predictions_2', attr_method + '_1', attr_method + '_2', attr_method, 'unmatched_atom_index_1' , 'unmatched_atom_index_2')

    plot_2D(['$r^2$(' + attr_method + ') = ' + str(f"{r2_whole:.2f}")], 'upper left', data_MMPs['delta_prediction'], data_MMPs['delta_sum_' + attr_method],
            '$\Delta$Predictions MMP', '$\Delta$Attributions MMP (whole Mol)', working_dir + 'PREDvsCONTRIBUTIONSwhole' + attr_method + '.png', 
            color, header,
            include_line=False, line_style='None')

    plot_2D(['$r^2$(' + attr_method + ') = ' + str(f"{r2_fragment:.2f}")], 'upper left', data_MMPs['delta_prediction'], data_MMPs['delta_sum_fragment_contributions_' + attr_method],
            '$\Delta$Predictions MMP', '$\Delta$Attributions MMP (Variable)', working_dir + 'PREDvsCONTRIBUTIONSfragment' + attr_method + '.png', 
            color, header,
            include_line=False, line_style='None')


    data_MMPs['delta_target'] = data_MMPs[Target_Column_Name + '_1'] - data_MMPs[Target_Column_Name + '_2']

    r2 = np.corrcoef(data_MMPs['delta_target'], data_MMPs['delta_sum_' + attr_method])[0,1]**2

    plot_2D(['$r^2$(' + attr_method + ') = ' + str(f"{r2:.2f}")], 'upper left', data_MMPs['delta_target'], data_MMPs['delta_sum_' + attr_method],
            '$\Delta$Reference MMP', '$\Delta$Attributions MMP (whole Mol)', working_dir + 'EXPvsCONTRIBUTIONSwhole' + attr_method + '.png', 
            color, header,
            include_line=False, line_style='None')
    
    return data_MMPs

def plot_const_histogram(data_MMPs, attr_method, color, working_dir, header=None):
    """
    Compute “constant” atom indices for each MMP pair and plot a histogram
    of their summed attributions.

    This function:
      - Finds atoms *not* in the unmatched fragment for both molecules.
      - Uses get_r2_and_summed_data_attributions_const to compute
        `delta_sum_const_contributions_{attr_method}`.
      - Calls plot_histogram_one_dataset to save a histogram of those deltas.

    Parameters:
        data_MMPs (pd.DataFrame):
            Must contain 'smiles_1', 'smiles_2', and
            'unmatched_atom_index_1/2' columns.
        attr_method (str):
            Base name of the attribution column (without '_fix').
        color (str):
            Matplotlib color code for the histogram bars.
        working_dir (str):
            Directory where the histogram PNG will be saved.
        header (str, optional):
            Plot title and part of the filename.

    Returns:
        pd.DataFrame: The original DataFrame augmented with:
            - 'const_indices_1', 'const_indices_2'
            - 'delta_sum_const_contributions_{attr_method}'
    """
    data_MMPs['const_indices_1'] = data_MMPs.apply(lambda row: get_unselected_atom_indices(row['smiles_1'], row['unmatched_atom_index_1']), axis=1)
    data_MMPs['const_indices_2'] = data_MMPs.apply(lambda row: get_unselected_atom_indices(row['smiles_2'], row['unmatched_atom_index_2']), axis=1)

    _, data_MMPs = get_r2_and_summed_data_attributions_const(data_MMPs, 'predictions_1', 'predictions_2', attr_method + '_1', attr_method + '_2', attr_method, 'const_indices_1' , 'const_indices_2')

    plot_histogram_one_dataset(data_MMPs, 'delta_sum_const_contributions_' + attr_method, '$\Delta$Attributions MMP (Constant)', color, attr_method, working_dir, header)
    
    return data_MMPs

def generate_heatmap(data, index, output_dir, smiles_column, attribution_column, ID_column, task_type):
    """
    Generate and save an atom‐level heatmap and colorbar for a single molecule.

    For the row at `data.loc[index]`, this function:
      1. Extracts the SMILES string and attribution list.
      2. Replaces NaN values with 0 and casts attributions to float.
      3. Determines vmax (the color scale max) based on
         - 70th percentile of absolute attributions for regression, or
         - a fixed 0.7 for classification.
      4. Draws the molecule with RDKit’s SimilarityMap with a blue–white–red scale.
      5. Saves the horizontal colorbar legend as `<attribution_column>_Legend.png`.
      6. Saves the rendered molecule heatmap as `<ID>_<attribution_column>.png`.
      7. Returns an IPython.display.Image of the heatmap file.

    Parameters:
        data (pd.DataFrame): Must contain columns
            `smiles_column` (SMILES string),
            `attribution_column` (list of floats per atom),
            and `id_column` (unique identifier).
        index (int): Row index in `data` to plot.
        output_dir (str): Directory where PNGs will be written.
        smiles_column (str): Name of the SMILES column in `data`.
        attribution_column (str): Name of the attribution‐list column in `data`.
        id_column (str): Name of the ID column for naming the output file.
        task_type (str): Either 'regression' or 'classification'; controls color scale.

    Returns:
        IPython.display.Image: Embedded image of the saved molecule heatmap.

    Notes:
        - Requires imports:
            from rdkit import Chem, Draw
            from rdkit.Chem.SimilarityMaps import GetSimilarityMapFromWeightsWithScale
            import matplotlib.pyplot as plt
            from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
            from matplotlib.colorbar import ColorbarBase
            from IPython.display import Image
            import numpy as np, math, os
        - If the molecule has fewer than 2 atoms, the function returns None.
    """
    smiles = data[smiles_column][index]
    attributions = data[attribution_column][index]
    attributions = [0 if isinstance(x, float) and math.isnan(x) else x for x in attributions]
    attributions = [float(x) for x in attributions]
    mol_id = data[ID_column][index]

    if task_type == 'regression':
        vmax = np.percentile(data[attribution_column].apply(lambda lst: max(map(abs, lst))),70)
    if task_type == 'classification':
        vmax = 0.7
    vmin = -vmax

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol)#to keep the explicit hydrogens
    
    if mol.GetNumAtoms() < 2: # to avoid errors later
        return

    draw2d = Draw.MolDraw2DCairo(300, 300)
    d = GetSimilarityMapFromWeightsWithScale(mol, attributions, draw2d, '#10384f', '#9C0D38', vmin, vmax)
    d.FinishDrawing()

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

    filename = f"{mol_id}_{attribution_column}.png"
    output_path = os.path.join(output_dir, filename)
    with open(output_path, "wb") as f:
        f.write(d.GetDrawingText())

    return Image(filename=output_path)


def GetSimilarityMapFromWeightsWithScale(mol, weights, draw2d, hex_color1, hex_color2, overall_min:float, overall_max:float, colorMap=None, *args,**kwargs, ):
    """
    Create a similarity‐map drawing with a dynamic blue–white–red color scale
    based on clipped weights and global min/max.

    Parameters:
        mol (rdkit.Chem.Mol):
            The molecule to render.
        weights (iterable of float):
            Atom‐level scores to map onto colors.
        draw2d (MolDraw2D):
            An RDKit 2D drawing canvas (e.g. MolDraw2DCairo).
        hex_color1 (str):
            Hex code for the negative‐end color (e.g. blue).
        hex_color2 (str):
            Hex code for the positive‐end color (e.g. red).
        overall_min (float):
            Global minimum weight used for clipping and scaling.
        overall_max (float):
            Global maximum weight used for clipping and scaling.
        colorMap:
            Must be None (custom colormap is built internally).
        *args, **kwargs:
            Forwarded to SimilarityMaps.GetSimilarityMapFromWeights.

    Returns:
        The object returned by RDKit’s GetSimilarityMapFromWeights,
        with the drawing ready on the provided canvas.

    """
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
    Map a single weight to an RGB color using a diverging two‐slope colormap.

    Parameters:
        weight (float):
            The value to map (can be negative, zero, or positive).
        neg_color (str or tuple):
            Color at vmin (e.g. '#10384f' for negative).
        neut_color (str or tuple):
            Color at center (zero).
        pos_color (str or tuple):
            Color at vmax (e.g. '#9C0D38' for positive).
        vmin (float):
            Minimum of the color scale.
        vmax (float):
            Maximum of the color scale.

    Returns:
        rgb (tuple of float):
            The mapped color as an (R, G, B) tuple in [0,1].
        cmap (LinearSegmentedColormap):
            The colormap instance used for the mapping.
    """
    norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    colors = [neg_color, neut_color, pos_color]
    cmap = LinearSegmentedColormap.from_list("custom_cmap", colors)
    rgba = cmap(norm(weight))
    return tuple(rgba[:3]), cmap