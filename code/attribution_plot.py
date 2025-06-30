import numpy as np
from pathlib import Path
import shutil
from smal.io import random_fle
from rdkit.Chem import Draw
from pathlib import Path
from rdkit import rdBase
from rdkit import Chem

rdBase.DisableLog("rdApp.*")

import PIL
from PIL import ImageDraw
from PIL.Image import Image

# from .helper import *

# from .config import ConfigDict, SecretDict

import math
import random
from typing import List
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import io
from PIL import Image
from collections import defaultdict



def tile_images_with_annots(
    imgs: List,
    annots: List[str],
    tile_w: int,
    tile_h: int,
    n_tiles_w: int,
    n_tiles_h: int,
    sampling_strategy="random",
    with_replacement=False,
    annot_color=(0, 0, 0),
) -> Image:
    assert sampling_strategy in ["random", "deterministic"]
    imgs_left = list(range(len(imgs)))
    w_total = n_tiles_w * tile_w
    h_total = n_tiles_h * tile_h
    tile_image = PIL.Image.new("RGB", (w_total, h_total))
    draw = ImageDraw.Draw(tile_image)
    x = 0
    y = 0
    while y < h_total:
        while x < w_total:
            if sampling_strategy == "random" and imgs_left:
                if sampling_strategy == "random":
                    idx = random.choice(imgs_left)
                img = imgs[idx]

            elif sampling_strategy == "deterministic":
                ix = x / tile_w
                iy = y / tile_h
                idx = int(round(iy * n_tiles_w + ix % n_tiles_w))

                if idx < len(imgs_left):
                    img = imgs[idx]
                else:
                    img = PIL.Image.new(
                        "RGB", (16, 16), color="white"
                    )  # will be resized anyways
            else:
                img = PIL.Image.new(
                    "RGB", (16, 16), color="white"
                )  # will be resized anyways
            tile_image.paste(img.resize((tile_w, tile_h)), (x, y))

            if sampling_strategy == "deterministic":
                if idx < len(annots):
                    annot = annots[idx]
                    draw.text((x, y), annot, annot_color)

            if imgs_left and sampling_strategy == "random":
                if idx < len(annots):
                    annot = annots[idx]
                else:
                    annot = ""
                draw.text((x, y), annot, annot_color)

                if not with_replacement:
                    imgs_left = [i for i in imgs_left if i != idx]

            x += tile_w
        x = 0
        y += tile_h

    return tile_image




def GetSimilarityMapFromWeightsWithScale(mol, weights, draw2d, overall_min:float, overall_max:float, colorMap=None, *args,**kwargs, ):
    assert colorMap is None, "custom colorMap may not be specified"

    weights = list(weights)
    weights = [min(max(w,overall_min),overall_max) for w in weights]

    this_max = max(weights)
    this_min = min(weights)

    
    
    
    green = np.array([0,0,1,]) 
    red = np.array([1,0,0,])
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
    
    green = green * scale_min + white * (1 - scale_min)
    red = red * scale_max + white * (1 - scale_max)
    
    
    
    custom_colors = [green, '#ffffff', red]
    custom_cmap = LinearSegmentedColormap.from_list("custom", custom_colors, N=256)
    
    d = GetSimilarityMapFromWeights(mol,list(weights),draw2d,colorMap=custom_cmap,
                                    *args,**kwargs,
                                   )
    return d


def GridSimilarityMaps(mols, weights, annots=None, *args,**kwargs,):
    """

    Displays a grid of the similarity maps of the given molecules 
    <mols> with the given atom attributions <weights>, 
    as well as optional annotations to go along with it.

    Fixes the issue described in:
    https://github.com/rdkit/rdkit/issues/7937
    by scaling every instance's color space according to the overall
    minimum / maximum, this way making the color scheme consistent
    throughout all provided instances.
    This is useful in practice because we don't want every example
    to span its own reference frame, but instead want the reference
    frame (= color scheme) to be consistent throughout, so
    that we can compare atom attributions across different examples.
    
    Example Usage:
    ================
    weights = [0.3878029485291862, 0.18429023419510093, 0.004651627523322409,
     0.993537583585071, 0.16879923401475366, 0.36112685890584145,
     0.47997814021490615, 0.5133176109505179, -0.1555064997665855,
     0.048520776948712, -0.025468312563339257, 0.37974263117933693,
     0.2491675745533543, 0.09273801343718749, 0.2727059394215573, 0.1264606792172728, 
     -0.025468312563339257, 0.37974263117933693, 0.2491675745533543, 0.09273801343718749,
     0.2727059394215573, 0.1264606792172728, - 0.11780881485086714, -1]
    
    smi = 'CC(=O)OCN1C(=O)NC(c2ccccc2)(c2ccccc2)C1=O'
    mol = Chem.MolFromSmiles(smi)
    display(GridSimilarityMaps([mol for _ in range(32)],[np.array(weights) * (random.random()-0.5) * 2 for i in range(32)]))
    """
    
    if len(mols) != len(weights):
        raise ValueError("mols and weights must be same length!")

    if annots is None:
        annots = ["" for _ in mols]
        
    overall_min = np.array(weights).reshape(-1).min()
    overall_max = np.array(weights).reshape(-1).max()

    print("min",overall_min,"max",overall_max,)

    images = []
    for mol,weight in zip(mols,weights):
        draw2d = rdMolDraw2D.MolDraw2DCairo(300,300)
        d = GetSimilarityMapFromWeightsWithScale(mol,weight,draw2d,overall_min,overall_max,*args,**kwargs,)
        d.FinishDrawing()
        image_data = d.GetDrawingText()
        images.append(Image.open(io.BytesIO(image_data)))


    N = len(images)
    n_cols = 4
    n_rows = math.ceil(N / n_cols)
    tile_image = tile_images_with_annots(images, annots=annots, tile_w=300,tile_h=300,n_tiles_w=n_cols,n_tiles_h=n_rows,sampling_strategy='deterministic',)
    return tile_image

    




