from contextlib import redirect_stderr, redirect_stdout
from collections import defaultdict
import copy
import io
import tempfile
import contextlib


import json
from pathlib import Path
import random
import uuid
import os
import joblib


import random
import numpy as np
import pandas as pd
from rdkit.Chem import Draw
from rdkit import Chem

import PIL
from PIL import ImageDraw
from PIL.Image import Image


import hashlib

from standardizer.config import _DIR_DATA


def warn(msg):
    print("\033[31;1;4m", msg, "\033[0m")


class Silencer:
    """
    A useful tool for silencing stdout and stderr.
    Usage:
    >>> with Silencer() as s:
    ...         print("kasldjf")

    >>> print("I catched:",s.out.getvalue())
    I catched: kasldjf
    <BLANKLINE>

    Note that nothing was printed and that we can later
    access the stdout via the out field. Similarly,
    stderr will be redirected to the err field.
    """

    def __init__(self):
        self.out = io.StringIO()
        self.err = io.StringIO()

    def __enter__(self):
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
        self.rs = redirect_stdout(self.out)
        self.re = redirect_stderr(self.err)
        self.rs.__enter__()
        self.re.__enter__()
        return self

    def __exit__(self, exctype, excinst, exctb):
        from rdkit import RDLogger

        RDLogger.EnableLog("rdApp.*")
        self.rs.__exit__(exctype, excinst, exctb)
        self.re.__exit__(exctype, excinst, exctb)


def to_smi(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(mol)


def from_smi(smi: str) -> Chem.Mol:
    return Chem.MolFromSmiles(smi)


def from_smi_else_none(smi: str) -> Chem.Mol:
    try:
        return Chem.MolFromSmiles(smi)
    except:
        return None


_opsin_cache_file = _DIR_DATA / ".opsin_cache.json"


def iupac_to_smiles(iupac: str) -> str:
    """
    Forces the actual parsing of the IUPAC via
    the OPSIN tool
    :param iupac: The IUPAC string to parse
    :return: The resulting SMILES string.
    """
    if _opsin_cache_file.exists():
        cache = json.loads(_opsin_cache_file.read_text())
    else:
        cache = {}
    if iupac in cache:
        return cache[iupac]
    try:
        tmp_fle = f"/tmp/__opsin_tmp__{random.randint(1000,10000000)}.in"
        with open(tmp_fle, "w") as fout:
            fout.write(iupac)
        smi = os.popen("opsin " + tmp_fle).read()
    finally:
        os.remove(tmp_fle)
    cache[iupac] = smi
    with open(_opsin_cache_file, "wt") as fout:
        json.dump(obj=cache, fp=fout)
    return smi


try:
    import cirpy
except:
    pass

_cas_cache_file = _DIR_DATA / ".cas_cache.json"


def cas_to_smiles(cas: str) -> str:
    """

    >>> cas_to_smiles('108-95-2')
    'Oc1ccccc1'
    """
    if _cas_cache_file.exists():
        cache = json.loads(_cas_cache_file.read_text())
    else:
        cache = {}
    if cas in cache:
        return cache[cas]
    try:
        rslt = cirpy.resolve(cas, "smiles")
    except:
        rslt = None

    cache[cas] = rslt

    with open(_cas_cache_file, "wt") as fout:
        json.dump(obj=cache, fp=fout)
    return cache[cas]


@contextlib.contextmanager
def working_directory(path):
    """Changes working directory and returns to previous on exit."""
    prev_cwd = Path.cwd()
    os.chdir(str(path))
    try:
        yield
    finally:
        os.chdir(str(prev_cwd))


def temp_dir():
    """
    Creates a path to a temporary directory
    """
    if Path("/tmp").exists():
        tmp_dir = Path("/tmp") / "smal"
    else:
        tmp_dir = Path(tempfile.gettempdir()) / "smal"
    tmp_dir.mkdir(
        exist_ok=True,
    )
    return tmp_dir


def random_fle(
    suf: str,
):
    """
    Creates a random file with the given suffix suf
    """
    fle_name = f"{uuid.uuid4().hex}.{suf}"
    return temp_dir() / fle_name


def download_file(
    url: str,
    dest: str,
) -> bool:
    """
    Attempts to download the given file given file first using curl, then
    using wget and finally using pythons urllib.

    A warning is printed in case this download fails.
    """
    import urllib.request

    if os.system(f"curl {url} > {dest}") and os.system(f"wget {url} -O {dest}"):
        try:
            urllib.request.urlretrieve(f"{url}", f"{dest}")
        except:
            warn(f"failed to download file {dest} from url {url}")
            warn(
                f"please download file from url: {url} and move it into location {dest}"
            )


def obj_store(ident: str, obj, store_dir: Path) -> None:
    store_dir = Path(store_dir)
    hsh = hashlib.sha1(ident.encode("utf8")).hexdigest()
    store_dir.mkdir(exist_ok=True)
    tmp_fle = (store_dir / hsh).with_suffix(f".{random.randint(1,100000000000000)}tmp")
    joblib.dump(value=obj, filename=tmp_fle)
    tmp_fle.rename(tmp_fle.with_suffix(".pkl"))


def obj_load(ident: str, store_dir: Path):
    store_dir = Path(store_dir)
    store_dir.mkdir(exist_ok=True)
    hsh = hashlib.sha1(ident.encode("utf8")).hexdigest()
    return joblib.load(store_dir / (hsh + ".pkl"))


def mol_to_image(
    mol,
    width=300,
    height=300,
) -> Image:
    try:
        tmp_fle = str(random_fle("jpg"))
        Draw.MolToImageFile(mol, filename=tmp_fle, format="JPG", size=(width, height))
        img = PIL.Image.open(tmp_fle)
        return copy.deepcopy(img)
    finally:
        if Path(tmp_fle).exists():
            Path(tmp_fle).unlink()


def read_pix_format(df, sep=None, skip_cols=None):
    """
    Converts the given dataframe from the <sep>-separated
    pix-format, into a expanded, regular, dataframe.
    >>> df = pd.DataFrame([
    ... {"id":1,"a":"1|2|3","b":"x|y|z","c":"something",},
    ... {"id":2,"a":"4|5|6","b":"i|o|p","c":"another word",},
    ... ])
    >>> read_pix_format(df)
       id             c  a  b
    0   1     something  3  z
    1   1     something  3  z
    2   1     something  3  z
    3   2  another word  6  p
    4   2  another word  6  p
    5   2  another word  6  p
    """
    if skip_cols is None:
        skip_cols = []
    if sep is None:
        sep = "|"

    df_out = []
    for row in df.to_dict("records"):
        rows_not_expanded = {}
        rows_expanded = defaultdict(list)
        for key, val in row.items():
            if isinstance(val, str) and sep in val and key not in skip_cols:
                for innerval in val.split(sep):
                    rows_expanded[key].append(innerval)
            else:
                rows_not_expanded[key] = val

        N_rows_expanded = {len(val) for val in rows_expanded.values()}
        if len(N_rows_expanded) not in [0, 1]:
            err_string = []
            for N in N_rows_expanded:
                err_string.append(f"the following rows have {N} entries:")
                for key, ent in rows_expanded.items():
                    if len(ent) == N:
                        err_string.append(f"key={key}:")
                        err_string.append(sep.join(ent))
            err_string = "\n".join(err_string)
            raise ValueError(f"could not parse row: {err_string}")

        N_rows_expanded = list(N_rows_expanded)[0] if N_rows_expanded else 0

        row_dct = dict(rows_not_expanded)
        if N_rows_expanded:
            for i in range(N_rows_expanded):
                for k, v in rows_expanded.items():
                    row_dct[k] = v[i]

                df_out.append(row_dct)
        else:
            df_out.append(row_dct)

    return pd.DataFrame(df_out)


def bde_fingerprint(
    df: "pd.DataFrame", fp_fun, bde_col=None, pop_fun=None
) -> "pd.DataFrame":
    """
    Computes the bond dissociation energy fingerprint,
    given a bde result dataframe as outputted by either
    our xtb bde tool or the alfabet bde program.

    Returns a dataframe consisting of smiles and fingerprints
    applying the formula

       fp[parent] = sum_{all fragments} pop_fun(bde) * (fp[fragment] - fp[parent])

    where pop_fun is a function that --- given the bde energy ---
    computes a corresponding "population" of the given fragment.
    We could think of pop_fun of providing a


    """
    assert "smiles" in df.columns
    assert "fragment1" in df.columns
    assert "fragment2" in df.columns
    if bde_col is None:
        bde_col = "bde"
    assert bde_col in df.columns
    df = df.copy()
    df["f1"] = [
        row["fragment1"]
        if len(row["fragment1"]) > len(row["fragment2"])
        else row["fragment2"]
        for _, row in df.iterrows()
    ]
    df["f2"] = [
        row["fragment2"]
        if len(row["fragment1"]) > len(row["fragment2"])
        else row["fragment1"]
        for _, row in df.iterrows()
    ]
    df["molp"] = df["smiles"].apply(else_none(Chem.MolFromSmiles))
    df["molf1"] = df["f1"].apply(else_none(Chem.MolFromSmiles))
    df["molf2"] = df["f2"].apply(else_none(Chem.MolFromSmiles))
    df["fpp"] = df["molp"].apply(fp_fun)
    df["fp1"] = df["molf1"].apply(fp_fun)
    df["fp2"] = df["molf2"].apply(fp_fun)

    rslt = []
    for smi in df.smiles.unique():
        bde_fp = np.zeros(fp_fun("CCC").shape)
        g = df[df.smiles == smi]
        for _, row in g.iterrows():
            if not Chem.CanonSmiles(row["fragment2"]) != Chem.CanonSmiles("[H]"):
                continue
            bde_fp += row["bde"] * (row["fp1"] - row["fpp"])
        rslt.append(
            {
                "smiles": smi,
                "bde_fp": bde_fp,
            }
        )
    return pd.concat(rslt)


def else_none(
    fun,
):
    """
    A useful wrapper that attempts to call the
    given function <fun>, but returns None
    in case any error is encountered.
    """

    def wrapper(*args, **kwargs):
        try:
            return fun(*args, **kwargs)
        except:
            return None

    return wrapper
