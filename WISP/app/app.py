
from io import StringIO
from pathlib import Path
import asyncio
import json
import io
import base64
import tornado
import tornado.ioloop
import tornado.web
from tornado import web
import pandas as pd
import numpy as np

import importlib

STATIC_FILE_DIR = Path(__file__).parent


# Some Utilities
DEVEL = True

def compose(*functions):
    """
    Composes the given functions.
    """
    def inner(arg):
        result = arg
        for func in reversed(functions):
            result = func(result)
        return result
    return inner

def info(*args,**kwargs):
    print(*args,**kwargs)

def debug(*args,**kwargs):
    print(*args,**kwargs)

def log_function_call(func):
    """

    >>> class C:
    ...     @log_function_call
    ...     def run(self):
    ...         print("hello")
    >>> o = C()
    >>> o.run() #doctest:SKIP
    call   -> run@<doctest molfactual.log_function_call[0]>:2
    hello
    return <- run
    """

    def wrapper(*args, **kwargs):
        if DEVEL:
            print(
                "call   ->",
                func.__name__
                + "@"
                + func.__code__.co_filename
                + ":"
                + str(func.__code__.co_firstlineno),
            )
        rslt = func(*args, **kwargs)
        if DEVEL:
            print("return <-", func.__name__)
        return rslt

    return wrapper



# The BaseHandler is the blueprint for a Request Handler
class BaseHandler(tornado.web.RequestHandler):
    pass



# The MainHandler displays the website
SITE = Path(__file__).parent / "site.html"
class MainHandler(BaseHandler):
    @log_function_call
    def read(self):
        with open(SITE) as reader:
            return reader.read()

    @log_function_call
    def get(self):
        self.write(self.read())


def resilient_read_csv(csv):
    try:
        return pd.read_csv(csv)
    except:
        try:
            return pd.read_csv(csv,sep="\t")
        except:
            return pd.read_csv(csv,sep=";")
    return None

def resilient_read_smiles(df):
    if "smiles" in df.columns:
        return df["smiles"].tolist()
    elif "SMILES" in df.columns:
        return df["SMILES"].tolist()
    else:
        smi_cols = [col for col in df.columns if 'smiles' in col.lower()]
        if len(smi_cols) == 1:
            return df[smi_cols[0]].tolist()
    
    return []

def resilient_read_target(df):
    if "y" in df.columns:
        return df["y"].tolist()
    elif "ground_truth" in df.columns:
        return df["ground_truth"].tolist()
    elif "target" in df.columns:
        return df["target"].tolist()
    else:
        smi_cols = [col for col in df.columns if 'smiles' in col.lower()]
        target_cols = [col for col in df.columns if col not in smi_cols]
        if len(smi_cols) == 1 and len(target_cols) == 1:
            return df[target_cols[0]].tolist()



from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
from copy import deepcopy

def mol_to_image(mol, width=300, height=300) -> "Image":
    rand_fle = "tmp.jpg"
    Draw.MolToImageFile(mol, filename=rand_fle, format="JPG", size=(width, height))
    img = Image.open(rand_fle)
    return deepcopy(img)

class DatasetHandler(BaseHandler):
    @log_function_call
    def post(self):
        try:
            req = json.loads(self.request.body)

            print("req:",req)
            csv = StringIO(req["csv"][0])
            df = resilient_read_csv(csv)

            if df is not None:
                smiles = resilient_read_smiles(df)
                target = resilient_read_target(df)
                if not target or (smiles and len(target) != len(smiles)):
                    target = ["N/A" for _ in smiles]
                    

                images = {}
                # This here would write molecule images as example output,
                # just for testing / illustration purposes:
                if False:
                    images = []
                    for smi in smiles:
                        output = io.BytesIO()
                        img = mol_to_image(Chem.MolFromSmiles(smi))
                        img.save(output, format="png")
                        hex_data = output.getvalue()
                        img64 = base64.b64encode(hex_data).decode("utf-8")
                        images.append(img64)
                
                result_dir = Path(__file__).parent.parent.parent / "example_images"
                for img_fle in result_dir.glob("*.png"):
                    img = Image.open(img_fle)
                    output = io.BytesIO()
                    img.save(output, format="png")
                    hex_data = output.getvalue()
                    img64 = base64.b64encode(hex_data).decode("utf-8")
                    images[img_fle.name] = img64

            else:
                smiles = []
                target = []
                images = {}

            
            resp = json.dumps(
                {
                    "smiles": smiles,
                    "target": target,
                    "images": images,
                    "status": "success"
                }
            )
            self.write(resp)
        except:
            resp = json.dumps(
                {"status": "failure"}
            )
            self.write(resp)
            if DEVEL:
                raise


async def main():
    application = tornado.web.Application(
        [
            (r"/", MainHandler),
            (r"/DatasetHandler", DatasetHandler),
            (r"/static/(.*)", tornado.web.StaticFileHandler, {"path": STATIC_FILE_DIR}),
        ],
        autoreload=True,
        cookie_secret="__TODO:_GENERATE_YOUR_OWN_RANDOM_VALUE_HERE__",
    )

    try:
        http_server = tornado.httpserver.HTTPServer(
            application,
            ssl_options={
                "certfile": "cert/cert.pem",
                "keyfile": "cert/key.pem",
            },
        )
    except:
        info("no cert file found, defaulting to http")
        http_server = tornado.httpserver.HTTPServer(
            application,
        )

    http_server.listen(8912)
    await asyncio.Event().wait()


if __name__ == "__main__":
    asyncio.run(main())