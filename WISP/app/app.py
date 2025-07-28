
from io import StringIO
import os
from pathlib import Path
import asyncio
import json
import io
import base64
import random
import secrets
import time
import tornado
import tornado.ioloop
import tornado.web
from tornado import web
import pandas as pd
import numpy as np

import concurrent.futures

PROCESS_POOL = concurrent.futures.ProcessPoolExecutor(max_workers=4,)

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
    def get_current_user(self):
        return self.get_signed_cookie("user")

from captcha.image import ImageCaptcha
import random
import string
import time
class LoginHandler(tornado.web.RequestHandler):

    IMAGE_CAPTCHAS = []

    sleep_time = 1

    def get(self):
        time.sleep(self.sleep_time * len(self.IMAGE_CAPTCHAS))
        cap = ImageCaptcha()
        characters = string.ascii_uppercase + string.digits
        cap_text = ''.join(secrets.choice(characters) for _ in range(6))
        self.IMAGE_CAPTCHAS.append(cap_text)
        self.IMAGE_CAPTCHAS = self.IMAGE_CAPTCHAS[-10:]
        img_fle = Path(f"captcha_{uuid.uuid4()}.png")
        cap.write(cap_text,str(img_fle))

        hex_data = img_fle.read_bytes()
        img64 = base64.b64encode(hex_data).decode("utf-8")
        src = "data:image/png;base64," + img64
        self.write('<html><body><form action="/login" method="post">'
                    "<p>Please complete the captcha (may take several seconds)</p>"
                    f'<img src="{src}"></img>'
                   'Captcha: <input type="text" name="name">'
                   '<input type="submit" value="Sign in">'
                   '</form></body></html>')

    def post(self):
        trial = self.get_argument("name")
        if trial in self.IMAGE_CAPTCHAS:
            print("LOGIN SUCCESS")
            self.set_signed_cookie("user", str(uuid.uuid4()))
            self.redirect("/")
        else:
            print("LOGIN FAILURE")
            self.redirect("/login")


# The MainHandler displays the website
SITE = Path(__file__).parent / "site.html"
class MainHandler(BaseHandler):
    @log_function_call
    def read(self):
        with open(SITE) as reader:
            return reader.read()

    @tornado.web.authenticated
    @log_function_call
    def get(self):
        self.write(self.read())


def resilient_read_csv(csv,nrows=None,):
    try:
        return pd.read_csv(csv,nrows=nrows,)
    except:
        try:
            return pd.read_csv(csv,sep="\t",nrows=nrows)
        except:
            return pd.read_csv(csv,sep=";",nrows=nrows)
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

    return None



from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
from copy import deepcopy

def mol_to_image(mol, width=300, height=300) -> "Image":
    rand_fle = "tmp.jpg"
    Draw.MolToImageFile(mol, filename=rand_fle, format="JPG", size=(width, height))
    img = Image.open(rand_fle)
    return deepcopy(img)


class MMPOverview(BaseHandler):
    @log_function_call
    def post(self):
        req = json.loads(self.request.body)
        job_id = req["job_id"]
        here = Path(__file__).parent
        working_dir = here / "working_dir" / f"{job_id}"

        mmp_fle = working_dir / "MMPs_with_attributions.csv"
        df = pd.read_csv(mmp_fle)

        # df columns:
        # smiles_1,smiles_2,ID_1,ID_2,transformation,constant,constant_atom_count,Atom Attributions_1,Atom Attributions_2,target_1,target_2
        df["target_diff"] = df["target_2"] - df["target_1"]

        resp = json.dumps(
            {
                "mmp_overview_data": json.dumps(
                    [
                        {"MMP_rule":row["transformation"],"target_diff":row["target_diff"],"smiles_1":row["smiles_1"],"smiles_2":row["smiles_2"]}
                        for _,row in df.iterrows()
                     ]),
                "status": "success",
            }
        )
        self.write(resp)

class MoleculePage(BaseHandler):
    @log_function_call
    def post(self):
        req = json.loads(self.request.body)
        print("req:",req)

        job_id = req["job_id"]

        here = Path(__file__).parent
        working_dir = here / "working_dir" / f"{job_id}"

        meta_fle = working_dir / "metadata.json"
        df = json.loads(meta_fle.read_text())

        mol_images = []
        if meta_fle.exists():
            for smi in df.smiles:
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    img = mol_to_image(mol)
                    output = io.BytesIO()
                    img.save(output, format="png")
                    hex_data = output.getvalue()
                    img64 = base64.b64encode(hex_data).decode("utf-8")
                    mol_images.append(img64)

        resp = json.dumps(
            {
                "mol_images": mol_images,
                "status": "success",
            }
        )
        self.write(resp)

class HeatMaps(BaseHandler):
    @log_function_call
    def post(self):
        req = json.loads(self.request.body)
        print("req:",req)

        job_id = req["job_id"]

        here = Path(__file__).parent
        working_dir = here / "working_dir" / f"{job_id}"
        heat_maps_dir = working_dir / "HeatMaps"

        meta_fle = working_dir / "metadata.json"
        meta_df = json.loads(meta_fle.read_text())

        heat_maps = []
        if heat_maps_dir.exists():
            for png_fle in heat_maps_dir.glob("*.png"):
                hex_data = png_fle.read_bytes()
                img64 = base64.b64encode(hex_data).decode("utf-8")

                if "legend" in png_fle.name.lower():
                    legend = img64
                else:
                    heat_maps.append(img64)

        resp = json.dumps(
            {
                "legend": legend,
                "heatmaps": heat_maps,
                "status": "success",
            }
        )
        self.write(resp)

class WispOverviewPage(BaseHandler):
    @log_function_call
    def post(self):
        req = json.loads(self.request.body)
        print("req:",req)

        images = {}
        
        job_id = req["job_id"]

        here = Path(__file__).parent
        working_dir = here / "working_dir" / f"{job_id}"
        for img_fle in working_dir.glob("*.png"):
            img = Image.open(img_fle)
            output = io.BytesIO()
            img.save(output, format="png")
            hex_data = output.getvalue()
            img64 = base64.b64encode(hex_data).decode("utf-8")
            images[img_fle.name] = img64

        resp = json.dumps(
            {
                "images": images,
                "status": "success",
            }
        )
        self.write(resp)


class GuessColumnsHandler(BaseHandler):
    @log_function_call
    def post(self):

        smiles_col = None
        target_col = None
        try:
            req = json.loads(self.request.body)
            print("req:",req)
            csv = StringIO(req["csv"][0])
            df = resilient_read_csv(csv,nrows=30,)
            for col in df.columns:
                if "smi" in col.lower():
                    smiles_col = col
                
                if "float" in str(df[col].dtype):
                    target_col = col

        except:
            pass

        resp = json.dumps(
            {
                "smiles_col": smiles_col,
                "target_col": target_col,
                "status": "success",
            }
        )
        self.write(resp)


from wisp.WISP import WISP # :o
def run_wisp(args,metafle):
    print("START","run_wisp")
    #time.sleep(20)
    WISP(**args)
    print("END","run_wisp")
    metadat = json.loads(metafle.read_text())
    metadat["status"] = "done"
    metafle.write_text(json.dumps(metadat))
    return


class JobStatusHandler(BaseHandler):
    @log_function_call
    def post(self):
        req = json.loads(self.request.body)
        print("req:",req)

        job_id = req["job_id"]

        here = Path(__file__).parent
        working_dir = here / "working_dir" / f"{job_id}"
        metafle = working_dir / "metadata.json"
        metadat = json.loads(metafle.read_text())

        resp = json.dumps(
            {
                "job_status": metadat["status"],
            }
        )
        self.write(resp)


from datetime import datetime
import uuid
def generate_job_id():
    return "__".join(
        [
            datetime.today().strftime('%Y_%m_%d_%H_%M_%S'),
            uuid.uuid4().hex,
         ])

class JobSubmissionHandler(BaseHandler):
    @log_function_call
    async def post(self):

        job_id = generate_job_id()

        try:
            req = json.loads(self.request.body)
            print("req:",req)
            csv = StringIO(req["csv"][0])
            df = resilient_read_csv(csv)

            if df is not None:
                smiles = resilient_read_smiles(df)
                #target = resilient_read_target(df)
                target = df["measured log solubility in mols per litre"].tolist()
                if target is None or (smiles and len(target) != len(smiles)):
                    target = ["N/A" for _ in smiles]
                    

                df_new = pd.DataFrame(
                    {"smiles": smiles, "target": target,
                     "ID": [str(i) for i in range(len(smiles))]}
                    )
                id_col = "ID"
                smiles_col = "smiles"
                target_col = "target"

                here = Path(__file__).parent

                working_dir = here / "working_dir" / f"{job_id}"
                working_dir.mkdir(exist_ok=True,)

                metafle = working_dir / "metadata.json"
                metadat = {
                            "smiles": smiles,
                            "target": target,
                            "job_id": job_id,
                            "status": "running",
                        }
                metafle.write_text(json.dumps(metadat))

                input_fle = working_dir / f"input_{job_id}.csv"
                df_new.to_csv(input_fle)

                if True:
                    loop = asyncio.get_running_loop()
                    _ = loop.run_in_executor(
                        PROCESS_POOL, run_wisp, 
                        {
                            # fix for internal wisp processing problems
                            "working_dir":str(working_dir)+str(os.path.sep), 
                            "input_dir":str(input_fle),
                            "ID_Column_Name":id_col,
                            "Smiles_Column_Name":smiles_col,
                            "Target_Column_Name":target_col,
                            "use_GNN":False,
                            },  metafle
                    )


            resp = json.dumps(
                {
                    "job_id": job_id,
                    "status": "success",
                }
            )
            self.write(resp)
            print("FINISHED PROCESSING SUBMIT JOB")

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
            (r"/JobSubmission", JobSubmissionHandler),
            (r"/JobStatus", JobStatusHandler),
            (r"/GuessColumnsHandler", GuessColumnsHandler),
            (r"/HeatMaps", HeatMaps),
            (r"/WispOverviewPage", WispOverviewPage),
            (r"/MMPOverview", MMPOverview),
            (r"/static/(.*)", tornado.web.StaticFileHandler, {"path": STATIC_FILE_DIR}),
            (r"/login", LoginHandler),
        ],
        autoreload=True,
        cookie_secret="__TODO:_GENERATE_YOUR_OWN_RANDOM_VALUE_HERE__",
        login_url="/login",
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

    http_server.listen(8914)
    await asyncio.Event().wait()


if __name__ == "__main__":
    asyncio.run(main())