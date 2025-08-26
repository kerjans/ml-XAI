
from io import StringIO
import os
# Silence any TQDM output as it is writing to stderr!!!!
os.environ["TQDM_DISABLE"] = "True"

from pathlib import Path
import asyncio
import json
import io
import base64
import random
import secrets
import sys
import time
from standardizer.io import Silencer, else_none
import tornado
import tornado.ioloop
import tornado.web
from tornado import web
import pandas as pd
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

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


HERE = Path(__file__).parent 

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


        templ = (HERE / "dialog_box.html").read_text()
        dlg = f"""<form action="/login" method="post">
                    <p>Please complete the captcha (may take several seconds)</p>
                    <img src="{src}"></img>
                   <br>
                   Captcha: <input type="text" name="name">
                   <input type="submit" value="Sign in">
                   <br>
                   <br>
                   <p class="small-print">
                   By signing in you agree to our use of cookies and client-side data storage
                   which is only used for technical reasons (login handling, keeping track of past jobs). 
                   <br>
                   Also, you agree to our <a href="/terms_of_use">terms of use</a>.
                   </p>
                   </form>
                   """
        self.write(templ.replace("{{MESSAGE}}",dlg))
        img_fle.unlink()

    def post(self):
        trial = self.get_argument("name")
        if trial in self.IMAGE_CAPTCHAS:
            print("LOGIN SUCCESS")
            self.set_signed_cookie("user", str(uuid.uuid4()))
            self.redirect("/wisp")
        else:
            print("LOGIN FAILURE")
            self.redirect("/login")

TERMS_OF_USE = HERE / "terms_of_use.txt"
class TermsOfUseHandler(BaseHandler):
    @log_function_call
    def read(self):
        terms = TERMS_OF_USE.read_text()
        page = f"<html><body><pre>{terms}</pre></body></html>"
        return page

    @log_function_call
    def get(self):
        self.write(self.read())

# The LandingPageHandler displays the initial landing web page
LANDING_PAGE = HERE / "landing_page.svg"
class LandingPageHandler(BaseHandler):
    @log_function_call
    def read(self):
        temp = (HERE / "landing_page.html").read_text()
        svg = (HERE / "landing_page.svg").read_text()
        svg_mobile = (HERE / "landing_page_mobile.svg").read_text()
        return temp.replace("{SVG}","\n".join([svg,svg_mobile,]))

    @log_function_call
    def get(self):
        self.write(self.read())

# The MainHandler displays the wisp web app
SITE = HERE / "site.html"
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

    @tornado.web.authenticated
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

    @tornado.web.authenticated
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

    @tornado.web.authenticated
    @log_function_call
    def post(self):
        req = json.loads(self.request.body)
        print("req:",req)

        job_id = req["job_id"]
        page_size = req["page_size"]
        current_page = req["current_page"]

        here = Path(__file__).parent
        working_dir = here / "working_dir" / f"{job_id}"
        heat_maps_dir = working_dir / "HeatMaps"

        meta_fle = working_dir / "metadata.json"
        meta_df = json.loads(meta_fle.read_text())

        n_pages = 0
        heat_maps = []
        if heat_maps_dir.exists():
            fles = sorted(heat_maps_dir.glob("*.png"))
            legend_fles = [fle for fle in fles if "legend" in fle.name.lower()]
            n_entries = len(fles) - len(legend_fles)  
            for png_fle in fles[page_size*current_page:page_size*(current_page + 1)]+legend_fles:
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
                "n_entries": n_entries,
                "status": "success",
            }
        )
        self.write(resp)

class WispOverviewPage(BaseHandler):

    @tornado.web.authenticated
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



def styled_html_table(df, columns):
    """
    Generates a modern styled HTML table from a DataFrame and specified columns
    with alternating row colors.

    Parameters:
    df (pd.DataFrame): The input DataFrame.
    columns (list): List of columns to include in the HTML table.

    Returns:
    str: A string containing the HTML representation of the table.
    """
    # Start the HTML table
    html = '<table border="1" style="border-collapse: collapse; width: 100%;">\n'
    
    # Create the header row
    html += '  <tr style="background-color: #1b3f90ff; color: white;">\n'
    for column in columns:
        html += f'    <th>{column}</th>\n'
    html += '  </tr>\n'
    
    # Create the data rows with alternating colors
    row_count = 0
    for _, row in df.iterrows():
        row_color = '#ffffff' if row_count % 2 == 1 else "#5781e5ff"  # White and light blue
        html += f'  <tr style="background-color: {row_color};">\n'
        for column in columns:
            html += f'    <td style="padding: 8px; border: 1px solid #ddd;">{row[column]}</td>\n'
        html += '  </tr>\n'
        row_count += 1
    
    # End the HTML table
    html += '</table>'
    
    return html

class ModelPerfOverview(BaseHandler):

    @tornado.web.authenticated
    @log_function_call
    def post(self):
        req = json.loads(self.request.body)

        job_id = req["job_id"]

        here = Path(__file__).parent
        working_dir = here / "working_dir" / f"{job_id}"
        perf_overview = working_dir / "Grid-Search.csv"
        df = pd.read_csv(perf_overview)
        df = df[["Model_Type","Feature","r2","MAE","RMSE",]]
        for col in ["r2","MAE","RMSE"]:
            df[col] = df[col].apply(lambda val: "{:.2f}".format(val))

        # split out any args e.g. "RandomForestClassifier()" => "RandomForestClassifier"
        df["Model_Type"] = df["Model_Type"].apply(lambda n: n.split("(")[0])
        df = df.sort_values("MAE")

        resp = json.dumps(
            {
                "status": "success",
                "model_perf_overview": styled_html_table(df,["Model_Type","Feature","r2","MAE","RMSE",]),
            }
        )
        self.write(resp)


def _sniff_smiles_columns(df,):
    smi_else_none = else_none(Chem.MolFromSmiles)
    for col in sorted(df.columns):
        col_score = 0
        for smi in df[col].sample(10,random_state=123):
            mol = smi_else_none(smi)
            mol_score = mol is not None and bool(mol.GetNumAtoms())
            col_score += mol_score
        if col_score > 0:
            yield col

def _sniff_float_columns(df,):
    float_else_none = else_none(float)
    for col in sorted(df.columns):
        col_score = 0
        for smi in df[col].sample(10,random_state=123):
            val = float_else_none(smi)
            col_score += val is not None

        if col_score > 0:
            yield col


class GuessColumnsHandler(BaseHandler):

    @tornado.web.authenticated
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
                
                if "float" in str(df[col].dtype) or "int" in str(df[col].dtype):
                    target_col = col

        except:
            pass

        resp = json.dumps(
            {
                "all_columns": list(df.columns),
                "all_smiles_columns": list(_sniff_smiles_columns(df)),
                "all_float_columns": list(_sniff_float_columns(df)),
                "smiles_col": smiles_col,
                "target_col": target_col,
                "status": "success",
            }
        )
        self.write(resp)


from WISP import WISP # :o
from WISP import plotting_helper
from contextlib import redirect_stdout,redirect_stderr
def run_wisp(args,metafle):

    working_dir = metafle.parent
    log_out_fle = working_dir / "out.log"
    log_err_fle = working_dir / "err.log"
    log_out_fle.touch(exist_ok=True,)
    log_err_fle.touch(exist_ok=True,)


    with open(log_out_fle,"w") as fout:
        with open(log_err_fle,"w") as ferr:
            with redirect_stderr(ferr):
                with redirect_stdout(fout):
                    print("START","run_wisp")
                    sys.stdout.flush()
                    WISP.DISPLAY_PLOTS = False
                    plotting_helper.DISPLAY_PLOTS = False
                    WISP.WISP(**args)
                    print("END","run_wisp")
                    sys.stdout.flush()
                    metadat = json.loads(metafle.read_text())

    metadat["status"] = "done"
    metadat["log_out"] = log_out_fle.read_text()
    metadat["log_err"] = log_err_fle.read_text()
    metafle.write_text(json.dumps(metadat))
    return

class FeatureImportanceHandler(BaseHandler):

    @tornado.web.authenticated
    @log_function_call
    def post(self):
        req = json.loads(self.request.body)
        job_id = req["job_id"]

        here = Path(__file__).parent
        working_dir = here / "working_dir" / f"{job_id}"

        feature_imp_fle = working_dir / "feature_imp.svg"
        feature_imp_fle.touch(exist_ok=True,)
        svg_content = feature_imp_fle.read_text()

        resp = json.dumps(
            {
                "status": "success",
                "feature_importance_plot": svg_content
            }
        )
        self.write(resp)



class JobStatusHandler(BaseHandler):

    @tornado.web.authenticated
    @log_function_call
    def post(self):
        req = json.loads(self.request.body)
        print("req:",req)

        job_id = req["job_id"]

        here = Path(__file__).parent
        working_dir = here / "working_dir" / f"{job_id}"
        metafle = working_dir / "metadata.json"
        metadat = json.loads(metafle.read_text())

        log_out_fle = working_dir / "out.log"
        log_out_fle.touch(exist_ok=True,)
        log_out_txt = log_out_fle.read_text()
        log_err_fle = working_dir / "err.log"
        log_err_fle.touch(exist_ok=True,)
        log_err_txt = log_err_fle.read_text()

        resp = json.dumps(
            {
                "job_status": metadat["status"],
                "job_name": metadat["job_name"],
                "log_out": log_out_txt,
                "log_err": log_err_txt,
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

    @tornado.web.authenticated
    @log_function_call
    async def post(self):

        job_id = generate_job_id()

        try:
            req = json.loads(self.request.body)
            print("req:",req)
            csv = StringIO(req["csv"][0])
            df = resilient_read_csv(csv)

            if df is not None:
                smiles_col = req["smiles_col"]
                target_col = req["target_col"]
                job_name = req["job_name"]
                smiles = df[smiles_col].tolist()
                target = df[target_col].tolist()

                df_new = pd.DataFrame(
                    {"smiles": smiles, "target": target,
                     "ID": [str(i) for i in range(len(smiles))]}
                    )
                
                if len(df_new) > 5000:
                    print("WARNING: DOWNSAMPLE data to 5000 entries")
                    df_new = df_new[~df_new["smiles"].apply(else_none(Chem.MolFromSmiles)).isna()]
                    df_new = df_new.sample(5000,random_state=123,)

                print("WARNING: DOWNSAMPLE by 50%")
                df_new = df_new.sample(frac=.5)

                id_col = "ID"

                here = Path(__file__).parent

                working_dir = here / "working_dir" / f"{job_id}"
                working_dir.mkdir(exist_ok=True,)

                metafle = working_dir / "metadata.json"
                metadat = {
                            "smiles": smiles,
                            "target": target,
                            "job_name": job_name,
                            "job_id": job_id,
                            "status": "running",
                        }
                metafle.write_text(json.dumps(metadat))

                input_fle = working_dir / f"input_{job_id}.csv"

                df_new.to_csv(input_fle)

                args = {
                    # fix for internal wisp processing problems
                    "working_dir":str(working_dir)+str(os.path.sep), 
                    "input_dir":str(input_fle),
                    "ID_Column_Name":id_col,
                    "Smiles_Column_Name":"smiles",
                    "Target_Column_Name":"target",
                    "use_GNN":True,
                    "fast_run":True,
                    }
                PROCESS_PARALLEL = True
                if PROCESS_PARALLEL:
                    loop = asyncio.get_running_loop()
                    _ = loop.run_in_executor(
                        PROCESS_POOL, run_wisp, 
                        args, metafle
                    )
                else:
                    os.environ["_WISP_NO_PARALLEL"] = "True"
                    run_wisp(args,metafle)


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
            (r"/", LandingPageHandler),
            (r"/terms_of_use", TermsOfUseHandler),
            (r"/JobSubmission", JobSubmissionHandler),
            (r"/JobStatus", JobStatusHandler),
            (r"/GuessColumnsHandler", GuessColumnsHandler),
            (r"/HeatMaps", HeatMaps),
            (r"/WispOverviewPage", WispOverviewPage),
            (r"/MMPOverview", MMPOverview),
            (r"/FeatureImportance", FeatureImportanceHandler),
            (r"/ModelPerfOverview",ModelPerfOverview),
            (r"/wisp", MainHandler),
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

    http_server.listen(80)
    await asyncio.Event().wait()


if __name__ == "__main__":
    asyncio.run(main())
