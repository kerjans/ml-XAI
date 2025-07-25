// see: https://stackoverflow.com/questions/58676967/load-textfile-via-drag-and-drop-on-textarea
// and: https://www.w3schools.com/html/html5_draganddrop.asp

console.log("loading ui.js");
function allowDrop(ev) {
    ev.preventDefault();
}


CURRENT_DATA = [];
const drop = function (evt) {
    window.my_event = evt;
    evt.preventDefault();
    var file = evt.dataTransfer.files[0];
    var reader = new FileReader();

    reader.onload = function (e) {
        const txt = e.target.result;
        //if (evt.target.id == "modelSetting") {
        CURRENT_DATA = [txt];

        // TODO: trigger dataset preview from here : )

        fetch("GuessColumnsHandler", {

            method: "POST",
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify({
                "csv": CURRENT_DATA
            })

        }).then(res => res.json()).then(res => {
            console.log("Request complete! response:", res);
            const status = res["status"];

            if (status == "success") {
                document.getElementById("smiles-column").value = res["smiles_col"];
                document.getElementById("target-column").value = res["target_col"];
            }

            else {
                alert("failure could not parse input!");
            }
        });
        //}
    };
    reader.readAsText(file, "UTF-8");
}

IMAGES = [];

const submitJob = function () {
    const but = document.getElementById("submit-button");
    but.disabled = true;
    but.style.opacity = 0.5;
    fetch("JobSubmission", {

        method: "POST",
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({
            "csv": CURRENT_DATA
        })

    }).then(res => res.json()).then(res => {
        console.log("Request complete! response:", res);
        const status = res["status"];

        if (status == "success") {
            document.getElementById("job-id-input").value = res["job_id"];
            but.disabled = false;
            but.style.opacity = 1.0;
        }

        else {
            alert("failure could not parse input!");
        }
    });
};

const get_current_job_id = function () {
    return document.getElementById("job-id-input").value;
};

const retrieveResults = function (dataset) {
    fetch("WispOverviewPage", {

        method: "POST",
        headers: {
            'Content-Type': 'application/json'
        }

        ,
        body: JSON.stringify({
            "job_id": get_current_job_id()
        })
    }).then(res => res.json()).then(res => {
        console.log("Request complete! response:", res);
        const status = res["status"];

        if (status == "success") {

            IMAGES = res["images"];
            //MOLECULE_IMAGES = res["molecule_images"];

            refreshFirstPage();
            //refreshSecondPage();
        }

        else {
            alert("failure could not parse input!");
        }
    }).then(

        fetch("HeatMaps", {

            method: "POST",
            headers: {
                'Content-Type': 'application/json'
            }

            ,
            body: JSON.stringify({
                "job_id": get_current_job_id()
            })
        }).then(res => res.json()).then(res => {
            console.log("Request complete! response:", res);
            const status = res["status"];

            if (status == "success") {

                //IMAGES = res["images"];
                MOLECULE_IMAGES = res["heatmaps"];

                //refreshFirstPage();
                refreshSecondPage();
            }

            else {
                alert("failure could not parse input!");
            }
        }
        )
    );
};

const styleImage = function (img) {
    if (img && img.style) {
        //img.style.height = '300px';
        img.style.width = '300px';
    }
};

const refreshFirstPage = function () {
    var tab = document.createElement("div");
    tab.style.display = "table";
    // clear the current results
    document.getElementById("result-div").innerHTML = "";

    // First Row: Showing the parity plots
    var cp = document.createElement("div");
    cp.style.display = "table-row";

    // Top Left -- How well can we explain predictions?
    const explain_pred_div = document.createElement("div");
    const explain_pred_label = document.createElement("div");
    explain_pred_label.innerText = "How well can we explain predictions?";
    const explain_pred_img = document.createElement('img');

    const col = "PREDvsCONTRIBUTIONSfragmentAtom Attributions_Training_Set.png";
    explain_pred_img.src = "data:image/png;base64," + IMAGES[col];
    styleImage(explain_pred_img);

    explain_pred_div.appendChild(explain_pred_label);
    explain_pred_div.appendChild(explain_pred_img);

    explain_pred_div.style.display = "inline";
    explain_pred_div.style.display = "table-cell";
    explain_pred_div.style.width = "350px";
    cp.appendChild(explain_pred_div);

    // Top Right -- How well can we explain reality?
    const explain_exp_div = document.createElement("div");
    const explain_exp_label = document.createElement("div");
    explain_exp_label.innerText = "How well can we explain reality?";
    const explain_exp_img = document.createElement('img');

    const col2 = "EXPvsCONTRIBUTIONSwholeAtom Attributions_Test_Set.png";
    explain_exp_img.src = "data:image/png;base64," + IMAGES[col2];
    styleImage(explain_exp_img);

    explain_exp_div.appendChild(explain_exp_label);
    explain_exp_div.appendChild(explain_exp_img);

    explain_exp_div.style.display = "inline";
    explain_exp_div.style.display = "table-cell";
    explain_exp_div.style.width = "350px";
    cp.appendChild(explain_exp_div);

    // Second Row: Showing the examples
    const cp2 = document.createElement("div");
    cp2.style.display = "table-row";


    if (false) {
        // Top Left -- How well can we explain predictions?
        const examples_div = document.createElement("div");
        const examples_label = document.createElement("div");
        examples_label.innerText = "Molecular explanation examples:";
        const examples_img = document.createElement('img');

        const col3 = "positive_examples_Atom Attributions-Training.png";
        examples_img.src = "data:image/png;base64," + IMAGES[col3];
        styleImage(examples_img);
        examples_img.style.width = "600px";

        examples_div.appendChild(examples_label);
        const slider_div = document.createElement("div");
        slider_div.innerHTML = ` <input type="range" style="width: 200px;" id="saturation-slider" min="0" max="10000" value="100">`;
        const slider = slider_div.children[0];

        slider.addEventListener('input', function () {
            const saturationValue = slider.value;

            examples_img.style.filter = `saturate($ {
                            saturationValue
                        }

                        %)`;
        });
        examples_div.appendChild(slider_div);
        examples_div.appendChild(examples_img);

        examples_div.style.display = "inline";
        examples_div.style.display = "table-cell";
        examples_div.style.width = "350px";
        cp2.appendChild(examples_div);
    }

    tab.append(cp);

    document.getElementById("result-div").appendChild(tab);
    document.getElementById("result-div").appendChild(cp2);
};

// An example showing how one could display molecules
// in a grid. Not needed any more, just for reference
// pagination handling of images
PAGE_SIZE = 4
CURRENT_PAGE = 0;
MOLECULE_IMAGES = [];

const refreshSecondPage = function () {
    const cp = document.getElementById("result-div-2");
    cp.innerHTML = "";
    const startP = CURRENT_PAGE * PAGE_SIZE;

    for (var q = 0; q < PAGE_SIZE; q++) {
        const image = MOLECULE_IMAGES[q + startP];
        const imgElement = document.createElement('img');
        //cp.innerHTML += `<img id="pngImage" alt="Base64 Image" />`;
        //document.getElementById("pngImage").src = "data:image/png;base64," + image;
        imgElement.src = "data:image/png;base64," + image;
        styleImage(imgElement);
        cp.appendChild(imgElement);
    }
}


window.onload = () => {
    var coll = document.getElementsByClassName("collapsible");
    var i;
    for (i = 0; i < coll.length; i++) {
        var elt = coll[i];
        elt.addEventListener(
            "click",
            function () {
                elt.classList.toggle("active");
                var ext = document.getElementById("left-column-ext");
                var min = document.getElementById("left-column-min");
                if (ext.style.display === "flex") {
                    ext.style.display = "none";
                    min.style.display = "flex";
                } else {
                    ext.style.display = "flex";
                    min.style.display = "none";
                }
            }
        );
    };

    // initialize state properly
    coll[0].click();

    coll = document.getElementsByClassName("collapsiblex");
    i = 0;

    const contents = { "Jobs": "job-id-div", "Overview": "result-div", "Evaluation": "result-div-2" }
    for (i = 0; i < coll.length; i++) {
        const elt = coll[i];
        const clicked_on = elt.innerText;
        elt.addEventListener("click", function () {
            Array.from(coll).forEach(et => et.classList.remove("collapsiblex-active"));
            elt.classList.add("collapsiblex-active");
            for (const [key, value] of Object.entries(contents)) {
                const other_elt = document.getElementById(value);
                if (key === clicked_on) {
                    other_elt.style.display = "block";
                } else {
                    other_elt.style.display = "none";
                }
            }
        });
    }

    // initialize state properly
    coll[0].click();
};
