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

const renderJobs = function () {
    const elt = document.getElementById("job-list");
    elt.innerHTML = "";
    const jobs = JSON.parse(localStorage.getItem("jobs"));
    jobs.forEach(
        function (job) {
            const b = document.createElement("button");
            b.classList.add("primary-button");
            b.classList.add("job-buttons");
            b.innerText = job;
            b.onclick = function (args) {
                Array.from(document.getElementsByClassName("job-buttons")).forEach(elt => elt.classList.remove("active-job-button"));
                b.classList.add("active-job-button");
                retrieveResults(job);
            };
            const bd = document.createElement("div");
            bd.appendChild(b);

            fetch("JobStatus", {
                method: "POST",
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    "job_id": job
                })

            }).then(res => res.json()).then(res => {
                console.log("Request complete! response:", res);
                const status = res["job_status"];
                const status_label = document.createElement("p");
                status_label.style.display = "inline";
                status_label.style.marginLeft = "10px";
                status_label.innerText = status;
                bd.appendChild(status_label);
                elt.appendChild(bd);
            });
        }
    );
};


function uniq(a) {
    return a.sort().filter(function (item, pos, ary) {
        return !pos || item != ary[pos - 1];
    });
}

const addJob = function (job_id) {
    jobs = JSON.parse(localStorage.getItem("jobs"));
    if (jobs === null) {
        jobs = [];
    }
    jobs.push(job_id);
    jobs = uniq(jobs).reverse();
    localStorage.setItem("jobs", JSON.stringify(jobs));
};

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
            addJob(res["job_id"]);
            renderJobs();
            but.disabled = false;
            but.style.opacity = 1.0;
        }

        else {
            alert("failure could not parse input!");
        }
    });
};

const retrieveResults = function (job_id) {

    fetch("WispOverviewPage", {

        method: "POST",
        headers: {
            'Content-Type': 'application/json'
        }

        ,
        body: JSON.stringify({
            "job_id": job_id
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
                "job_id": job_id
            })
        }).then(res => res.json()).then(res => {
            console.log("Request complete! response:", res);
            const status = res["status"];

            if (status == "success") {

                //IMAGES = res["images"];
                MOLECULE_IMAGES = res["heatmaps"];
                LEGEND_IMAGE = res["legend"];

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
    explain_pred_img.classList.add("zoomable");

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

    explain_exp_img.classList.add("zoomable");

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
PAGE_SIZE = 16
COL_SIZE = 4
CURRENT_PAGE = 0;
MOLECULE_IMAGES = [];
LEGEND_IMAGE = null;

const refreshSecondPage = function () {
    const cp = document.getElementById("result-div-2");
    cp.innerHTML = "";
    const startP = CURRENT_PAGE * PAGE_SIZE;

    const div_gallery = document.createElement("div");

    const but_div = document.createElement("div");
    const but_left = document.createElement("button");
    const N_PAGES = Math.ceil(MOLECULE_IMAGES.length / PAGE_SIZE);
    but_left.id = "left-button";
    but_left.innerText = "<";
    const but_right = document.createElement("button");
    but_right.id = "right-button";
    but_right.innerText = ">";

    const page_label = document.createElement("p");
    page_label.innerText = `${CURRENT_PAGE + 1} / ${N_PAGES}`;
    page_label.style.color = "black";
    page_label.style.display = "inline";

    but_left.classList.add("nav-button");
    but_right.classList.add("nav-button");

    but_left.onclick = function (evt) {
        CURRENT_PAGE = Math.max(0, CURRENT_PAGE - 1);
        refreshSecondPage();
    };
    but_right.onclick = function (evt) {
        CURRENT_PAGE = Math.min(N_PAGES - 1, CURRENT_PAGE + 1);
        refreshSecondPage();
    };

    but_div.style.position = "absolute";
    but_div.style.top = "-2%";
    but_div.style.right = "-2%";
    but_div.appendChild(but_left);
    but_div.appendChild(page_label);
    but_div.appendChild(but_right);
    div_gallery.appendChild(but_div);

    const imgLeg = document.createElement('img');
    imgLeg.src = "data:image/png;base64," + LEGEND_IMAGE;
    imgLeg.style.width = "150px";
    imgLeg.classList.add("zoomable");

    div_gallery.appendChild(imgLeg)

    const elt_ul = document.createElement("ul");

    elt_ul.classList.add("gallery");
    for (var q = 0; q < PAGE_SIZE; q++) {
        const image = MOLECULE_IMAGES[q + startP];
        const imgElement = document.createElement('img');
        //cp.innerHTML += `<img id="pngImage" alt="Base64 Image" />`;
        //document.getElementById("pngImage").src = "data:image/png;base64," + image;
        if (image) {
            imgElement.src = "data:image/png;base64," + image;
        }
        //styleImage(imgElement);
        imgElement.style.width = "150px";

        const elt_li = document.createElement("li");
        elt_li.appendChild(imgElement)
        elt_ul.appendChild(elt_li)
    }
    div_gallery.appendChild(elt_ul);
    cp.appendChild(div_gallery);
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

    const contents = { "Jobs": "job-id-div", "Evaluation": "result-div", "Atom Contributions": "result-div-2", "MMP Overview": "result-div-3" }
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



const renderMMPOverview = function (data) {

    const svg = d3.select("#mmp-overview");
    const width = +svg.attr("width");
    const height = +svg.attr("height");
    const margin = { top: 20, right: 30, bottom: 30, left: 120 };


    if (!data)
        // Toy data: array of objects with target_diff and MMP_rule
        data = [
            { target_diff: 0.2, MMP_rule: 'Rule A' },
            { target_diff: 0.4, MMP_rule: 'Rule A' },
            { target_diff: 0.35, MMP_rule: 'Rule A' },
            { target_diff: 0.5, MMP_rule: 'Rule A' },
            { target_diff: 0.1, MMP_rule: 'Rule B' },
            { target_diff: 0.12, MMP_rule: 'Rule B' },
            { target_diff: 0.15, MMP_rule: 'Rule B' },
            { target_diff: 0.18, MMP_rule: 'Rule B' },
            { target_diff: 0.3, MMP_rule: 'Rule C' },
            { target_diff: 0.28, MMP_rule: 'Rule C' },
            { target_diff: 0.32, MMP_rule: 'Rule C' },
            { target_diff: 0.34, MMP_rule: 'Rule C' }
        ];

    // Group by MMP_rule
    const groups = Array.from(
        d3.group(data, d => d.MMP_rule),
        ([key, values]) => ({ MMP_rule: key, values: values.map(d => d.target_diff) })
    );

    // Sort alphabetically
    groups.sort((a, b) => d3.ascending(a.MMP_rule, b.MMP_rule));

    const allDiffs = data.map(d => d.target_diff);

    const x = d3.scaleLinear()
        .domain(d3.extent(allDiffs))
        .range([margin.left, width - margin.right]);

    const y = d3.scaleBand()
        .domain(groups.map(g => g.MMP_rule))
        .range([margin.top, height - margin.bottom])
        .paddingInner(0.3);

    const kde = kernelDensityEstimator(kernelEpanechnikov(0.02), x.ticks(100));

    const maxDensity = d3.max(groups, g => d3.max(kde(g.values), d => d[1]));

    const yScale = d3.scaleLinear()
        .domain([0, maxDensity])
        .range([0, y.bandwidth()]);

    // X Axis
    svg.append("g")
        .attr("transform", `translate(0,${height - margin.bottom})`)
        .call(d3.axisBottom(x).ticks(6))
        .append("text")
        .attr("x", width / 2)
        .attr("y", 25)
        .attr("fill", "black")
        .attr("text-anchor", "middle")
        .text("target_diff");

    // Draw the curves
    groups.forEach(g => {
        const density = kde(g.values);

        const area = d3.area()
            .x(d => x(d[0]))
            .y0(y(g.MMP_rule) + y.bandwidth() / 2)
            .y1(d => y(g.MMP_rule) + y.bandwidth() / 2 - yScale(d[1]))
            .curve(d3.curveBasis);

        svg.append("path")
            .datum(density)
            .attr("class", "area")
            .attr("d", area);

        svg.append("text")
            .attr("x", margin.left - 10)
            .attr("y", y(g.MMP_rule) + y.bandwidth() / 2)
            .attr("text-anchor", "end")
            .attr("class", "group-label")
            .text(g.MMP_rule);
    });

    // KDE functions
    function kernelDensityEstimator(kernel, X) {
        return function (V) {
            return X.map(x => [x, d3.mean(V, v => kernel(x - v))]);
        };
    }

    function kernelEpanechnikov(k) {
        return function (v) {
            return Math.abs(v /= k) <= 1 ? 0.75 * (1 - v * v) / k : 0;
        };
    }
};

window.addEventListener("load", (event) => {
    renderJobs();
    renderMMPOverview();
});