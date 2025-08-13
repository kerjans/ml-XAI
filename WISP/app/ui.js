






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

                const any_errors = res["log_err"].length > 0;
                var status = null;
                if (any_errors) {
                    status = "failed";
                }
                else {
                    status = res["job_status"];
                }
                const status_label = document.createElement("p");
                status_label.style.display = "inline";
                status_label.style.marginLeft = "10px";
                status_label.innerText = status;
                bd.appendChild(status_label);

                const closeButton = document.createElement('button');
                closeButton.innerText = "close";
                closeButton.classList.add("closeButton");

                closeButton.addEventListener('click', () => {
                    overlay.style.display = 'none';
                    document.getElementById("overlayText").innerText = "";
                });

                const overlay = document.getElementById('overlay');

                const openLog = document.createElement('button');
                openLog.classList.add("logButton");

                openLog.innerText = "log";
                openLog.addEventListener('click', () => {
                    const p = document.createElement("p");
                    p.innerText = res["log_out"];
                    overlay.style.display = 'flex';
                    document.getElementById("overlayText").appendChild(p);
                    document.getElementById("overlayText").appendChild(closeButton);

                });

                const openErr = document.createElement('button');

                if (any_errors) {
                    openErr.classList.add("errButton-active");
                }
                else {
                    openErr.classList.add("errButton-inactive");
                }

                openErr.innerText = "err";
                openErr.addEventListener('click', () => {
                    const p = document.createElement("p");
                    p.innerText = res["log_err"];
                    overlay.style.display = 'flex';
                    document.getElementById("overlayText").appendChild(p);
                    document.getElementById("overlayText").appendChild(closeButton);

                });


                bd.appendChild(openLog);
                bd.appendChild(openErr);

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
            "csv": CURRENT_DATA,
            "smiles_col": document.getElementById("smiles-column").value,
            "target_col": document.getElementById("target-column").value
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

CURRENT_JOB_ID = 0;
MMP_OVERVIEW_DATA = [];
MODEL_PERF_OVERVIEW = "";
const retrieveResults = function (job_id) {
    CURRENT_JOB_ID = job_id;
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
            MODEL_PERF_OVERVIEW = res["model_perf_overview"];
            IMAGES = res["images"];
            //MOLECULE_IMAGES = res["molecule_images"];

            refreshFirstPage();
            //refreshSecondPage();
        }

        else {
            alert("failure could not parse input!");
        }
    });


    fetch("FeatureImportance", {
        method: "POST",
        headers: {
            "Content-Type": "application/json"
        },
        body: JSON.stringify({
            "job_id": job_id
        })
    }
    ).then(
        res => res.json()
    ).then(
        res => {
            console.log("Request complete! response:", res);
            const status = res["status"];

            if (status == "success") {
                const elementId = "result-div-4";
                d3.select(`#${elementId}`).html(res["feature_importance_plot"]);
            }

        }
    );
    fetch("MMPOverview", {
        method: "POST",
        headers: {
            "Content-Type": "application/json"
        },
        body: JSON.stringify({
            "job_id": job_id
        })
    }).then(res => res.json()).then(res => {
        console.log("Request complete! response:", res);
        const status = res["status"];

        if (status == "success") {
            //IMAGES = res["images"];
            MMP_OVERVIEW_DATA = res["mmp_overview_data"];
            if (typeof MMP_OVERVIEW_DATA === 'string' || MMP_OVERVIEW_DATA instanceof String)
                MMP_OVERVIEW_DATA = JSON.parse(MMP_OVERVIEW_DATA);

            //refreshFirstPage();
            //refreshSecondPage();
            renderMMPOverview(MMP_OVERVIEW_DATA);
        }

        else {
            alert("failure could not parse input!");
        }
    });
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

    // const col = "PREDvsCONTRIBUTIONSfragmentAtom Attributions_Training_Set.png";
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

    // const col2 = "EXPvsCONTRIBUTIONSwholeAtom Attributions_Test_Set.png";
    const col2 = "EXPvsCONTRIBUTIONSwholeAtom Attributions_Test_Set.png";
    explain_exp_img.src = "data:image/png;base64," + IMAGES[col2];
    styleImage(explain_exp_img);

    explain_exp_div.appendChild(explain_exp_label);
    explain_exp_div.appendChild(explain_exp_img);

    explain_exp_div.style.display = "inline";
    explain_exp_div.style.display = "table-cell";
    explain_exp_div.style.width = "350px";
    cp.appendChild(explain_exp_div);

    // Second Row: Showing the predictive
    // performance on test data
    const cp2 = document.createElement("div");
    cp2.style.display = "table-row";


    const pred_test_div = document.createElement("div");
    const pred_test_label = document.createElement("div");
    pred_test_label.innerText = "Predictive performance on test set";
    const pred_test_img = document.createElement('img');

    pred_test_img.src = "data:image/png;base64," + IMAGES["20-80-split-true-pred.png"];
    styleImage(pred_test_img);

    pred_test_div.appendChild(pred_test_label);
    pred_test_div.appendChild(pred_test_img);

    pred_test_div.style.display = "inline";
    pred_test_div.style.display = "table-cell";
    pred_test_div.style.width = "350px";
    cp2.appendChild(pred_test_div);

    const tab_model_overview = document.createElement("div");
    tab_model_overview.innerHTML = MODEL_PERF_OVERVIEW;
    cp2.appendChild(tab_model_overview);

    tab.append(cp);
    tab.append(cp2);

    document.getElementById("result-div").appendChild(tab);
};

// An example showing how one could display molecules
// in a grid. Not needed any more, just for reference
// pagination handling of images
PAGE_SIZE = 16
COL_SIZE = 4
CURRENT_PAGE = 0;
N_PAGES = 0;
MOLECULE_IMAGES = [];
LEGEND_IMAGE = null;


const refreshSecondPage = function () {
    const job_id = CURRENT_JOB_ID;

    fetch("HeatMaps", {

        method: "POST",
        headers: {
            'Content-Type': 'application/json'
        },
        body: JSON.stringify({
            "job_id": job_id,
            "current_page": CURRENT_PAGE,
            "page_size": PAGE_SIZE
        })
    }).then(res => res.json()).then(res => {
        console.log("Request complete! response:", res);
        const status = res["status"];

        if (status == "success") {

            //IMAGES = res["images"];
            MOLECULE_IMAGES = res["heatmaps"];
            LEGEND_IMAGE = res["legend"];
            N_PAGES = Math.ceil(res["n_entries"] / PAGE_SIZE);

            //refreshFirstPage();
            refreshSecondPageElements();
        }

        else {
            alert("failure could not parse input!");
        }
    }
    );
}

const refreshSecondPageElements = function () {
    const cp = document.getElementById("result-div-2");
    cp.innerHTML = "";
    const startP = CURRENT_PAGE * PAGE_SIZE;

    const div_gallery = document.createElement("div");

    const but_div = document.createElement("div");
    const but_left = document.createElement("button");
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
        const image = MOLECULE_IMAGES[q];
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

    const contents = { "Jobs": "job-id-div", "Evaluation": "result-div", "Atom Contributions": "result-div-2", "MMP Overview": "result-div-3", "Feature Importance": "result-div-4" }
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
                    if (key == "Atom Contributions") {
                        refreshSecondPage();
                    }
                } else {
                    other_elt.style.display = "none";
                }
            }
        });
    }

    // initialize state properly
    coll[0].click();
};

RDKit = null;

initRDKitModule().then(function (instance) {
    RDKit = instance;
    console.log("RDKit loaded.");
});

const displaySmiles = function (smiles) {
    const mol = RDKit.get_mol(smiles);
    return mol.get_svg();
};


const renderMMPOverview = function (data) {

    const svg = d3.select("#mmp-overview");
    svg.selectAll("*").remove(); // Clear the canvas

    //const margin = { top: 20, right: 30, bottom: 30, left: 120 };
    const margin = { top: 40, right: 60, bottom: 60, left: 200 };



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
    var groups = Array.from(
        d3.group(data, d => d.MMP_rule),
        ([key, values]) => ({
            MMP_rule: key,
            values: values.map(d => d.target_diff),
            smiles_1: values.map(d => d.smiles_1),
            smiles_2: values.map(d => d.smiles_2),
        })
    );


    // Sort alphabetically
    // groups.sort((a, b) => d3.ascending(a.MMP_rule, b.MMP_rule));
    // Sort by value
    groups.sort((a, b) => d3.mean(a.values) - d3.mean(b.values));

    // Keep only the top/bottom 10 elements
    const top10 = groups.slice(-10); // highest means
    const bottom10 = groups.slice(0, 10); // lowest means

    // Merge and re-sort however you want (e.g., descending by mean)
    const filter_top = false;
    if (filter_top)
        groups = [...bottom10, ...top10];

    const rowHeight = 20;
    const height = groups.length * rowHeight + margin.top + margin.bottom;

    // Update the SVG height
    svg.attr("height", height);

    const width = +svg.attr("width");
    //const height = +svg.attr("height");

    const allDiffs = data.map(d => d.target_diff);

    const x = d3.scaleLinear()
        .domain(d3.extent(allDiffs))
        .range([margin.left, width - margin.right]);

    const y = d3.scaleBand()
        .domain(groups.map(g => g.MMP_rule))
        .range([margin.top, height - margin.bottom])
        .paddingInner(0.3);

    // This one here gave poor results:
    // const kde = kernelDensityEstimator(kernelEpanechnikov(0.02), x.ticks(100));
    const kde = kernelDensityEstimator(kernelGaussian(0.25), x.ticks(100));

    const maxDensity = d3.max(groups, g => d3.max(kde(g.values), d => d[1]));

    const AMPLIFICATION_FACTOR = 1.0;
    const yScale = d3.scaleLinear()
        .domain([0, maxDensity])
        .range([0, y.bandwidth() * AMPLIFICATION_FACTOR]);

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
    groups.forEach((g, i) => {
        const density = kde(g.values);

        const area = d3.area()
            .x(d => x(d[0]))
            .y0(y(g.MMP_rule) + y.bandwidth() / 2)
            .y1(d => y(g.MMP_rule) + y.bandwidth() / 2 - yScale(d[1]))
            .curve(d3.curveBasis);

        svg.append("path")
            .datum(density)
            .attr("class", "area")
            .attr("fill", i % 2 === 0 ? "#2452baff" : "#a2bff4ff") // alternate shades
            .attr("d", area);

        svg.append("text")
            .attr("x", margin.left - 10)
            .attr("y", y(g.MMP_rule) + y.bandwidth() / 2)
            .attr("text-anchor", "end")
            .attr("class", "group-label")
            .text(g.MMP_rule);

        function minDistance(value, array) {
            return array.reduce((minDist, current) => {
                const dist = Math.abs(value - current);
                return Math.min(minDist, dist);
            }, Infinity);
        }

        const jitter_with_memo_and_span = function (val, memo, span_sub, span_add) {
            const jit_min = -span_sub;
            const jit_max = +span_add;
            const thresh = 5;
            const cur_dist = minDistance(val, memo);

            if (cur_dist < thresh) {
                const jit_amnt = Math.random() * (jit_max - jit_min) + jit_min;
                const jittered_val = val + jit_amnt;
                memo.push(jittered_val);
                return jittered_val;
            } else {
                memo.push(val);
                return val;
            }
        };

        var _jitterx_memo = [];
        var _jittery_memo = [];
        const jitterx = function (val) {
            return jitter_with_memo_and_span(val, _jitterx_memo, 10, 10);
        }
        const jittery = function (val) {
            return jitter_with_memo_and_span(val, _jittery_memo, 10, 0);
        }


        // Scatterplot the single datapoints
        g.values.forEach((value, i) => {
            svg.append("circle")
                .attr("cx", jitterx(x(value)))
                .attr("cy", jittery(y(g.MMP_rule)) + y.bandwidth() / 2) // vertically center
                .attr("r", 3)
                .attr("fill", "black")
                .attr("opacity", 0.6)
                .on("mouseover", function (event) {
                    d3.select(this)
                        .transition()
                        .duration(150)
                        .attr("r", 6)
                        .attr("fill", "red");

                    const toolt = d3.select("#tooltip")
                        .style("opacity", 1)
                        .html(`<p style="margin:0px; padding:0px;">target_diff: ${value.toFixed(3)}</p> ${displaySmiles(g.smiles_1[i])} <p class="mmp-arrow"> => </p> ${displaySmiles(g.smiles_2[i])}`)
                        .style("left", `${event.pageX + 10}px`)
                        .style("top", `${event.pageY - 20}px`);
                    toolt.selectAll("svg").style("width", "150px").style("height", "150px").style("margin", "5px");
                    toolt.selectAll(".mmp-arrow").style("position", "absolute").style("top", "50%").style("left", "48%");
                })
                .on("mouseout", function () {
                    d3.select(this)
                        .transition()
                        .duration(150)
                        .attr("r", 3)
                        .attr("fill", "black");

                    d3.select("#tooltip")
                        .style("opacity", 0);
                });;
        });
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
    function kernelGaussian(bandwidth) {
        return function (v) {
            return (1 / (Math.sqrt(2 * Math.PI) * bandwidth)) * Math.exp(-0.5 * (v / bandwidth) ** 2);
        };
    }
};

window.addEventListener("load", (event) => {
    renderJobs();
    renderMMPOverview();

});




