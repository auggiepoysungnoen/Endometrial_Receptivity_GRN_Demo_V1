# Endometrial GRN Explorer · demo_v1

Interactive map of the **stromal-cell gene regulatory network** of the human
endometrium, comparing **fertile** vs **recurrent implantation failure (RIF)**
at **LH+7** (window of implantation). Networks inferred with CellOracle.

**Live demo:** https://auggiepoysungnoen.github.io/Endometrial_Receptivity_GRN_Demo_V1/

> **demo_v1** is a review release, not for citation. Links are model
> predictions, not experiments; literature annotations are collected
> automatically from NCBI / CollecTRI, and the 36 gene notes are drafts
> pending lab review (see `local_host_development/curation/`).

## What's here

| Path | Contents |
|---|---|
| `local_host_development/site/` | The website (static HTML/JS; published by GitHub Pages) |
| `local_host_development/script/` | Build scripts (territory layout, literature annotations), validation suite, local server |
| `local_host_development/curation/` | Draft implantation notes + review workflow |
| `local_host_development/workflow/` | Design, functionality, progress log, validation report |
| `METHODS_AND_DATA_GUIDE.md` | Methods and data dictionary |

The exported analysis data (`data/`) and paper figures (`figures/`) are **not**
in this repository; the build scripts need them and are kept for provenance.

## Run locally

```bash
cd local_host_development/script && ./serve_local.sh
# open http://127.0.0.1:8000/local_host_development/site/
```

## Validation (demo_v1)

`script/validate_site.py`: 44/44 data checks against the original exports ·
`script/tests/browser_smoke.js`: 40/40 in-browser checks. Report:
`local_host_development/workflow/validation_report.txt`.

---

# Nandini_Website — handoff package for the interactive GRN website

This folder is a self-contained package for building a GitHub-deployed
website with an interactive version of the gene-regulatory-network (GRN)
figures in `figures/05_network_maps/`. It has everything needed to build
that site **without rerunning any of the analysis pipeline** (no data
cleaning, no CellOracle/dynGENIE3 refitting) — all of it is already-computed
output, exported or copied into plain CSV/JSON.

**Start here:** [`METHODS_AND_DATA_GUIDE.md`](METHODS_AND_DATA_GUIDE.md) —
the full write-up: what this project is, how the two network-inference
methods (CellOracle vs. dynGENIE3) work and differ, exactly how the static
figures pick which genes/edges to draw (and why they cap it), and a
complete data dictionary for every file in `data/`.

## What's in this folder

- **`code/`**
  - `05_network_maps.py` — the actual script that generated the static PNGs
    in `figures/05_network_maps/`, copied unmodified. This is the ground
    truth for exactly how node size/color, edge width/color, and the
    curated-gene cap work — read it alongside the guide.
  - `_style.py` — the shared color palette / style config that script
    imports (node/edge colors referenced in the guide come from here).
  - `export_full_network_data.py` — the new script that generated
    everything in `data/`. Re-run it (inside the `celloracle_grn` conda env,
    from the repo root) any time the upstream CellOracle/dynGENIE3 results
    change and you want to refresh this folder's data.

- **`data/`** — one subfolder per compartment (`epithelial/`, `stromal/`,
  `tcell/`, `unk/`, `activatednk/`), each with the **full, uncapped**
  network data (not just the ~15-26 genes shown in the static PNG), plus a
  shared `tf_list.csv`. Full column-by-column description in the guide.

- **`figures/05_network_maps/`** — copies of the 14 published static PNGs
  (WOI-trajectory, Fertile-vs-RIF, and dynGENIE3 network maps for all 5
  compartments) — visual reference for what the current print figure looks
  like and what the interactive version is meant to improve on.

## Key thing to know before building

The static PNGs deliberately show only a small curated subset of each
compartment's network (a readability constraint for print, not a data
limitation). The data in this folder is **not** capped that way — it's the
full network CellOracle and dynGENIE3 actually fit (~2,000 edges per
compartment per timepoint/condition group, among ~2,000-2,500 modeled
genes). The interactive site can and should show more than the static
figure does; the guide explains exactly what's available and how the
original figure chose its subset, so the site can offer "the paper's
highlighted genes" as a sensible default view with an option to expand to
the full network.
