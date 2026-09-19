"""
Export the FULL (uncapped) network data behind figures/05_network_maps.py,
for the interactive website (Nandini_Website/).

Why this script exists: 05_network_maps.py only ever draws a small curated
subset of each compartment's network (union of top-15 WOI-trajectory and
top-15 RIF-dysregulated driver genes, ~15-25 genes) because a full network
(~2,000 edges among ~2,000-2,500 modeled genes) is unreadable as a static
PNG. The full network was never exported to a plain file, though — it only
exists inside the pickled results/celloracle/*.celloracle.links objects,
which require the `celloracle` Python package to open. This script loads
each compartment's .celloracle.links once and dumps everything to flat
CSV/JSON so a web developer can build an interactive graph that shows AS
MANY nodes/edges as they want, not just the curated subset in the paper
figures.

Run once, from the pipeline root, inside the `celloracle_grn` conda env:
    conda activate celloracle_grn
    python3 Nandini_Website/code/export_full_network_data.py

Reads only already-computed results (results/celloracle/*.celloracle.links,
results/celloracle/base_GRN_hg38.parquet, results/celloracle/comparisons/*,
results/dyngenie3/*_importances.csv). Does not touch the raw or processed
AnnData objects and does not re-fit anything, so it's cheap to re-run if the
upstream CellOracle results are ever regenerated.

Output layout, all under Nandini_Website/data/:
  tf_list.csv                          Every candidate TF in the base GRN (gene, is_TF=True).
                                        Any gene NOT in this list is "target-only" (can be
                                        regulated, cannot itself be a regulator in this framework).
  <compartment>/
    network_snapshots.csv/.json        Full per-cluster edge list (~2,000 edges x 6 clusters),
                                        i.e. exactly what feeds the curated subgraph in
                                        05_network_maps.py, but for EVERY gene pair, not just
                                        the curated hub set.
    node_centrality.csv/.json          Full per-cluster, per-gene centrality table (every gene
                                        that appears as a node in that cluster's network), with
                                        an is_TF flag added.
    curated_hub_genes.csv              Exactly the ~15-25 genes shown in the static PNG for this
                                        compartment, with a is_TF flag -- so the website can
                                        reproduce "what the paper figure shows" as a default view
                                        before a user expands to the full network.
    consecutive_timepoint_edge_diffs.csv / driver_diffs.csv
                                        Copied from results/celloracle/comparisons/ -- gained/
                                        lost/retained edges and centrality-shift ranking for each
                                        of the 4 consecutive Fertile timepoint pairs, uncapped.
    fertile_vs_rif_edge_diff.csv / driver_diff.csv
                                        Same, for the Fertile-LH+7 vs RIF-LH+7 comparison.
    dyngenie3_importances_full.csv     Copied as-is from results/dyngenie3/ (every candidate
                                        regulator x target pair dynGENIE3 scored -- can be 1M+
                                        rows; see README for why you'll likely want to filter
                                        this before shipping it to a browser).
    dyngenie3_importances_top100.csv   Convenience cut: top 100 targets per regulator by
                                        importance score (still the FULL gene universe, just
                                        pre-filtered to a browser-friendly size). Not present for
                                        Tcell (no dynGENIE3 model -- see README).
"""
import json
import os
import shutil

import celloracle as co
import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.chdir(REPO_ROOT)

OUT_DIR = "Nandini_Website/data"
COMPARTMENTS = ["epithelial", "stromal", "tcell", "unk", "activatednk"]
COMPARTMENTS_NO_DYNGENIE3 = {"tcell"}
FERTILE_TIMEPOINT_ORDER = ["LH+3", "LH+5", "LH+7", "LH+9", "LH+11"]
TOP_N_DRIVERS = 15  # must match src/figures/05_network_maps.py's TOP_N_DRIVERS

os.makedirs(OUT_DIR, exist_ok=True)


def get_tf_set():
    base_grn = pd.read_parquet("results/celloracle/base_GRN_hg38.parquet")
    return set(c for c in base_grn.columns if c not in ("peak_id", "gene_short_name"))


def get_curated_hubs(compartment):
    woi = pd.read_csv(f"results/celloracle/comparisons/{compartment}_top_WOI_regulators_integrated.csv", index_col=0)
    rif = pd.read_csv(f"results/celloracle/comparisons/{compartment}_top_RIF_regulators_integrated.csv", index_col=0)
    woi_set = set(woi.head(TOP_N_DRIVERS).index)
    rif_set = set(rif.head(TOP_N_DRIVERS).index)
    return sorted(woi_set | rif_set), woi_set, rif_set


def df_to_csv_and_json(df, path_no_ext):
    df.to_csv(f"{path_no_ext}.csv", index=False)
    df.to_json(f"{path_no_ext}.json", orient="records")


def main():
    tf_set = get_tf_set()
    tf_df = pd.DataFrame(sorted(tf_set), columns=["gene"])
    tf_df["is_TF"] = True
    tf_df.to_csv(f"{OUT_DIR}/tf_list.csv", index=False)
    print(f"Wrote {OUT_DIR}/tf_list.csv ({len(tf_df)} candidate TFs)")

    for compartment in COMPARTMENTS:
        print(f"\n{'='*60}\n{compartment}\n{'='*60}")
        comp_dir = f"{OUT_DIR}/{compartment}"
        os.makedirs(comp_dir, exist_ok=True)

        # --- full per-cluster edge snapshots (the genuinely new export) ---
        links = co.load_hdf5(file_path=f"results/celloracle/{compartment}.celloracle.links")
        all_edges = []
        for cluster, df in links.filtered_links.items():
            d = df[["source", "target", "coef_mean", "coef_abs", "p", "-logp"]].copy()
            d.insert(0, "cluster", cluster)
            all_edges.append(d)
        edges_df = pd.concat(all_edges, ignore_index=True)
        edges_df["source_is_TF"] = edges_df["source"].isin(tf_set)
        edges_df["target_is_TF"] = edges_df["target"].isin(tf_set)
        df_to_csv_and_json(edges_df, f"{comp_dir}/network_snapshots")
        print(f"  network_snapshots: {len(edges_df)} edges across {edges_df['cluster'].nunique()} clusters")

        # --- full per-cluster node centrality table ---
        links.get_network_score()
        nodes_df = links.merged_score.reset_index().rename(columns={"index": "gene"})
        nodes_df["is_TF"] = nodes_df["gene"].isin(tf_set)
        df_to_csv_and_json(nodes_df, f"{comp_dir}/node_centrality")
        print(f"  node_centrality: {len(nodes_df)} (gene, cluster) rows")

        # --- curated hub set actually shown in the static PNG ---
        hubs, woi_set, rif_set = get_curated_hubs(compartment)
        hub_df = pd.DataFrame({"gene": hubs})
        hub_df["is_TF"] = hub_df["gene"].isin(tf_set)
        hub_df["in_top_WOI_drivers"] = hub_df["gene"].isin(woi_set)
        hub_df["in_top_RIF_drivers"] = hub_df["gene"].isin(rif_set)
        hub_df.to_csv(f"{comp_dir}/curated_hub_genes.csv", index=False)
        print(f"  curated_hub_genes: {len(hub_df)} genes ({hub_df['is_TF'].sum()} TFs) -- "
              f"this is what the static PNG shows")

        # --- copy already-existing comparison CSVs (uncapped, no rerun needed) ---
        comp_src = "results/celloracle/comparisons"
        for tp_a, tp_b in zip(FERTILE_TIMEPOINT_ORDER[:-1], FERTILE_TIMEPOINT_ORDER[1:]):
            ca, cb = f"Fertile_{tp_a}", f"Fertile_{tp_b}"
            for kind, dest in [("edges", "consecutive_timepoint_edge_diffs"), ("drivers", "consecutive_timepoint_driver_diffs")]:
                src = f"{comp_src}/{compartment}_{ca}_vs_{cb}_{kind}.csv"
                if os.path.exists(src):
                    dst_dir = f"{comp_dir}/{dest}"
                    os.makedirs(dst_dir, exist_ok=True)
                    shutil.copy2(src, f"{dst_dir}/{ca}_vs_{cb}.csv")
        for kind, dest in [("edges", "fertile_vs_rif_edge_diff.csv"), ("drivers", "fertile_vs_rif_driver_diff.csv")]:
            src = f"{comp_src}/{compartment}_Fertile_LH+7_vs_RIF_LH+7_{kind}.csv"
            if os.path.exists(src):
                shutil.copy2(src, f"{comp_dir}/{dest}")
        print("  copied consecutive-timepoint and Fertile-vs-RIF edge/driver diff CSVs")

        # --- dynGENIE3 (skip for Tcell -- no pseudotime model) ---
        if compartment in COMPARTMENTS_NO_DYNGENIE3:
            print("  no dynGENIE3 results for this compartment (pseudotime validation gate failed)")
            continue
        dyn_src = f"results/dyngenie3/{compartment}_importances.csv"
        shutil.copy2(dyn_src, f"{comp_dir}/dyngenie3_importances_full.csv")
        dyn = pd.read_csv(dyn_src)
        dyn_top = (
            dyn.sort_values("importance", ascending=False)
            .groupby("regulator", group_keys=False)
            .head(100)
            .sort_values(["regulator", "importance"], ascending=[True, False])
        )
        dyn_top.to_csv(f"{comp_dir}/dyngenie3_importances_top100.csv", index=False)
        print(f"  dynGENIE3: {len(dyn)} full rows copied, {len(dyn_top)} rows in the top-100-per-regulator cut")

    print("\nExport complete.")


if __name__ == "__main__":
    main()
