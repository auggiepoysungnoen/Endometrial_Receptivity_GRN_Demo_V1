"""
Figure section 5: Comprehensive (but curated) GRN network diagrams.

Design (see conversation / data/progress.txt for the reasoning): a literal
render of a full CellOracle network (2,000 edges among ~2,500 genes) is an
unreadable hairball. Instead of restricting to one hub's ego-network (the
previous approach), this shows a small, curated set of genes we've already
identified as important — the union of the top-20 WOI-trajectory driver genes
and top-20 RIF-dysregulated driver genes per compartment (TOP_N_DRIVERS; the
same top-20 used in Section 4's charts and results/summary_report.txt, so the
network figures agree with the rest of the paper) — and draws every
regulatory edge that exists BETWEEN genes in that set. Every gene shown is
one we already have independent evidence matters.

The SAME gene set and a SINGLE fixed node layout (computed once, from the
union of edges across every view below) are reused across all three figure
types in this section, so gene positions are visually consistent whether
you're looking at a CellOracle snapshot or the dynGENIE3 view — you can find
"ZEB1" in the same spot in every panel.

Three figures per compartment:
  A) woi_trajectory_network_{compartment}.png — 5-panel small multiples, one
     per Fertile timepoint (LH+3->LH+11), CellOracle. Node size = that
     timepoint's eigenvector centrality; edge width = |regression
     coefficient|. Fixed layout across panels turns this into a flipbook of
     one network rewiring over time, not five unrelated diagrams.
  B) fertile_vs_rif_network_{compartment}.png — 2-panel (Fertile-LH+7 vs
     RIF-LH+7), same encoding, edge color additionally shows gained/lost/
     retained status (extends the old single-hub ego-network to the full
     curated set).
  C) dyngenie3_network_{compartment}.png — 1 panel, continuous-dynamics
     counterpart. Node size = aggregate dynGENIE3 TF importance (summed
     across ALL of that gene's target edges in the full importance table,
     not just the curated subset). Edge width = that edge's dynGENIE3
     importance score. Capped to the top ~40 edges among the curated set (raw
     count there is 150-220 — dynGENIE3 considers many more candidate
     regulator-target pairs than CellOracle's base-GRN-constrained fitting,
     so it needs its own density cap to match).

Node fill distinguishes actual transcription factors (candidate regulators in
the base GRN — genes that can have outgoing edges) from target-only genes
(genes that showed up as important by centrality/importance but aren't
themselves candidate regulators in this framework) — this is a real and
informative distinction: about a third of the curated genes in each
compartment are TFs, the rest are convergence points regulated by them.
"""
import sys
sys.path.insert(0, "figures")
from _style import (apply_style, save_fig, DIVERGING_LOW, DIVERGING_HIGH,
                     INK_SECONDARY, INK_MUTED, TARGET_NODE_FILL, MUTED_NODE, TF_NODE_FILL)

import numpy as np
import pandas as pd
import celloracle as co
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

apply_style()
OUT_DIR = "figures/05_network_maps"
FERTILE_TIMEPOINT_ORDER = ["LH+3", "LH+5", "LH+7", "LH+9", "LH+11"]
TOP_N_DRIVERS = 15  # back up from 10 per user request
DYNGENIE3_EDGE_CAP = 40

NODE_SIZE_RANGE = (120, 1400)
EDGE_WIDTH_RANGE = (0.8, 6.0)
RETAINED_COLOR = "#c3c2b7"


def scale(values, out_range):
    values = np.asarray(values, dtype=float)
    lo, hi = values.min(), values.max()
    if hi <= lo:
        return np.full_like(values, out_range[0])
    return out_range[0] + (values - lo) / (hi - lo) * (out_range[1] - out_range[0])


def get_curated_hubs(compartment):
    woi = pd.read_csv(f"results/celloracle/comparisons/{compartment}_top_WOI_regulators_integrated.csv", index_col=0)
    rif = pd.read_csv(f"results/celloracle/comparisons/{compartment}_top_RIF_regulators_integrated.csv", index_col=0)
    return sorted(set(woi.head(TOP_N_DRIVERS).index) | set(rif.head(TOP_N_DRIVERS).index))


def get_tf_set():
    base_grn = pd.read_parquet("results/celloracle/base_GRN_hg38.parquet")
    return set(c for c in base_grn.columns if c not in ("peak_id", "gene_short_name"))


def celloracle_edges_within(links, cluster, hubs):
    df = links.filtered_links[cluster]
    sub = df[df["source"].isin(hubs) & df["target"].isin(hubs)].copy()
    sub["abs_coef"] = sub["coef_mean"].abs()
    return sub[["source", "target", "coef_mean", "abs_coef"]]


def dyngenie3_edges_within(compartment, hubs, cap):
    dyn = pd.read_csv(f"results/dyngenie3/{compartment}_importances.csv")
    within = dyn[dyn["regulator"].isin(hubs) & dyn["target"].isin(hubs)].copy()
    within = within.sort_values("importance", ascending=False).head(cap)
    within = within.rename(columns={"regulator": "source"})
    return within, dyn


def compute_fixed_layout(all_edges_union, hubs, tf_set):
    """Two-layer circular layout: TFs (regulators) on an inner ring, target-only
    genes on an outer ring, with a wide, guaranteed empty gap between the two
    rings. History: spring_layout collapsed TFs into an unreadable overlapping
    cluster; a compact centered grid (tried next, per an earlier "TF in the
    middle" request) read fine for a handful of TFs but crowded badly once a
    compartment had 13-14 TFs (uNK/ActivatedNK). Back to two concentric rings
    (per user follow-up), this time with an explicit minimum inner-to-outer
    gap so the two layers always read as visually distinct "layers with
    space between them," not just two touching circles. Targets are ordered
    to sit near whichever TF(s) regulate them most strongly, to keep edges
    reasonably short and reduce crossings.
    """
    tf_nodes = sorted(h for h in hubs if h in tf_set)
    target_nodes = sorted(h for h in hubs if h not in tf_set)

    # order target nodes by their primary (highest-weight) regulator's position
    # among tf_nodes, so they cluster near the TF that most influences them
    edge_weight = {}
    for u, v in all_edges_union:
        edge_weight.setdefault(v, []).append(u)
    tf_order_idx = {tf: i for i, tf in enumerate(tf_nodes)}

    def primary_tf_index(target):
        regs = [r for r in edge_weight.get(target, []) if r in tf_order_idx]
        return tf_order_idx[regs[0]] if regs else len(tf_nodes) / 2

    target_nodes = sorted(target_nodes, key=primary_tf_index)

    pos = {}
    # Inner ring (TFs): radius scales with TF count so per-node arc length
    # (and label spacing) stays roughly constant regardless of how many TFs
    # a given compartment has (3 for Epithelial, 14 for uNK).
    inner_r = max(0.95, 0.95 * len(tf_nodes) / 7)
    for i, tf in enumerate(tf_nodes):
        angle = 2 * np.pi * i / max(len(tf_nodes), 1)
        pos[tf] = np.array([inner_r * np.cos(angle), inner_r * np.sin(angle)])

    # Outer ring (targets): scales with target count, and always sits at
    # least GAP beyond the inner ring — a real empty band between the two
    # layers, not just "whatever the count-based formula happens to give."
    GAP = 1.6
    outer_r = max(2.4, 2.4 * len(target_nodes) / 22, inner_r + GAP)
    for i, tgt in enumerate(target_nodes):
        angle = 2 * np.pi * i / max(len(target_nodes), 1)
        pos[tgt] = np.array([outer_r * np.cos(angle), outer_r * np.sin(angle)])
    return pos


def draw_network(ax, pos, hubs, tf_set, edges_df, node_sizes, title, tf_color, edge_color_fn=None,
                  edge_color_default=RETAINED_COLOR, weight_col="abs_coef"):
    G = nx.DiGraph()
    G.add_nodes_from(hubs)
    for _, row in edges_df.iterrows():
        G.add_edge(row["source"], row["target"], weight=row[weight_col],
                    color=edge_color_fn(row) if edge_color_fn else edge_color_default)

    if len(edges_df):
        widths_lookup = dict(zip(zip(edges_df["source"], edges_df["target"]),
                                  scale(edges_df[weight_col], EDGE_WIDTH_RANGE)))
    else:
        widths_lookup = {}
    edge_widths = [widths_lookup.get((u, v), 1.0) for u, v in G.edges()]
    edge_colors = [G[u][v]["color"] for u, v in G.edges()]

    nx.draw_networkx_edges(ax=ax, G=G, pos=pos, edge_color=edge_colors, width=edge_widths,
                            alpha=0.8, arrows=True, arrowsize=8, arrowstyle="-|>",
                            connectionstyle="arc3,rad=0.05", node_size=[node_sizes.get(n, 200) for n in G.nodes()])

    # Nodes with no edge in THIS panel's network are grayed out, not hidden —
    # they're still real, independently-important curated genes (top
    # driver-shift/RIF-dysregulation genes over the full ~2,500-gene
    # network), they just aren't wired to another curated gene at this
    # specific timepoint/condition. Distinguishing them avoids implying every
    # node is equally embedded in the small subgraph shown here.
    degree = dict(G.degree())
    node_colors = []
    node_edgecolors = []
    for n in G.nodes():
        is_isolated = degree.get(n, 0) == 0
        node_colors.append(MUTED_NODE if is_isolated else (tf_color if n in tf_set else TARGET_NODE_FILL))
        node_edgecolors.append("white")
    nx.draw_networkx_nodes(ax=ax, G=G, pos=pos, node_size=[node_sizes.get(n, 200) for n in G.nodes()],
                            node_color=node_colors, edgecolors=node_edgecolors, linewidths=1.3)

    # Radial label offset with quadrant-based alignment: every node sits on a
    # circle centered at the origin (inner ring for TFs, outer for targets),
    # so pushing the label further out along its own angle and choosing
    # ha/va by quadrant keeps it clear of its own marker without needing a
    # collision-avoidance library (see figures/05 history for why that
    # approach was dropped — it can re-introduce the exact overlap it's
    # meant to fix on a dense ring layout).
    for n, (x, y) in pos.items():
        is_tf = n in tf_set
        r = np.hypot(x, y)
        angle = np.arctan2(y, x)
        # Base offset clears an average-sized marker; scaled up further for
        # large-centrality nodes so the label doesn't collide with an
        # enlarged marker (found via figure QC: HPSE2/VCAN/ID1 labels were
        # overlapping their own oversized nodes at high centrality/importance).
        marker_frac = np.sqrt(node_sizes.get(n, NODE_SIZE_RANGE[0]) / NODE_SIZE_RANGE[1])
        offset = (0.16 if is_tf else 0.13) + 0.16 * marker_frac
        tx = x + offset * np.cos(angle) if r > 1e-6 else x
        ty = y + offset * np.sin(angle) if r > 1e-6 else y
        ha = "left" if np.cos(angle) > 0.35 else ("right" if np.cos(angle) < -0.35 else "center")
        va = "bottom" if np.sin(angle) > 0.35 else ("top" if np.sin(angle) < -0.35 else "center")
        is_isolated = degree.get(n, 0) == 0
        ax.text(tx, ty, n, fontsize=6.5 if is_isolated else 8.5,
                fontweight="bold" if (is_tf and not is_isolated) else "normal",
                ha=ha, va=va, color=INK_MUTED if is_isolated else INK_SECONDARY,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.75))

    ax.set_title(title, fontsize=13, color=INK_SECONDARY)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    xs, ys = zip(*pos.values())
    pad = 0.55
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_ylim(min(ys) - pad, max(ys) + pad)
    # Layout is computed as a perfect circle (compute_fixed_layout); without
    # this, matplotlib stretches x/y independently to fill each subplot's
    # (non-square) box, squashing the ring into an oval. 'box' shrinks the
    # axes box to match instead of distorting the data — the fix is extra
    # background space on the shorter side, not a re-stretched circle.
    ax.set_aspect("equal", adjustable="box")




if __name__ == "__main__":
    tf_set = get_tf_set()

    # Display name per compartment — NOT just compartment.capitalize(), which
    # mangles "unk" -> "Unk" and "activatednk" -> "Activatednk".
    DISPLAY_NAME = {
        "epithelial": "Epithelial",
        "stromal": "Stromal",
        "tcell": "Tcell",
        "unk": "uNK",
        "activatednk": "ActivatedNK",
    }
    # Tcell's pseudotime failed Step 5's validation gate, so it has no dynGENIE3
    # results (results/dyngenie3/tcell_importances.csv doesn't exist) — same
    # compartment excluded in src/07_integrate_and_report.py's
    # COMPARTMENTS_NO_DYNGENIE3. Section C (the dynGENIE3 network figure) is
    # skipped for it; Sections A/B (CellOracle-only) are unaffected.
    COMPARTMENTS_NO_DYNGENIE3 = {"tcell"}

    for compartment in ["epithelial", "stromal", "tcell", "unk", "activatednk"]:
        print(f"\n{'='*60}\n{compartment}\n{'='*60}")
        has_dyngenie3 = compartment not in COMPARTMENTS_NO_DYNGENIE3
        hubs = get_curated_hubs(compartment)
        print(f"Curated gene set: {len(hubs)} genes ({sum(h in tf_set for h in hubs)} TFs)")

        links = co.load_hdf5(file_path=f"results/celloracle/{compartment}.celloracle.links")
        links.get_network_score()
        merged_score = links.merged_score

        co_edges_by_cluster = {}
        for tp in FERTILE_TIMEPOINT_ORDER:
            cluster = f"Fertile_{tp}"
            co_edges_by_cluster[cluster] = celloracle_edges_within(links, cluster, hubs)
        co_edges_by_cluster["RIF_LH+7"] = celloracle_edges_within(links, "RIF_LH+7", hubs)

        if has_dyngenie3:
            dyn_within, dyn_full = dyngenie3_edges_within(compartment, hubs, DYNGENIE3_EDGE_CAP)
            dyn_tf_importance = dyn_full.groupby("regulator")["importance"].sum()
        else:
            print("No dynGENIE3 results for this compartment (pseudotime validation gate "
                  "failed) — skipping the dynGENIE3 network figure.")
            dyn_within = pd.DataFrame(columns=["source", "target", "importance"])

        # Fixed layout: union of edges across every view in this section
        union_edges = set()
        for df in co_edges_by_cluster.values():
            union_edges.update(zip(df["source"], df["target"]))
        union_edges.update(zip(dyn_within["source"], dyn_within["target"]))
        pos = compute_fixed_layout(union_edges, hubs, tf_set)

        title_prefix = DISPLAY_NAME[compartment]

        # --- A) WOI trajectory small multiples ---
        # Wider than a naive 5x square grid — per-panel width is what actually
        # controls label crowding (font size is fixed in points, so more physical
        # width per panel = more breathing room around each node at the same
        # data-space radius). Bumped twice: once for isolated-node labels, again
        # after the compartments with many TFs (ActivatedNK 13, uNK 14) were
        # still crowded here despite being clean in the wider single/2-panel
        # figures — this figure's per-panel width is the binding constraint.
        fig, axes = plt.subplots(1, 5, figsize=(42, 8.5))
        for ax, tp in zip(axes, FERTILE_TIMEPOINT_ORDER):
            cluster = f"Fertile_{tp}"
            cent = merged_score[merged_score["cluster"] == cluster]["eigenvector_centrality"]
            node_sizes = dict(zip(cent.index, scale(cent.reindex(hubs).fillna(0), NODE_SIZE_RANGE)))
            draw_network(ax, pos, hubs, tf_set, co_edges_by_cluster[cluster], node_sizes,
                         f"{title_prefix} {tp}", tf_color=TF_NODE_FILL)
        legend_elems = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=TF_NODE_FILL,
                   markeredgecolor="white", markersize=11, label="Transcription factor"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=TARGET_NODE_FILL,
                   markeredgecolor="white", markersize=11, label="Target-only gene"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=MUTED_NODE,
                   markeredgecolor="white", markersize=11, label="No edge in this panel"),
        ]
        fig.legend(handles=legend_elems, loc="lower center", bbox_to_anchor=(0.5, -0.08), ncol=3, frameon=False, fontsize=11)
        fig.suptitle(f"{title_prefix}: WOI trajectory network (CellOracle, curated gene set)", y=1.02, fontsize=18)
        save_fig(fig, f"{OUT_DIR}/woi_trajectory_network_{compartment}.png")
        plt.close(fig)

        # --- B) Fertile vs RIF ---
        edges_fertile = co_edges_by_cluster["Fertile_LH+7"]
        edges_rif = co_edges_by_cluster["RIF_LH+7"]
        fertile_keys = set(zip(edges_fertile["source"], edges_fertile["target"]))
        rif_keys = set(zip(edges_rif["source"], edges_rif["target"]))

        fig, axes = plt.subplots(1, 2, figsize=(15, 7.5))
        cent_f = merged_score[merged_score["cluster"] == "Fertile_LH+7"]["eigenvector_centrality"]
        cent_r = merged_score[merged_score["cluster"] == "RIF_LH+7"]["eigenvector_centrality"]
        draw_network(axes[0], pos, hubs, tf_set, edges_fertile,
                     dict(zip(cent_f.index, scale(cent_f.reindex(hubs).fillna(0), NODE_SIZE_RANGE))),
                     f"{title_prefix} A - Fertile (LH+7)", tf_color=TF_NODE_FILL,
                     edge_color_fn=lambda row: RETAINED_COLOR if (row["source"], row["target"]) in rif_keys else DIVERGING_LOW)
        draw_network(axes[1], pos, hubs, tf_set, edges_rif,
                     dict(zip(cent_r.index, scale(cent_r.reindex(hubs).fillna(0), NODE_SIZE_RANGE))),
                     f"{title_prefix} B - RIF (LH+7)", tf_color=TF_NODE_FILL,
                     edge_color_fn=lambda row: RETAINED_COLOR if (row["source"], row["target"]) in fertile_keys else DIVERGING_HIGH)
        legend_elems = [
            Line2D([0], [0], color=DIVERGING_HIGH, lw=3, label="Gained in RIF"),
            Line2D([0], [0], color=DIVERGING_LOW, lw=3, label="Lost from Fertile"),
            Line2D([0], [0], color=RETAINED_COLOR, lw=3, label="Retained"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=TF_NODE_FILL,
                   markeredgecolor="white", markersize=11, label="Transcription factor"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=TARGET_NODE_FILL,
                   markeredgecolor="white", markersize=11, label="Target-only gene"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=MUTED_NODE,
                   markeredgecolor="white", markersize=11, label="No edge in this panel"),
        ]
        fig.legend(handles=legend_elems, loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=6, frameon=False, fontsize=11)
        fig.suptitle(f"{title_prefix}: Fertile vs RIF network (CellOracle, curated gene set)", y=1.03, fontsize=18)
        save_fig(fig, f"{OUT_DIR}/fertile_vs_rif_network_{compartment}.png")
        plt.close(fig)

        if not has_dyngenie3:
            continue

        # --- C) dynGENIE3 continuous-dynamics network ---
        fig, ax = plt.subplots(figsize=(9.5, 8.5))
        node_sizes = dict(zip(hubs, scale([dyn_tf_importance.get(h, 0) for h in hubs], NODE_SIZE_RANGE)))
        draw_network(ax, pos, hubs, tf_set, dyn_within, node_sizes,
                     f"{title_prefix}: dynGENIE3 continuous-dynamics network",
                     tf_color=TF_NODE_FILL, weight_col="importance")
        legend_elems = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=TF_NODE_FILL,
                   markeredgecolor="white", markersize=11, label="Transcription factor"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=TARGET_NODE_FILL,
                   markeredgecolor="white", markersize=11, label="Target-only gene"),
            Line2D([0], [0], marker="o", color="w", markerfacecolor=MUTED_NODE,
                   markeredgecolor="white", markersize=11, label="No edge in this panel"),
        ]
        ax.legend(handles=legend_elems, loc="lower center", bbox_to_anchor=(0.5, -0.15), ncol=3, frameon=False, fontsize=10.5)
        save_fig(fig, f"{OUT_DIR}/dyngenie3_network_{compartment}.png")
        plt.close(fig)

    print("\nSection 05 complete.")
