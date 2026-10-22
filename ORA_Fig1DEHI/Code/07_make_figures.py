"""
Generate the Figure 7D-F / S5C-D panels: two-panel (Whole gene set vs
IRF5-bound) comparison plots of the curated pathway shortlist ("For plot"
sheet in ORA_Clustered.xlsx), one figure pair per quadrant (Q1 iHPC, Q3 TLR7
Mo). Reproduces every plot variant found in this folder from one parametrized
plotting function:

  {tag}_lollipop[.png/.pdf]                    curated order,   indep. x-axis
  {tag}_lollipop_indepX / _sharedX              curated order,   indep/shared x
  {tag}_lollipop_FEsort_indepX / _sharedX       sorted by FE,    indep/shared x  (final version)
  {tag}_FoldEnrichment[_indepX]                 sorted by FE,    shared/indep x, blue titles
  {tag}_GeneRatio[_indepX]                      sorted by ratio, shared/indep x, blue titles
  {tag}_ORA_plot                                dot plot (no stems), Reds colormap
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from config import QUADRANTS, ORA_RAW_DIR, ORA_CLUSTERED_XLSX, ROOT, SIG_ADJP, key_name

NOT_SIG_COLOR = "black"


def load_panel_data(qkey, analysis, terms):
    raw = pd.read_csv(ORA_RAW_DIR / f"{key_name(qkey, analysis)}.csv")
    rows = []
    for term in terms:
        match = raw[raw["Pathway / Term"].str.strip().str.lower() == term.strip().lower()]
        if match.empty:
            rows.append({"term": term, "sig": False})
            continue
        row = match.sort_values("Adj. P-value").iloc[0]
        sig = row["Adj. P-value"] < SIG_ADJP
        k, _K = (int(x) for x in row["Overlap"].split("/"))
        rows.append({
            "term": term, "sig": sig,
            "fold_enrichment": row["Fold Enrichment"] if sig else 0.0,
            "gene_count": k, "adj_p": row["Adj. P-value"] if sig else np.nan,
            "neglog10_p": -np.log10(row["Adj. P-value"]) if sig else 0.0,
        })
    return pd.DataFrame(rows)


def query_list_size(qkey, analysis):
    from config import GENE_LIST_DIR
    return len((GENE_LIST_DIR / f"{key_name(qkey, analysis)}.txt").read_text().split())


def two_panel_plot(qkey, out_stub, metric="fold_enrichment", sort_by_metric=False,
                    shared_x=False, style="lollipop", title_color="black",
                    cmap_name="viridis"):
    q = QUADRANTS[qkey]
    terms = [t for t in pd.read_excel(ORA_CLUSTERED_XLSX, sheet_name="For plot")[qkey].dropna().tolist()]

    all_df = load_panel_data(qkey, "all", terms)
    irf5_df = load_panel_data(qkey, "irf5", terms)
    n_all = query_list_size(qkey, "all")
    n_irf5 = query_list_size(qkey, "irf5")
    for df, n in ((all_df, n_all), (irf5_df, n_irf5)):
        df["gene_ratio"] = df["gene_count"].fillna(0) / n

    order = terms
    if sort_by_metric:
        order = all_df.sort_values(metric, ascending=False)["term"].tolist()
    all_df = all_df.set_index("term").loc[order].reset_index()
    irf5_df = irf5_df.set_index("term").loc[order].reset_index()

    xcol = "fold_enrichment" if metric == "fold_enrichment" else "gene_ratio"
    xlabel = "Fold Enrichment" if metric == "fold_enrichment" else "Gene Ratio"
    vmax = max(all_df["neglog10_p"].max(), irf5_df["neglog10_p"].max(), 1e-9)
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=0, vmax=vmax)

    fig, axes = plt.subplots(1, 2, figsize=(20, max(4, 0.5 * len(order) + 1)))
    fig.subplots_adjust(wspace=0.9)
    panel_titles = [f"{qkey} {q['label']} — All Genes", f"{qkey} {q['label']} — IRF5-Bound"]
    for ax, df, title in zip(axes, (all_df, irf5_df), panel_titles):
        y = np.arange(len(df))[::-1]
        sizes = 20 + 8 * df["gene_count"].fillna(0)
        colors = [cmap(norm(v)) if s else NOT_SIG_COLOR for v, s in zip(df["neglog10_p"], df["sig"])]
        if style == "lollipop":
            for yi, xi, sig in zip(y, df[xcol], df["sig"]):
                if sig:
                    ax.plot([0, xi], [yi, yi], color="gray", lw=1.2, zorder=1)
        ax.scatter(df[xcol].where(df["sig"], 0), y, s=sizes, c=colors, zorder=2, edgecolor="none")
        ax.set_yticks(y)
        ax.set_yticklabels(df["term"])
        ax.set_xlabel(xlabel)
        ax.set_title(title, color=title_color, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    if shared_x:
        xmax = max(ax.get_xlim()[1] for ax in axes)
        for ax in axes:
            ax.set_xlim(0, xmax)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=axes, shrink=0.5, pad=0.02)
    cbar.set_label("-log10(adj. p-value)")

    fig.savefig(ROOT / f"{out_stub}.png", dpi=300, bbox_inches="tight")
    fig.savefig(ROOT / f"{out_stub}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_stub}.png/.pdf")


def main():
    for qkey, q in QUADRANTS.items():
        tag = q["file_tag"]
        two_panel_plot(qkey, f"{tag}_lollipop", metric="fold_enrichment", sort_by_metric=False,
                        shared_x=False, style="lollipop", title_color="black")
        two_panel_plot(qkey, f"{tag}_lollipop_indepX", metric="fold_enrichment", sort_by_metric=False,
                        shared_x=False, style="lollipop", title_color="black")
        two_panel_plot(qkey, f"{tag}_lollipop_sharedX", metric="fold_enrichment", sort_by_metric=False,
                        shared_x=True, style="lollipop", title_color="black")
        two_panel_plot(qkey, f"{tag}_lollipop_FEsort_indepX", metric="fold_enrichment", sort_by_metric=True,
                        shared_x=False, style="lollipop", title_color="black")
        two_panel_plot(qkey, f"{tag}_lollipop_FEsort_sharedX", metric="fold_enrichment", sort_by_metric=True,
                        shared_x=True, style="lollipop", title_color="black")
        two_panel_plot(qkey, f"{tag}_FoldEnrichment", metric="fold_enrichment", sort_by_metric=True,
                        shared_x=True, style="lollipop", title_color="tab:blue")
        two_panel_plot(qkey, f"{tag}_FoldEnrichment_indepX", metric="fold_enrichment", sort_by_metric=True,
                        shared_x=False, style="lollipop", title_color="tab:blue")
        two_panel_plot(qkey, f"{tag}_GeneRatio", metric="gene_ratio", sort_by_metric=True,
                        shared_x=True, style="lollipop", title_color="tab:blue")
        two_panel_plot(qkey, f"{tag}_GeneRatio_indepX", metric="gene_ratio", sort_by_metric=True,
                        shared_x=False, style="lollipop", title_color="tab:blue")
        two_panel_plot(qkey, f"{tag}_ORA_plot", metric="gene_ratio", sort_by_metric=False,
                        shared_x=False, style="dot", title_color="tab:blue", cmap_name="Reds")


if __name__ == "__main__":
    main()
