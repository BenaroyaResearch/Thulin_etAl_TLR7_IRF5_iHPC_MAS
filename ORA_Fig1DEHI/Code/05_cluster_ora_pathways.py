"""
Reduce redundancy among significant ORA terms by greedily clustering terms
whose overlapping-gene sets are highly similar (Jaccard similarity on the
"Genes" hit list, since near-duplicate GO/pathway terms tend to be enriched
by the same handful of genes). The most significant term in a cluster becomes
its "Cluster Representative". Writes ORA_Clustered.xlsx with one detail sheet
and one summary sheet per quadrant/analysis, plus an auto-drafted "For plot"
sheet (top representatives per quadrant) that mirrors the manually curated
shortlist used for the final figures -- edit that sheet by hand to match the
exact published panel before re-running 06/07.
"""
import pandas as pd
from config import (
    QUADRANTS, ANALYSES, ORA_RAW_DIR, ORA_CLUSTERED_XLSX,
    SIG_ADJP, CLUSTER_SIMILARITY_THRESHOLD, key_name, sheet_name,
)


def overlap_coefficient(set_a, set_b):
    if not set_a or not set_b:
        return 0.0
    return len(set_a & set_b) / min(len(set_a), len(set_b))


def cluster_terms(df, threshold=CLUSTER_SIMILARITY_THRESHOLD):
    """Complete-linkage clustering: a term joins a cluster only if it clears
    the similarity threshold against every existing member (not just the
    most-significant one), which avoids runaway chaining through loosely
    related terms."""
    df = df.sort_values("Adj. P-value", ascending=True).reset_index(drop=True)
    df["gene_set"] = df["Genes"].apply(lambda g: set(g.split(";")) if isinstance(g, str) else set())

    clusters = []  # list of dicts: members (row indices, most-significant first)
    cluster_of_row = [None] * len(df)
    for i, row in df.iterrows():
        joined = None
        for c_id, cluster in enumerate(clusters):
            if all(overlap_coefficient(row["gene_set"], df.loc[m, "gene_set"]) >= threshold
                   for m in cluster["members"]):
                joined = c_id
                break
        if joined is not None:
            clusters[joined]["members"].append(i)
            cluster_of_row[i] = joined
        else:
            clusters.append({"members": [i]})
            cluster_of_row[i] = len(clusters) - 1

    df["Cluster ID"] = [cluster_of_row[i] + 1 for i in range(len(df))]
    df["Cluster Representative"] = df["Cluster ID"].apply(lambda cid: df.loc[clusters[cid - 1]["members"][0], "Pathway / Term"])
    df["Cluster Size"] = df["Cluster ID"].apply(lambda cid: len(clusters[cid - 1]["members"]))
    return df, clusters


def build_detail_sheet(df):
    cols = ["Cluster ID", "Cluster Representative", "Cluster Size", "Database",
            "Pathway / Term", "Overlap", "Fold Enrichment", "P-value"]
    return df[cols].sort_values(["Cluster ID"]).reset_index(drop=True)


def build_summary_sheet(df, clusters):
    rows = []
    for c_id, cluster in enumerate(clusters, start=1):
        members = cluster["members"]
        rep_row = df.loc[members[0]]
        other_terms = [df.loc[m, "Pathway / Term"] for m in members[1:]]
        rows.append({
            "Cluster ID": c_id,
            "Cluster Representative": rep_row["Pathway / Term"],
            "Cluster Size": len(members),
            "Database": rep_row["Database"],
            "Adj. P-value (Rep)": rep_row["Adj. P-value"],
            "Fold Enrichment (Rep)": rep_row["Fold Enrichment"],
            "Odds Ratio (Rep)": rep_row["Odds Ratio"],
            "Member Pathways": "; ".join(other_terms) if other_terms else "—",
        })
    return pd.DataFrame(rows).sort_values("Adj. P-value (Rep)").reset_index(drop=True)


def draft_for_plot_column(summary_df, top_n=15):
    return summary_df.sort_values("Adj. P-value (Rep)")["Cluster Representative"].head(top_n).tolist()


def main():
    detail_sheets, summary_sheets, for_plot_cols = {}, {}, {}
    for qkey in QUADRANTS:
        key = key_name(qkey, "all")
        raw = pd.read_csv(ORA_RAW_DIR / f"{key}.csv")
        sig = raw[raw["Adj. P-value"] < SIG_ADJP].reset_index(drop=True)
        clustered, clusters = cluster_terms(sig)
        detail_sheets[sheet_name(qkey, "all")] = build_detail_sheet(clustered)
        summary = build_summary_sheet(clustered, clusters)
        summary_sheets[f"{sheet_name(qkey, 'all')} Summary"] = summary
        for_plot_cols[qkey] = draft_for_plot_column(summary)
        print(f"{key}: {len(sig)} significant terms -> {len(clusters)} clusters")

        key_irf5 = key_name(qkey, "irf5")
        raw_irf5 = pd.read_csv(ORA_RAW_DIR / f"{key_irf5}.csv")
        sig_irf5 = raw_irf5[raw_irf5["Adj. P-value"] < SIG_ADJP].reset_index(drop=True)
        clustered_irf5, clusters_irf5 = cluster_terms(sig_irf5)
        detail_sheets[sheet_name(qkey, "irf5")] = build_detail_sheet(clustered_irf5)
        summary_sheets[f"{sheet_name(qkey, 'irf5')} Summary"] = build_summary_sheet(clustered_irf5, clusters_irf5)
        print(f"{key_irf5}: {len(sig_irf5)} significant terms -> {len(clusters_irf5)} clusters")

    max_len = max(len(v) for v in for_plot_cols.values())
    for_plot_df = pd.DataFrame({
        qkey: for_plot_cols[qkey] + [None] * (max_len - len(for_plot_cols[qkey]))
        for qkey in QUADRANTS
    })

    with pd.ExcelWriter(ORA_CLUSTERED_XLSX, engine="openpyxl") as writer:
        for name, df in detail_sheets.items():
            df.to_excel(writer, sheet_name=name, index=False)
        for_plot_df.to_excel(writer, sheet_name="For plot", index=False)
        for name, df in summary_sheets.items():
            df.to_excel(writer, sheet_name=name[:31], index=False)
    print(f"Wrote {ORA_CLUSTERED_XLSX}")
    print("NOTE: 'For plot' was auto-drafted from top cluster representatives.")
    print("Edit it by hand to match the exact curated shortlist used in the published figure.")


if __name__ == "__main__":
    main()
