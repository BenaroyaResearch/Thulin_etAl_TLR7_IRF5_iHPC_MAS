"""
Over-representation analysis (ORA) of each gene list against Enrichr gene-set
libraries, via gseapy. Writes ORA_Results.xlsx (one sheet per quadrant/analysis,
significant terms only) and caches the full raw results (incl. member genes,
needed by the clustering and comparison-summary steps) under intermediate/ora_raw.
"""
import pandas as pd
import gseapy as gp
from config import (
    QUADRANTS, ANALYSES, GENE_LIST_DIR, ORA_RAW_DIR, ORA_RESULTS_XLSX,
    GENE_SET_LIBRARIES, ENRICHR_ORGANISM, SIG_ADJP, BACKGROUND_SIZE,
    key_name, sheet_name,
)


def fold_enrichment(overlap_str, n_query, background=BACKGROUND_SIZE):
    k, K = (int(x) for x in overlap_str.split("/"))
    return round((k / n_query) / (K / background), 3)


def run_ora_for_gene_list(key):
    genes = (GENE_LIST_DIR / f"{key}.txt").read_text().split()
    n_query = len(genes)
    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=GENE_SET_LIBRARIES,
        organism=ENRICHR_ORGANISM,
        outdir=None,
    )
    df = enr.results.copy()
    df["Fold Enrichment"] = df["Overlap"].apply(lambda o: fold_enrichment(o, n_query))
    df = df.rename(columns={
        "Gene_set": "Database",
        "Term": "Pathway / Term",
        "P-value": "P-value",
        "Adjusted P-value": "Adj. P-value",
        "Odds Ratio": "Odds Ratio",
        "Combined Score": "Combined Score",
        "Genes": "Genes",
    })
    df = df.sort_values("Adj. P-value", ascending=True).reset_index(drop=True)

    raw_path = ORA_RAW_DIR / f"{key}.csv"
    df.to_csv(raw_path, index=False)

    display_cols = [
        "Database", "Pathway / Term", "Overlap", "Fold Enrichment",
        "P-value", "Adj. P-value", "Odds Ratio", "Combined Score",
    ]
    sig = df[df["Adj. P-value"] < SIG_ADJP][display_cols].reset_index(drop=True)
    return sig


def main():
    with pd.ExcelWriter(ORA_RESULTS_XLSX, engine="openpyxl") as writer:
        for qkey in QUADRANTS:
            for analysis in ANALYSES:
                key = key_name(qkey, analysis)
                print(f"Running Enrichr ORA for {key} ...")
                sig = run_ora_for_gene_list(key)
                sig.to_excel(writer, sheet_name=sheet_name(qkey, analysis), index=False)
                print(f"  {len(sig)} significant terms (Adj. P < {SIG_ADJP})")
    print(f"Wrote {ORA_RESULTS_XLSX}")


if __name__ == "__main__":
    main()
