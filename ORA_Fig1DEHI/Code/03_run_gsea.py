"""
Pre-ranked GSEA (gseapy) for each quadrant/analysis gene universe, restricted to
the pathways already found significant by ORA (02_run_ora.py). The GSEA_Results
workbook carries over each term's ORA Overlap/Fold Enrichment and appends the
GSEA enrichment score (ES), normalized ES (NES), and GSEA p-values.

Note: rnaseq_log2FC in the source spreadsheet is rounded to 2 decimals, so a
large fraction of genes are tied in rank. GSEA's ES/NES therefore depend on how
ties are broken; exact values may differ slightly from a prior run that used a
different (unrecorded) tie-break. Pathway membership and ORA statistics are
unaffected and reproduce exactly.
"""
import pandas as pd
import gseapy as gp
from config import (
    QUADRANTS, ANALYSES, RANKED_DIR, ORA_RAW_DIR, GSEA_RESULTS_XLSX,
    GENE_SET_LIBRARIES, GSEA_MIN_SIZE, GSEA_MAX_SIZE, GSEA_PERMUTATIONS, GSEA_SEED,
    SIG_ADJP, key_name, sheet_name,
)


def run_gsea_for_key(key):
    rnk = pd.read_csv(RANKED_DIR / f"{key}.rnk", sep="\t", header=None, names=["gene", "score"])
    prerank = gp.prerank(
        rnk=rnk,
        gene_sets=GENE_SET_LIBRARIES,
        min_size=GSEA_MIN_SIZE,
        max_size=GSEA_MAX_SIZE,
        permutation_num=GSEA_PERMUTATIONS,
        seed=GSEA_SEED,
        threads=4,
        outdir=None,
        no_plot=True,
    )
    gsea_df = prerank.res2d.copy()
    gsea_df[["Database", "Pathway / Term"]] = gsea_df["Term"].str.split("__", n=1, expand=True)
    gsea_df = gsea_df.rename(columns={"NOM p-val": "GSEA P-value", "FDR q-val": "GSEA Adj. P-value"})
    gsea_df = gsea_df[["Database", "Pathway / Term", "ES", "NES", "GSEA P-value", "GSEA Adj. P-value"]]

    ora_sig = pd.read_csv(ORA_RAW_DIR / f"{key}.csv")
    ora_sig = ora_sig[ora_sig["Adj. P-value"] < SIG_ADJP]
    ora_sig = ora_sig.rename(columns={"P-value": "ORA P-value", "Adj. P-value": "ORA Adj. P-value"})

    merged = ora_sig.merge(gsea_df, on=["Database", "Pathway / Term"], how="inner")
    merged = merged.rename(columns={"ORA P-value": "P-value", "ORA Adj. P-value": "Adj. P-value"})
    merged = merged[[
        "Database", "Pathway / Term", "Overlap", "Fold Enrichment",
        "ES", "NES", "P-value", "Adj. P-value",
    ]]
    return merged.sort_values("Adj. P-value", ascending=True).reset_index(drop=True)


def main():
    with pd.ExcelWriter(GSEA_RESULTS_XLSX, engine="openpyxl") as writer:
        for qkey in QUADRANTS:
            for analysis in ANALYSES:
                key = key_name(qkey, analysis)
                print(f"Running prerank GSEA for {key} ...")
                merged = run_gsea_for_key(key)
                merged.to_excel(writer, sheet_name=sheet_name(qkey, analysis), index=False)
                print(f"  {len(merged)} ORA-significant terms with matching GSEA stats")
    print(f"Wrote {GSEA_RESULTS_XLSX}")


if __name__ == "__main__":
    main()
