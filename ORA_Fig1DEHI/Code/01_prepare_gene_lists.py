"""
Build the four gene lists / ranked lists (Q1 all, Q1 IRF5-bound, Q3 all, Q3 IRF5-bound)
used downstream for ORA and GSEA, from the RNA-ATAC quadrant calls in
"Q1 and Q3-IRF.xlsx".
"""
import pandas as pd
from config import (
    SOURCE_XLSX, ALL_GENES_SHEET, IRF5_BOUND_SHEET,
    QUADRANTS, ANALYSES, GENE_LIST_DIR, RANKED_DIR, key_name,
)


def load_quadrant_genes(sheet, quadrant_code):
    df = pd.read_excel(SOURCE_XLSX, sheet_name=sheet)
    df = df[df["quadrant"] == quadrant_code].copy()
    df["gene_upper"] = df["gene"].str.upper()
    # A gene can appear more than once (multiple ATAC peaks); keep the row
    # with the largest-magnitude RNA-seq log2FC as that gene's representative value.
    df["abs_log2fc"] = df["rnaseq_log2FC"].abs()
    df = df.sort_values("abs_log2fc", ascending=False).drop_duplicates("gene_upper")
    return df[["gene_upper", "rnaseq_log2FC"]].rename(columns={"rnaseq_log2FC": "log2FC"})


def write_gene_list(df, key):
    path = GENE_LIST_DIR / f"{key}.txt"
    path.write_text("\n".join(df["gene_upper"]) + "\n")
    return path


def write_ranked_list(df, key):
    ranked = df.sort_values("log2FC", ascending=False)
    path = RANKED_DIR / f"{key}.rnk"
    ranked.to_csv(path, sep="\t", header=False, index=False)
    return path


def main():
    for qkey, qinfo in QUADRANTS.items():
        for analysis in ANALYSES:
            sheet = ALL_GENES_SHEET if analysis == "all" else IRF5_BOUND_SHEET
            df = load_quadrant_genes(sheet, qinfo["code"])
            key = key_name(qkey, analysis)
            gene_path = write_gene_list(df, key)
            rnk_path = write_ranked_list(df, key)
            print(f"{key}: {len(df)} unique genes -> {gene_path.name}, {rnk_path.name}")


if __name__ == "__main__":
    main()
