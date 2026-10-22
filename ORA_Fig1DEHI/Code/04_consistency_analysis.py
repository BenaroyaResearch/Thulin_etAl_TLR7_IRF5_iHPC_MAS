"""
For each quadrant, compare the set of significant pathways found in the
All-genes ORA vs the IRF5-bound-subset ORA, and append a "Consistency" sheet
(shared/unique counts, Jaccard index) to both ORA_Results.xlsx and
GSEA_Results.xlsx.
"""
import openpyxl
import pandas as pd
from config import (
    QUADRANTS, ORA_RAW_DIR, ORA_RESULTS_XLSX, GSEA_RESULTS_XLSX,
    SIG_ADJP, key_name,
)


def significant_terms(key):
    df = pd.read_csv(ORA_RAW_DIR / f"{key}.csv")
    df = df[df["Adj. P-value"] < SIG_ADJP]
    return set(zip(df["Database"], df["Pathway / Term"]))


def consistency_rows(qkey):
    q = QUADRANTS[qkey]
    all_set = significant_terms(key_name(qkey, "all"))
    irf5_set = significant_terms(key_name(qkey, "irf5"))
    shared = all_set & irf5_set
    only_all = all_set - irf5_set
    only_irf5 = irf5_set - all_set
    union = all_set | irf5_set
    jaccard = len(shared) / len(union) if union else 0.0
    pct_shared_in_irf5 = (len(shared) / len(all_set) * 100) if all_set else 0.0

    return [
        ("Metric", "Value"),
        (f"Significant pathways – All {q['label']} genes", len(all_set)),
        (f"Significant pathways – IRF5-bound {q['label']} genes", len(irf5_set)),
        ("Shared significant pathways", len(shared)),
        ("Only in All-gene analysis", len(only_all)),
        ("Only in IRF5-bound analysis", len(only_irf5)),
        ("Jaccard similarity index", round(jaccard, 4)),
        ("% of total pathways shared in IRF5 analysis", f"{pct_shared_in_irf5:.1f}%"),
    ]


def append_sheet(xlsx_path, sheet_title, rows):
    wb = openpyxl.load_workbook(xlsx_path)
    if sheet_title in wb.sheetnames:
        del wb[sheet_title]
    ws = wb.create_sheet(sheet_title)
    for row in rows:
        ws.append(row)
    wb.save(xlsx_path)


def main():
    for qkey, q in QUADRANTS.items():
        rows = consistency_rows(qkey)
        sheet_title = f"{qkey} {q['label']} - Consistency"
        for xlsx_path in (ORA_RESULTS_XLSX, GSEA_RESULTS_XLSX):
            append_sheet(xlsx_path, sheet_title, rows)
        print(f"{sheet_title}: {rows[3][1]} shared, Jaccard={rows[6][1]}")


if __name__ == "__main__":
    main()
