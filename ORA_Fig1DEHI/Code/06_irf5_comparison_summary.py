"""
For the curated pathway shortlist (the "For plot" sheet in ORA_Clustered.xlsx),
build a side-by-side comparison of ORA stats for the whole up-regulated gene
set vs. its IRF5-bound subset. Writes ORA_IRF5_Comparison_Summary.xlsx.
"""
import openpyxl
import pandas as pd
from config import QUADRANTS, ORA_RAW_DIR, ORA_CLUSTERED_XLSX, IRF5_COMPARISON_XLSX, SIG_ADJP, key_name

NOT_ENRICHED = "Not significant / not enriched"

QUADRANT_TITLE = {
    "Q1": ("Q1 - iHPC", "Quadrant 1: Pathways Up in iHPC", "Q1", "Up in iHPC"),
    "Q3": ("Q3 - TLR7 Monocytes", "Quadrant 3: Pathways Up in TLR7 Monocytes", "Q3", "Up in TLR7 Monocytes"),
}


def load_curated_list(qkey):
    for_plot = pd.read_excel(ORA_CLUSTERED_XLSX, sheet_name="For plot")
    return [t for t in for_plot[qkey].dropna().tolist()]


def lookup(term, raw_df):
    match = raw_df[raw_df["Pathway / Term"].str.strip().str.lower() == term.strip().lower()]
    if match.empty:
        return None
    row = match.sort_values("Adj. P-value").iloc[0]
    if row["Adj. P-value"] >= SIG_ADJP:
        return None
    return row


def build_rows(qkey):
    all_raw = pd.read_csv(ORA_RAW_DIR / f"{key_name(qkey, 'all')}.csv")
    irf5_raw = pd.read_csv(ORA_RAW_DIR / f"{key_name(qkey, 'irf5')}.csv")
    terms = load_curated_list(qkey)

    rows = []
    for term in terms:
        all_row = lookup(term, all_raw)
        irf5_row = lookup(term, irf5_raw)
        if all_row is None:
            row_out = [term, NOT_ENRICHED, None, None, None]
        else:
            k = int(all_row["Overlap"].split("/")[0])
            row_out = [term, all_row["Fold Enrichment"], k, all_row["Adj. P-value"],
                       "; ".join(sorted(all_row["Genes"].split(";")))]
        if irf5_row is None:
            row_out += [NOT_ENRICHED, None, None, None]
        else:
            k = int(irf5_row["Overlap"].split("/")[0])
            row_out += [irf5_row["Fold Enrichment"], k, irf5_row["Adj. P-value"],
                        "; ".join(sorted(irf5_row["Genes"].split(";")))]
        rows.append(row_out)
    return rows


def write_sheet(wb, qkey):
    sheet_title, banner, _, direction = QUADRANT_TITLE[qkey]
    ws = wb.create_sheet(sheet_title)
    ws.append([banner])
    ws.append([f'Over-representation analysis (ORA) — whole up-regulated gene set vs. the IRF5-bound subset '
               f'of that gene set. Pathways shown are the curated {qkey} list (column '
               f'{"A" if qkey == "Q1" else "B"}, "For plot" sheet).'])
    ws.append([])
    ws.append(["Pathway", f"Whole Gene Set ({direction})", None, None, None,
               f"IRF5-Bound Genes ({direction})", None, None, None])
    ws.append([None, "Fold Enrichment", "Gene Count", "Adj. P-value", "Genes in Overlap",
               "Fold Enrichment", "Gene Count", "Adj. P-value", "Genes in Overlap"])
    for row in build_rows(qkey):
        ws.append(row)


def main():
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    for qkey in QUADRANTS:
        write_sheet(wb, qkey)
    wb.save(IRF5_COMPARISON_XLSX)
    print(f"Wrote {IRF5_COMPARISON_XLSX}")


if __name__ == "__main__":
    main()
