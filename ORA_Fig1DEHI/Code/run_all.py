"""
This code was generated using Claude Code, Anthropic (Opus 5.5) and carefully by a bioinformatician.

GOAL: This code performs ORA on genes that are unregulated and have increased in iHPCs (Fig. 1 C, quadrant I, green dots) and TLR7.1 Ly6CHI monocytes (Fig. 1 C, quadrant III, blue dots), respectively. ORA analysis is also performed on IRF5-binding iHPCs (Fig. 1 G, quadrant I, green dots) and IRF5-binding TLR7.1 Ly6CHI monocytes (Fig. 1 G, quadrant III, blue dots).

After performing ORA this code subsequently graphs the top unique paths for these cell types. 

Run the full gene set enrichment pipeline end to end, in order.
Requires internet access (queries the Enrichr web API via gseapy).

Regenerates, in the parent folder:
  ORA_Results.xlsx, GSEA_Results.xlsx, ORA_Clustered.xlsx,
  ORA_IRF5_Comparison_Summary.xlsx, and all Q1_iHPC_*/Q3_TLR7Mo_* figure files.

IMPORTANT: this OVERWRITES the existing output files in this folder. If you
want to keep the current versions, copy them elsewhere first.

After step 5, ORA_Clustered.xlsx's "For plot" sheet is auto-drafted (top
cluster representatives) -- open it and edit the curated pathway shortlist by
hand to match the intended figure panel before running steps 6-7.
"""
import importlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

STEPS = [
    "01_prepare_gene_lists",
    "02_run_ora",
    "03_run_gsea",
    "04_consistency_analysis",
    "05_cluster_ora_pathways",
    "06_irf5_comparison_summary",
    "07_make_figures",
]


def main():
    for step in STEPS:
        print(f"\n=== {step} ===")
        mod = importlib.import_module(step)
        mod.main()


if __name__ == "__main__":
    main()
