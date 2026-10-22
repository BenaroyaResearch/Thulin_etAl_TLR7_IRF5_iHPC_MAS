# Gene tracks: Fig. 2D and Fig. S3E

`gene_tracks_pub_corrected.py` draws the ATAC-seq gene tracks: Ly6C<sup>HI</sup> monocyte and
iHPC coverage, the gene model, and the genomic axis, with the promoter (TSS −2000 / +200 bp)
shaded. The script is self-contained and needs no other project code.

## Run

From this `Code/` folder:

```
# Fig. 2D
python3 gene_tracks_pub_corrected.py C1qb Batf Tnfaip3 Pecam1
# Fig. S3E
python3 gene_tracks_pub_corrected.py Fcgr4 Il10 Ets2 Il12b Il18bp Nfkbia
```

It writes `<Gene>_gene_track_pub.pdf` and `.png` to `../Generated figures_corrected/`. Use
`--outdir` to write somewhere else. If any gene fails, the script exits with an error.

## Inputs

| File | Use |
|---|---|
| `../Data/TLR7-Ly6CMon_average_mm10.bw` | Ly6C<sup>HI</sup> monocyte ATAC coverage |
| `../Data/TLR7-iHPC_average_mm10.bw` | iHPC ATAC coverage |
| `cache/<Gene>_mm10.json` | Gene models (UCSC refGene, mm10, retrieved July–Sept 2026) |

The view window for each gene is set in `WINDOWS` at the top of the script.

## Figure settings

- Linear coverage: the mean bigWig signal in each of 1,000 bins across the window (`--nbins`).
- Both tracks of a gene share one y-axis, from 0 to 1.1 × the larger maximum. The maximum
  differs between genes.
- A 0.3-pt frame surrounds the coverage tracks (`--box-lw`; `--no-box` turns it off).
- Text is Arial and stays editable in the PDF.

## Tested with

Python 3.9.6, matplotlib 3.9.4, numpy 2.0.2, pyBigWig 0.3.18 (macOS).

```
pip install matplotlib numpy pyBigWig
```
