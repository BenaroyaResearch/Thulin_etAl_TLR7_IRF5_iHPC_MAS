#!/usr/bin/env python3
"""
This code was generated using Claude Code, Anthropic (Opus 5.5) and carefully by a bioinformatician.

Overview: This code generates the gene tracks shown in Fig. 2D and Fig. S3E. All genes were
selected from the heatmap (Fig. 2C) and are genes that were both R848-induced in an
IRF5-dependent manner in BM monocytes in vitro (Fig. 2B) and that were upregulated with
increased accessibility in iHPCs compared to Ly6Chi monocytes (Fig. 1C, quadrant 1, green dots).

bigWig files for TLR7.1 Ly6Chi Monocytes and iHPCs can be found in the superseries GSE351192 at GEO.

Compact, publication-ready genome tracks for a select set of genes, including:
    1. Ly6C-HI Mo ATAC accessibility (bigWig coverage)
    2. iHPCs ATAC accessibility (bigWig coverage)
    3. Gene model (exons / introns / direction), gene name in italic
    4. Genomic coordinate axis
A shaded band marks the proximal promoter (TSS -2000 / +200 bp).

Coverage is the mean bigWig signal in each of --nbins (default 1,000) bins across the view
window. Both tracks of a gene share one linear y-axis scale (0 to 1.1 x the larger maximum);
the scale differs between genes.

View window is set explicitly per gene in WINDOWS below (bp, on the gene's own chromosome)
rather than derived automatically.

Gene models (UCSC refGene, mm10) are read from the local cache in Code/cache/.

All text is Arial and left as editable text in the saved PDF (pdf.fonttype 42) /
SVG (svg.fonttype none).

Usage:
    # Fig. 2D
    python3 gene_tracks_pub_corrected.py C1qb Batf Tnfaip3 Pecam1
    # Fig. S3E
    python3 gene_tracks_pub_corrected.py Fcgr4 Il10 Ets2 Il12b Il18bp Nfkbia
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
import pyBigWig

# --------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))         # "Code" folder
ROOT = os.path.dirname(HERE)                               # project folder
DATA_DIR = os.path.join(ROOT, "Data")
CACHE = os.path.join(HERE, "cache")                        # cached UCSC refGene gene models
FIGDIR = os.path.join(ROOT, "Generated figures_corrected")
GENOME = "mm10"

# One entry per coverage track, drawn top to bottom.
TRACKS = [
    {
        "bw":    os.path.join(DATA_DIR, "TLR7-Ly6CMon_average_mm10.bw"),
        "label": "Ly6C$^{\\mathrm{HI}}$ Mo",
        "color": "#4f2e89",     # purple
    },
    {
        "bw":    os.path.join(DATA_DIR, "TLR7-iHPC_average_mm10.bw"),
        "label": "iHPCs",
        "color": "#ef6faa",     # pink
    },
]

# --------------------------------------------------------------------------
# View windows (bp). None = the gene's 3' end (txEnd on + strand, txStart on
# - strand), padded outward by GENE_BOUND_PAD_* so the end isn't flush with
# the frame. The right-hand comment is the window the script actually draws.
# --------------------------------------------------------------------------
WINDOWS = {
    # Fig. 2D
    "C1qb":    (None,         136_896_000),   # chr4:136,879,325-136,896,000
    "Batf":    (85_684_500,   None),          # chr12:85,684,500-85,710,316
    "Tnfaip3": (None,         19_017_000),    # chr10:19,000,104-19,017,000
    "Pecam1":  (None,         106_735_000),   # chr11:106,650,173-106,735,000
    # Fig. S3E
    "Fcgr4":   (171_015_000,  None),          # chr1:171,015,000-171,030,561
    "Il10":    (131_015_000,  None),          # chr1:131,015,000-131,025,770
    "Ets2":    (95_699_000,   None),          # chr16:95,699,000-95,722,151
    "Il12b":   (44_397_000,   None),          # chr11:44,397,000-44,415,561
    "Il18bp":  (102_013_000,  102_028_000),   # chr7:102,013,000-102,028,000
    "Nfkbia":  (55_486_000,   55_496_000),    # chr12:55,486,000-55,496,000
}
GENE_BOUND_PAD_FRAC = 0.05   # pad a gene-end bound by 5% of the window...
GENE_BOUND_PAD_MIN = 800     # ...but by at least 800 bp

# --------------------------------------------------------------------------
# Style
# --------------------------------------------------------------------------
PROM_UPSTREAM = 2000         # promoter band: TSS - 2000 bp ...
PROM_DOWNSTREAM = 200        # ... to TSS + 200 bp
PROM_COLOR = "#FFF3CD"       # pale yellow
GENE_COLOR = "#2B2B2B"

YMAX_HEADROOM = 1.10         # y-axis max = 1.1 x the larger track maximum
FRAME_LW = 0.3               # line width of the frame around the coverage tracks

# Figure geometry (inches). Row heights are relative units, scaled by INCH_PER_UNIT.
FIG_W = 3.67                 # figure width
COV_H, GENE_H, AXIS_H = 0.85, 0.34, 0.40   # coverage-track / gene-model / axis row heights
INCH_PER_UNIT = 0.42
TOP_MARGIN = 0.20            # room for the title
BOT_MARGIN = 0.32            # room for the Mb tick labels + chromosome axis title
GS_LEFT, GS_RIGHT = 0.22, 0.975   # plot area, as a fraction of figure width (left = track labels)

# Scale bar, in the axis row's own 0-1 y-range (kept inside the row so it
# can't run into the gene model above).
SCALEBAR_Y = 0.55            # height of the bar
SCALEBAR_LABEL_Y = 0.68      # bottom of its "N kb" label
SCALEBAR_RIGHT_INSET = 0.04  # gap between the bar and the right edge, as a fraction of the window
TICK_EDGE_MARGIN = 0.05      # no x tick within 5% of the window edges (its label would overhang)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,      # editable text in Illustrator / Inkscape
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    # The Ly6C^HI label uses mathtext for its superscript; point mathtext at
    # Arial so it matches the rest of the text.
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
})


# --------------------------------------------------------------------------
# Data helpers
# --------------------------------------------------------------------------
def fetch_gene_model(gene):
    """Return the cached gene model: chrom, strand, txStart/End, cdsStart/End, exons."""
    cache_file = os.path.join(CACHE, f"{gene}_{GENOME}.json")
    if not os.path.exists(cache_file):
        raise FileNotFoundError(f"No cached gene model for '{gene}' ({cache_file}).")
    with open(cache_file) as fh:
        return json.load(fh)


def coverage(bw_path, chrom, start, end, nbins):
    """Mean bigWig signal in each of nbins equal bins across [start, end); no data = 0."""
    bw = pyBigWig.open(bw_path)
    vals = bw.stats(chrom, start, end, type="mean", nBins=nbins)
    bw.close()
    return np.array([v if v is not None else 0.0 for v in vals], dtype=float)


def end_of_gene(gm):
    """The transcript's 3' end."""
    return gm["txEnd"] if gm["strand"] == "+" else gm["txStart"]


def promoter_span(gm):
    """(start, end) of the promoter band around the TSS, strand-aware."""
    if gm["strand"] == "-":
        tss = gm["txEnd"]
        return tss - PROM_DOWNSTREAM, tss + PROM_UPSTREAM
    tss = gm["txStart"]
    return tss - PROM_UPSTREAM, tss + PROM_DOWNSTREAM


def resolve_window(gene, gm):
    """Turn a WINDOWS entry into (start, end) in bp."""
    s_bp, e_bp = WINDOWS[gene]
    raw_start = s_bp if s_bp is not None else end_of_gene(gm)
    raw_end = e_bp if e_bp is not None else end_of_gene(gm)
    raw_start, raw_end = min(raw_start, raw_end), max(raw_start, raw_end)

    # A bound set to the gene's own 3' end would sit exactly on the frame edge,
    # which looks cramped and can flip the direction arrow in draw_gene (its
    # tip gets clamped past its tail). Nudge that side outward.
    pad = max(GENE_BOUND_PAD_MIN, round(GENE_BOUND_PAD_FRAC * (raw_end - raw_start)))
    start = raw_start - pad if s_bp is None else raw_start
    end = raw_end + pad if e_bp is None else raw_end
    return start, end


# --------------------------------------------------------------------------
# Drawing helpers
# --------------------------------------------------------------------------
def shade_promoter(axes_list, prom_start, prom_end, view_start, view_end):
    ps = max(prom_start, view_start)
    pe = min(prom_end, view_end)
    if pe <= ps:
        return
    for ax in axes_list:
        ax.axvspan(ps, pe, color=PROM_COLOR, alpha=1.0, zorder=0, linewidth=0)


def draw_gene(ax, gm, start, end):
    """Intron line, exon boxes (UTR thin, CDS thick), direction arrow, italic gene name."""
    y0 = 0.5
    ax.plot([gm["txStart"], gm["txEnd"]], [y0, y0],
            color=GENE_COLOR, lw=0.8, zorder=1, solid_capstyle="butt")
    utr_h, cds_h = 0.28, 0.50
    for es, ee in gm["exons"]:
        ax.add_patch(Rectangle((es, y0 - utr_h / 2), ee - es, utr_h,
                               color=GENE_COLOR, lw=0, zorder=2))
        cs, ce = max(es, gm["cdsStart"]), min(ee, gm["cdsEnd"])
        if ce > cs:
            ax.add_patch(Rectangle((cs, y0 - cds_h / 2), ce - cs, cds_h,
                                   color=GENE_COLOR, lw=0, zorder=3))

    span = end - start
    shaft = span * 0.025
    sign = 1 if gm["strand"] == "+" else -1
    anchor = end_of_gene(gm)
    arrow = dict(arrowstyle="-|>", color="black", lw=1.2,
                 mutation_scale=14, shrinkA=0, shrinkB=0)
    if start <= anchor <= end:
        # 3' end is on-screen: draw the arrowhead right at it.
        x_tail = anchor
        x_tip = anchor + sign * shaft
        x_tip = min(x_tip, end - span * 0.01) if sign > 0 else max(x_tip, start + span * 0.01)
        ax.annotate("", xy=(x_tip, y0), xytext=(x_tail, y0), arrowprops=arrow, zorder=5)
    else:
        # 3' end is off-screen: clamping it into view would invert the apparent
        # direction, so draw a small direction marker on the visible part instead.
        vis_lo, vis_hi = max(gm["txStart"], start), min(gm["txEnd"], end)
        if vis_hi > vis_lo:
            mid = (vis_lo + vis_hi) / 2
            half = min(shaft, (vis_hi - vis_lo) * 0.3)
            x_tail, x_tip = (mid - half, mid + half) if sign > 0 else (mid + half, mid - half)
            ax.annotate("", xy=(x_tip, y0), xytext=(x_tail, y0), arrowprops=arrow, zorder=5)

    ax.set_xlim(start, end)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.text(-0.02, 0.5, gm["gene"], transform=ax.transAxes,
            ha="right", va="center", fontsize=9, fontstyle="italic")


def nice_scale(span):
    """Scale-bar length: the smallest round value >= 1/5 of the window."""
    target = span / 5.0
    for s in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000]:
        if s >= target:
            return s
    return 100000


def scale_label(bp):
    return f"{bp // 1000} kb" if bp >= 1000 else f"{bp} bp"


# --------------------------------------------------------------------------
# Main plot
# --------------------------------------------------------------------------
def plot_pub(gene, ymax=None, nbins=1000, box=True, box_lw=FRAME_LW, outdir=FIGDIR):
    gm = fetch_gene_model(gene)
    chrom = gm["chrom"]

    if gene in WINDOWS:
        view_start, view_end = resolve_window(gene, gm)
    else:
        # Fallback for genes without a set window: gene body + 10% padding,
        # widened if needed so the promoter band isn't flush with the frame.
        span0 = gm["txEnd"] - gm["txStart"]
        pad = max(1500, int(0.10 * span0))
        prom_margin = 300
        ps, pe = promoter_span(gm)
        view_start = min(gm["txStart"] - pad, ps - prom_margin)
        view_end = max(gm["txEnd"] + pad, pe + prom_margin)
    span_bp = view_end - view_start

    # ---- Coverage ---------------------------------------------------------
    covs = [coverage(t["bw"], chrom, view_start, view_end, nbins) for t in TRACKS]
    x = np.linspace(view_start, view_end, nbins)
    shared_ymax = (ymax if ymax else max(c.max() for c in covs) * YMAX_HEADROOM) or 1.0

    # ---- Figure geometry --------------------------------------------------
    gene_row = len(TRACKS)
    axis_row = gene_row + 1
    height_ratios = [COV_H] * len(TRACKS) + [GENE_H, AXIS_H]
    total_units = sum(height_ratios)
    fig_h = total_units * INCH_PER_UNIT + TOP_MARGIN + BOT_MARGIN
    gs_top, gs_bottom = 1 - TOP_MARGIN / fig_h, BOT_MARGIN / fig_h

    fig = plt.figure(figsize=(FIG_W, fig_h))
    gs = fig.add_gridspec(len(height_ratios), 1, height_ratios=height_ratios,
                          hspace=0, left=GS_LEFT, right=GS_RIGHT,
                          top=gs_top, bottom=gs_bottom)

    # ---- Coverage tracks --------------------------------------------------
    all_axes = []
    track_axes = []
    for i, t in enumerate(TRACKS):
        ax_t = fig.add_subplot(gs[i])
        track_axes.append(ax_t)
        all_axes.append(ax_t)
        ax_t.fill_between(x, covs[i], 0, color=t["color"], lw=0, step="mid")
        ax_t.set_xlim(view_start, view_end)
        ax_t.set_ylim(0, shared_ymax)
        ax_t.set_yticks([0, shared_ymax])
        ylabs = ax_t.set_yticklabels(["0", f"{shared_ymax:.0f}"])
        # The tracks touch, so nudge each label into its own track to keep one
        # track's "0" from colliding with the next track's max.
        ylabs[0].set_verticalalignment("bottom")
        ylabs[1].set_verticalalignment("top")
        ax_t.tick_params(axis="y", length=2, labelsize=6, pad=1)
        ax_t.spines[["top", "right", "bottom"]].set_visible(False)
        ax_t.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
        ax_t.text(-0.02, 0.5, t["label"], transform=ax_t.transAxes, ha="right",
                  va="center", fontsize=8, color="black")

    # ---- Gene model -------------------------------------------------------
    ax_g = fig.add_subplot(gs[gene_row], sharex=track_axes[0])
    all_axes.append(ax_g)
    draw_gene(ax_g, gm, view_start, view_end)

    # ---- Coordinate axis + scale bar --------------------------------------
    ax_x = fig.add_subplot(gs[axis_row], sharex=track_axes[0])
    all_axes.append(ax_x)
    ax_x.set_xlim(view_start, view_end)
    ax_x.set_ylim(0, 1)
    ax_x.spines[["top", "left", "right"]].set_visible(False)
    ax_x.yaxis.set_visible(False)
    # Round tick spacing (1, 2 or 5 x 10^n bp); drop ticks right at the edges,
    # whose labels would run past the frame.
    ticks = MaxNLocator(nbins=4, steps=[1, 2, 5, 10]).tick_values(view_start, view_end)
    edge = span_bp * TICK_EDGE_MARGIN
    ticks = [tk for tk in ticks if view_start + edge <= tk <= view_end - edge]
    ax_x.set_xticks(ticks)
    ax_x.set_xticklabels([f"{tk / 1e6:.3f}" for tk in ticks])
    ax_x.tick_params(axis="x", length=2.5, pad=1, labelsize=6)
    ax_x.set_xlabel(f"{chrom} (Mb, {GENOME})", fontsize=7, labelpad=2)

    sb = nice_scale(span_bp)
    x1 = view_end - span_bp * SCALEBAR_RIGHT_INSET
    x0 = x1 - sb
    ax_x.plot([x0, x1], [SCALEBAR_Y, SCALEBAR_Y], color="black", lw=1.0,
              solid_capstyle="butt", clip_on=False)
    ax_x.text((x0 + x1) / 2, SCALEBAR_LABEL_Y, scale_label(sb), ha="center", va="bottom",
              fontsize=6)

    # ---- Promoter shading (continuous band across every row) --------------
    ps, pe = promoter_span(gm)
    shade_promoter(all_axes, ps, pe, view_start, view_end)

    fig.suptitle(gene, fontsize=10, fontstyle="italic",
                 y=1 - (TOP_MARGIN * 0.30) / fig_h)

    if box:
        # Frame only the coverage tracks, not the gene model / axis rows.
        tracks_frac = sum(height_ratios[:len(TRACKS)]) / total_units
        box_bottom = gs_top - tracks_frac * (gs_top - gs_bottom)
        fig.add_artist(Rectangle(
            (GS_LEFT, box_bottom), GS_RIGHT - GS_LEFT, gs_top - box_bottom,
            transform=fig.transFigure, fill=False, edgecolor="black",
            linewidth=box_lw, zorder=20, clip_on=False))

    os.makedirs(outdir, exist_ok=True)
    base = os.path.join(outdir, f"{gene}_gene_track_pub")
    fig.savefig(base + ".pdf", dpi=300)
    fig.savefig(base + ".png", dpi=300)
    plt.close(fig)
    print(f"[{gene}] window {chrom}:{view_start:,}-{view_end:,} "
          f"({span_bp:,} bp) | promoter {ps:,}-{pe:,}")
    print(f"    wrote {base}.pdf")
    print(f"    wrote {base}.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("genes", nargs="+", help="gene symbol(s), e.g. C1qb Batf")
    ap.add_argument("--ymax", type=float, default=None,
                    help="fixed y-axis max (default: per gene, from the data)")
    ap.add_argument("--nbins", type=int, default=1000, help="coverage bins across the window")
    ap.add_argument("--box", action=argparse.BooleanOptionalAction, default=True,
                    help="frame the coverage tracks (default: on; --no-box to turn off)")
    ap.add_argument("--box-lw", type=float, default=FRAME_LW, help="frame line width")
    ap.add_argument("--outdir", default=FIGDIR, help="output folder")
    args = ap.parse_args()

    failed = []
    for gene in args.genes:
        try:
            plot_pub(gene, ymax=args.ymax, nbins=args.nbins, box=args.box,
                     box_lw=args.box_lw, outdir=args.outdir)
        except Exception as ex:
            print(f"[{gene}] ERROR: {ex}", file=sys.stderr)
            failed.append(gene)
    if failed:
        sys.exit(f"Failed: {', '.join(failed)}")
