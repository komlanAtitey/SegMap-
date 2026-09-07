"""Command-line argument parsing."""
from .common import *  # noqa: F401,F403


def parse_args() -> argparse.Namespace:  # define parse_args
    parser = argparse.ArgumentParser(
    # create parser
    description="Xenium + Cellpose (FULL-RES tiled fidelity) upgraded with speedups and robust merging.")
    parser.add_argument(
    "--xenium_dir",
    type=str,
    required=True,
     help="Xenium output folder (e.g. outs_0058080).")  # xenium dir arg
    parser.add_argument(
    "--morphology",
    type=str,
    default=None,
     # morphology arg
     help="Path to morphology.ome.tif (default: <xenium_dir>/morphology.ome.tif).")
    parser.add_argument(
    "--transcripts",
    type=str,
    default=None,
     # transcripts arg
     help="Path to transcripts.parquet (default: <xenium_dir>/transcripts.parquet).")
    parser.add_argument(
    "--outdir",
    type=str,
    default=None,
     # outdir arg
     help="Output directory (default: <xenium_dir>/cellpose_pipeline_outputs).")
    parser.add_argument(
    "--segmentation_mode",
    type=str,
    default="tile",
    choices=["tile"],
     help="Segmentation mode (only 'tile' supported).")  # seg mode arg
    parser.add_argument(
    "--tile_size",
    type=int,
    default=2048,
     help="Tile size in full-res pixels (default: 2048).")  # tile size arg
    parser.add_argument(
    "--tile_overlap",
    type=int,
    default=128,
     help="Tile overlap in full-res pixels (default: 128).")  # overlap arg
    parser.add_argument(
    "--tile_batch",
    type=int,
    default=8,
     # tile batch arg
     help="Tile batch size (CPU only). GPU ignores this. (default: 8).")
    parser.add_argument(
    "--skip_tile_std",
    type=float,
    default=0.01,
     # skip threshold arg
     help="Skip tiles with std below this after normalization (default: 0.01).")
    parser.add_argument(
    "--cpu_workers",
    type=int,
    default=1,
     help="CPU workers for tiling (GPU uses 1).")  # cpu workers arg
    parser.add_argument(
    "--gpu",
    action="store_true",
     help="Use GPU for Cellpose if available.")  # gpu flag arg
    parser.add_argument(
    "--force_cpu_on_gpu_error",
    action="store_true",
     help="If set, a GPU error on a tile falls back to CPU for that tile.")  # fallback flag arg
    parser.add_argument(
    "--n_clusters",
    type=int,
    default=12,
     help="Requested number of clusters (default: 12).")  # clusters arg
    parser.add_argument(
    "--max_cells_umap_fit",
    type=int,
    default=15000,
     help="Max cells to fit UMAP (default: 15000).")  # umap cap arg
    parser.add_argument(
    "--do_scanpy_markers",
    action="store_true",
     help="If set, compute Scanpy markers (slow).")  # scanpy markers arg
    parser.add_argument(
    "--seed",
    type=int,
    default=0,
     help="Random seed.")  # seed arg
    parser.add_argument(
        "--resolution",
        type=float,
        default=0.3,
        help="Microns per pixel"
    )

    parser.add_argument(
        "--tune_resolution",
        dest="tune_resolution",
        action="store_true",
        default=True,
        help="Sweep Leiden resolutions around --resolution and keep the best by "
             "subsampled silhouette (improves cluster separation). Enabled by default."
    )
    parser.add_argument(
        "--no_tune_resolution",
        dest="tune_resolution",
        action="store_false",
        help="Disable resolution tuning and use --resolution as-is."
    )

    parser.add_argument(
        "--spatial_weight",
        type=float,
        default=0.0,
        help="Weight of spatial-coordinate + local-density terms in the final "
             "(X_bayes) clustering space. 0.0 = transcriptomically homogeneous "
             "cell-type clusters (recommended for a clean cell-type UMAP). "
             "Increase toward ~0.3-1.0 to recover spatial-domain behavior."
    )

    parser.add_argument(
        "--no_preview",
        action="store_true",
        help="Skip the 6 Cellpose preview patches / mosaic (visualization only; "
             "does not affect segmentation or any results). Small speedup."
    )
    parser.add_argument(
        "--skip_seg_exports",
        action="store_true",
        help="Skip the full-resolution segmentation PNG exports (~5%% of runtime, "
             "visualization only; does not affect any results)."
    )

    parser.add_argument(
        "--reuse_mask",
        type=str,
        default="",
        help="Skip Cellpose and reuse a saved nucleus mask for fast downstream "
             "iteration (~1 min vs ~45 min). Pass 'auto' to load "
             "<outdir>/data/full_mask_nuclei.tif from a previous run, or an "
             "explicit path. Expansion/QC/clustering still run, so you can tune "
             "those without re-segmenting. Re-run without this flag if you change "
             "segmentation params (diameter, flow, cellprob, min_size, tile_size)."
    )

    parser.add_argument(
        "--expand_pixels",
        type=int,
        default=0,
        help="Expand each nucleus mask outward by this many pixels before "
             "assigning transcripts (0=off; nucleus-only). Captures cytoplasmic "
             "transcripts for richer, more separable cell profiles. At Xenium "
             "~0.2125 um/px, ~10-25 px ~= 2-5 um. This is the main lever for real "
             "cluster separability with nucleus-based segmentation."
    )

    # --- Overlap resolution / non-overlapping tessellation (post-Cellpose) ---
    parser.add_argument(
        "--resolve_overlaps",
        type=lambda v: str(v).lower() not in ("0", "false", "no", "off"),
        default=True,
        help="Resolve overlapping expanded territories into a NON-OVERLAPPING, "
             "nucleus-seeded cell tessellation before transcript assignment "
             "(default: on). Contested pixels are awarded to the nearest nucleus; "
             "geometry is repaired and implausible objects dropped. Set to "
             "false/0/off to fall back to plain label dilation (expand_labels)."
    )
    parser.add_argument(
        "--resolve_mode",
        type=str,
        default="expand_labels",
        choices=["expand_labels", "watershed"],
        help="How overlapping expansions are resolved. 'expand_labels' (default) "
             "reproduces classic label expansion EXACTLY and only removes true "
             "overlaps by awarding contested pixels to the nearer nucleus -- "
             "preserves tissue morphology and leaves acellular space empty. "
             "'watershed' fills contact seams into a gap-free tessellation and is "
             "required for edge-aware refinement (--edge_weight > 0). Use "
             "'watershed' only if you specifically want gap-free territories."
    )
    parser.add_argument(
        "--edge_weight",
        type=float,
        default=0.0,
        help="Edge-aware refinement strength for overlap resolution. 0.0 = "
             "distance-only tessellation (default); ~0.3-0.7 blends the DAPI/"
             "morphology intensity gradient into the watershed so cell boundaries "
             "snap to real image edges (membranes/nuclear rims) rather than to "
             "purely geometric medial lines. Requires --expand_pixels > 0."
    )
    parser.add_argument(
        "--min_cell_area",
        type=int,
        default=15,
        help="Minimum plausible cell area (px) after overlap resolution; smaller "
             "objects are dropped as fragments/artefacts."
    )
    parser.add_argument(
        "--max_cell_area",
        type=int,
        default=0,
        help="Maximum plausible cell area (px) after overlap resolution "
             "(0 = no upper bound). Removes implausibly large merged territories."
    )
    parser.add_argument(
        "--min_solidity",
        type=float,
        default=0.0,
        help="Minimum cell solidity (area / convex-hull area) after overlap "
             "resolution (0 disables the check; ~0.8-0.9 removes fragmented or "
             "non-convex implausible shapes)."
    )

    parser.add_argument(
        "--umap_init",
        type=str,
        default="spectral",
        choices=["spectral", "paga"],
        help="UMAP initialization. 'spectral' (default) is the standard init; "
             "'paga' initializes from a cluster-connectivity graph so distinct "
             "clusters separate more while related clusters stay adjacent "
             "(preserves biological relationships + internal structure)."
    )

    parser.add_argument(
        "--umap_neighbors",
        type=int,
        default=30,
        help="Neighbors for the shared kNN graph used by UMAP (and Leiden). "
             "Matches Seurat RunUMAP n.neighbors."
    )
    parser.add_argument(
        "--umap_metric",
        type=str,
        default="cosine",
        help="Distance metric for the shared kNN graph (UMAP + Leiden). "
             "Matches Seurat RunUMAP metric. e.g. cosine, euclidean."
    )
    parser.add_argument(
        "--umap_min_dist",
        type=float,
        default=0.2,
        help="UMAP min_dist (layout only). Lower packs clusters tighter and "
             "opens space between them."
    )
    parser.add_argument(
        "--umap_spread",
        type=float,
        default=1.5,
        help="UMAP spread (layout only). Higher scales the embedding up."
    )

    parser.add_argument(
        "--merge_threshold",
        type=float,
        default=0.0,
        help="Consolidate over-split clusters whose separation (between-centroid "
             "distance / summed within-cluster radius) is below this value "
             "(0=off; try ~1.0). Fewer, larger clusters -> stronger DE significance."
    )
    parser.add_argument(
        "--merge_min_clusters",
        type=int,
        default=2,
        help="Lower bound on cluster count for --merge_threshold consolidation."
    )

    parser.add_argument(
        "--normalization",
        type=str,
        default="lognorm",
        choices=["lognorm", "pearson"],
        help="Normalization before PCA. 'lognorm' = CP10k+log1p+scale (default); "
             "'pearson' = analytic Pearson residuals (better separability / marker "
             "contrast for sparse count panels)."
    )
    parser.add_argument(
        "--auto_n_pcs",
        action="store_true",
        help="Pick the number of PCs from the variance-ratio elbow instead of a "
             "fixed --n_pcs (drops noise PCs on targeted panels)."
    )
    parser.add_argument(
        "--subtract_background",
        action="store_true",
        help="Subtract ambient background estimated from Xenium negative-control "
             "probes/codewords before clustering (denoises profiles)."
    )
    parser.add_argument(
        "--remove_doublets",
        action="store_true",
        help="Drop predicted doublets via Scrublet before clustering (removes "
             "between-cluster cells that blunt marker metrics)."
    )

    parser.add_argument(
        "--n_pcs",
        type=int,
        default=30,
        help="Number of principal components for the neighbor graph / UMAP / "
             "clustering. For targeted Xenium panels (a few hundred genes), 50 "
             "includes mostly-noise PCs that fragment the UMAP into spurious "
             "sub-groups; ~20-30 is usually cleaner."
    )

    parser.add_argument(
        "--min_genes_per_cell",
        type=int,
        default=30,
        help="QC: drop cells with fewer than this many distinct genes detected. "
             "30 is strict for nucleus-only segmentation (which captures few "
             "transcripts/cell) and can discard most cells; lower to ~5-10 to keep them."
    )
    parser.add_argument(
        "--min_cells_per_gene",
        type=int,
        default=10,
        help="QC: drop genes detected in fewer than this many cells."
    )
    parser.add_argument(
        "--min_counts_per_cell",
        type=int,
        default=0,
        help="QC: drop cells with fewer than this many total transcripts (0=off). "
             "A value like 20-30 removes very sparse nuclei that smear the UMAP, "
             "and is gentler than raising --min_genes_per_cell alone."
    )

    parser.add_argument(
        "--cellpose_model",
        type=str,
        default="nuclei",
        help="Cellpose built-in model. 'nuclei' is correct for DAPI nuclei "
             "segmentation (improves recall vs the generic cyto model). "
             "Other options: 'cyto3', 'cyto2', 'cyto'."
    )
    parser.add_argument(
        "--min_size",
        type=int,
        default=15,
        help="Minimum object size in pixels passed to Cellpose. Lower it "
             "(e.g. 8) to keep small/dim nuclei that would otherwise be dropped."
    )

    parser.add_argument(
        "--top_gene_cap",
        type=int,
        default=5050,
        help="Cap on the number of most-abundant genes kept when building the "
             "cell-by-gene matrix (top genes by total counts)."
    )

    parser.add_argument(
        "--spatial_denoise_iters",
        type=int,
        default=2,
        help="Spatial denoising passes applied to cell-type labels (0 disables). "
             "Uses spatial neighbors to fix low-confidence, isolated misassignments "
             "without overriding transcriptomic identity."
    )
    parser.add_argument(
        "--spatial_denoise_agree",
        type=float,
        default=0.65,
        help="Fraction of a cell's spatial neighbors that must agree on one label "
             "before a low-confidence cell is reassigned during spatial denoising."
    )

    parser.add_argument(
        "--dapi_zproject",
        type=str,
        default="none",
        choices=["none", "max", "mean"],
        help="Collapse the DAPI z-stack across focal planes before segmentation. "
             "'none' picks the single best-focus plane (default); 'max' takes the "
             "brightest value per pixel (recommended for nuclei — recovers cells "
             "in focus at different depths); 'mean' averages planes."
    )

    parser.add_argument(
        "--cellpose_diameter",
        type=float,
        default=18,
        help="Expected object diameter in pixels"
    )

    parser.add_argument(
        "--flow_threshold",
        type=float,
        default=0.4,
        help="Cellpose flow threshold"
    )

    parser.add_argument(
        "--cellprob_threshold",
        type=float,
        default=0.2,
        help="Cellpose cell probability threshold"
    )
    return parser.parse_args()  # return parsed args
