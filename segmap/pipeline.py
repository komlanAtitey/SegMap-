"""End-to-end pipeline orchestration (main)."""
from .common import *  # noqa: F401,F403
from .assignment import assign_transcripts_to_cells_fullres, parquet_schema_columns, save_cell_by_gene_csv
from .cli import parse_args
from .clustering import cluster_separability_report, merge_indistinct_clusters, run_leiden_best, silhouette_samples_fast, silhouette_score_fast, spatial_denoise_labels
from .export_seurat import export_seurat_bundle
from .io_ome import load_one_plane_ome, pick_top_patches_by_signal, save_mosaic, save_patch_polygon_image
from .markers import auto_annotate_clusters, compute_scanpy_markers_all_cells, fast_cluster_markers_foldchange
from .bayes_markers import run_bayes_markers, couple_pj_to_qi  # Bayesian marker model (Fig. 1D) + p_j->q_i coupling (Fig. 1E)
from .preprocessing import exact_kmeans_clustering, safe_filter_adata, subtract_negcontrol_background, remove_doublets
from .refinement import save_segmap_stage, segmap_refinement_publication
from .segmentation import CellposeRunner, save_segmentation_exports_fullres, tile_cellpose_fullres_upgraded
from .utils import downsample_mean, normalize01, pick_first_present
from .visualization import save_cluster_tables_csv, save_spatial_plot, save_umap_plot


def _load_morphology_plane(MORPHOLOGY_IMAGE):
    """Load the DAPI/morphology plane as a 2-D float array for use as the
    edge-aware cue in overlap resolution. Returns None on failure."""
    try:
        img, _plane, _axes = load_one_plane_ome(MORPHOLOGY_IMAGE, auto_pick=True,
                                                score_ds=16, zproject="none")
        return np.asarray(img, dtype=np.float32)
    except Exception:
        return None


def _segment_full(args, MORPHOLOGY_IMAGE, FIG_DIR, GPU, MODEL_TYPE, MIN_SIZE,
                  TILE_SIZE, TILE_OVERLAP, CPU_WORKERS, SKIP_TILE_STD,
                  FORCE_CPU_FALLBACK, TILE_BATCH, BOUNDARY_WIDTH, PASTEL_ALPHA,
                  PATCH_SIZE_FULL, N_PATCHES, PATCH_STRIDE_FULL, PATCH_PREVIEW_DS,
                  PATCH_DS, _timer):
    """Load DAPI, render preview patches, and run full-res tiled Cellpose.
    Returns the pre-expansion nucleus label mask. Skipped entirely in reuse mode."""
    dapi, chosen_plane, axes = load_one_plane_ome(MORPHOLOGY_IMAGE, auto_pick=True, score_ds=16, zproject=str(getattr(args, "dapi_zproject", "none")))
    dapi = dapi.astype(np.float32, copy=False)
    print(f"Loaded DAPI plane index={chosen_plane} axes={axes} shape={dapi.shape}", flush=True)

    coords = pick_top_patches_by_signal(dapi, PATCH_SIZE_FULL, N_PATCHES,
                                        PATCH_STRIDE_FULL, preview_ds=PATCH_PREVIEW_DS)
    print("Selected patch coords (full-res x,y):", coords, flush=True)
    patch_pngs = []
    for i, (x0, y0) in enumerate(coords, start=1):
        if getattr(args, "no_preview", False):
            break  # preview patches are visualization-only; skip for speed
        try:
            patch_full = dapi[y0:y0 + PATCH_SIZE_FULL, x0:x0 + PATCH_SIZE_FULL]
            patch_ds = normalize01(downsample_mean(patch_full, PATCH_DS))
            runner = CellposeRunner(
                gpu=GPU, model_type=MODEL_TYPE, diameter=args.cellpose_diameter,
                flow_threshold=args.flow_threshold,
                cellprob_threshold=args.cellprob_threshold, min_size=MIN_SIZE,
            )
            masks = runner.eval(patch_ds)
            out_png = os.path.join(FIG_DIR, f"cellpose_polygons_patch{i:02d}.png")
            save_patch_polygon_image(masks, out_png, boundary_width=BOUNDARY_WIDTH,
                                     pastel_alpha=PASTEL_ALPHA)
            patch_pngs.append(out_png)
        except Exception as e:
            print(f"Preview patch segmentation failed for patch {i}: {e}", flush=True)

    if len(patch_pngs) > 0:
        mosaic_png = os.path.join(FIG_DIR, f"cellpose_polygons_mosaic_{len(patch_pngs)}patches.png")
        save_mosaic(patch_pngs, mosaic_png, cols=3)

    print("[SEG] Running FULL-RES tiled Cellpose ...", flush=True)
    t_seg0 = time.perf_counter()
    _timer.lap("Load DAPI + preview patches")
    mask_full = tile_cellpose_fullres_upgraded(
        dapi, tile_size=TILE_SIZE, overlap=TILE_OVERLAP, gpu=GPU,
        flow_threshold=args.flow_threshold,
        cellprob_threshold=args.cellprob_threshold,
        cellpose_diameter=args.cellpose_diameter, cpu_workers=CPU_WORKERS,
        skip_threshold=SKIP_TILE_STD, force_cpu_on_gpu_error=FORCE_CPU_FALLBACK,
        tile_batch=TILE_BATCH, model_type=MODEL_TYPE, min_size=MIN_SIZE,
    )
    print(f"[SEG] Full-res tiled segmentation done in {time.perf_counter() - t_seg0:.2f}s; "
          f"cells={int(mask_full.max())}", flush=True)
    _timer.lap("Cellpose segmentation (tiled)")
    return mask_full, axes, chosen_plane


def main() -> None:  # main entry
    args = parse_args()  # parse CLI
    FLOW_THRESHOLD = args.flow_threshold
    CELLPROB_THRESHOLD = args.cellprob_threshold
    RESOLUTION = args.resolution
    CELLPOSE_DIAMETER = args.cellpose_diameter
    SPATIAL_WEIGHT_BAYES = float(args.spatial_weight)  # spatial/density weight in X_bayes
    TOP_GENE_CAP = int(args.top_gene_cap)  # cap on genes kept in the cell-by-gene matrix
    MODEL_TYPE = str(args.cellpose_model)   # cellpose model (nuclei for DAPI)
    MIN_SIZE = int(args.min_size)           # min object size (px); lower keeps smaller nuclei
    N_PCS = int(args.n_pcs)                  # PCA components (fewer for targeted panels)
    
    t0 = time.perf_counter()  # start timer
    _timer = StageTimer()  # per-stage timing

    if not CELLPOSE_AVAILABLE:  # require cellpose
        raise RuntimeError("Cellpose is required; install with `pip install cellpose`.")  # error

    XENIUM_DIR = args.xenium_dir  # xenium dir
    MORPHOLOGY_IMAGE = args.morphology or os.path.join(XENIUM_DIR, "morphology.ome.tif")  # morphology path
    TRANSCRIPTS_FILE = args.transcripts or os.path.join(XENIUM_DIR, "transcripts.parquet")  # transcripts path
    OUTPUT_DIR = args.outdir or os.path.join(XENIUM_DIR, "cellpose_pipeline_outputs")  # output dir
    FIG_DIR = os.path.join(OUTPUT_DIR, "figures")  # group 1: all figures (PNGs)
    DATA_DIR = os.path.join(OUTPUT_DIR, "data")     # group 2: all numerical data + Seurat bundle

    os.makedirs(FIG_DIR, exist_ok=True)  # make dirs
    os.makedirs(DATA_DIR, exist_ok=True)  # make dirs

    SEED = int(args.seed)  # seed
    N_CLUSTERS = int(args.n_clusters)  # clusters
    MAX_CELLS_UMAP_FIT = int(args.max_cells_umap_fit)  # umap cap
    GPU = bool(args.gpu)  # gpu flag
    TILE_SIZE = int(args.tile_size)  # tile size
    TILE_OVERLAP = int(args.tile_overlap)  # tile overlap
    TILE_BATCH = int(args.tile_batch)  # tile batch
    SKIP_TILE_STD = float(args.skip_tile_std)  # skip threshold
    CPU_WORKERS = int(args.cpu_workers)  # cpu workers
    FORCE_CPU_FALLBACK = bool(args.force_cpu_on_gpu_error)  # fallback
    DO_SCANPY_MARKERS = bool(args.do_scanpy_markers)  # scanpy markers

    print("XENIUM_DIR:", XENIUM_DIR, flush=True)  # log
    print("MORPHOLOGY_IMAGE:", MORPHOLOGY_IMAGE, flush=True)  # log
    print("TRANSCRIPTS_FILE:", TRANSCRIPTS_FILE, flush=True)  # log
    print("GPU:", GPU, "tile_size:", TILE_SIZE, "tile_overlap:", TILE_OVERLAP, "tile_batch:", TILE_BATCH, "skip_tile_std:", SKIP_TILE_STD, flush=True)  # log
    print("Requested clusters:", N_CLUSTERS, "max_cells_umap_fit:", MAX_CELLS_UMAP_FIT, flush=True)  # log
    print("Scanpy available:", SCANPY_AVAILABLE, "UMAP available:", UMAPLEARN_AVAILABLE, flush=True)  # log

    #@@@@@@@@@@@ marker_sets = load_marker_sets(MARKERS_YAML)  # load markers
    #@@@@@@@@@@@if marker_sets:  # if loaded
    #@@@@@@@@@@@    print(f"Loaded marker YAML with {len(marker_sets)} cell #@@@@@@@@@@@ types.", flush=True)  # log
        #@@@@@@@@@@@ else:  # none
           #@@@@@@@@@@@ print("WARNING: No marker YAML loaded; predicted cell types will be 'Unknown'.", flush=True)  # warn

    # ---- Reuse mode: load a previously-saved NUCLEUS mask and skip the entire
    # expensive front end (DAPI load + preview + Cellpose, ~90% of runtime). Lets
    # you iterate on expansion/QC/clustering in ~1 min instead of re-segmenting.
    nuclei_mask_path = os.path.join(DATA_DIR, "full_mask_nuclei.tif")
    _reuse = str(args.reuse_mask).strip()
    reuse_path = (nuclei_mask_path if _reuse.lower() == "auto" else _reuse) if _reuse else ""
    do_reuse = bool(reuse_path) and os.path.exists(reuse_path)
    if _reuse and not do_reuse:
        print(f"[REUSE] requested but mask not found at '{reuse_path}'; "
              f"running full segmentation instead.", flush=True)

    axes, chosen_plane = None, -1  # metadata; populated when segmenting fresh
    if do_reuse:
        print(f"[REUSE] Loading nucleus mask from {reuse_path}; skipping DAPI load, "
              f"preview, and Cellpose.", flush=True)
        mask_full = tiff.imread(reuse_path).astype(np.int32, copy=False)
        print(f"[REUSE] Loaded mask shape={mask_full.shape}, nuclei={int(mask_full.max())}", flush=True)
        _timer.lap("Load cached nucleus mask")
    else:
        mask_full, axes, chosen_plane = _segment_full(
            args, MORPHOLOGY_IMAGE, FIG_DIR, GPU, MODEL_TYPE, MIN_SIZE,
            TILE_SIZE, TILE_OVERLAP, CPU_WORKERS, SKIP_TILE_STD,
            FORCE_CPU_FALLBACK, TILE_BATCH, BOUNDARY_WIDTH, PASTEL_ALPHA,
            PATCH_SIZE_FULL, N_PATCHES, PATCH_STRIDE_FULL, PATCH_PREVIEW_DS,
            PATCH_DS, _timer,
        )
        # Save the pre-expansion nucleus mask so future --reuse_mask runs can
        # re-apply a different --expand_pixels without re-segmenting.
        tiff.imwrite(nuclei_mask_path, mask_full.astype(np.uint32), compression="zlib")
        print(f"[SEG] Saved nucleus mask: {nuclei_mask_path}", flush=True)

    H, W = mask_full.shape  # image dims (mask is full-res H x W; used in metadata)

    # ---- Nucleus expansion with OVERLAP RESOLUTION (applies whether segmented
    # fresh or reused). Instead of plain label dilation (which can produce
    # overlapping/ambiguous territories where neighbouring expansions collide),
    # grow each nucleus into a NON-OVERLAPPING tessellation: contested pixels go
    # to the nearest nucleus, optionally snapping to image edges (edge-aware).
    # QC statistics are recorded in the run log and a JSON sidecar. ----
    _resolve_qc = None
    if int(args.expand_pixels) > 0:
        t_exp = time.perf_counter()
        n_before = int(mask_full.max())
        _use_resolve = bool(getattr(args, "resolve_overlaps", True))
        if _use_resolve:
            from .resolve_overlaps import resolve_from_nucleus_mask
            _edge_w = float(getattr(args, "edge_weight", 0.0))
            _edge_img = None
            if _edge_w > 0:
                # use the DAPI/morphology plane as the edge cue when available
                try:
                    _edge_img = _load_morphology_plane(MORPHOLOGY_IMAGE)  # H x W float
                    if _edge_img is not None and _edge_img.shape != mask_full.shape:
                        _edge_img = None  # shape mismatch -> fall back to distance-only
                except Exception as _e:
                    print(f"[SEG] edge image unavailable ({_e}); distance-only tessellation.", flush=True)
                    _edge_img = None
            # mode: 'expand_labels' (default) reproduces classic expansion exactly
            # and only removes true overlaps -> preserves tissue morphology.
            # 'watershed' fills contact seams / enables edge-aware; opt-in only.
            _mode = str(getattr(args, "resolve_mode", "expand_labels"))
            if _edge_w > 0 and _mode == "expand_labels":
                _mode = "watershed"  # edge-aware requires the watershed surface
            _res = resolve_from_nucleus_mask(
                mask_full, expand_pixels=int(args.expand_pixels),
                edge_image=_edge_img, edge_weight=_edge_w,
                min_area=int(getattr(args, "min_cell_area", 15)),
                max_area=(int(args.max_cell_area) if getattr(args, "max_cell_area", 0) else None),
                min_solidity=float(getattr(args, "min_solidity", 0.0)),
                mode=_mode,
                verbose=True,
            )
            mask_full = _res["labels"].astype(np.int32, copy=False)
            _resolve_qc = _res["qc"]
            print(f"[SEG] Overlap-resolved expansion by {int(args.expand_pixels)} px "
                  f"({n_before} nuclei -> {int(mask_full.max())} cells) in "
                  f"{time.perf_counter()-t_exp:.1f}s", flush=True)
            print(f"[SEG][QC] input_overlap={_resolve_qc['overlap_fraction_input']:.4f} "
                  f"residual_overlap={_resolve_qc['residual_overlap_fraction']:.4f} "
                  f"cells_dropped={_resolve_qc['cells_dropped_implausible']} "
                  f"median_area_px={_resolve_qc['median_cell_area_px']} "
                  f"median_solidity={_resolve_qc['median_solidity']} "
                  f"edge_aware={_resolve_qc['edge_aware']}", flush=True)
            # persist QC sidecar for the run log / downstream reporting
            try:
                with open(os.path.join(DATA_DIR, "overlap_resolution_qc.json"), "w") as _f:
                    json.dump(_resolve_qc, _f, indent=2)
            except Exception as _e:
                print(f"[SEG][QC] could not write QC json: {_e}", flush=True)
        else:
            from skimage.segmentation import expand_labels
            mask_full = expand_labels(mask_full, distance=int(args.expand_pixels)).astype(np.int32, copy=False)
            print(f"[SEG] Expanded nuclei by {int(args.expand_pixels)} px "
                  f"({n_before} labels) in {time.perf_counter()-t_exp:.1f}s "
                  f"[plain dilation; overlap resolution disabled]", flush=True)
        _timer.lap("Nucleus expansion + overlap resolution")

    mask_path = os.path.join(DATA_DIR, "full_mask_cellpose_fullres_tiled.tif")  # mask path
    tiff.imwrite(mask_path, mask_full.astype(np.uint32), compression="zlib")  # save mask
    print("Saved full-res mask:", mask_path, flush=True)  # log
    if not do_reuse and not getattr(args, "skip_seg_exports", False):
        save_segmentation_exports_fullres(mask_full, FIG_DIR, boundary_width=BOUNDARY_WIDTH)  # exports
        _timer.lap("Segmentation PNG exports")

    cols_schema = parquet_schema_columns(TRANSCRIPTS_FILE)  # schema

    if cols_schema is None:  # fallback read
        tx = pd.read_parquet(TRANSCRIPTS_FILE)  # read parquet
        cols_schema = list(tx.columns)  # cols
        x_col = pick_first_present(cols_schema, X_CANDIDATES)  # x col
        y_col = pick_first_present(cols_schema, Y_CANDIDATES)  # y col
        g_col = pick_first_present(cols_schema, GENE_CANDIDATES)  # gene col
        if x_col is None or y_col is None or g_col is None:  # guard
            raise ValueError(f"Could not detect x/y/gene columns. Available: {cols_schema[:80]}")  # error
        tx = tx[[x_col, y_col, g_col]].copy()  # subset
    else:  # fast read
        x_col = pick_first_present(cols_schema, X_CANDIDATES)  # x col
        y_col = pick_first_present(cols_schema, Y_CANDIDATES)  # y col
        g_col = pick_first_present(cols_schema, GENE_CANDIDATES)  # gene col
        if x_col is None or y_col is None or g_col is None:  # guard
            raise ValueError(f"Could not detect x/y/gene columns from parquet schema. Available: {cols_schema[:80]}")  # error
        tx = pd.read_parquet(TRANSCRIPTS_FILE, columns=[x_col, y_col, g_col])  # read needed

    tx[g_col] = tx[g_col].astype(str)  # ensure gene str
    tx = tx.dropna(subset=[x_col, y_col, g_col])  # drop na

    _timer.lap("Segmentation PNG exports")
    tx_assigned, assign_summary = assign_transcripts_to_cells_fullres(mask_full, tx, x_col, y_col)  # assignment
    print("Assigned transcripts:", int(tx_assigned.shape[0]), flush=True)  # log
    print("Assigned cells (unique):", int(tx_assigned["cell_id"].nunique()), flush=True)  # log
    print("Assignment mapping used:", assign_summary["label"], flush=True)  # log

    assign_json = os.path.join(DATA_DIR, "transcript_assignment_summary.json")  # json path
    #@@@@@@with open(os.path.join(OUTPUT_DIR, "cluster_labels_unsupervised.json"), "w") as f:
         #@@@@@@@@@json.dump(cluster_labels, f, indent=2)
    print("Saved:", assign_json, flush=True)  # log

    if int(tx_assigned.shape[0]) == 0:  # none assigned
        raise ValueError("No transcripts assigned to any cell AFTER mapping attempts. This indicates transcripts and morphology are in different coordinate frames.")  # error

    gene_counts = tx_assigned[g_col].value_counts()  # gene counts
    top_genes = set(gene_counts.head(int(TOP_GENE_CAP)).index)  # top genes
    tx_assigned = tx_assigned[tx_assigned[g_col].isin(top_genes)].copy()  # filter
    print("Genes in matrix (post-cap):", int(tx_assigned[g_col].nunique()), flush=True)  # log

    cell_cat = pd.Categorical(tx_assigned["cell_id"])  # cell categories
    gene_cat = pd.Categorical(tx_assigned[g_col])  # gene categories
    row = cell_cat.codes.astype(np.int32, copy=False)  # row indices
    col = gene_cat.codes.astype(np.int32, copy=False)  # col indices
    data = np.ones(int(tx_assigned.shape[0]), dtype=np.int32)  # data ones
    X = sp.coo_matrix((data, (row, col)), shape=(len(cell_cat.categories), len(gene_cat.categories))).tocsr()  # csr matrix

    obs = pd.DataFrame(index=cell_cat.categories.astype(np.int32))  # obs df
    obs.index.name = "cell_id"  # index name

    props = measure.regionprops_table(mask_full, properties=("label", "centroid", "area"))  # props
    cells_df = pd.DataFrame(props).rename(columns={"label": "cell_id", "centroid-1": "x_centroid_px", "centroid-0": "y_centroid_px", "area": "area_px"})  # df
    cells_df["cell_id"] = cells_df["cell_id"].astype(np.int32)  # int
    obs = obs.join(cells_df.set_index("cell_id"), how="left")  # join

    var = pd.DataFrame(index=gene_cat.categories.astype(str))  # var df
    var.index.name = "gene"  # name

    adata = ad.AnnData(X=X, obs=obs, var=var)  # create AnnData

    # Negative-control background subtraction (before filtering, while control
    # features are still present in the matrix).
    if bool(args.subtract_background):
        adata = subtract_negcontrol_background(adata)

    print(f"Cells before filtering: {adata.n_obs}, genes: {adata.n_vars}", flush=True)  # log

    _timer.lap("Transcript assignment + matrix build")
    adata = safe_filter_adata(adata, min_cells_floor=int(MIN_CELLS_FLOOR),
                              min_genes=int(args.min_genes_per_cell),
                              min_cells=int(args.min_cells_per_gene),
                              min_counts=int(args.min_counts_per_cell))  # filter
    
        # --- export cell-by-gene CSV for UMAP (new) ---
    try:
        csv_umap_path = os.path.join(DATA_DIR, "cell_by_gene_for_umap.csv")
        # If you prefer gzip, set compress=True
        save_cell_by_gene_csv(adata, csv_umap_path, use_layer=None, index=True, compress=False)
    except Exception as e:
        print(f"WARNING: failed to save cell-by-gene CSV for UMAP: {type(e).__name__}: {e}", flush=True)
        
    print(f"Cells after filtering: {adata.n_obs}, genes: {adata.n_vars}", flush=True)  # log

    if int(adata.n_obs) == 0:  # guard
        raise ValueError("0 cells after safe filtering; upstream assignment produced almost no usable signal.")  # error

    h5ad_path = os.path.join(DATA_DIR, "xenium_cellpose_fullres_tiled.h5ad")  # h5ad path
    adata.write(h5ad_path)  # write h5ad
    print("Saved AnnData:", h5ad_path, flush=True)  # log
   
    # Doublet removal (after QC, before clustering): drops between-cluster cells.
    if bool(args.remove_doublets):
        adata = remove_doublets(adata, seed=SEED)
        if int(adata.n_obs) == 0:
            raise ValueError("0 cells after doublet removal.")

    # ===== RUN PCA FIRST (CRITICAL) =====
    _timer.lap("QC filter + intermediate write")
    adata = exact_kmeans_clustering(
        adata,
        n_clusters=N_CLUSTERS,
        n_pcs=N_PCS,
        seed=SEED,
        key="cluster",
        normalization=str(args.normalization),
        auto_n_pcs=bool(args.auto_n_pcs)
    )

    print("[CHECK] PCA exists:", "X_pca_fast" in adata.obsm, flush=True)
    print("[CHECK] PCA shape:", adata.obsm["X_pca_fast"].shape, flush=True)
    print("[CLUSTER] Using resolution:", RESOLUTION, flush=True)
    # ========================= CLUSTERING (LEIDEN) =========================
    if SCANPY_AVAILABLE:

        print("[CLUSTER] Running Leiden clustering...", flush=True)
        
        from sklearn.preprocessing import StandardScaler

        # ===== Normalize PCA =====
        X_pca = adata.obsm["X_pca_fast"].copy()

        X_pca = StandardScaler().fit_transform(X_pca)

        # ===== Normalize spatial =====
        X_spatial = adata.obs[
            ["x_centroid_px", "y_centroid_px"]
        ].values.astype(np.float32)

        X_spatial = StandardScaler().fit_transform(X_spatial)

        # ===== Tunable spatial contribution =====
        SPATIAL_WEIGHT = 0.0 #@@ 0.15

        # ===== Combined representation =====
   #@@     X_combined = np.hstack([
   #@@          X_pca,
   #@@         SPATIAL_WEIGHT * X_spatial
   #@@    ])
        X_combined = X_pca.copy()
        
        adata.obsm["X_combined"] = X_combined

        print(
            "[CHECK] X_combined shape:",
            X_combined.shape,
            flush=True
        )

        # Build graph (shared by Leiden + UMAP); n_neighbors/metric match the
        # requested UMAP settings.
        sc.pp.neighbors(
            adata,
            use_rep="X_combined",
            n_neighbors=int(args.umap_neighbors),
            metric=str(args.umap_metric)
        )

        # Leiden clustering (igraph backend + optional resolution tuning)
        chosen_res = run_leiden_best(
            adata,
            base_resolution=RESOLUTION,
            tune=bool(args.tune_resolution),
            key="cluster_leiden",
            seed=SEED,
            rep="X_pca_fast",
            min_clusters=max(3, int(round(0.6 * N_CLUSTERS))),  # avoid resolution collapse
        )
        print(f"[CLUSTER] Leiden done at resolution={chosen_res}", flush=True)

        adata.obs["cluster"] = adata.obs["cluster_leiden"].astype(int)

        # Consolidate over-split clusters (fewer, larger, purer clusters -> stronger
        # DE significance + better interpretability). Off by default.
        if float(args.merge_threshold) > 0 and "X_pca_fast" in adata.obsm:
            n0 = int(adata.obs["cluster"].nunique())
            merged = merge_indistinct_clusters(
                adata.obs["cluster"].to_numpy(),
                adata.obsm["X_pca_fast"],
                threshold=float(args.merge_threshold),
                min_clusters=max(2, int(args.merge_min_clusters)),
            )
            adata.obs["cluster"] = merged.astype(int)
            print(f"[MERGE] consolidated {n0} -> {int(adata.obs['cluster'].nunique())} "
                  f"clusters (threshold={float(args.merge_threshold)}).", flush=True)

        # 🔥 REQUIRED for Scanpy markers
        print(
            "[CHECK] clusters:",
            adata.obs["cluster"].value_counts().sort_values(ascending=False),
            flush=True
        )
        adata.obs["cluster"] = adata.obs["cluster"].astype("category")

        print(
            "[CHECK] clusters:",
            adata.obs["cluster"].nunique(),
            flush=True
        )

        print(
            "[CHECK] dtype:",
            adata.obs["cluster"].dtype,
            flush=True
        )

        # ==========================================================
        # Silhouette diagnostics
        # ==========================================================
        from sklearn.metrics import (
            silhouette_score,
            silhouette_samples
        )

        try:

            cluster_labels = (
                adata.obs["cluster"]
                .astype(int)
                .to_numpy()
            )

            score = silhouette_score_fast(
                adata.obsm["X_pca_fast"],
                cluster_labels,
                seed=SEED,
            )

            sil = silhouette_samples_fast(
                adata.obsm["X_pca_fast"],
                cluster_labels,
                seed=SEED,
            )

            adata.obs["silhouette"] = sil

            print(
                f"[CLUSTER] Mean silhouette = {score:.4f}",
                flush=True
            )

            print(
                f"[CLUSTER] Median silhouette = {np.nanmedian(sil):.4f}",
                flush=True
            )

            print(
                f"[CLUSTER] Min silhouette = {np.nanmin(sil):.4f}",
                flush=True
            )

            print(
                f"[CLUSTER] Max silhouette = {np.nanmax(sil):.4f}",
                flush=True
            )

        except Exception as e:

            print(
                f"[CLUSTER] Silhouette computation failed: {e}",
                flush=True
            )

        # ==========================================================
        # UMAP
        # ==========================================================
        # min_dist=0.05 + spread=2.0 produces tight clumps joined by long
        # filaments that overlap; min_dist~0.3 / spread~1.0 with spectral init
        # gives well-separated, compact groups.
        # min_dist controls within-cluster packing: lower -> tighter clusters and
        # more empty space between them. spread scales the overall layout.
        # Cluster separation WITH preserved relationships: initialize UMAP from
        # PAGA, a cluster-connectivity graph. Related clusters stay adjacent,
        # distinct clusters separate more cleanly than under a spectral/random
        # init, and within-cluster structure is preserved. The PAGA connectivities
        # (adata.uns['paga']) are themselves an interpretable readout of which cell
        # types are related. Falls back to spectral if PAGA is unavailable.
        # min_dist/spread remain layout-only knobs on top of this.
        init = "spectral"
        if str(getattr(args, "umap_init", "spectral")) == "paga" \
                and "cluster" in adata.obs and adata.obs["cluster"].nunique() >= 2:
            try:
                sc.tl.paga(adata, groups="cluster")
                sc.pl.paga(adata, plot=False)  # computes adata.uns['paga']['pos'] used by init
                init = "paga"
                print("[UMAP] PAGA-initialized embedding "
                      "(separation + preserved relationships).", flush=True)
            except Exception as e:
                print(f"[UMAP] PAGA init unavailable ({e}); using spectral.", flush=True)
                init = "spectral"

        sc.tl.umap(
            adata,
            min_dist=float(args.umap_min_dist),
            spread=float(args.umap_spread),
            init_pos=init,
        )

        adata.obsm["X_umap_fast"] = adata.obsm["X_umap"]

        print(
            "UMAP computed using Leiden graph.",
            flush=True
        )

    else:
        print("[WARNING] Scanpy not available → fallback to KMeans", flush=True)

        exact_kmeans_clustering(
            adata,
            n_clusters=N_CLUSTERS,
            n_pcs=N_PCS,
            seed=SEED,
            key="cluster",
            normalization=str(args.normalization),
            auto_n_pcs=bool(args.auto_n_pcs)
        )

    # ========================= FAST MARKERS =========================
    markers_fast = fast_cluster_markers_foldchange(
        adata,
        cluster_key="cluster",
        top_n=10
    )

    markers_fast_csv = os.path.join(DATA_DIR, "cluster_markers_all_fast.csv")
    markers_fast.to_csv(markers_fast_csv, index=False)

    print("Saved FAST markers:", markers_fast_csv, flush=True)

#@@@@@@@@ komlan
#@@@@@@@@         if umap_res is not None:  # if success
#@@@@@@@@             Y_umap, idx_fit = umap_res  # unpack
#@@@@@@@@             adata.obsm["X_umap_fast"] = Y_umap  # store
#@@@@@@@@             print("UMAP computed for all cells.", flush=True)  # log
#@@@@@@@@         else:  # none
#@@@@@@@@             print("WARNING: UMAP not available; install with `pip install umap-learn`.", flush=True)  # warn

    markers_fast = fast_cluster_markers_foldchange(adata, cluster_key="cluster", top_n=10)  # markers
    markers_fast_csv = os.path.join(DATA_DIR, "cluster_markers_all_fast.csv")  # path
    markers_fast.to_csv(markers_fast_csv, index=False)  # save
    print("Saved FAST markers:", markers_fast_csv, flush=True)  # log
    
    # ------------------------- SegMap INITIAL outputs -------------------------
    save_segmap_stage(
        adata,
        DATA_DIR,
        DATA_DIR,
        "initial"
    )
    #@@@@@@@@ komlan
    # ===== INITIAL markers (ALL CELLS) =====
    compute_scanpy_markers_all_cells(
        adata,
        groupby="cluster",
        output_path=os.path.join(DATA_DIR, "cluster_markers_all_cells_initial.csv")
    )
    #@@@@@@@@ komlan
    # ============================================================
    # 🔹 Save initial clusters (ADD HERE)
    # ============================================================
    adata.obs["cluster_initial"] = adata.obs["cluster"].copy()
    print("[CHECK] initial clusters saved:", adata.obs["cluster_initial"].nunique(), flush=True)
    
    # ------------------------- SegMap refinement -------------------------
    _timer.lap("PCA/HVG + Leiden + initial UMAP + markers")
    adata = segmap_refinement_publication(
        adata,
        N_CLUSTERS,
        DATA_DIR
    )

    # ===== Ensure cluster_final exists =====
    if "cluster_final" not in adata.obs:
        raise RuntimeError("cluster_final not found after refinement.")

    # ============================================================
    # Bayesian marker model (Fig. 1D): per-gene spatial-vs-null
    # model selection -> posterior marker prob (p_j) + FSV, and
    # cluster-specific Bayesian markers. Non-fatal add-on.
    # ============================================================
    if bool(getattr(args, "do_bayes_markers", True)):
        _bayes_fit = run_bayes_markers(
            adata,
            cluster_key="cluster_final",
            output_dir=DATA_DIR,
            coord_keys=("x_centroid_px", "y_centroid_px"),
            k=int(getattr(args, "bayes_k", 15)),
            m=int(getattr(args, "bayes_m", 100)),
            alpha=float(getattr(args, "bayes_fdr", 0.05)),
            max_cells=(int(args.bayes_max_cells)
                       if getattr(args, "bayes_max_cells", 0) else None),
            seed=SEED,
        )
        # Close the D->E loop: gene-relevance p_j gates the cell-state feature
        # space, so q_i(phi) and its entropy H_i inherit marker uncertainty.
        if _bayes_fit is not None:
            try:
                couple_pj_to_qi(adata, _bayes_fit, cluster_key="cluster_final",
                                coord_keys=("x_centroid_px", "y_centroid_px"))
                adata.obs[["posterior_entropy_coupled",
                           "posterior_confidence_coupled"]].to_csv(
                    os.path.join(DATA_DIR, "posterior_entropy_coupled.csv"))
                print("[BAYES] p_j->q_i coupling done "
                      "(marker-weighted cell-state entropy written).", flush=True)
            except Exception as _e:
                print(f"[BAYES] coupling skipped (non-fatal): {_e}", flush=True)
        
    # ============================================================
    # 🔥 STEP 2 — BUILD ADVANCED BAYESIAN REPRESENTATION
    # ============================================================

    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors
    from scipy.special import softmax
    from scipy.sparse import csr_matrix

    print("\n[STEP2] Building advanced Bayesian latent space...", flush=True)

    # ============================================================
    # 1. PCA NORMALIZATION + WHITENING
    # ============================================================

    X_pca = adata.obsm["X_pca_fast"].copy()

    X_pca = StandardScaler().fit_transform(X_pca)

    # additional whitening
    pca_whiten = PCA(
        n_components=min(50, X_pca.shape[1]),
        whiten=True,
        random_state=SEED
    )

    X_pca = pca_whiten.fit_transform(X_pca)

    print("[STEP2] PCA whitened:", X_pca.shape, flush=True)

    # ============================================================
    # 2. SPATIAL NORMALIZATION
    # ============================================================

    X_spatial = adata.obs[
        ["x_centroid_px", "y_centroid_px"]
    ].values.astype(np.float32)

    X_spatial = StandardScaler().fit_transform(X_spatial)

    # ============================================================
    # 3. SOFT BAYESIAN POSTERIOR
    # ============================================================

    adata.obs["cluster_final"] = (
        adata.obs["cluster_final"]
        .astype("category")
    )

    cluster_codes = (
        adata.obs["cluster_final"]
        .cat.codes
        .values
    )

    K = int(adata.obs["cluster_final"].nunique())

    print("[STEP2] K =", K, flush=True)

    # soft posterior instead of hard one-hot
    posterior_soft = softmax(
        adata.obsm["X_pca_fast"][:, :K],
        axis=1
    )

    # entropy sharpening
    posterior_soft = posterior_soft ** 1.5
    posterior_soft = (
        posterior_soft /
        posterior_soft.sum(1, keepdims=True)
    )

    # ============================================================
    # 4. LOCAL GRAPH DENSITY
    # ============================================================

    nn_density = NearestNeighbors(
        n_neighbors=15,
        metric="euclidean"
    )

    nn_density.fit(X_pca)

    distances, _ = nn_density.kneighbors(X_pca)

    density = 1.0 / (
        distances.mean(1) + 1e-8
    )

    density = (
        density - density.mean()
    ) / (
        density.std() + 1e-8
    )

    density = density.reshape(-1, 1)

    # ============================================================
    # 5. CONTRASTIVE FEATURE AMPLIFICATION
    # ============================================================

    _timer.lap("SegMap Bayesian refinement")
    X_bayes = np.hstack([
        1.5 * X_pca,
        0.8 * posterior_soft,
        (0.35 * SPATIAL_WEIGHT_BAYES) * X_spatial,
        (0.25 * SPATIAL_WEIGHT_BAYES) * density
    ])

    # ============================================================
    # 6. FEATURE STABILIZATION
    # ============================================================

    rng = np.random.default_rng(SEED)

    X_bayes += (
        0.01 * rng.normal(size=X_bayes.shape)
    )

    X_bayes = StandardScaler().fit_transform(X_bayes)

    adata.obsm["X_bayes"] = X_bayes

    print(
        "[STEP2] X_bayes:",
        X_bayes.shape,
        flush=True
    )
    # ============================================================
    # 7. GRAPH CONSTRUCTION (COSINE)
    # ============================================================
#@@@@@@@@@@@
    print(
        "[STEP2] Computing Bayesian graph...",
        flush=True
    )

    sc.pp.neighbors(
        adata,
        use_rep="X_bayes",
        n_neighbors=30,
        metric="cosine"
    )

    # ============================================================
    # 8. LEIDEN REFINEMENT
    # ============================================================

    sc.tl.leiden(
        adata,
        resolution=RESOLUTION * 1.25,
        key_added="cluster_final_refined",
        flavor="igraph",
        n_iterations=4
    )

    # The X_bayes Leiden is a SPATIAL-DOMAIN / niche label. Keep it separate so it
    # never overrides transcriptomic cell identity (that override was what made the
    # "final" clusters look worse than the "initial" ones).
    adata.obs["spatial_domain"] = (
        adata.obs["cluster_final_refined"].astype("category")
    )

    # Cell types stay transcriptomic (cluster_initial). Spatial context is used only
    # to DENOISE low-confidence, isolated misassignments (confidence-gated majority
    # vote over spatial neighbors) -> clean cell types that consider spatial domains.
    _coords_xy = adata.obs[["x_centroid_px", "y_centroid_px"]].to_numpy()
    _conf = (adata.obs["posterior_confidence"].to_numpy()
             if "posterior_confidence" in adata.obs else None)
    _celltypes = spatial_denoise_labels(
        adata.obs["cluster_initial"].astype(int).to_numpy(),
        _coords_xy,
        confidence=_conf,
        k=15,
        agree_frac=float(args.spatial_denoise_agree),
        n_iter=int(args.spatial_denoise_iters),
        seed=SEED,
    )
    adata.obs["cluster_final"] = pd.Categorical(
        pd.Series(_celltypes.astype(str), index=adata.obs_names)
    )

    print(
        "[STEP2] cell-type clusters (denoised):",
        adata.obs["cluster_final"].nunique(),
        "| spatial domains:",
        adata.obs["spatial_domain"].nunique(),
        flush=True
    )

    # ============================================================
    # 9. HIGH-SEPARATION UMAP
    # ============================================================

    print(
        "[STEP2] Computing final UMAP...",
        flush=True
    )

    # Display cell types on the transcriptomic embedding they were derived from,
    # so the final UMAP is as clean as the initial one (and we skip a redundant
    # UMAP recompute on the hybrid space).
    adata.obsm["X_umap_final"] = adata.obsm["X_umap_fast"]

    # ============================================================
    # 10. OPTIONAL GRAPH CLEANING
    # ============================================================

    # remove tiny disconnected islands
    cluster_sizes = adata.obs["cluster_final"].value_counts()

    small_clusters = cluster_sizes[cluster_sizes < 25].index

    mask_small = adata.obs["cluster_final"].isin(small_clusters)

    # IMPORTANT: allow new category
    adata.obs["cluster_final"] = adata.obs["cluster_final"].cat.add_categories(["noise"])

    adata.obs.loc[mask_small, "cluster_final"] = "noise"

    adata.obs["cluster_final"] = adata.obs["cluster_final"].astype("category")

    print("[STEP2] Bayesian representation complete.", flush=True)

    # ============================================================
    # FINAL BAYESIAN FEATURE SPACE
    # ============================================================

    adata.obs["cluster_final"] = (
        adata.obs["cluster_final"]
        .astype("category")
    )

    adata.obsm["X_bayes"] = X_bayes

    # ============================================================
    # 🔥 L2 NORMALIZATION OF BAYESIAN SPACE
    # ============================================================

    from sklearn.preprocessing import normalize

    adata.obsm["X_bayes"] = normalize(
        adata.obsm["X_bayes"],
        norm="l2"
    )
    
    # ============================================================
    # 🔥 NORMALIZED + BALANCED BAYESIAN FEATURE SPACE
    # ============================================================
    print("[CHECK] X_bayes std:", np.std(X_bayes), flush=True)
    print("[CHECK] std spatial:", np.std(X_spatial), flush=True)
    
    from sklearn.metrics import silhouette_score

    # --------------------------------------------------
    # BEFORE Bayesian refinement
    # --------------------------------------------------
    score_before = silhouette_score_fast(
        np.asarray(adata.obsm["X_bayes"]),
        adata.obs["cluster_initial"].astype("category").cat.codes.to_numpy(),
        seed=SEED,
    )

    print(
        f"[STEP2] Silhouette before Bayesian Leiden = {score_before:.4f}",
        flush=True
    )


    sil = silhouette_samples_fast(
        np.asarray(adata.obsm["X_bayes"]),
        adata.obs["cluster_final"].cat.codes.to_numpy(),
        seed=SEED,
    )
    adata.obs["silhouette"] = sil
    for cl in adata.obs["cluster_final"].cat.categories:
        vals = sil[
            (adata.obs["cluster_final"] == cl).to_numpy()
        ]
        vals = vals[np.isfinite(vals)]  # drop subsample NaNs
        if vals.size == 0:
            continue
        print(
            f"Cluster {cl}: "
            f"mean={vals.mean():.3f} "
            f"median={np.median(vals):.3f}",
            flush=True
        )
    
    # --------------------------------------------------
    # AFTER Bayesian refinement
    # --------------------------------------------------
    score_after = silhouette_score_fast(
        np.asarray(adata.obsm["X_bayes"]),
        adata.obs["cluster_final"].astype("category").cat.codes.to_numpy(),
        seed=SEED,
    )

    print(
        f"[STEP2] Silhouette after Bayesian Leiden = {score_after:.4f}",
        flush=True
    )

    print(
        f"[STEP2] Silhouette improvement = "
        f"{score_after - score_before:.4f}",
        flush=True
    )

    # --------------------------------------------------
    # Diagnostic checks
    # --------------------------------------------------
    changed = (
        adata.obs["cluster_initial"].astype(str)
        !=
        adata.obs["cluster_final"].astype(str)
    ).sum()

    print(
        f"[CHECK] cells changed: {changed}",
        flush=True
    )

    print(
        "[CHECK] initial clusters:",
        adata.obs["cluster_initial"].nunique(),
        flush=True
    )

    print(
        "[CHECK] final clusters:",
        adata.obs["cluster_final"].nunique(),
        flush=True
    )

    # ------------------------- Unsupervised annotation -------------------------
    cluster_labels = auto_annotate_clusters(
        adata,
        cluster_key="cluster_final"
    )

    adata.obs["cell_type_unsupervised"] = (
        adata.obs["cluster_final"].map(cluster_labels)
    )

    # 🔹 unify cluster key
    adata.obs["cluster"] = adata.obs["cluster_final"]
    
    # ============================================================
    # 🔍 CHECK how many cells changed (ADD HERE)
    # ============================================================
    changed = (
        adata.obs["cluster_initial"].astype(str) !=
        adata.obs["cluster_final"].astype(str)
    ).sum()
    changed_pct = changed / adata.n_obs
    print("[CHECK] % cells changed:", changed_pct, flush=True)
    print("[CHECK] cells changed:", changed, flush=True)
    
    # ============================================================
    # 🔍 CHECK cluster changes (ADD HERE)
    # ============================================================
    print("[CHECK] initial clusters:", adata.obs["cluster_initial"].nunique(), flush=True)
    print("[CHECK] final clusters:", adata.obs["cluster_final"].nunique(), flush=True)

    
    # ============================================================
    # 🔍 CHECK cluster changes (ADD HERE)
    # ============================================================
    print("[CHECK] initial clusters:", adata.obs["cluster"].nunique(), flush=True)
    print("[CHECK] final clusters:", adata.obs["cluster_final"].nunique(), flush=True)
    
    # ===== FINAL markers (ALL CELLS) =====
    compute_scanpy_markers_all_cells(
        adata,
        groupby="cluster_final",
        output_path=os.path.join(DATA_DIR, "cluster_markers_all_cells_final.csv")
    )

    # ------------------------- Validate annotation -------------------------
    assert "cell_type_unsupervised" in adata.obs.columns, \
        "cell_type_unsupervised missing before downstream steps!"

    # ------------------------- Save cluster labels -------------------------
    labels_json = os.path.join(DATA_DIR, "cluster_labels_unsupervised.json")

    print("DEBUG: cluster_labels exists?", "cluster_labels" in locals(), flush=True)
    assert "cluster_labels" in locals(), "cluster_labels NOT defined before saving!"

    with open(labels_json, "w") as f:
        json.dump(cluster_labels, f, indent=2)

    print("Saved:", labels_json, flush=True)

    # ========================= NOW SAFE TO PLOT =========================

    # ------------------------- CSV outputs -------------------------
    _timer.lap("X_bayes clustering + final UMAP + final markers")
    save_cluster_tables_csv(adata, DATA_DIR)
    cluster_separability_report(
        adata, cluster_key="cluster_final", rep="X_pca_fast",
        out_csv=os.path.join(DATA_DIR, "cluster_separability.csv"), seed=SEED)

    # ------------------------- UMAP plots -------------------------
    if "X_umap_fast" in adata.obsm:
        save_umap_plot(
            adata,
            os.path.join(FIG_DIR, "umap_clusters.png"),
            "cluster",
            "UMAP (clusters)"
        )

        save_umap_plot(
            adata,
            os.path.join(FIG_DIR, "umap_celltypes_unsupervised.png"),
            "cell_type_unsupervised",
            "UMAP (unsupervised cell types)"
        )

    # ------------------------- Spatial plots -------------------------
    save_spatial_plot(
        adata,
        os.path.join(FIG_DIR, "spatial_clusters.png"),
        "cluster",
        "Spatial (clusters)"
    )

    save_spatial_plot(
        adata,
        os.path.join(FIG_DIR, "spatial_celltypes_unsupervised.png"),
        "cell_type_unsupervised",
        "Spatial (unsupervised cell types)"
    )

    # ------------------------- Spatial-domain (niche) plot -------------------------
    if "spatial_domain" in adata.obs:
        save_spatial_plot(
            adata,
            os.path.join(FIG_DIR, "spatial_domains.png"),
            "spatial_domain",
            "Spatial domains (niches)"
        )

    # ------------------------- FINAL SegMap outputs -------------------------
    save_segmap_stage(
        adata,
        DATA_DIR,
        DATA_DIR,
        "final"
    )

    _timer.lap("Plots + CSV exports")
    # ------------------------- Save final AnnData -------------------------
    adata_final_path = os.path.join(DATA_DIR, "xenium_cellpose_fullres_tiled_final.h5ad")
    adata.write(adata_final_path)
    print("Saved final AnnData:", adata_final_path, flush=True)

    # ------------------------- Seurat-ready export -------------------------
    export_seurat_bundle(adata, os.path.join(DATA_DIR, "seurat"))
    _timer.lap("Write final AnnData")

    # ------------------------- Summary -------------------------
    dt = time.perf_counter() - t0

    summary = {
        "xenium_dir": XENIUM_DIR,
        "morphology_image": MORPHOLOGY_IMAGE,
        "transcripts_file": TRANSCRIPTS_FILE,
        "axes": axes,
        "chosen_plane": int(chosen_plane),
        "image_H": int(H),
        "image_W": int(W),
        "n_cells_mask": int(mask_full.max()),
        "n_cells_matrix": int(adata.n_obs),
        "n_genes_matrix": int(adata.n_vars),
        "n_clusters": int(adata.obs["cluster"].nunique()),
        "assignment": assign_summary,
        "overlap_resolution_qc": _resolve_qc,
        "runtime_seconds": float(dt),
    }

    summary_json = os.path.join(DATA_DIR, "pipeline_summary.json")
    with open(summary_json, "w") as f:
        json.dump(summary, f, indent=2)

    print("Saved:", summary_json, flush=True)
    print(f"[DONE] Total runtime: {dt/60:.2f} minutes", flush=True)
    _timer.summary()
