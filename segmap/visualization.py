"""UMAP / spatial plots and cluster CSV exports."""
from .common import *  # noqa: F401,F403


def save_umap_plot(
    adata: ad.AnnData,
    out_png: str,
    color_key: str,
    title: str,
    embed_key: Optional[str] = None,
) -> None:
    # Use the embedding that MATCHES the labels being shown. The final clusters
    # live in the X_bayes space, whose embedding is X_umap_final; plotting them on
    # the transcriptomic-only X_umap_fast makes coherent clusters look scrambled.
    if embed_key is None:
        embed_key = "X_umap_final" if "X_umap_final" in adata.obsm else "X_umap_fast"
    if embed_key not in adata.obsm:
        return
    Y = np.asarray(adata.obsm[embed_key])
    labels = adata.obs[color_key].astype(str).to_numpy()

    # Legend in numeric order when labels are numeric (avoids 0,1,10,11,2,3,...)
    uniq = list(np.unique(labels))
    def _key(s):
        try:
            return (0, float(s))
        except Exception:
            return (1, s)
    uniq = sorted(uniq, key=_key)

    # High-contrast categorical palette (60 distinct colors) instead of tab20,
    # which recycles perceptually-similar hues once you pass ~10 clusters.
    base_colors = []
    for nm in ("tab20", "tab20b", "tab20c"):
        base_colors += list(plt.get_cmap(nm).colors)
    color_of = {}
    ci = 0
    for lab in uniq:
        if lab == "noise":
            color_of[lab] = (0.82, 0.82, 0.82)  # grey for noise
        else:
            color_of[lab] = base_colors[ci % len(base_colors)]
            ci += 1

    # Draw large clusters first so small ones land on top (visible, not buried).
    sizes = {lab: int((labels == lab).sum()) for lab in uniq}
    draw_order = sorted(uniq, key=lambda l: sizes[l], reverse=True)

    fig, ax = plt.subplots(figsize=(11, 9), dpi=200)
    for z, lab in enumerate(draw_order):
        mask = labels == lab
        ax.scatter(
            Y[mask, 0], Y[mask, 1],
            s=4, alpha=0.6,
            color=color_of[lab],
            linewidths=0,
            rasterized=True,           # keeps the saved file light despite many points
            zorder=z + 2,              # small clusters (drawn later) sit on top
            label=f"{lab} (n={sizes[lab]})",
        )

    # Centroid labels so identity is readable regardless of color collisions.
    # Only for short ids (cluster numbers); long cell-type names would clutter.
    annotate_centroids = all(len(str(l)) <= 4 for l in uniq if l != "noise")
    for lab in uniq:
        if lab == "noise" or not annotate_centroids:
            continue
        m = labels == lab
        if m.sum() == 0:
            continue
        cx, cy = np.median(Y[m, 0]), np.median(Y[m, 1])
        ax.text(cx, cy, str(lab), fontsize=9, fontweight="bold",
                ha="center", va="center", zorder=1000,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.6))

    ax.set_title(title)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")

    # Legend in numeric order (not draw order)
    handles = {h.get_label(): h for h in ax.collections}
    ordered = [f"{lab} (n={sizes[lab]})" for lab in uniq]
    ax.legend(
        [handles[o] for o in ordered if o in handles],
        [o for o in ordered if o in handles],
        bbox_to_anchor=(1.02, 1), loc="upper left",
        fontsize=8, markerscale=2.5, frameon=False,
    )

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()
    print("Saved:", out_png, "(embedding:", embed_key + ")", flush=True)


def save_spatial_plot(
    adata: ad.AnnData,
    out_png: str,
    color_key: str,
    title: str
) -> None:

    if (
        ("x_centroid_px" not in adata.obs.columns) or
        ("y_centroid_px" not in adata.obs.columns)
    ):
        return
    x = adata.obs["x_centroid_px"].to_numpy(dtype=np.float32)
    y = adata.obs["y_centroid_px"].to_numpy(dtype=np.float32)
    labels = adata.obs[color_key].astype(str).to_numpy()
    unique_labels = np.unique(labels)
    cmap = plt.get_cmap("tab20", len(unique_labels))

    fig, ax = plt.subplots(figsize=(10, 8), dpi=250)
    for idx, lab in enumerate(unique_labels):

        mask = labels == lab
        n_cells = int(mask.sum())

        ax.scatter(
            x[mask],
            y[mask],
            s=5,
            alpha=0.85,
            color=cmap(idx),
            label=f"{lab} (n={n_cells})"
        )
    ax.invert_yaxis()

    ax.set_title(title)
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (px)")

    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        fontsize=8,
        markerscale=2,
        frameon=False
    )
    plt.tight_layout()
    plt.savefig(
        out_png,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()
    print("Saved:", out_png, flush=True)  # log


def save_cluster_tables_csv(adata: ad.AnnData, fig_dir: str) -> None:  # save CSV tables for cluster plots
    os.makedirs(fig_dir, exist_ok=True)  # ensure output directory exists

    spatial_cols = ["cluster"]
    if "cell_type_unsupervised" in adata.obs.columns:
        spatial_cols.append("cell_type_unsupervised")
    
    for c in ["x_centroid_px", "y_centroid_px"]:  # centroid coordinate columns
        if c in adata.obs.columns:  # if present
            spatial_cols.append(c)  # include column

    spatial_df = adata.obs[spatial_cols].copy()  # build spatial table
    spatial_df.insert(0, "cell_id", adata.obs_names.astype(str))  # add cell_id column
    spatial_csv = os.path.join(fig_dir, "spatial_clusters.csv")  # output path
    spatial_df.to_csv(spatial_csv, index=False)  # write CSV
    print("Saved:", spatial_csv, flush=True)  # log

    if "X_umap_fast" in adata.obsm:  # only if UMAP coordinates exist
        Y = np.asarray(adata.obsm["X_umap_fast"])  # get UMAP coords array
        umap_cols = ["cluster"]  # always include cluster
        
        if "cell_type_unsupervised" in adata.obs.columns:
            umap_cols.append("cell_type_unsupervised")

        umap_df = adata.obs[umap_cols].copy()  # build base table
        umap_df.insert(0, "cell_id", adata.obs_names.astype(str))  # add cell_id
        umap_df.insert(1, "umap1", Y[:, 0].astype(np.float32))  # add UMAP1
        umap_df.insert(2, "umap2", Y[:, 1].astype(np.float32))  # add UMAP2

        umap_csv = os.path.join(fig_dir, "umap_clusters.csv")  # output path
        umap_df.to_csv(umap_csv, index=False)  # write CSV
        print("Saved:", umap_csv, flush=True)  # log
    else:  # if UMAP not present
        print("NOTE: UMAP not present; skipping fig/umap_clusters.csv", flush=True)  # log
