# =============================================================================
# ADD to pipeline.py — Tier-1 iterative refinement (data-agnostic, off by default,
# side-analysis). Place AFTER cluster_final + tx_assigned + cells_df exist
# (same location as the single-pass boundary-refine hook). Wrap in the isolation
# guard so a future edit that mutates a shared object fails loudly.
# =============================================================================
if bool(getattr(args, "run_iterative_refine", False)):
    try:
        from .iterative_refine import run_iterative_refine
        from .refine_guard import run_refinement_guarded  # isolation assertion
        import numpy as _np

        _genes = list(adata.var_names)
        _cent = (cells_df.set_index("cell_id")
                 .reindex(range(1, int(mask_full.max()) + 1))
                 [["x_centroid_px", "y_centroid_px"]].to_numpy(float))
        _uxy = tx_assigned[[x_col, y_col]].to_numpy(float)
        _gmap = {g: k for k, g in enumerate(_genes)}
        _fi = tx_assigned[g_col].map(_gmap).to_numpy()
        _ok = ~_np.isnan(_fi); _fi = _np.where(_ok, _fi, -1).astype(int)
        _a0 = (tx_assigned["cell_id"].to_numpy() - 1).astype(int)
        _a0[(_a0 < 0) | (_a0 >= _cent.shape[0]) | (~_ok)] = -1

        # Tier-1-only by default: no validators, no biology. (Tier 2/3 plug-ins can
        # be added later by passing validators=[...] here.)
        run_refinement_guarded(
            {"adata": adata, "tx_assigned": tx_assigned},
            run_iterative_refine,
            _a0, _uxy, _fi, _cent, _cent.shape[0], len(_genes),
            os.path.join(DATA_DIR, "iterative_refine"),
            update_kwargs=dict(tau=float(getattr(args, "refine_tau", 0.2)),
                               lam=float(getattr(args, "refine_lambda", 1.0)),
                               inner_iter=int(getattr(args, "refine_max_iter", 5))),
            max_iter=int(getattr(args, "iter_max_iter", 5)),
            tol=float(getattr(args, "iter_tol", 1e-3)),
        )
    except Exception as _e:
        print(f"[ITER] skipped (non-fatal): {_e}", flush=True)

# --- ADD to cli.py parse_args() ---
#   parser.add_argument("--run_iterative_refine",
#       type=lambda v: str(v).lower() not in ("0","false","no","off"), default=False,
#       help="Run data-agnostic iterative refinement (Tier-1 diagnostics) as an "
#            "optional side-analysis; writes data/iterative_refine/. Off by default.")
#   parser.add_argument("--iter_max_iter", type=int, default=5,
#       help="Max outer iterations for iterative refinement.")
#   parser.add_argument("--iter_tol", type=float, default=1e-3,
#       help="Convergence tol: stop when fraction reassigned between iterations < tol.")
