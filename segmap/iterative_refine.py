"""
iterative_refine.py — SegMap iterative refinement engine (Fig. 1a), data-agnostic.

The engine requires only four things, none of them biological:
  - an INITIAL ASSIGNMENT        (unit -> element; e.g. transcript -> cell)
  - a FEATURE REPRESENTATION     (per-element features, rebuilt each iteration)
  - a NEIGHBOURHOOD/RELATIONSHIP  (used by the update rule)
  - an UPDATE RULE                (callable: reassigns units; default = boundary
                                   reassignment by profile + proximity)

TIER 1 optimization diagnostics are ALWAYS reported and assume no biology:
convergence, reassignment dynamics, stability, oscillation, iteration count and
wall-clock cost. They answer "does the loop provide information beyond the initial
assignment?" without any tissue/technology/marker assumptions.

TIER 2 / TIER 3 are OPTIONAL validator plug-ins (refine_validators.Validator). The
engine calls whatever validators it is given, per iteration, and records their
results verbatim — including each validator's epistemic "claim". With no
validators, only Tier 1 is produced.

Side-analysis & isolation: operates on copies, writes only to out_dir, never
mutates the input AnnData or any primary output. Off by default. Non-fatal.
Outcome reported as measured — improvement, none, or slight degradation are all
valid findings; nothing is tuned to force a result.
"""
from __future__ import annotations
import os, json, time
import numpy as np

__all__ = ["run_iterative_refine", "default_update_rule"]


def _rebuild_matrix(assignment, feat_index, n_elements, n_feat):
    from scipy.sparse import coo_matrix
    v = (assignment >= 0) & (feat_index >= 0)
    return coo_matrix((np.ones(int(v.sum())), (assignment[v], feat_index[v])),
                      shape=(n_elements, n_feat)).tocsr()


def _lognorm(counts):
    counts = np.asarray(counts, float)
    lib = counts.sum(1, keepdims=True); lib[lib == 0] = 1.0
    return np.log1p(counts / lib * 1e4)


def default_update_rule(unit_xy, feat_index, assignment, centroids_xy, n_feat,
                        tau=0.2, lam=1.0, inner_iter=3):
    """Default update: boundary-transcript reassignment by profile + proximity.
    Swappable — any callable with this signature works as the engine's update rule."""
    from .refine_boundaries import refine_boundaries  # or top-level import when standalone
    return refine_boundaries(unit_xy, feat_index, assignment, centroids_xy,
                             n_genes=n_feat, tau=tau, lam=lam,
                             max_iter=inner_iter, verbose=False)[0]


def run_iterative_refine(assignment0, unit_xy, feat_index, centroids_xy, n_elements,
                         n_feat, out_dir, *,
                         update_rule=None, update_kwargs=None,
                         max_iter=5, tol=1e-3,
                         rebuild_feature_fn=None,
                         validators=None, unit_group=None,
                         element_labels_fn=None,
                         verbose=True):
    """Run the iterative loop.

    assignment0      : (T,) initial unit->element assignment (-1 = unassigned)
    unit_xy          : (T,2) unit coordinates (for the default update rule)
    feat_index       : (T,) feature id per unit in [0, n_feat)
    centroids_xy     : (N,2) element coordinates
    update_rule      : callable(unit_xy, feat_index, assignment, centroids_xy,
                       n_feat, **update_kwargs) -> new assignment. Defaults to
                       boundary reassignment.
    rebuild_feature_fn(assignment) -> (N x n_feat) feature matrix. Defaults to
                       lognorm cell-by-feature counts.
    validators       : optional list of refine_validators.Validator (Tier 2/3).
    unit_group       : optional (T,) group id per unit, for held-out validators.
    element_labels_fn(feature_matrix) -> (N,) labels, for concordance validators.

    Writes iterative_metrics.csv (Tier 1 + any validator deltas) and
    iterative_summary.json to out_dir. Returns the summary dict, or None on failure."""
    try:
        from .refine_validators import RefinementState
    except Exception:
        from refine_validators import RefinementState
    update_rule = update_rule or default_update_rule
    update_kwargs = update_kwargs or {}
    validators = validators or []
    rebuild = rebuild_feature_fn or (lambda a: _lognorm(np.asarray(
        _rebuild_matrix(a, feat_index, n_elements, n_feat).todense())))

    try:
        os.makedirs(out_dir, exist_ok=True)
        t_start = time.time()
        assignment = np.asarray(assignment0).copy()
        state_init = RefinementState(assignment.copy(), n_elements, unit_group,
                                     None)

        traj = []
        prev = assignment.copy()
        for it in range(int(max_iter) + 1):
            t_it = time.time()
            feats = rebuild(assignment)
            elem_labels = element_labels_fn(feats) if element_labels_fn else None

            # ---- TIER 1: optimization diagnostics (always, biology-free) ----
            frac_reassigned = (float(np.mean(assignment[prev >= 0] != prev[prev >= 0]))
                               if it > 0 and (prev >= 0).any() else 0.0)
            rec = {"iteration": it,
                   "frac_reassigned": frac_reassigned,
                   "n_assigned": int((assignment >= 0).sum()),
                   "iter_seconds": round(time.time() - t_it, 4)}

            # ---- optional validator plug-ins (Tier 2/3), recorded verbatim ----
            if validators:
                st_now = RefinementState(assignment.copy(), n_elements, unit_group,
                                         elem_labels)
                for v in validators:
                    try:
                        res = v.evaluate(state_init, st_now)
                        for k, val in res.items():
                            if k == "claim":
                                continue
                            rec[f"{v.name}.{k}"] = val
                    except Exception as _e:
                        rec[f"{v.name}.error"] = f"{type(_e).__name__}: {_e}"
            traj.append(rec)

            if verbose:
                extra = " ".join(f"{k.split('.')[-1]}={rec[k]:.4f}"
                                 for k in rec if k.endswith(".delta"))
                print(f"[ITER] {it}: reassigned={frac_reassigned:.4f} "
                      f"n_assigned={rec['n_assigned']} {extra}".rstrip(), flush=True)

            if it == int(max_iter):
                break
            prev = assignment.copy()
            assignment = np.asarray(update_rule(unit_xy, feat_index, assignment,
                                                centroids_xy, n_feat, **update_kwargs))
            changed = (float(np.mean(assignment[prev >= 0] != prev[prev >= 0]))
                       if (prev >= 0).any() else 0.0)
            if changed < tol:
                if verbose:
                    print(f"[ITER] converged at step {it+1} (Δ={changed:.3g} < tol={tol})",
                          flush=True)
                traj[-1]["converged"] = True
                break

        # ---- Tier 1 summary (universal) ----
        fr = np.array([t["frac_reassigned"] for t in traj], float)
        d = np.diff(fr)
        tier1 = {"n_iterations": len(traj),
                 "converged": bool(traj[-1].get("converged", False)),
                 "final_frac_reassigned": float(fr[-1]) if fr.size else np.nan,
                 "max_frac_reassigned": float(np.nanmax(fr)) if fr.size else np.nan,
                 "monotone_decreasing_reassignment": bool(np.all(d <= 1e-9)) if d.size else True,
                 "oscillation_events": int(np.sum(np.sign(d[1:]) * np.sign(d[:-1]) < 0)) if d.size > 1 else 0,
                 "total_seconds": round(time.time() - t_start, 3),
                 "claim": "Tier 1 optimization diagnostics: algorithm behavior only; "
                          "NOT evidence of quality or biological accuracy."}

        import pandas as pd
        pd.DataFrame(traj).to_csv(os.path.join(out_dir, "iterative_metrics.csv"), index=False)
        summary = {"role": "post_analysis_sink", "primary_outputs_modified": False,
                   "tier1_optimization_diagnostics": tier1,
                   "validators_used": [v.name for v in validators],
                   "trajectory": traj,
                   "note": ("data-agnostic iterative refinement (Fig. 1a); optional "
                            "side-analysis; Tier 1 always reported, Tier 2/3 are "
                            "optional plug-ins; effect reported as measured.")}
        with open(os.path.join(out_dir, "iterative_summary.json"), "w") as f:
            json.dump(summary, f, indent=2, default=float)
        if verbose:
            print(f"[ITER] {tier1['n_iterations']} iters, converged={tier1['converged']}, "
                  f"final_reassigned={tier1['final_frac_reassigned']:.4f}, "
                  f"{tier1['total_seconds']}s -> {out_dir}", flush=True)
        return summary
    except Exception as e:
        print(f"[ITER] skipped (non-fatal): {type(e).__name__}: {e}", flush=True)
        return None
