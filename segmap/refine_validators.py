"""
refine_validators.py — generic, optional validator interface for iterative
refinement. Tier 1 (optimization diagnostics) lives in the refinement engine and
is ALWAYS reported. Tiers 2 and 3 are OPTIONAL PLUG-INS that conform to this
interface; the engine calls whatever validators it is handed and nothing more.

No biology, no gene panels, no Baysor, no tissue assumptions are hard-coded here.
A validator is any object with:

    name : str
    evaluate(state_before, state_after) -> dict
        returns a metrics dict; MUST include a "claim" key stating the epistemic
        level of the result (what it does and does NOT establish).

`state_*` is a RefinementState (below): a minimal, modality-agnostic container of
the current assignment + feature representation. Validators read from it; they
never mutate it and never touch the primary workflow.

Two reference plug-ins are provided (HeldoutFeatureValidator, ConcordanceValidator)
but are entirely optional — the engine runs with an empty validator list.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np


@dataclass
class RefinementState:
    """Modality-agnostic snapshot the engine passes to validators.
    - assignment  : (T,) int, unit -> element id (e.g. transcript -> cell), -1 = none
    - unit_group  : (T,) group id per unit for group-based validators (optional)
    - element_labels : (N,) label per element (e.g. per-cell cluster), optional
    - n_elements  : int
    Nothing here is gene/tissue/technology specific."""
    assignment: np.ndarray
    n_elements: int
    unit_group: np.ndarray = None
    element_labels: np.ndarray = None


class Validator:
    """Interface. Subclass and implement evaluate()."""
    name = "validator"
    def evaluate(self, before: RefinementState, after: RefinementState) -> dict:
        raise NotImplementedError


# --------------------------------------------------------------------------- #
# OPTIONAL PLUG-IN — Tier 2: held-out feature generalization (data-agnostic).
# The caller decides what a "feature group" is (gene set, pathway, spatial module,
# image feature, network community, ...). This class only counts group concentration.
# --------------------------------------------------------------------------- #
class HeldoutFeatureValidator(Validator):
    name = "tier2_heldout_feature"
    def __init__(self, holdout_groups):
        self.holdout = np.asarray([str(g) for g in holdout_groups])

    def _coherence(self, st: RefinementState):
        if st.unit_group is None:
            return np.nan
        keep = np.isin(st.unit_group, self.holdout)
        valid = (st.assignment >= 0) & keep
        if not valid.any():
            return np.nan
        elem = st.assignment[valid]; grp = st.unit_group[valid]
        uniq = np.asarray(sorted(set(grp.tolist())))
        gi = {g: k for k, g in enumerate(uniq)}
        M = np.zeros((st.n_elements, uniq.size))
        np.add.at(M, (elem, np.array([gi[g] for g in grp])), 1.0)
        tot = M.sum(1); sc = tot > 0
        return float(np.mean(M[sc].max(1) / tot[sc])) if sc.any() else np.nan

    def evaluate(self, before, after):
        c0, c1 = self._coherence(before), self._coherence(after)
        return {"holdout_coherence_before": c0, "holdout_coherence_after": c1,
                "delta": (c1 - c0) if (np.isfinite(c0) and np.isfinite(c1)) else np.nan,
                "claim": ("TIER 2 internal validation: generalization to features "
                          "excluded from the refinement objective; NOT biological "
                          "accuracy (held-out features share the same assumptions).")}


# --------------------------------------------------------------------------- #
# OPTIONAL PLUG-IN — Tier 3: concordance with an external labelling.
# --------------------------------------------------------------------------- #
class ConcordanceValidator(Validator):
    name = "tier3_external_concordance"
    def __init__(self, external_labels):
        self.ext = np.asarray(external_labels)

    def evaluate(self, before, after):
        from sklearn.metrics import adjusted_rand_score
        ext = self.ext
        if before.element_labels is None or after.element_labels is None:
            return {"claim": "TIER 3: no element labels provided", "delta": np.nan}
        ok = np.arange(min(len(ext), before.n_elements))
        a0 = float(adjusted_rand_score(ext[ok], np.asarray(before.element_labels)[ok]))
        a1 = float(adjusted_rand_score(ext[ok], np.asarray(after.element_labels)[ok]))
        return {"concordance_before_ARI": a0, "concordance_after_ARI": a1,
                "delta": a1 - a0,
                "claim": ("TIER 3 external validation: concordance with an orthogonal "
                          "method; NOT accuracy (neither method is ground truth).")}


def make_prespecified_split(group_ids, holdout_fraction=0.5, seed=0):
    """Deterministic, seeded refine/holdout partition of group ids — fixed before
    analysis so the reserved set is never chosen from results."""
    g = sorted(set(map(str, group_ids)))
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(g))
    n_hold = max(1, int(round(holdout_fraction * len(g))))
    hold = {g[i] for i in perm[:n_hold]}
    refine = [x for x in g if x not in hold]
    holdout = [x for x in g if x in hold]
    return refine, holdout
