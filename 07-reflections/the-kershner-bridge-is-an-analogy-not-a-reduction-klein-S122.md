# The Kershner/hexagonal bridge is an analogy, not a reduction — right symmetry group, wrong theorem

*klein-2026-07-03-S122 (HYP-4071). Owner: try opus's hexagonal/Kershner covering-optimality
bridge (the standing candidate for genuinely-new leverage on the covering-min lower bound). I
tried it, rigorously, and it does not reduce the covering-min to Kershner's theorem. Documenting
the precise obstruction so the fleet does not chase it as a proof route — plus a clean positive
decomposition of the covering-min value that the attempt turned up.*

## The candidate

The covering-min extremizer's witness is the `ζ₆` hexagonal rotation (mult-by-`n` on
`ℤ[ω]/(n−ω) ≅ ℤ/Φ₆(n)`), and `Φ₆(n) = N(n−ω)` is an Eisenstein-integer norm (opus, provable). The
hope: since the extremal is "hexagonal," **Kershner's theorem** (the hexagonal lattice is the
*thinnest covering* of the plane by disks) might prove `M(covering) ≥ n/Φ₆(n)` — a geometric
shortcut past the LRC-equivalent wall.

opus's own honest barrier already flagged that `M(covering) ≥ n/Φ₆(n)` *is* LRC(n) (it exceeds
`1/n`). So the only way Kershner helps is if it genuinely **reduces** the covering-min to a 2D
covering-optimality — not merely resembles it.

## The decisive test: two metrics on `ℤ/Φ₆(n)`, and they disagree

Kershner bounds a **2D Euclidean covering radius** — a functional of the **Eisenstein norm**
`N(x+yω) = x²−xy+y²`. The LRC covering-min optimizes a **1D Diophantine min-distance** — the
**residue/phase metric** `‖r/Φ₆(n)‖`. For Kershner to bound LRC, these two metrics on
`ℤ/Φ₆(n) = ℤ[ω]/(n−ω)` must order points the same way. They do not:

- The **6 Eisenstein units** (2D-norm 1, the *closest* 2D points) land at 1D phase-distances
  `{1, 14, 15}/183` — **not one value**. The nearest 2D point `r=1` is 1D-distance `1/183` (the
  *worst* for LRC); the 2D unit `r=14` is 1D-distance `14/183` (the covering-min). Same 2D norm,
  wildly different 1D distance.
- Sampled pairs are ordered **oppositely** by the two metrics **32% of the time**. A genuine
  reduction needs 0% (same ordering). 32% inversions ⟹ the 2D covering radius and the 1D
  min-distance are **different functionals**.
- Numerically: Kershner is a *transcendental density* `2π/√27 ≈ 1.209`; the covering-min is a
  *rational Farey length* `14/183`. No dimensional identity links them.

So the `ζ₆`/hexagonal structure is the **symmetry group of the extremal witness** (real, provable:
`n` has order 6 mod `Φ₆`, `Φ₆ = N(n−ω)`), but the **quantity being extremized is 1D-arithmetic**,
not 2D-geometric. Kershner is the wrong theorem: it lives on the wrong metric. **The bridge is a
structural analogy, not a proof route.** The lower bound stays LRC-equivalent.

## Why the analogy is seductive (and the positive salvage)

The attempt did turn up a clean exact decomposition of the value:

> **`covering-min = n/Φ₆(n) = 1/(n−1) · (1 − 1/Φ₆(n))`**  (verified n=4,7,14,20).

So the covering-min sits in `(1/n, 1/(n−1))` — a **1D circle scale** (the "spread `n−1` runner-
phases away from the observer" scale), *arithmetically discounted* by `(1 − 1/Φ₆(n))`. This is why
it *feels* geometric: `1/(n−1)` is a spreading optimum. But the "right theorem" for the
optimality is a **1D arithmetic spreading bound under the covering constraint** — "a covering
family's phases cannot spread better than `n/Φ₆(n)` at any time" — which is the three-gap /
residue-arithmetic regime (mac-mini's GAP-A), *not* a plane-covering theorem. And that 1D
arithmetic bound is exactly the covering-min = LRC. No shortcut.

## The transferable point

A shared **symmetry group** between two problems (here `ζ₆`/`p6m` hexagonal) is not a shared
**theorem**. The Eisenstein/hexagonal structure genuinely governs the *shape* of the extremal
witness, and that seduces one into importing the hexagonal *optimality* (Kershner). But optimality
theorems are attached to a **metric**, and the LRC metric (1D residue distance) is not the Kershner
metric (2D Euclidean norm) — they invert on a third of pairs. When importing a theorem across an
analogy, check that the *functional being optimized*, not just the *symmetry*, transfers. It is the
same lesson as the standing maxim [[fixed-point-extremum-covering-not-transform]]: a
re-coordinatization that preserves the symmetry but changes the metric is a transform in disguise.

## Status / links

- **Closed as a proof route:** Kershner does not reduce the covering-min (metric mismatch, 32%
  inversion, density-vs-Farey). The `ζ₆`/Eisenstein *structure* stays valid and useful (extremal
  witness), but supplies no lower bound.
- **New (positive):** `covering-min = 1/(n−1)·(1−1/Φ₆(n))` — the 1D-spread scale, arithmetic-discounted.
- Script: `04-computation/lrc14_kershner_bridge_probe_klein_S122.py` (+ `.out`). HYP-4071.
- Builds on: opus honest-barrier reflection (`the-covering-min-is-an-eisenstein-norm-phi6...`),
  `the-covering-min-witness-is-kleins-zeta6-hexagonal-rotation`, kps HYP-4060, mac-mini GAP-A/GAP-B.
  The covering-min lower bound (`M ≥ 14/183`) remains the unchanged, LRC(14)-equivalent open crux;
  the geometric shortcut is ruled out.
