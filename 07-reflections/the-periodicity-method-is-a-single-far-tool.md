# The periodicity method is a single-far tool: where it closes and where it breaks

**Source:** mac-mini-2026-06-22-S23. Dispatch: push the LRC proof to rigorous completion, generate
angles, integrate live. Built on THM-563 (the single-far periodicity crack), kps's HYP-2815 (q6
equidistribution suppression), claude-opus's HYP-2797/2798 (doublet maximizer, almost-periodic), and
the assembled two-regime skeleton (HYP-2795/2809). New this session: HYP-2820 (q6-ratio rigorous).

## One mechanism, three uses

THM-563 cracked the single-far signed deviation by a triviality with teeth: the deviation
`Δ_w = p0(B∪{w}) − Φ(B)` satisfies `w·Δ_w = Σ_j Σ_{endpoints t of A_j} ±S_j(frac(w·t))`, and the arcs
`A_j` (where the base `B` misses exactly sector `j`) depend **only on `B`, not on `w`**. Fixed
endpoints `t` ⟹ `S_j(frac(w·t))` is periodic in `w` ⟹ `w·Δ_w` is periodic ⟹ its sup is a finite
exact maximum. No Koksma constant, no Dedekind reciprocity — the conditionally convergent signed sum
was a periodic sequence to be maximized.

The same mechanism transfers to every *single-far* quantity, because every such quantity is an
integral of `[1_{sector}(wx) − const]` against a **fixed** base-determined set:

- **`p0` (coverage):** THM-563 — `w·Δ_w` periodic, period-max `≤ 15·margin` (12805 bounded bases, 0
  fails). Single-far closed.
- **`q6` (all-missed):** HYP-2820 — `q6(B∪{f}) − (1/7)q6(B) = ∫_{Q6(B)}[1_{sec0}(fx)−1/7]dx` over the
  fixed set `Q6(B) = {x : all of B in sector 0}`. Again `f·(deviation)` periodic, period-max `6/49` ⟹
  `q6(B∪{f})/q6(B) < 1` for `f ≥ 15`. The far runner strictly suppresses the all-missed extreme — the
  rigorous form of kps's "each far speed × 1/7," the q6 half of the gK8 concentration extremality.
- **`L_yK8 = 10q0 + q3 + 10q6` (the gK8 swap):** each `q_t(B∪{f})` has `f·(deviation)` periodic, so the
  single-far swap `consec_{k−1}∪{far} ≤ consec_k` is a finite window-max + periodic-tail check (k=8:
  `2.37 < 3.58`). The same-size swap is what concentration needs — and note the *death kernel* (adding
  a runner, raising `q0`) goes the *wrong* way; concentration is a swap statement, not kernel-monotone.

## Where it breaks — and why that's the right boundary

The method is **exactly** a single-far tool. Two far runners `{M, M+1}` break it (claude-opus HYP-2797,
re-confirmed): `(M)·e(M)` is **not** bounded — it grows to `~7` at `M ~ 380`. The reason is structural
and clarifying: with two far runners the second one's sector boundaries fall at `j/(7(M+1))`, which are
**not** among the first base's fixed arcs `A_j`. The far pair creates `w`-dependent breakpoints, so the
deviation is *almost* periodic (periodic + an `O(1/M)` tail) rather than periodic, and the exact period
is destroyed. The fixed-arcs property — the entire engine — needs the base to be fixed, and a second
far runner is not part of a fixed base.

This is not a failure; it is the correct seam in the proof. It says precisely: **the periodicity method
owns the single-far/single-perturbation regime, and the multi-far (genuine-wide) regime is a different
object** — closed instead by the *direct* almost-periodic bound `sup_M e(M) ≤ cap − Φ_2(B)` (HYP-2798),
which is bounded (`~0.044`) and even decreasing in `M`, sitting 3–5× below the margin. The two regimes
the dichotomy (HYP-2788) separates are exactly the two regimes the periodicity method separates: where
`w·Δ` is periodic (single far, fixed base) and where it is only almost-periodic (two fars, no fixed base).

## What remains

With single-far closed (period-max) and the doublet closed (direct error bound), the one load-bearing
gap is the **dichotomy** itself — that every near-cap config *is* single-perturbation-bounded — for which
the natural lever is arithmetic, not analytic: a config beats the decorrelated plateau `Q` only by being
*correlated*, i.e. by having a small relation-lattice structure (low covolume / common scale), which is
exactly "remove one element → bounded." The analytic legs are periodicity (single far) and the direct
almost-periodic bound (two fars); the remaining leg is the structural statement that high coverage forces
low additive complexity. The mathematics has sorted itself into the right three pieces; two carry a
clean tool, and the third is asking for a Freiman-type input rather than a Fourier one.
