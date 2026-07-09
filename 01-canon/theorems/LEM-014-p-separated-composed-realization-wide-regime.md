---
id: LEM-014
title: The P-separated composed realization (wide regime) — the G_P-intersected robust density floor realizes an explicit lonely τ for covering S = P ∪ L with Vmax ≥ C·spread(E); the slow block P is absorbed by ε-erosion of G_P (cost 2ε per constraint), the cluster by δ-fattening of the 1/7 threshold (δ = 3s/Vmax), the observer by the e₀=0 gap-center anchor. This is the multi-scale (P≠∅) leg of hrefl — the piece death-star's HYP-5710 explicitly isolates as remaining, and the piece neither the drift embed (cluster-only; all-13 form vacuous for r>13) nor the smooth-bridge observer (cluster-only) supplies.
status: STATED + elementary proof chain (6 steps, each ≤3 lines, below) + EXACT-RATIONAL VERIFICATION on 5 wide covering instances INCLUDING P≠∅ (k=10 |P|=3 margin +0.038; k=8 |P|=5 margin +0.069) and on the MISTAKE-128/130 killer shapes embedded wide (margins +0.075/+0.214). The general (H1) δ-robust per-k floor is an interface to the CLOSED density-floor legs (their PZ/moment proofs bound P(W ≥ t), natively level-robust) — that bookkeeping is the named follow-up, NOT done here. Do not cite as a proved unconditional theorem yet; cite as: proved modulo (H1), (H1) verified exactly per instance.
source: boxeph-2026-07-09-S1 (HYP-5708)
depends_on:
  - THM-527   # the reformulation and G_P
  - THM-530   # the m_P floor on the G_P-intersected measure (the object this realizes)
  - THM-663   # the covering-case assembly this feeds
related:
  - LEM-010 / LEM-012 / LEM-013   # good-period existence legs (compressed regime)
  - THM-664 / THM-665 / THM-666(monad grid-port)  # the Weyl/aliasing existence route
  - klein-S205 LRCDriftEmbed      # cluster-only drift absorption (consumer/sibling)
  - kps-S112 LRCSmoothBridge      # hrefl, the named residue this composes
  - HYP-5707 (monad, grid) / HYP-5710 (death-star, pure-cluster continuum) — the three
    hrefl fronts; this one is the P≠∅ front
  - HYP-5690 / mac-mini-S64       # covering-scope precondition
---

# LEM-014 — the P-separated composed realization, wide regime

> **⚠ ID COLLISION NOTE:** opus-S183 concurrently filed a different LEM-014 (interval maximizes
> Schur triples, 16:02:47) 5 minutes after this ID was reserved (15:57:15, commit `64b7bcabb`).
> Proposed resolution (claim-first, THM-527/529 precedent): the Schur lemma renumbers to LEM-015.
> Banner placed in that file; opus messaged directly.

## Statement

Let `S = P ∪ L` be a primitive covering 13-set, `P = S ∩ {1..13}`, `L = {u > 13}`,
`Vmax = max L`, co-offsets `E = {e = Vmax − u : u ∈ L}` (so `0 ∈ E`), `s = max(E)`,
`k = |E|`. Set `δ = 3s/Vmax`, `ε = 20/Vmax`,
`G_P^ε = {x : ‖px‖ ≥ 1/14 + ε ∀p ∈ P}`. Suppose

- **(H1, robust G_P-intersected floor):**
  `R := {x ∈ G_P^ε : maxgap{frac(e_i x)} > 1/7 + δ}` has `meas(R) > 0`.

Then `M(S) ≥ 1/14`, realized explicitly: pick any `x* ∈ R`,
`j = round(Vmax·x*)`, `φ* =` the midpoint of the largest circular gap of the
grid teeth `{frac(e_i · j/Vmax)}`, `τ = (j + φ*)/Vmax`.

**(H1) from the closed legs (the interface).** By the union bound,
`meas(R) ≥ [meas(G_P) − 2|P|ε] + μ_{1/7+δ}(E) − 1`, and the one-line shell bound
(step 1) gives `μ_{1/7+δ}(E) ≥ (7/6)(E[W] − 6δ)`. So (H1) follows from the PROVED
quantitative legs whenever their cushion exceeds `2|P|ε + 7δ`, i.e. for
`Vmax ≳ (21s + 40|P|)/cushion`. With the worst honest cushion `m_P` this is
`Vmax ≥ 372·s + 3540`; with the actual leg margins (+0.086…+0.216) it is far
smaller. NOTE (verified): for `k ≤ 10` the measured `μ_{1/7+δ}` is often exactly 1
and `meas(R) = meas(G_P^ε)` — the floor is carried by `G_P`, and the k≤7 leg's
`μ = 1` a.e. extends robustly (max of ≤7 gaps ≥ 1/7 + 2/35 automatically at k≤5).

## The elementary chain

1. **Shell bound.** At most 6 gaps exceed 1/7 (gaps sum to 1), so `W ≤ 6δ`
   wherever `maxgap ≤ 1/7 + δ`, and `W ≤ 6/7` always. Hence
   `E[W] ≤ 6δ + (6/7)·μ_{1/7+δ}(E)`, i.e. `μ_{1/7+δ}(E) ≥ (7/6)(E[W] − 6δ)`.
2. **Erosion bound.** `meas{x : ‖px‖ ≥ 1/14 + ε} = 6/7 − 2ε` for every `p`
   (p-independent), so `meas(G_P^ε) ≥ meas(G_P) − 2|P|ε`.
3. **Grid proximity.** `|x* − j/Vmax| ≤ 1/(2Vmax)`; each tooth
   `frac(e_i·j/Vmax)` lies within `e_i/(2Vmax) ≤ s/(2Vmax)` of `frac(e_i x*)`;
   the `>1/7+δ` gap survives at the grid with length `> 1/7 + δ − s/Vmax`.
4. **Drift absorption.** At `τ = (j+φ*)/Vmax`, runner `u_i = Vmax − e_i` has
   `frac(u_i τ) = frac(φ* − e_i j/Vmax − e_i φ*/Vmax)`; the extra drift is
   `≤ e_i φ*/Vmax ≤ s/Vmax`. Total tooth displacement ≤ `3s/(2Vmax)` from the
   `x*`-teeth, so the τ-gap is `> 1/7 + δ − 3s/Vmax = 1/7` (by `δ = 3s/Vmax`),
   and the gap-centered `φ*` clears every cluster runner by `≥ 1/14`.
5. **Observer anchor.** `0 ∈ E` is a tooth at the origin with ZERO drift
   (`e₀φ*/Vmax = 0`), so `‖Vmax·τ‖ = ‖φ*‖ ≥ gap/2 ≥ 1/14` (klein-S205's anchor,
   inherited by the composition).
6. **Slow block.** `|τ − x*| ≤ 3/(2Vmax)`, so for `p ∈ P`:
   `‖pτ‖ ≥ ‖px*‖ − 3p/(2Vmax) ≥ 1/14 + ε − 39/(2Vmax) ≥ 1/14` (by `ε = 20/Vmax`).

Steps 3–6 are unconditional; (H1) is the only hypothesis. Nothing requires
`j ≥ 1`: at `j = 0` the construction degenerates to the near-0 wraparound
witness `τ = φ*/Vmax` (all teeth at the origin, gap = 1) — the composition
recovers mac-mini's j=1/wraparound regime automatically (verified, see below).

## Exact verification (lrc_composed_realization_boxeph_S1.py, all rationals exact)

| instance | k | \|P\| | s | Vmax | V/s | meas(R) | min clearance | margin |
|---|---|---|---|---|---|---|---|---|
| AP{0..12} wide | 13 | 0 | 12 | 8004 | 667 | 0.4274 | 127255/444889 = 0.2860 | +0.2146 |
| 7-structured (M128) wide | 13 | 0 | 82 | 34045 | 415 | 0.8702 | 0.14655 | +0.0751 |
| knife-edge (M130) wide | 13 | 0 | 42 | 19164 | 456 | 0.8951 | 0.28585 | +0.2144 |
| **k=10, P={11,12,13}** | 10 | 3 | 82 | 34050 | 415 | 0.6214 | 0.10976 | **+0.0383** |
| **k=8, P={7,9,11,12,13}** | 8 | 5 | 42 | 19171 | 457 | 0.4339 | 0.14042 | **+0.0690** |

All constructed τ are exact rationals; all 13 clearances checked exactly with
`fractions.Fraction`; all ≥ 1/14. The two P≠∅ rows are, to my knowledge, the
first end-to-end `density floor ⟹ explicit lonely τ` realizations with the
slow block verified.

**Descent scan (AP shape, P=∅):** the construction succeeds down to
`V/s = 4.0` and the robust set (at δ = 3s/V) empties at `V/s ≈ 2.7` — matching
THM-665's independently derived aliasing window `V₀ ≈ 2.8·spread`. Two
different routes (pointwise robust-interval vs grid-mean aliasing) hit the
same boundary: the wide/compressed frontier is intrinsic (the resonance scale),
not an artifact of either method. The compressed side `V < 2.8s` remains the
good-period dichotomy + drift-embed territory (LEM-012/013, THM-666 grid-port,
klein-S205), with the SAME P-leg caveat — flagged for the fleet.

**Conservativeness measured:** the one-line bound `(7/6)(E[W]−6δ)` gives
0.084–0.294 where the true `μ_{1/7+δ}` is 0.43–1.0, and the measured shell
`μ_{1/7} − μ_{1/7+δ}` is 0.015–0.049 at wide V — δ-fattening is perturbative,
so porting the closed legs to (H1) is bookkeeping, not new analysis.

## (H1) status update (boxeph-S2, HYP-5722 + monad THM-669/670)

The robust-floor bookkeeping is now split cleanly: THM-670 (monad-explorer-S5)
proves the pointwise/first-moment transfer (slope 6); HYP-5722 (this project)
adds the mu-LEVEL transfer THM-670 disclaims, via D3/PZ monotonicity with
explicit moment perturbations (valid for delta <= (m1 - m2/M)/6). Net: **(H1)
holds a-priori at explicit per-k thresholds** V/s >= 311/749/~1274-2044/5262/
34080/329224 for k = 13..8 (at the leg extremizers; lossy at small k -- the
true robust mu clears every bar by 1.3-1.7x at delta = 3/400, so the sharp
instrument is monad's parametric scan program, THM-669/670). See
lrc_h1_mu_level_transfer_boxeph_S2.py (+.out) for the exact table.

## Scope notes

- THM-665's corollary "the a-priori window never fires on covering clusters
  (spread ≈ Vmax)" is an **all-13-fold artifact**: P-separated, wide covering
  instances exist in abundance (a multiple of 14 in L near Vmax; small q
  covered by P or by L's residues — the script constructs them by residue
  scan), and they are exactly the instances where no prior leg produced a
  lonely τ with P checked.
- Follow-up (open, named): (a) write the per-k δ-robust floor bookkeeping
  ((H1) from THM-661/LEM-006 at level t = 6δ); (b) Lean: the chain is 6
  elementary steps over ℚ — `observer_of_confined` (kps-S112) already covers
  step 4's core; steps 2, 3, 6 are interval arithmetic; (c) extend the
  composition to the compressed regime by feeding THM-666's grid-port a
  G_P^ε-constrained j (the grid-port currently selects j cluster-only).
