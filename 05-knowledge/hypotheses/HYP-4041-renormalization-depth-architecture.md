---
id: HYP-4041
title: THE RENORMALIZATION-DEPTH ARCHITECTURE for LRC(14) -- the three closers are COMPLEMENTARY and together cover all covering-compressed families. spread13 (global ratio<=13, PROVED) + bounded-denominator census (kps lonely14_of_ratio, for bounded-magnitude non-cluster families) + renormalization/fine-t (opus HYP-3901, for large near-equal clusters). KEY RECONCILIATION: the families that break kps's empirical q<=35 census bound are exactly the ALIGNED large-cluster families -- and those are LOOSE (M ~ 3x the danger radius, lonely at FINE t), so they belong to the renormalization route, NOT the census. The census's small-q failure on them is a PINNED-PHI artifact (small rationals force {Nt} into danger), not tightness. Recursion depth ~log(max-speed) = the discrepancy cost (arXiv:2607.00876 / opus HYP-4013 / my HYP-4040).
status: SYNTHESIS + numerically validated reconciliation; not a proof (the renormalization step = opus HYP-3901, still a program). The census-breakers-are-loose fact is VERIFIED (M=0.24=3.4x at N=1000, q_witness=41).
source: mac-mini-2026-07-03-S23
related:
  - HYP-4040   # my S22: census bound GROWS ~log(magnitude); no uniform small-q bound (the reason magnitude must be reduced)
  - HYP-3901   # opus deep-cluster renormalization (the magnitude-reduction step)
  - HYP-4013   # opus S48: star-discrepancy density floor + arXiv:2607.00876 bridge (the depth = discrepancy cost)
  - HYP-3984   # kps S28 bounded-denominator route / spread13 / lonely14_of_ratio (the base census)
  - MISTAKE-095 # my S21 band {15..33} -> {15..~50}; the aligned drifts are the census-breakers here
results:
  - 04-computation/reconcile_kps_denominator_bound_macmini_20260703.py
  - 04-computation/renorm_absorbs_cluster_macmini_20260703.py
  - 05-knowledge/results/renorm_absorbs_cluster_macmini_20260703.out
---

# HYP-4041 -- the renormalization-depth architecture

## The reconciliation that started it
kps-S28's bounded-denominator route says every hard covering-compressed LRC(14) family is lonely at some
`p/q` with `q <= 35`, independent of magnitude, and reduces LRC(14) to a bounded finite search. My HYP-4040
proves the witness denominator is UNBOUNDED (`q(S_X) > X` for lcm families, `Theta(log max-speed)`). These
looked contradictory. Resolving them gives the architecture.

The census bound `q <= 35` is **under-sampled** (the same MISTAKE-095 I made): ALIGNED covering-compressed
families -- far runners `q_i * round(N/q_i)` (near-equal `~N`, each `≡ 0 mod q_i`) -- reach witness `q = 44, 45`
at every magnitude (`73` such families with `q > 35` at `N = 1000`; `reconcile_kps_denominator_bound...out`).
So the census bound genuinely grows.

BUT these census-breakers are **LOOSE**, not tight: for the `N=1000`, `q_witness=41` family, `M = 0.243 =`
**`3.4x` the danger radius** `1/14`, achieved at a FINE `t* = 0.49975` (near `1/2`). Their small-`q` failure is
a **pinned-`phi` artifact**: at `t = a/q` (small `q`) the cluster phase `phi = {N a / q}` is pinned, and the
alignment parks it in danger; at a FINE `t`, `phi` is free and the near-equal cluster `{N + c_i}` (short arc
`{c_i t}`) slides into the safe region. So they are **renormalization families, not census families**.

## The mechanism (why large clusters are loose and fine-t-lonely)
For a far runner `N + c_i`, `||(N + c_i) t|| = ||phi + c_i t||` with `phi = {N t}`. As `t` ranges over a window
of width `1/N`, `phi` sweeps a full period while `c_i t` barely moves. So the cluster is `{phi + c_i t}` = a
rigid translate (by the free `phi`) of the short arc `{c_i t}` (length = spread`*t`). If that arc is shorter
than the safe region `[1/14, 13/14]` (length `6/7`), some `phi` places the whole cluster safely. The slow rest
`R` is placed by the coarse `t`; the fast cluster is absorbed by `phi`. This is the deep-cluster
renormalization ([[HYP-3901]]) seen at the level of a single witness.

## The architecture (three complementary closers)
For a covering family `S` (WLOG `gcd = 1`):
1. **Global ratio `<= 13`** => `spread13_lonely` (kps, PROVED): `t = 1/(min+max)` puts every runner in the full
   band `[1/14, 13/14]`. Includes the tight `{AP, GW}` locus at the boundary.
2. **A large near-equal top cluster** (max `>= B`, top-two within a small ratio) => **renormalize** it: absorb
   the cluster by the fast phase `phi`, reduce to the smaller family `R` (+ the difference core). Magnitude
   DROPS. Recurse. These families are LOOSE (verified `M ~ 3x`), so the reduction has margin.
3. **Bounded magnitude, no large cluster** (max `< B`) => **bounded-denominator census** (kps
   `lonely14_of_ratio`): a finite search over `q <= Q0(B)` (with the corrected `Q0 ~ 50`, not `35`).

Every covering-compressed family lands in one bucket. The census-breakers (aligned large clusters) go to (2),
not (3) -- so the census only needs to cover bounded-magnitude non-cluster families, where a fixed bound is
legitimate. The recursion in (2) has depth `~log(max-speed / B)` -- the **discrepancy cost** shared with
arXiv:2607.00876 (binary-tree continual counting) via opus's HYP-4013 star-discrepancy floor and my HYP-4040.

## What is proved vs open
- PROVED: spread13 (kps); the census-breakers-are-loose fact (this file, numeric, `M ~ 3x`); the no-uniform-
  census-bound lower bound (HYP-4040).
- OPEN (the one hinge): the renormalization step (2) -- "peel a large near-equal cluster, magnitude drops,
  loneliness preserved with margin" -- is opus's HYP-3901 program. Making it a theorem (a clean scale-
  separation lemma: `R` lonely with slack `delta` + cluster spread`*t < 6/7` => `S` lonely) closes the
  architecture. This is the highest-leverage next target, and it is where kps's census + opus's renormalization
  meet.

## Consequence for the proof plan
Do NOT build the census on a fixed `q <= 35` bound -- it is false for large clusters. Instead: peel large
clusters by renormalization first (they are loose, so this is safe), then run the census on the bounded base.
The proof is `spread13` + a renormalization descent + a bounded-magnitude census -- three finite/elementary
pieces, with the renormalization descent as the sole remaining analytic obligation.
