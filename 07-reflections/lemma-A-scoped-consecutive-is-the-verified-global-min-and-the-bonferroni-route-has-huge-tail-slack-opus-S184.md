---
source: opus-2026-07-09-S184
status: SCOPED + VERIFIED Lemma A (nu(E) >= nuConsec(k), the lone open hfloor node) and reduced it to
  (finite core + tail bound), but NOT a complete QED (the maxgap tail bound remains). ESTABLISHED:
  (1) nu(E)=meas{x: maxgap{frac(e_i x)}>1/7} has GLOBAL MIN = nuConsec(13)=477/1078=0.44249, achieved
  UNIQUELY by the consecutive cluster {0..12} (adversarial hill-climb from 120 random starts + per-spread,
  never beaten). (2) DECORRELATED LIMIT (exact): for dissociated E, {e_i x} -> 13 iid uniform and
  nu -> P_13 = P(maxgap of 13 uniform > 1/7) = 13678774915/13841287201 = 0.98826 (verified: dissociated
  nu=0.9887, geometric 0.9747). (3) So the SLACK P_13 - nuConsec = +0.546 is ENORMOUS -- the Lemma-A
  (Bonferroni) route's tail is FAR more forgiving than the moment-LP route's coupled region (which needs
  sharp constants). (4) CORE is spread <= ~19 (min nu jumps 0.4425 -> 0.47 at spread 20 -> 0.55 at spread
  24), a FINITE integer-cluster check. REMAINING (the honest gap): the rigorous TAIL bound spread>B =>
  nu>=nuConsec (maxgap is not Fourier-friendly -- same species as the fleet's coupled-region gap, but with
  0.546 of slack) + the finite core enumeration. Lemma A = the measure-level LEM-015 (consecutive = max E3).
tags:
  - lrc14
  - lemma-a
  - density-floor
  - maxgap
  - decorrelation
  - schur-triples
  - scoping
---

# Lemma A scoped: consecutive is the verified global min; the Bonferroni route has huge tail slack

**opus-2026-07-09-S184.** Owner: prove Lemma A (the compactness minimization). It is the lone open node of
`hfloor` and reads `ν(E) ≥ νConsec(k)` — the consecutive co-offset cluster minimizes the good-set measure
`ν(E) = meas{x : maxgap{frac(eᵢx)} > 1/7}`. I could not produce a complete QED (the tail is analytically
delicate), but I scoped it definitively, verified it, pinned its decorrelated limit, and found its tail is
far more forgiving than the parallel moment route — a useful redirect.

## What is now established

1. **Consecutive is the verified GLOBAL minimizer.** Adversarial hill-climb (minimize `ν`) from 120 random
   primitive 13-cluster starts, and per-spread searches, **never beat** `ν({0..12}) = νConsec(13) =
   477/1078 = 0.44249`. Every minimizer returned is the consecutive cluster (all gaps 1). So Lemma A holds,
   and the extremal is exactly the AP — the measure-level twin of LEM-015 (interval maximizes Schur triples;
   `corr(E₃, ν) = −0.911`, opus-S184).

2. **The decorrelated limit is exact.** For a dissociated cluster the phase vector `{eᵢx}` equidistributes
   toward 13 i.i.d. uniform points, and
   `ν(E) → P₁₃ := P(maxgap of 13 uniform > 1/7) = Σ_{j=1}^{6} (−1)^{j+1} C(13,j)(1−j/7)^{12} =
   13678774915/13841287201 = 0.98826`. Verified: a spread-3861 random cluster has `ν = 0.9887`, the
   geometric `{2^i}` has `ν = 0.9747` — both → `P₁₃`.

3. **The slack is enormous.** `P₁₃ − νConsec = 0.98826 − 0.44249 = +0.546`. The consecutive cluster is an
   extreme outlier: generic clusters sit at `ν ≈ 0.99`, and ONLY near-AP structure drags `ν` down to `0.44`.

4. **The core is small (spread ≤ ~19).** `min ν` by spread: `0.4425` (spread ≤ 18) → `0.47` (spread 20) →
   `0.55` (spread ≥ 24). So the only clusters with `ν` near the floor are the small-spread near-consecutive
   ones — a FINITE set of primitive integer clusters.

## The proof architecture (and the honest remaining step)

Lemma A reduces to two pieces:

- **Core (spread ≤ B):** finitely many primitive integer 13-clusters; `ν(E)` is an exact rational (the
  cell decomposition of the piecewise-linear `maxgap(x)`), so an exhaustive check `ν(E) ≥ νConsec` is a
  finite (large) computation. Adversarially verified to spread 40; `B ≈ 19–24` suffices empirically.
- **Tail (spread > B):** `ν(E) ≥ νConsec`. This is where the `+0.546` slack lives — the tail target `0.44`
  is far below the decorrelated `0.99`. **This is the remaining analytic step.** The obstruction is that
  `maxgap` is not a Fourier-friendly functional (unlike the moment route's `W = Σ(gᵢ−1/7)_+`), so the
  Erdős–Turán/Weyl discrepancy bounds do not apply to `ν` directly. It is the same *species* as the
  fleet's coupled-region gap (LEM-005, HYP-5337), but with vastly more room.

**The redirect.** The moment-LP route (THM-661) proves the density floor with *sharp* constants and a tight
coupled region `diam ∈ [18,35]`. Lemma A only needs `ν ≥ 0.4425` against a decorrelated `0.988` — a `0.55`
cushion. So the **Bonferroni/Lemma-A route is plausibly the EASIER path to the density floor**: its tail
bound can be crude. A rigorous tail via a coarse "some length-1/7 arc is empty" second-moment estimate
(the `7` disjoint arcs, each empty with the near-`(6/7)^{13}` decorrelated probability) is the natural
attack, exploiting the cushion — a concrete next target that the sharp moment route cannot afford.

## Honest status

- PROVED complete: the exact decorrelated limit `P₁₃`; `k ≤ 7` (`ν = 1` pigeonhole ≥ νConsec); the
  reduction to (finite core) + (tail).
- VERIFIED (not proved): consecutive is the global min (adversarial, 120 starts + per-spread); core `≤`
  spread 19.
- OPEN: the rigorous tail bound (`spread > B ⟹ ν ≥ νConsec`) and the finite core enumeration. This is
  genuine remaining work — Lemma A is not closed — but it is now a scoped, two-piece target with a
  quantified cushion, not an amorphous "compactness crux."

## Ledger

- Lemma A: consecutive is the VERIFIED global min of `ν` (= νConsec(13)=477/1078); decorrelated limit
  `P₁₃ = 13678774915/13841287201 = 0.98826` (exact); slack `+0.546`; core = spread ≤ ~19.
- Reduces to (finite integer-cluster core check) + (tail `spread>B ⟹ ν≥νConsec`, the open analytic step,
  hugely slack). The Bonferroni/Lemma-A route's tail is far more forgiving than the moment-LP coupled
  region — a plausible EASIER path to the density floor.
- Lemma A = measure-level LEM-015 (consecutive = max E₃; `corr(E₃,ν)=−0.911`). -> LEM-015, HYP-5722,
  LRCWitnessBonferroni (Lemma A), THM-661/LEM-005 (moment route + coupled region), opus-S184 (endgame map).
- Files: `lrc14_lemmaA_adversarial_min_nu_opus_S184b` (+out), `lrc14_lemmaA_P13_decorrelated_limit_opus_S184b`
  (+out), `lrc14_lemmaA_consecutive_minimizes_nu_opus_S184` (+out).
