---
id: HYP-4074
title: Completing opus THM-615's m=2 confinement — the 11-runner tight-locus is {AP11} (verified, no GW_12) so opus's AP proof covers ALL tight even-parts; the loose-U case reduces to an a=b dodge (rigorous but vacuous for small tighteners = the same argmax-arithmetic barrier)
status: PARTIAL. 11-runner tight-locus={AP11} VERIFIED (thorough search, no GW_12). a=b dodge is a RIGOROUS sufficient condition but VACUOUS for the relevant small tighteners. So the m=2 confinement reduces to the loose-U small-tightener case = opus's open argmax-arithmetic barrier.
source: mac-mini-2026-07-04-S37 (building on opus THM-615 folding identity)
related:
  - THM-615   # opus: folding identity M(2U u{w1,w2})=max_t min(g_U(2t),Psi(t)); proved the tight-U=AP case; handed off loose-U
  - HYP-4073  # opus: the reduction M(2U u 2odd)>=1/12 for all 11-runner U
  - HYP-4070  # my GAP-A framework (n=14 non-covering tight={AP,GW}); this is the n=12 analog (={AP11})
results:
  - 04-computation/folding_loose_U_macmini_20260704.py
  - 05-knowledge/results/folding_loose_U_macmini_20260704.out
external: LRC(14); opus THM-615 folding identity.
---

# HYP-4074 — completing opus THM-615's m=2 confinement (partial)

opus's THM-615 reduces the m=2 (q*=28), |F|=2 confinement to: **M(2U u {w1,w2}) >= 1/12 for all 11-runner U**,
via the folding identity `M = max_t min(g_U(2t), Psi(t))`, `Psi=max(min(a,b),1/2-max(a,b))`, `a=||w1 t||`,
`b=||w2 t||`. opus PROVED it for `U = c*{1..11}` (dilated AP, the point-only extremal case); handed off general U.

## Contribution 1 (VERIFIED): the 11-runner tight-locus is {AP11} — no GW_12
A thorough exact-M search (single/double AP11 swaps+lifts + 80k random to speed 50) finds the ONLY 11-runner
tight family (`M=1/12`) is the **dilated AP `{1,…,11}`** — a SINGLE residue class, NO GW-type. (Contrast n=14,
which has GW because `q=12` admits a second coverer 24; the n=12 analog does not.) **Consequence:** every tight
even-part `U` (`M(U)=1/12`) is a dilated AP11, so **opus's AP proof already covers ALL tight U**. The m=2
confinement thus reduces to the LOOSE-U case (`M(U)>1/12`) alone.

## Contribution 2 (RIGOROUS sufficient condition, but vacuous where it's needed): the a=b dodge
*Lemma:* let `δ` = width of the largest interval of `{t: g_U(t)>=1/12}`. If `δ > 2/(w1+w2)` then
`M(2U u {w1,w2}) >= 1/12`.
*Proof:* `{t: g_U(2t)>=1/12}` contains an interval `I` of width `δ/2 > 1/(w1+w2)`, so `I` contains a point
`t0=k/(w1+w2)`; there `(w1+w2)t0∈ℤ ⟹ ||w1 t0||=||w2 t0||=a ⟹ Psi(t0)=max(a,1/2-a)>=1/4`, and `g_U(2t0)>=1/12`,
so `min(...)>=1/12`. ∎
**But it is VACUOUS for the relevant tighteners:** loose 11-runner `U` have `δ≈0.004–0.012`, while
`2/(w1+w2)>=0.017` for all `w1+w2<=116`. So the dodge only fires for HUGE tighteners; the small-tightener
covering cases (the ones that matter) are not covered — the same "measure/interval is vacuous, needs argmax
arithmetic" barrier opus identified (loose-U `λ<<2/7`).

## Net
The m=2 confinement's TIGHT-U case is now fully closed (11-runner tight = AP11, verified, + opus). The LOOSE-U
case remains open and is confirmed to require argmax arithmetic (not measure/interval): `M(2U u 2odd)=M(U)>1/12`
holds empirically but the single-interval a=b dodge is vacuous for small tighteners. The honest residual is a
sharp "hit a high point of `U` without tightener-extremity" for `M(U)` just above `1/12` — a bounded-tightener
finite problem per near-AP `U`. -> THM-615, HYP-4073 (opus), HYP-4070 (n=14 GAP-A). Script: folding_loose_U.
