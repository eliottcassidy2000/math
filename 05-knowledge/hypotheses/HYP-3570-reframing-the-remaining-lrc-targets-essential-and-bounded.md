---
id: HYP-3570
title: The remaining LRC(14) proof targets, reframed under the obstruction/relational/moment lens, collapse into ONE statement -- the danger relation's self-composition D.D^T (= the 2nd moment = the R-equivariant obstruction class) is ESSENTIAL (non-coboundary; the sigma-odd/witness/counting half) AND BOUNDED (CV(N_R)^2 < threshold; the sigma-even/floor/Lebesgue half) under the Gamma_0(N) change-of-base; the four scattered obligations (A witness-density, B witness-floor, C gK8-concentration, D doublet-R-tail) and the gatekeeper (R'>0) are facets of these two halves, and the metagraph is the PROVED template (CV(H)^2~2/n clean; b_1^-/b_1->1/2 robust; SC=trace(R)>0 essential)
status: SYNTHESIS / reorganization of the proof targets (not a new proof). Grounded in proved metagraph facts (THM-589, HYP-3565, THM-587) and the existing obligation structure (THM-527/579, OPEN-Q-108). The concrete next step (Gamma_0(N) 2nd moment) is named, not done.
source: mac-mini-2026-06-29-S26
related:
  - HYP-3564  # relations not things: the proof target = the danger relation is essential (non-coboundary)
  - HYP-3562  # the measure of the obstruction: essential = topological = set-independent
  - HYP-3553  # the Gamma_0(N) congruence change-of-base (the clean reference)
  - THM-589   # the metagraph rehearsal: CV(H)^2~2/n (the BOUNDED template, proved)
  - HYP-3565  # b_1^-/b_1 -> 1/2 (the ESSENTIAL template: the obstruction is robustly half)
  - THM-579   # the gatekeeper CV(N_R)^2 < m_Q/(1-m_Q)
  - THM-527   # the witness-density obligation (A)
  - OPEN-Q-108 # the covering floor / critical path
---

# HYP-3570 -- the remaining LRC targets: ESSENTIAL and BOUNDED

## The one statement
All remaining LRC(14) work is two properties of ONE object -- the danger relation `D(v,t)=[||vt||<1/14]`
composed with itself, `D.D^T`, which is simultaneously the **second moment**, the **R-equivariant
obstruction class**, and the **pair-correlation** (the only invariant, klein THM-588):
> **`D.D^T` is ESSENTIAL (non-coboundary) AND BOUNDED (`CV(N_R)^2 < m_Q/(1-m_Q)`) under the `Gamma_0(N)`
> change-of-base.**
ESSENTIAL gives existence (the floor is `>0`, a lonely point exists); BOUNDED gives the uniform gatekeeper
(THM-579). These are the `sigma`-odd and `sigma`-even halves (S23's two measures).

## The four obligations are facets of the two halves
| old target | reframed as | half |
|---|---|---|
| **A** THM-527 witness density `G2>0 => M>=1/14` | the obstruction is NONZERO (counting measure, `phi(n)/2`) | ESSENTIAL (`sigma`-odd) |
| **B** witness floor `k=8..13` | the counting measure survives where Lebesgue vanishes (the extremal) | ESSENTIAL (`sigma`-odd) |
| **C** gK8 concentration `max_E <= 10 cap` | the 2nd moment controls the max (Chen-Stein/Poisson) | BOUNDED (`sigma`-even) |
| **D** doublet R-tail `R(M)=O(1/M)` | the pair-correlation TAIL of `D.D^T` | BOUNDED (`sigma`-even) |
| **floor** `R'>0` / `CV(N_R)^2` | ESSENTIAL x BOUNDED (existence + concentration) | both |
So the witness side (A,B) is the obstruction's ESSENTIALITY (it does not factor; the bilinear `vt` forbids
separation, HYP-3564); the floor/cap side (C,D,gatekeeper) is its BOUNDEDNESS (the relation's
self-composition is finite-variance under the right base).

## Why the reframe changes the attack
- **ESSENTIAL is topological, not analytic.** The witness-density obligation A ("the linchpin") is, reframed,
  the non-vanishing of an equivariant Euler class -- a SET-INDEPENDENT topological invariant (S23), not a
  per-config inequality. Certified by non-separability of `D` (the bilinear `vt`, HYP-3564) and the Borsuk-
  Ulam counting measure (`phi(n)/2`), which survives at the extremal where the floor measure vanishes (S23).
  This dissolves the "prove `G2>0` for every config" framing into "the class is essential."
- **BOUNDED is a change-of-base, not a sharper estimate.** klein-S4: `CV(N_R)^2` is set-dependent and
  UNBOUNDED per-set. klein-S5: that is the symptom of the `Z_14`-collapse (degenerate at cusps), vs the
  metagraph's `S_n`-collapse (clean). So the move is NOT to estimate harder but to CHANGE THE BASE: the
  `Gamma_0(N)` congruence (HYP-3553) makes the 2nd moment depend only on `N` (`psi/phi/J2`), set-
  independent -- "the floor manufactures the transitive symmetry it lacks" (klein-S5). Obligations C and D
  are then bounded by the SAME congruence second moment (Han-Lee), not four separate analyses.

## The metagraph is the PROVED template (the rehearsal)
Every half has a proved finite mirror on the metagraph:
- BOUNDED: `CV(H)^2 = W(n)/n! - 1 ~ 2/n` (THM-589) -- the second moment IS bounded under the `S_n`-collapse.
  The LRC target is to reproduce this under `Gamma_0(N)`.
- ESSENTIAL: `b_1^-/b_1 -> 1/2` (HYP-3565) -- the R-odd obstruction is HALF the cycle space, robustly
  nonzero; and `SC = trace(R) > 0` (THM-587, Lefschetz) forces SC tournaments to exist WITHOUT construction.
  The LRC target is the same essentiality for the danger-relation obstruction.

## The concrete next step
Compute the `Gamma_0(14)` congruence second moment (Han-Lee, arXiv:2507.05905) and show it bounds
`CV(N_R)^2 < cap_r/(1-cap_r)` SET-INDEPENDENTLY (depending only on `N=14` via `psi(14)=24, phi(14)=6,
J2(14)=144`). That single computation discharges the gatekeeper and obligations C/D at once; the witness
side A/B is then the essentiality (the class is non-zero), which the non-separability + counting measure
already certify. The proof becomes: ESSENTIAL (topological, done in spirit) x BOUNDED (one congruence
2nd moment), not five separate analytic estimates.

## What it buys
A reorganization that replaces five scattered analytic obligations with two structural properties of one
relation, each with a PROVED metagraph template and a single concrete next computation (`Gamma_0(14)`).
The witness/floor split becomes the essential/bounded (counting/Lebesgue) split of the obstruction; the
covering becomes the change-of-base; the proof becomes "the danger relation does not factor and has
bounded self-composition under the right symmetry."
