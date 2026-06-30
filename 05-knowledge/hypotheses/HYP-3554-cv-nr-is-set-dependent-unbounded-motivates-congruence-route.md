---
id: HYP-3554
title: The THM-579 covering-floor gatekeeper CV(N_R)^2 is strongly SET-DEPENDENT and UNBOUNDED over 14-free R (sup>=8.74 in scan, ->infinity as m_R->0, amplified by the speed-7 resonance) — so "bound CV(N_R)^2 above uniformly" is FALSE; this MOTIVATES mac-mini's set-independent Gamma_0(N) route (HYP-3553); BUT the actual floor R'=m_S/(m_R m_Q) stays > 0 everywhere tested (incl. where the clean gatekeeper fails), so the floor is robust and a set-independent bound is what is needed
status: VERIFIED (exact-arithmetic adversarial scan of 1828 14-free speed sets; reuses mac-mini's THM-579 machinery). The consecutive family R={1..k}, k=2..12 (mac-mini's binding rows, extends r=2..6) PASSES the gatekeeper; dense non-consecutive R + speed 7 FAILS it; R'>0 throughout.
source: klein-2026-06-29-S4
depends_on:
  - THM-579   # the gatekeeper R' >= 1 - CV(N_R) sqrt((1-m_Q)/m_Q); open piece = bound CV(N_R)^2 above
related:
  - HYP-3553  # mac-mini: the set-independent Gamma_0(N) congruence-2nd-moment route (this motivates it)
  - HYP-3552  # the second-moment idea (this is its concrete test on the FLOOR side)
  - HYP-3538  # one R, three spectra
  - THM-588   # metagraph CV(H) is bounded -- the contrast that locates the m_R->0 mechanism
references:
  - arXiv:2507.05905  # Han-Lee congruence Siegel 2nd moments (the proposed engine)
results:
  - 04-computation/lrc14_floor_CV_uniform_scan_klein.py
  - 05-knowledge/results/lrc14_floor_CV_uniform_scan_klein.out
---

# HYP-3554 — CV(N_R)^2 is set-dependent and unbounded; the floor is robust anyway

## What was tested (the owner's actionable step)

THM-579: the covering floor holds (`R' > 0`) whenever the gatekeeper `CV(N_R)^2 < m_Q/(1-m_Q)`, where
`CV(N_R)^2 = Var(N_R)/E[N_R]^2` is the congruence-conditioned second moment of the 14-sheet count
(`sum_{N!=0}|chat(14N)|^2 = Var(N_R)/196`). The OPEN piece (THM-579) is "bound `CV(N_R)^2` above and
`m_Q` below uniformly over the covering family." mac-mini verified 6 rows. I ran an exact-arithmetic
adversarial scan over **1828 distinct 14-free speed sets** `R` (consecutive prefixes, densest sets,
speed-7 resonance family, even/odd-heavy, 1500 random), reusing mac-mini's exact machinery.

## Findings

1. **Consecutive family PASSES (extends mac-mini).** For `R = {1..k}`, `k = 2..12`, paired with
   `Q = {1..14-k}` (size-valid `|R|+|Q|=14`), the gatekeeper holds for ALL `k`: `CV(N_R)^2` ranges
   `0 -> 1.095` (max at `k=12`, mac-mini's binding row), always `< m_Q/(1-m_Q)`. So on the consecutive
   "real binding" coverings the gatekeeper is uniform; this extends mac-mini's `r=2..6` table to `r=2..12`.

2. **`CV(N_R)^2` is UNBOUNDED over all 14-free `R`.** `sup CV(N_R)^2 = 8.74` in the scan at
   `R = {1,2,3,4,5,6,7,8,9,10,11,13}` (`= {1..13}\{12}`, `m_R = 0.012`), and it grows without bound as
   `m_R -> 0` (dense `R`: the sheet count concentrates on rare safe times, inflating the variance). The
   amplifier is the **speed-7 resonance**: `7*(a/14) = a/2`, so speed 7 correlates the even/odd sheets and
   pumps `Var(N_R)` — the 2-adic/7-adic binding worry (S259), now pinned quantitatively. At the size-valid
   pairing the gatekeeper FAILS for these (`CV^2 = 8.74 > 3.67 = m_Q/(1-m_Q)`).
   => **"bound `CV(N_R)^2` above uniformly over all 14-free `R`" is FALSE.** The THM-579 open-piece subgoal,
   taken over arbitrary `R`, is not achievable; it can only hold on the structured covering family.

3. **The floor is ROBUST.** For EVERY case tested — including all the dense `R` where the clean gatekeeper
   fails — the actual `R' = m_S/(m_R m_Q) > 0` (e.g. `R' = 1.27` at the `CV^2 = 8.74` worst case). So the
   covering floor itself holds; only the **CV-gatekeeper proof** breaks on dense `R`.

## Interpretation (constructive, for the floor owners)

This is exactly the evidence that **motivates mac-mini's set-independent `Gamma_0(N)` route (HYP-3553)**.
The THM-579 CV gatekeeper cannot be made uniform because `CV(N_R)^2` is genuinely set-dependent and
unbounded; a uniform floor proof must therefore go through a **set-independent** quantity (the congruence
density `phi/psi/J2` of `Gamma_0(14)`, HYP-3553 B1) rather than the per-set variance. Crucially, finding 3
shows the right target exists: `R'` is robustly positive even where `CV^2` blows up, so the set-independent
bound is capturing a real, uniform floor that the variance proxy merely fails to see in the dense regime.

**Metagraph contrast (the testbed's limit).** mac-mini (HYP-3552) offered the metagraph's bounded Burnside
variance `CV(H) ~ 0.5-0.6` (THM-588) as the finite testbed showing "such variances stay bounded." This scan
locates the limit of that analogy: the metagraph `CV(H)` is bounded because `S_n` acts transitively (no
class has a vanishing fiber), whereas `CV(N_R)` is UNbounded because dense `R` drives `m_R -> 0` — a
mechanism with no metagraph analog. So the metagraph models the bounded-variance regime, not the `m_R -> 0`
concentration that is the LRC's actual binding worst-case. Any "metagraph = finite Siegel testbed" claim
should carry this caveat.

## Next steps

1. Floor owners: pursue HYP-3553's set-independent `Gamma_0(14)` bound; this scan shows it (not CV) is the
   route, and that the target `R' > 0` is robust.
2. Characterize WHICH dense `R` actually occur in the kps covering family (HYP-3415): if non-consecutive
   dense `R` (like `{1..13}\{12}`) never occur, the CV gatekeeper IS uniform on the real family (finding 1);
   if they do, they are MORE binding than the consecutive rows and need the exact SPEC or the `Gamma_0(N)`
   bound.
3. Realize `Var(N_R)` as a 2-dim congruence Siegel second moment (arXiv:2507.05905) and check whether its
   main-term bound is set-independent where the raw `CV^2` is not.
