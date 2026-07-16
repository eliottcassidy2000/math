# HYP-7018 — Compact-core flatness: the resonant-mode law in the balanced regime

**Status:** CLAIMED / IN PROGRESS (death-star-2026-07-16-S19; owner directive "prove
compact-core flatness"). Verify-first until results land.

## The claim being built

klein's THM-883-resonant-mode (cont.4) proves S(ta) = t(1−e(a/7))m̂_s(a) + O(M_slow) by
conditioning on the FASTEST runner — unavailable for compact cores (no scale separation:
O(M_slow) = O(M) swamps the main term). This session extends the law into the balanced
regime by EQUIDISTRIBUTION instead of separation:

1. **The independent-limit miss-measures (exact rationals, hand-derived, to be verified):**
   with the stationary runner pinning section 0 and five balanced others uniform,
   A* = 360/16807 (others miss exactly {s,c}, c ∉ {0,s}; A_0 = 0), B* = 120/16807
   (others miss exactly {s}); note A* = 3B*. Then
   **m̂*_s(a) = −(A*+B*)e(as/7) − A*, |m̂*| ≥ B* > 0.**
2. **The law:** for incoherent compact cores (no small internal relations), each owner e
   carries resonant modes |S(ea + slow-shift)| → e·|1−e(a/7)|·|m̂^{(e)}_s(a)| with
   m̂^{(e)} → m̂* at the core's equidistribution rate. Since m̂* ≠ 0: **absolute-constant
   compact-core flatness is FALSE asymptotically** — the same law as the far bank, slope
   ε ≈ |1−ω^a||m̂*| ∈ [~0.014, ~0.10].
3. **The program-range theorem (the true "flatness"):** the slope is so small that resonant
   modes stay below the flat fluctuation floor until e ~ 10³–10⁴; on the bounded-Vmax range
   the assembly needs (V ≤ ~500 per the sweeps), compact-core flatness holds with explicit
   constants — provable by the law + finite per-instance decidability (THM-881 w-freeness).
4. **Self-similarity:** consecutive-type compact cores have slow DIFFERENCE systems — their
   m̂^{(e)} is governed by the difference core (the renormalization/difference flow of
   HYP-3901 appears exactly here; dilate covariance S_{cE}(cn) = c·S_E(n) is the exact
   scaling face).

## Session verification plan

Resonant-frequency probes (no full Z_P scans needed): compute |S(n)| directly at the
predicted n = ea + (small combos) for compact cores to c ~ 10³⁺; measure approach to the
law; exact A*, B* combinatorics referee; the dilate-covariance exact check.

-> klein THM-883-resonant-mode (cont.4), HYP-7017 (S18 census), HYP-6994 (refuted),
THM-881/880/729, HYP-3901; death-star-S19.
