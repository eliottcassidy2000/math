# HYP-3902: the tower floor holds — exact rearrangement constants F_j > 1/36 for all cluster sizes j = 7..11 (range ≤ 14), worst-case compact position, AP pattern extremal

**Status:** CONFIRMED (exact rational computation over all primitive patterns; F1/F2 structure lemmas proved)
**Instance:** opus-2026-07-01-S33
**Script:** `04-computation/lrc14_tower_floor_rearrangement_opus_20260701_S33.py` (+ `.out`)
**Consumes:** HYP-3900 (peel), HYP-3901 (renormalization), THM-592 (kink taxonomy), (S1) union bound.
**Consumed by:** the assembled proof `03-artifacts/drafts/lrc14-11core-floor-assembled-proof-opus-20260701-S33.md`.

## The functional

For a deep cluster {N + c_i} with primitive difference pattern c (size j, range ≤ 14) over a compact
part B (|B| = 11 − j), the depth-1 floor needs inf over patterns and positions of ∫_{L_B} D_c(t)dt.
Worst-case position at fixed measure m = meas(L_B) ≥ (j−4)/7 is the increasing-rearrangement integral

    Q_c(m) = ∫_0^m D_c^*(u) du = ∫_0^∞ (m − ψ_c(s))_+ ds .

D_c is piecewise linear with kinks of exactly two types (THM-592(iii) taxonomy in t): opposite-endpoint
δt ≡ ±1/7 (convex) and same-endpoint δt ≡ 0 (concave peaks); grid (7m+e)/(7δ), e ∈ {−1,0,1}. Omitting
e = 0 under-estimates Q (conservative; MISTAKE-092, caught by the slope canary).

## The constants (exact, worst-40 per size re-verified in rationals; ~6300 patterns scanned)

    F_7  = 559/11025      = 0.050703  at (0,1,2,3,4,5,6)          1.825 x (1/36)
    F_8  = 184019/3246495 = 0.056682  at (0..7)                   2.041 x
    F_9  = 244547/3522610 = 0.069422  at (0..8)                   2.499 x
    F_10 = 56333/617400   = 0.091242  at (0,2,3,4,5,6,7,8,10,12)  3.285 x
    F_11 = 63941/432180   = 0.147950  at (0,2,4,5,6,7,8,9,10,12,14) 5.326 x

All clear 1/36; monotone in j (more cluster = fewer compact = safer, matching kps HYP-3950's
monotone-in-outlier-count law). The j=7 argmin is the CONSECUTIVE (AP) pattern = the renormalization
fixed point (HYP-3901), and F_7 = 0.0507 > pentagon 0.0323: the compact census extremizer stays the
global binding case with a proved 57% gap at the limit.

## Structure lemmas (proved in the assembled doc, verified here)

- **F1 (continuous-Fraenkel rigidity, j=7):** D_c(t) = 0 iff the 7 danger arcs tile iff centers form a
  translate of (1/7)Z; primitive patterns have zeros only among {k/7 : k=1..6} (verified: consecutive
  pattern zeros exactly that set).
- **F2 (slope at a tiling):** one-sided slopes at a zero = half the cyclic total variation of the
  tiling order ≥ range(c) ≥ 6 (verified exactly: 6 both sides at t = 1/7).
- j ≥ 8: interval zero-sets exist (consecutive j=8: measure 44/735) — absorbed by Q, no case work.

## Honest scope

The floor is proved AT THE N→∞ LIMIT for primitive pattern range ≤ 14, worst-case position, via the
freeze lemma F3 (proved at O(√(Δ/N)); observed O(1/N)). Remaining to full closure of Case 3: the F3
rate sharpening (shrinks the finite height band from N* ~ 10^8 to ~10^3), and the large-range pattern
recursion (R2) — both named, neither conceptual. See the assembled doc §6 ledger.

-> HYP-3900, HYP-3901, THM-592, THM-593, MISTAKE-092, kps HYP-3950, mac-mini HYP-3850 (Mirsky-Newman
   alternative for finite N), OPEN-Q-108.
