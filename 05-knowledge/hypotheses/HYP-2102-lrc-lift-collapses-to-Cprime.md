---
id: HYP-2102
status: REDUCTION (C' => LRC, rigorous) + rigorous partial (large multiples) + verified n=4..14; small-multiple residual open
source: opus-2026-06-02-S571
related:
  - HYP-2097
  - HYP-2095
  - THM-369
---

# HYP-2102: the lift lemma collapses to C' ("a multiple of n forces positive measure")

**SPLIT (simplification):** for speed set S test t=1/n. (a) NO multiple of n => every ||v_i/n||>=1/n => t=1/n is a witness, M>=1/n (THM-369, free). (b) A multiple of n (v=nw) kills every j/n; the only hard case. So the lift lemma = the hard case = C'.
**C' (S564, as a tight-set property):** n|v for some v in S => M(S)>1/n (loose).
**REDUCTION (headline, rigorous):** C' => LRC(n). Every config has no mult of n (=>1/n witness) or has one (=>M>1/n by C'); either way M>=1/n. So proving "a multiple of n can't be tight" proves ALL of LRC; S564 saw C' only as a symptom of tight sets, but its converse direction IS the whole conjecture.
**VERIFIED n=4..14, ALL multiplier sizes** (v=n exactly=hardest w=1, small w<=3, large w>=4): every mult-of-n config LOOSE, 0 tight-with-multiple, min safe-measure always >0 (~0.02-0.04). (`lrc_lift_Cprime_residual_s571.py`, `lrc_lift_lemma_measure_bound_s571.py`.)
**CRUDE BOUND FAILS:** mu(safe S)>=mu(safe S')-2/n is too lossy -- min mu of (n-2)-runner configs at 1/n is far below 2/n (n=6, S'={1,3,4,5}: 0.05<<1/3). Real mechanism: v's thin evenly-spaced arcs can't COVER safe(S').
**RIGOROUS PARTIAL (large multiples dodgeable, uses PROVEN LRC(n-1)):** if v=nw>(n-1)*max(other speeds) then M(S)>1/n. Proof: S'=S\{v} (n-2 runners) has t0 with min||v_i t0||>=1/(n-1) by LRC(n-1); on interval I half-width delta=1/(n(n-1)V') around t0 all S'>1/n; v's arcs have radius rho=1/(n^2 w); delta>rho <=> nw>(n-1)V' => I wider than an arc => contains a v-safe sub-interval => positive measure.
**RESIDUAL (open, sharp):** small multiples v<=(n-1)*max(others) (down to v=n). Loose empirically; needs equidistribution -- the AP of arcs (period 1/(nw), total 2/n) meets safe(S') in ~(2/n)mu<mu, can't align to cover. A 3-distance/Weyl statement about one AP vs a fixed union of intervals.
**SPEEDUPS:** no-mult configs free (1/n witness); only the mult-of-n slice needs a fast positive-measure certificate (one safe interval), not exact M; safe_measure via endpoint enumeration (kn+-1)/(nv), grid-free; residual is small-multiples band only.

**See:** `07-reflections/lrc-lift-lemma-reduces-to-Cprime-multiple-of-n-loose-s571.md`, `04-computation/lrc_lift_Cprime_residual_s571.py` (+.out), `lrc_lift_lemma_measure_bound_s571.py` (+.out); HYP-2097 (64-class), HYP-2095 (lift lemma), S564 (C'), THM-369.
