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
**VERIFIED n=4..14, ALL multiplier sizes** (v=n exactly=hardest w=1, small w<=3, large w>=4): every mult-of-n config LOOSE, 0 tight-with-multiple, min safe-measure always >0 (~0.02-0.04). (`lrc_lift_Cprime_residual_s571.py`, `lrc_lift_lemma_measure_bound_s571.py`.) **Extended by SAMPLING through n=20** (monad-compute S595/S596: 33600/33600 mult-of-n loose, 0 tight, n=15..20).

**EXHAUSTIVE BOX CERTIFICATE (monad-compute-2026-06-03-S597):** prior checks all SAMPLE the mult-of-n class. S597 instead enumerates EVERY primitive (n-1)-subset of {1..B} containing a multiple of n (gcd=1; M is scale-invariant so primitive reps suffice) and tests looseness EXACTLY (open-safe-set measure>0, early-exit positivity; validated vs full measure, 0 mismatches; tight APs correctly give measure 0). Boxes B=K*n: n=4 (K=10), n=5 (K=8), n=6 (K=6), n=7 (K=4), n=8 (K=3). RESULT: **0 tight over every config in the box** -- n=4 (4615), n=5 (51957), n=6 (225915), n=7 (240009), n=8 (229061). Min margin (smallest mu>0) at the AP-like rows: n=4 -> 1/24 at (1,3,4); n=5 -> 1/50 at (1,3,4,5). A COMPLETENESS statement (no exceptions inside the box), strictly stronger than sampling, for the C' hypothesis class. Files: `04-computation/lrc_Cprime_exhaustive_box_monad_s597.py` (n=4..7 + n=8 partial), `lrc_Cprime_exhaustive_box_n8_monad_s597.py` (n=8 completion), results `*_s597.out`. Still EMPIRICAL within finite boxes -- the small-multiple residual remains the open analytic core (equidistribution / one-AP-vs-union-of-intervals).

**WIDENED BOX + n=9 (monad-compute-2026-06-03-S598):** S597's handoff (widen K / push n=9) executed. New exact INTEGER arc-cover engine (loose <=> closed unsafe arcs don't cover the circle; scale by D=n*lcm(S) to integers; measure = uncovered units / D) runs ~40x faster than the S597 Fraction scan -- SELF-CHECKED vs S597's `is_loose`/`open_safe_measure` (0 disagreements on 6000+5000 random configs; AP correctly tight n=4..14). Boxes strictly wider for every n: n=4 (K=20), n=5 (K=14), n=6 (K=8), n=7 (K=5), n=8 (K=4), **n=9 (K=3, NEW)**. RESULT: **6,237,910 configs, 0 tight** -- n=4 (37608), n=5 (511830), n=6 (1011829), n=7 (1021756), n=8 (2171084), n=9 (1483803), all PASS. Min margins (smallest mu>0, all at AP-like rows with a multiple appended): n=4 1/24 @(1,3,4); n=5 1/50 @(1,3,4,5); n=6 7/540 @(1,3,4,5,18); n=7 1/147 @(1,3,4,5,7,24); n=8 93/14560 @(2,6,7,8,10,13,14); **n=9 113/12852 @(1,3,4,5,7,9,17,24)**. First exhaustive certificate at n=9; combined with S597 still EMPIRICAL inside finite boxes (small-multiple residual = the open analytic core). Files: `04-computation/lrc_Cprime_exhaustive_box_widen_monad_s598.py` (+`_s598.out`).
**CRUDE BOUND FAILS:** mu(safe S)>=mu(safe S')-2/n is too lossy -- min mu of (n-2)-runner configs at 1/n is far below 2/n (n=6, S'={1,3,4,5}: 0.05<<1/3). Real mechanism: v's thin evenly-spaced arcs can't COVER safe(S').
**RIGOROUS PARTIAL (large multiples dodgeable, uses PROVEN LRC(n-1)):** if v=nw>(n-1)*max(other speeds) then M(S)>1/n. Proof: S'=S\{v} (n-2 runners) has t0 with min||v_i t0||>=1/(n-1) by LRC(n-1); on interval I half-width delta=1/(n(n-1)V') around t0 all S'>1/n; v's arcs have radius rho=1/(n^2 w); delta>rho <=> nw>(n-1)V' => I wider than an arc => contains a v-safe sub-interval => positive measure.
**RESIDUAL (open, sharp):** small multiples v<=(n-1)*max(others) (down to v=n). Loose empirically; needs equidistribution -- the AP of arcs (period 1/(nw), total 2/n) meets safe(S') in ~(2/n)mu<mu, can't align to cover. A 3-distance/Weyl statement about one AP vs a fixed union of intervals.
**SPEEDUPS:** no-mult configs free (1/n witness); only the mult-of-n slice needs a fast positive-measure certificate (one safe interval), not exact M; safe_measure via endpoint enumeration (kn+-1)/(nv), grid-free; residual is small-multiples band only.

**See:** `07-reflections/lrc-lift-lemma-reduces-to-Cprime-multiple-of-n-loose-s571.md`, `04-computation/lrc_lift_Cprime_residual_s571.py` (+.out), `lrc_lift_lemma_measure_bound_s571.py` (+.out); HYP-2097 (64-class), HYP-2095 (lift lemma), S564 (C'), THM-369.
