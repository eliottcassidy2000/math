---
id: THM-1071
title: THE SOFT-WEYL DIVIDING LINE + EXACT CROSS-VALIDATION + the counting lemma refuted + r=4 finite horn complete; THM-1097 now supplies the separate all-scale r=4 bridge
status: (I) VERIFIED exactly. (II) MEASURED. (III) REFUTED. (IV) the one-frequency bound is PROVED. (V) PROVED FINITE-EXACT at KB=400 on all 220 cores with an independent tail referee. The later THM-1097 independently proves uniform r=4 by a sharp periodic-discrepancy component theorem, not by the bounded scans recorded here
source: kind-pasteur-2026-07-18-S128 (cont.61; owner: run the r=4 finite horn, but also pull and investigate angles other agents are working concurrently)
depends_on:
  - THM-1051, THM-1061   # my r=2 / r=3 closures and the horn-scaling law
related:
  - THM-1042             # klein: the component-length obstruction — (I) cross-validates it, (II) synthesises with it
  - THM-1037             # death-star: soft Weyl proves the position lemma — (IV) explains why it works there
  - HYP-7505, HYP-7402   # boxeph: soft Weyl closes density off-resonance / is wrong for one-line rigidity
script: 04-computation/klein_synthesis_kps_S128c61.py, kill_fraction_kps_S128c61.py, weyl_dividing_line_kps_S128c61.py, r4_measure_horn_kps_S128c61.py, r4_feasibility_kps_S128c61.py, r4_finite_horn_v2_kps_S128c61.py (+ .out)
---

# THM-1071 — the soft-Weyl dividing line, and a synthesis with three concurrent threads

## (I) Exact cross-validation with klein THM-1042

klein-S327 computed L_max(B) for B = {1,…,k} to show additive certificates can never absorb
a consecutive speed. I computed S(P) independently for the measure horn. **All nine rows
agree:**

| B | μ (mine) | μ (klein) | L_max | 1/L_max (mine) | 1/L_max (klein) |
|---|---|---|---|---|---|
| {1..7} | 0.33469 | 0.33469 | 3/49 | 16.3 | 16.3 |
| {1..9} | 0.18107 | 0.18107 | 2/63 | 31.5 | 31.5 |
| {1..11} | 0.05633 | 0.05633 | 1/77 | 77.0 | 77.0 |

Two agents, separate implementations, five-decimal agreement. This also validates my
THM-1051 table: its core-drop-12 row (largest component 0.0129870, measure 0.0563337) is
exactly klein's {1..11} row.

**And it reconciles their negative with my positive.** klein's criterion is that a base B
admits a speed w only if w > 1/L_max(B) — so a *consecutive* speed (k+1 < 1/L_max) can never
be absorbed. My measure horn never tries: cores are handled exactly, and only the
**killers** go through the absorption step, and killers exceed 143 ≫ 77. klein's theorem is
the statement that my method's hypothesis is necessary; my theorems are the case where it
holds. Same phenomenon, opposite sides.

## (II) Threshold complementarity — neither form dominates

klein tracks (μ, N) and loses 2dμ + 2dN/w per speed, giving a threshold **N/(6μ)**. I track
the largest component and get **1/(3L)**. Since L ≥ μ/N always, one might expect klein's to
win by up to 2×. Measured:

| setting | removed | μ′ | N | N/(6μ′) | 1/(3L) | winner |
|---|---|---|---|---|---|---|
| r=2 worst | 873 | 0.0905 | 102 | **187.9** | 339.5 | klein |
| r=2 drop12 | 873 | 0.0482 | 62 | **214.6** | 339.5 | klein |
| r=3 worst | 864, 897 | 0.0745 | 194 | 434.3 | **430.4** | mine |
| r=4 sample | 860,880,897 | 0.1083 | 346 | 532.6 | **354.1** | mine |

klein's wins while the set is intact; mine wins once removal has fragmented it (N grows,
their bound degrades; L tracks the best surviving piece). **Use min of the two.** At r=2
that drops the split point from 874 to ~188, shrinking THM-1051's finite horn substantially
— a free improvement to my own theorem, from klein's formulation.

## (III) The counting lemma is refuted

A family is uncertified exactly when the killers' kill-sets **cover** the core's safe (q,a)
set, so Σᵢ|kill(kᵢ)| < |bits(P)| would certify with no enumeration at all. Independence
would give kill-fraction ≈ 2⌈q/14⌉/q ≈ 1/7, and r ≤ 6 would follow. Measured worst
kill-fractions:

| core size | 11 | 10 | 9 |
|---|---|---|---|
| worst frac | 0.750 | 0.636 | 0.596 |
| r · frac | 1.50 | 1.91 | 2.38 |

All above 1 — **the counting route is dead.** The mechanism is visible in the data:
divisibility inflates the kill set (494 = 2·13·19 kills 50%, 780 = 2²·3·5·13 kills 45%), but
not only divisibility — 559 = 13·43 has *no* divisor in [15,40] and still kills 47.7%.
Conditioning on the core being safe makes killers more likely to be unsafe. Positive
correlation, not independence, and that is the same resonance that defeats every independence
heuristic in this problem.

## (IV) The dividing line: soft Weyl is a one-frequency tool

Four agents produced apparently contradictory verdicts on soft Weyl in the same week:
death-star proves the position lemma with it (THM-1037), boxeph closes density off-resonance
with it (HYP-7505) but calls it the wrong tool for one-line rigidity (HYP-7402), and my
THM-1061(V) shows the relation-support ladder diverges. **All four are correct, about two
different objects.**

**One frequency.** For G = ⋃ᵢ[aᵢ,bᵢ] with C components, integration by parts gives

> ∫_G e(kt) dt = Σᵢ [e(kbᵢ) − e(kaᵢ)]/(2πik),  so |∫_G e(kt) dt| ≤ **C/(π|k|)**.

No cancellation is needed — the constant is just the component count. Verified numerically
(ratios 0.02–0.36 against the bound). Its geometric twin is my measure horn's
μ/7 + 2C/(7k), which is **asymptotically sharp**:

| k | 157 | 900 | 20000 |
|---|---|---|---|
| bound / truth (drop 1) | 2.78 | 1.27 | **1.01** |
| bound / truth (drop 12) | 7.41 | 1.56 | **1.02** |

This one estimate is death-star's soft Weyl, boxeph's off-resonance argument, klein's
2dN/w loss term, and my boundary term — **the same inequality in four notations.**

**The product.** Expanding ∏ over all 13 speeds gives the relation-support ladder, whose
terms grow (w₂ = +1.12, w₃ = −5.23, w₄ = +12.06; partial sums 2.12, −3.11, +8.95 against a
true value of 0). The one-speed error 2C/(7k) decays like 1/k; the product's does not decay
at all.

> **Soft Weyl works whenever exactly one oscillating factor is tested against a fixed set,
> and fails whenever the 13-fold product is expanded.** That is the whole dividing line, and
> it predicts which future attempts will work.

## (V) r=4 finite horn — complete; all-scale bridge later supplied by THM-1097

Measure-horn samples put the r=4 threshold near 355. The finite horn at KB = 400 was made
feasible by a **sound** pruning: a quadruple can only be uncertified if its kill-sets cover
bits(P), which requires Σ frac ≥ 1, so quadruples failing that are certified automatically
(this cuts ~2·10⁹ raw quadruples to ~1.5·10⁸). The run reached

> **160 of 220 cores, 119,489,369 quadruples tested, ZERO uncertified**

at its core-160 progress line, but the frozen process continued.  Its terminal
summary is

> **220 of 220 cores, 142,475,077 necessary-condition quadruples, ZERO uncertified.**

A separately structured standard-library C++ referee treats each tail core as
an exact hypergraph set-cover instance and proves the stronger statement that
no cover of size at most four exists for cores 160–219.  Optimized,
unoptimized, and sanitizer builds are byte-identical.  The later
`143,112,134` headline double-counted core 160; see THM-1081 for the accounting
repair.

This completes the finite horn below `KB=400`.  By itself it does **not**
close r=4 globally.  THM-1097 now supplies the formerly missing analytic
three-comb tail plus exact guarded complement, so the stratum is uniform for
a separate rigorous reason.

## Named next
- Use THM-1097's sharp three-comb theorem as the all-scale r=4 bridge; do not
  extrapolate the historical threshold samples.
- Apply min(N/(6μ), 1/(3L)) retroactively to THM-1051 and THM-1061 to shrink both finite
  horns — at r=2 the split point falls from 874 to ~188.
- The positive correlation in (III) is the real object: it is why every independence-based
  bound in this problem fails, and quantifying it is more valuable than another census.
