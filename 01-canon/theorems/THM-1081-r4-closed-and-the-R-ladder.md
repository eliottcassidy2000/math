---
id: THM-1081
title: THE r=4 CLUSTERED STRATUM IS NOW UNIFORMLY CLOSED BY THM-1097; its earlier finite horn is independently complete, with the canonical all-220 total corrected to 142,475,077 and a DPLL tail referee
status: PROVED uniformly for r=4 by THM-1097's sharp periodic discrepancy, analytic guard surfaces, and 39,778,595-triple exact endpoint bank. Independently, the below-400 finite horn is PROVED FINITE-EXACT on all 220 cores and its 60-core tail has a separate hypergraph-DPLL referee. The old 0.98453 R value remains bounded-window telemetry, not the uniform proof
source: kind-pasteur-2026-07-18-S128 (cont.62; owner: finish the remaining 60 cores and scan the r=4 threshold properly)
depends_on:
  - THM-1071    # the r=4 partial run this completes, and the pruning it introduced
  - THM-1051, THM-1061   # the r=2 / r=3 closures whose finite horns (V) shows were redundant
  - THM-1097    # later sharp three-comb theorem supplying the genuine all-scale bridge
related:
  - THM-1042    # klein: the component-length obstruction; R < 1 is its quantitative converse
script: 04-computation/r4_finite_horn_v2_kps_S128c61.py, r4_finite_horn_tail_kps_S128c62.py, r4_finite_horn_tail_referee_codex_S66.cpp, R_exhaustive_kps_S128c62.py, R_scan_r3r4_kps_S128c62.py, R_scan_r4_kps_S128c62.py (+ .out)
---

# THM-1081 — r=4 uniformly closed; finite horn independently complete

> **Uniform supersession (codex-S73; THM-1097).**  The finite-horn accounting
> below remains exact.  Its former all-scale gap is now closed by a different
> proof: the sharp discrepancy `|J∩D_k|<=|J|/7+6/(49k)`, an analytic
> three-comb tail, and a guarded endpoint bank covering every omitted triple.
> No bounded `R` scan is used for the uniform conclusion.

## (I) The r=4 finite horn is complete

| run | scope | exact work | failures |
|---|---|---|---|
| cont.61 canonical full output | all 220 cores | 142,475,077 necessary-condition quadruples | 0 |
| cont.62 tail replay | cores 160–219 | 23,622,765 quadruples | 0 |
| independent C++ DPLL | cores 160–219 | 39,592 memo states; no cover of size at most four | 0 |

The pruning is sound: a quadruple can only be uncertified if its four kill-sets *cover* the
core's safe (q,a) set, which requires Σ frac ≥ 1; quadruples failing that certify
automatically. So the check is exhaustive despite testing only the tail.

The old total `143,112,134=119,489,369+23,622,765` is not a valid work
total: the progress line at core 160 was printed after processing core 160,
while the tail replay starts at core 160, so that core's 637,057 quadruples
were counted twice.  The canonical full output already continued through all
220 cores and gives the correct total `142,475,077`.  The conclusion is
unchanged.  The independent referee changes the problem to exact hypergraph
set cover and proves the stronger no-cover-of-size-at-most-four predicate on
all sixty tail cores.

## (II) Scanning the threshold properly — the right quantity is R

I had been scanning the absolute threshold T. That is the wrong object. The measure horn
removes the r−1 smaller killers and needs the largest to exceed T, and the largest always
exceeds every removed killer, so what matters is

> **R := T / k_max-removed**,  and **R < 1 ⟹ the measure horn certifies with no finite horn.**

For r=2 this is cheap enough to settle exhaustively — all 12 cores, every
k₁ ∈ (13·max P, 4000):

> **max R = 0.51852**, attained at core-drop 6 with k₁ = 160; **zero** cases where the
> removal swallows S(P).

## (III) The worst case is at small killers — my earlier scans looked in the wrong place

| r | worst-case killers | 13·max P | max R |
|---|---|---|---|
| 2 | 160 | 156 | 0.51852 |
| 3 | (153, 159) | 143 | **119/158 = 0.75316** |
| 4 | (150, 156, 158) | 143 | 0.98453 |

Every worst case sits **just above 13·max P**, i.e. at the smallest admissible killers,
tightly clustered. My cont.60 and cont.61 scans sampled the *top* of the killer range on the
reasoning that dense bad-sets chop the safe set finest — that reasoning was wrong, and the
ratios those scans reported (0.389–0.434, and the "generic 7/18") are the behaviour of the
easy end, not the worst case. The 7/18 value is still the correct *asymptotic* constant;
it is simply not the maximum.

## (IV) The R-ladder

> **0.51852 (r=2 scan) → 119/158 = 0.75316 (r=3 theorem) → 0.98453 (r=4 window)**

Rising steeply, with **1.5% of margin left at r=4** and a straightforward extrapolation
above 1 at r=5.

## (V) What this does to the earlier theorems

- **THM-1051 (r=2)** remains uniform by its explicit three-way split: exact
  both-small horn, exact mixed scan, and a simultaneous two-comb bound on a
  fixed core-safe interval.  The bounded R scan is not a replacement proof.
- **THM-1061 (r=3)** is now uniform for a different rigorous reason:
  THM-1094 proves the exact two-comb component inequality at all scales.  Its
  finite horn is therefore redundant, but the earlier sampled R argument was
  not a proof (MISTAKE-163).
- **r=4 is now uniform.**  The old `R=0.98453` remains a bounded-window
  measurement and is not promoted.  THM-1097 proves the sharper component
  target `L>1/(7k3)` at every scale; a fourth killer cannot cover it.  The
  below-400 horn is therefore an independent verification.
- **r=5 does break the older `1/(3L)` measure target on exact rows,** but
  MISTAKE-164 shows that the below-235 finite horn does not cover its
  all-scale complement.  Uniform `r=5` therefore remains open; a bounded
  ladder value is telemetry, not a reduction theorem.

The structural reading is now precise: finite scans do not determine an
all-scale rate, but the sharp periodic discrepancy does.  THM-1094 closes two
removals and THM-1097 closes three.  Its asymptotic toothpick quotient has a
real method wall at four removals, so uniform `r=5` needs additional overlap
or endpoint self-similarity rather than another sampled ray.

## Named next
- r=5: prove a four-comb endpoint/self-similarity tail.  THM-1101's below-235
  horn is bounded evidence only after MISTAKE-164 withdraws its sampled
  all-scale inference.
- Retire the split points from THM-1051 and THM-1061 in favour of "R < 1, measure horn
  alone" — it is a strictly simpler statement of both theorems.
