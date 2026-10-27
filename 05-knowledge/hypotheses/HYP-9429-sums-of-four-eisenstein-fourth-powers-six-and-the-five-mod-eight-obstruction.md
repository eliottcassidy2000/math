---
id: HYP-9429
title: "Sums of four fourth powers of Eisenstein integers: a rational integer n = 5 (mod 8) is never one (PROVED, THM-4826); among 0..50 every other n is one (bases at most 20) except 6, which has no representation with bases at most 100 although it is representable locally everywhere; pre-registered: (1) 6 has no representation with bases at most 150, (2) at most three n in 51..200 with n not 5 (mod 8) lack a representation with bases at most 60"
status: OPEN (pre-registered 2026-10-27 before the runs of (1) and (2); results to be appended below the line)
session: opus-2026-10-27-S66
source: 05-knowledge/results/two_grids_ring_or_ideal_class_13_37_20261027.md
related:
  - THM-4826 (the n = 5 (mod 8) obstruction; Eisenstein fourth powers mod 16 are 0, mu_3, 9 mu_3)
  - session78 FOURTH_POWERS.md (the rational-integer window: four fourth powers below 81)
---

# HYP-9429 -- four Eisenstein fourth powers: the 5 (mod 8) obstruction, and 6

**Setting.** A rational integer n is a *sum of four Eisenstein fourth powers* if n = α⁴ + β⁴ + γ⁴ + δ⁴ with α, β, γ, δ ∈
Z[ω] (zeros allowed). This is the owner's four-fourth-powers question read on the Eisenstein grid. Unlike the rational
case, the terms are complex and can cancel, so large bases can represent small n.

**What is known (THM-4826 (6), PROVED).** Every Eisenstein fourth power is ≡ 0, u or 9u (mod 16) with u ∈ {1, ω, ω²}.
Hence a rational sum of four of them is never ≡ 5 (mod 8). In 0..50 this excludes exactly 5, 13, 21, 29, 37, 45, so the
owner's 13 and 37 stay excluded on the Eisenstein grid. Local tests modulo 16, 32, 64, 48, 80, 112, 144, 9, 27, 5, 25, 7
and 49 find no further obstruction in 0..50.

**What was found (FINITE, a lower bound).** With bases |a|, |b| ≤ 20 (a + bω), every n in 0..50 with n ≢ 5 (mod 8) is a
sum of four Eisenstein fourth powers, except 6. Some need large cancelling bases: 47 uses (−20, −10), (−13, −17),
(−8, −15), (−8, 6).

**Six.** If n ≡ 6 (mod 8), the residue count forces exactly two bases ≡ ω and two ≡ ω² (mod 2). An exhaustive search
restricted to that pattern, with bases up to 100, finds no representation of 6. The same search finds 30, which is also
≡ 6 (mod 8). 6 is representable modulo every modulus tested, and by Hensel's lemma at 2 (a unit term, precision 2⁵)
it is representable in Z₂[ω]. So no local obstruction explains it.

**Heuristic.** The equation is in the critical dimension: eight real variables, two real equations of degree 4. The
expected number of representations with bases ≤ B therefore grows like c(n)·log B, as for sums of three cubes. Zero
representations up to B = 100 suggests a small c(6), not impossibility.

**Pre-registered predictions** (filed before the runs; the runs use the class-restricted meet-in-the-middle search of
the source note §5):
1. 6 has no representation with all bases |a|, |b| ≤ 150.
2. At most three integers n in 51..200 with n ≢ 5 (mod 8) have no representation with all bases |a|, |b| ≤ 60.

**Falsifier.** A representation of 6 with bases ≤ 150 refutes (1). Four or more n in (2) refute (2).

---

*Results (appended after the runs, 2026-10-27; program `04-computation/experiments/eisenstein_four_fourth_powers_20261027.c`, validated
on the B = 20 window and on 22 and 38, which share 6's residue class mod 16):*

1. **CONFIRMED.** The exact box search (modes `box`, |a|, |b| ≤ 150) finds no representation of 6: 0 hits among
   32,062,500 × 31,640,625 pair-sum combinations, in 14 s.
2. **REFUTED.** With bases ≤ 60 (mode `all`), 27 of the 150 integers in 51..200 have no representation:
   - 19 are the obstructed n ≡ 5 (mod 8): 53, 61, …, 197;
   - 8 are not obstructed: 54, 60, 75, 78, 102, 124, 140, 186.

   Eight is more than the predicted three. The prediction underestimated how sparse representations are at this size.
   - **No local obstruction.** None of the eight, nor 6, is obstructed modulo 2⁸, 3⁵, 5³ or 7², and composite moduli
     add nothing by the Chinese remainder theorem.
   - **Pattern.** Six of the nine misses (6, 54, 60, 75, 78, 102) are multiples of 3, against about one third
     expected. This is not significant at this size.
   - **Reading.** The misses are consistent with the critical-dimension heuristic (counts growing like c(n)·log B).
     They are not explained by any congruence.

**Pre-registered (filed after the results above, before its run):**

3. At least 3 of the 8 values 54, 60, 75, 78, 102, 124, 140, 186 acquire a representation with all bases |a|, |b| ≤ 120
   (mode `all`).
   - **What it tests.** The log-growth reading predicts that most misses fill in as B grows. A hidden global
     obstruction predicts that none do.
   - **Falsifier.** Zero, one or two of the eight acquire a representation.

**Result of (3): REFUTED.** Mode `all` with bases ≤ 120 (421,704,361 pair sums, 102 s) represents only 2 of the 8:
- 102 = (113 + 99ω)⁴ + (51 − 41ω)⁴ + (104 + 31ω)⁴ + (106 + 13ω)⁴;
- 186 = (37 − 44ω)⁴ + (16 + 40ω)⁴ + (54 + 36ω)⁴ + (77 + 14ω)⁴.

54, 60, 75, 78, 124 and 140 remain unrepresented, and so does 6. The misses persist longer than the rough log-growth
rate behind (3) assumed. Neither reading is settled: a slower fill-in and a hidden global (non-congruence) obstruction
both remain possible.

**Post-hoc check (not pre-registered;
`04-computation/experiments/eisenstein_fourth_power_obstructions_20261027.py`).** The 2-adic × 3-adic local densities
were computed exactly by four-fold convolution modulo 128 and 243.
- **Stability.** By Hensel's lemma the densities are stable from 2⁵ for n ≢ 0 (mod 16) and from 3³ for n ≢ 0 (mod 9);
  both are checked on 1..200. For the misses they are printed at every modulus. 54 = 2·27 settles from 3⁴.
- **The misses.** 6 of the 7 persistent misses have below-median density (the median over n ≢ 5 (mod 8), n ≤ 200, is
  2.667). Their densities and ranks among 175 (ties broken by n) are:

  | n | density | rank |
  |---|---|---|
  | 124 | 0.444 | 11 |
  | 140 | 0.444 | 12 |
  | 60 | 0.593 | 19 |
  | 6 | 1.778 | 42 |
  | 78 | 1.778 | 50 |
  | 75 | 2.370 | 66 |
  | 54 | 3.556 | 109 |

  The two that filled in rank 53 (102) and 117 (186).
- **Why 124 and 140 are low.** n ≡ 4 (mod 8) forces all four bases ≡ 1 (mod 2), which gives the lowest density class
  (0.444).
- **But low density does not decide it.** Every other member of that class up to 200 is represented: 4, 20, 28, 44,
  52, 68, 76, 92, 100, 116, 148, 164, 172, 188, 196.

So low local density raises the chance of a miss, as the heuristic says, but does not determine it. The status of the
seven is OPEN.

**Status after the runs.** (1) CONFIRMED, (2) REFUTED, (3) REFUTED.
- **Open question.** Are 6, 54, 60, 75, 78, 124, 140 sums of four fourth powers of Eisenstein integers at all?
- **Settling it.** A representation settles each one. An integral Brauer–Manin computation, or another global
  argument, would be needed to rule one out.
