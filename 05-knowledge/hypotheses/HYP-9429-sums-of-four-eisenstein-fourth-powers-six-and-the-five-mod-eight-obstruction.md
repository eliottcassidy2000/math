---
id: HYP-9429
title: "Sums of four fourth powers of Eisenstein integers: a rational integer n = 5 (mod 8) is never one, and every other n is one in every completion of Q(omega) (PROVED, THM-4826 (6)-(7)); is every rational n != 5 (mod 8) a global sum (integral Hasse principle)? 6 and 54 are, but only with large bases (6 needs bases up to 541); 60, 75, 78, 124, 140 have no representation with bases at most 650 (OPEN); pre-registered tests (1) CONFIRMED, (2) and (3) REFUTED"
status: >
  OPEN for 60, 75, 78, 124, 140 and for the general question. Pre-registered 2026-10-27: (1) and (2) in 06ca4c5d48 before
  their runs, (3) in d0a9ad6895 before its run; (1) CONFIRMED, (2) REFUTED, (3) REFUTED. Independent audit
  (04-computation/experiments/two_grids_class_group_20261027_audit_S66/AUDIT.md): pre-registration integrity verified;
  6 and 54 are representable; the other five have no representation with bases <= 650; corrections applied
  (MISTAKE-690).
session: opus-2026-10-27-S66
source: 05-knowledge/results/two_grids_ring_or_ideal_class_13_37_20261027.md
related:
  - THM-4826 (the n = 5 (mod 8) obstruction (6) and the local criterion (7); Eisenstein fourth powers mod 16 are 0, mu_3, 9 mu_3)
  - session78 FOURTH_POWERS.md (the rational-integer window: four fourth powers below 81)
---

# HYP-9429 -- four Eisenstein fourth powers: the 5 (mod 8) obstruction, and the integral Hasse principle

**Setting.** A rational integer n is a *sum of four Eisenstein fourth powers* if n = α⁴ + β⁴ + γ⁴ + δ⁴ with α, β, γ, δ ∈
Z[ω] (zeros allowed). This is the owner's four-fourth-powers question read on the Eisenstein grid. Unlike the rational
case, the terms are complex and can cancel, so large bases can represent small n.

**What is known (THM-4826 (6), PROVED).** Every Eisenstein fourth power is ≡ 0, u or 9u (mod 16) with u ∈ {1, ω, ω²}.
Hence a rational sum of four of them is never ≡ 5 (mod 8). In 0..50 this excludes exactly 5, 13, 21, 29, 37, 45, so the
owner's 13 and 37 stay excluded on the Eisenstein grid. Local tests modulo 16, 32, 64, 48, 80, 112, 144, 9, 27, 5, 25, 7
and 49 find no further obstruction in 0..50.
- *[Added after the audit.]* THM-4826 (7) now proves that n ≡ 5 (mod 8) is the *only* local obstruction: every other
  rational n is a sum of four fourth powers in every completion of Q(ω).

**What was found (FINITE, a lower bound).** With bases |a|, |b| ≤ 20 (a + bω), every n in 0..50 with n ≢ 5 (mod 8) is a
sum of four Eisenstein fourth powers, except 6.
- *[Corrected after the audit: the first draft said "some need large cancelling bases: 47 uses (−20, −10), …".]* The
  terms do cancel, but small bases suffice: every represented n ≤ 50 has a representation with |a|, |b| ≤ 5, for example
  47 = (4 + 2ω)⁴ + (2 − 2ω)⁴ + (4 + ω)⁴ + (3 − ω)⁴.

**Six.** If n ≡ 6 (mod 8), the residue count forces exactly two bases ≡ ω and two ≡ ω² (mod 2). An exhaustive search
restricted to that pattern, with bases up to 100, finds no representation of 6. The same search finds 30, which is also
≡ 6 (mod 8).
- *[Corrected after the audit: the first draft argued only at 2.]* 6 is representable in every completion of Q(ω):
  - by Hensel's lemma in Z₂[ω] (a unit term, precision 2⁵) and in Z₃[ω] (a unit term mod (1 − ω));
  - by a direct check at the residue fields F₂₅, F₇, F₁₃ and F₁₉;
  - by Weil's bound above them.

  So no local obstruction explains it.

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
   - **No local obstruction.** *[Corrected after the audit: the first draft argued from tests modulo 2⁸, 3⁵, 5³ and 7²
     plus the Chinese remainder theorem, which says nothing about primes ≥ 11.]* None of the eight, nor 6, is
     obstructed at any place, by THM-4826 (7): Hensel at 2 and 3, direct checks at F₂₅, F₇, F₁₃ and F₁₉, and Weil's bound
     above. The tests modulo 2⁸, 3⁵, 5³ and 7² agree.
   - **Pattern.** *[Corrected after the audit: the first draft said "six of the nine".]* Seven of the nine misses
     (6, 54, 60, 75, 78, 102, 186) are multiples of 3, against 45 of 132 (0.34) among the candidates.
     - The naive hypergeometric tail is P(≥ 7) = 0.0075, and 0.045 for 5 of the 7 persistent misses.
     - This is post hoc, not a test.
     - Multiples of 3 have the larger 3-adic density (32/27 against 8/9), and 54 = 9·6 is tied to 6 by scaling with
       1 + 2ω = √−3.
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
rate behind (3) assumed.

**Protocol note (from the audit).**
- The registration names a "class-restricted" search, but (2) and (3) were run with the unrestricted exhaustive mode
  `all`. That is a stronger search than the one named.
- The source note was first committed after the runs.
- The audit verified from git history that each prediction precedes the commit that reports its run.

**The audit's searches (independent code, three methods; `two_grids_class_group_20261027_audit_S66/`).**
- **6 is representable**, with least box |a|, |b| ≤ 541:

  6 = (12 − 149ω)⁴ + (150 + 541ω)⁴ + (123 + 427ω)⁴ + (297 − 309ω)⁴.

  There is nothing with bases ≤ 540. With bases ≤ 800 this representation and its complex conjugate are the only ones,
  up to signs and order. The bases lie in the classes that THM-4826 (6) forces for n ≡ 6 (mod 16).
- **54 is representable**, by scaling 6 with (1 + 2ω)⁴ = 9:

  54 = (310 + 173ω)⁴ + (−932 − 241ω)⁴ + (−731 − 181ω)⁴ + (915 + 903ω)⁴.

  There is nothing with bases ≤ 800.
- Both identities were re-verified exactly by the author, by integer arithmetic in (1, ω) and symbolically.
- **60, 75, 78, 124, 140**: no representation with bases ≤ 650 (3.4·10¹⁰ to 5.6·10¹⁰ pair sums per target).
- **Calibration.** Class-mates have few representations at these sizes, and the counts barely grow from B = 200 to
  B = 500. In the class of 6 at B = 200: 22: 9, 38: 3, 86: 1, 134: 1, 166: 1.
  - The non-obstructed misses in 51..200 fall 17, 15, 8, 6, 6, 6 at B = 20, 30, 60, 120, 300, 500.
  - Persistence to large B is therefore weak evidence of non-representability. 6 and 54 support the
    sparse-representation reading.

**Post-hoc check (not pre-registered;
`04-computation/experiments/eisenstein_fourth_power_obstructions_20261027.py`).** The 2-adic × 3-adic local densities
were computed by four-fold convolution modulo 256 and 243.
- **Exactness.** The values are exact for every n ≤ 200. *[After the audit:]* the first draft used modulus 128, which is
  not exact for the nine multiples of 16 (16, 32, 48, 64, 128, 144, 160, 176, 192), because
  σ₂(16m) = σ₂ᵘ(16m) + σ₂(m). With modulus 256 the median and every rank are unchanged.
- **Stability.** By Hensel's lemma the densities are stable from 2⁵ for n ≢ 0 (mod 16) and from 3³ for n ≢ 0 (mod 9).
  54 = 2·27 settles from 3⁴.
- **The five persistent misses all have below-median density** (the median over n ≢ 5 (mod 8), n ≤ 200, is 2.667).
  Densities and ranks among 175 (ties broken by n):

  | n | density | rank | status |
  |---|---|---|---|
  | 124 | 0.444 | 11 | none with bases ≤ 650 |
  | 140 | 0.444 | 12 | none with bases ≤ 650 |
  | 60 | 0.593 | 19 | none with bases ≤ 650 |
  | 78 | 1.778 | 50 | none with bases ≤ 650 |
  | 75 | 2.370 | 66 | none with bases ≤ 650 |
  | 6 | 1.778 | 42 | representable (bases up to 541) |
  | 54 | 3.556 | 109 | representable (bases up to 932) |

  The two that filled in at B ≤ 120 rank 53 (102) and 117 (186).
- **Why 124 and 140 are low.** n ≡ 4 (mod 8) forces all four bases ≡ 1 (mod 2), which gives the lowest density class
  (0.444).
- **But low density does not decide it.** Every other member of that class up to 200 is represented: 4, 20, 28, 44,
  52, 68, 76, 92, 100, 116, 148, 164, 172, 188, 196.

**Status after the runs and the audit.** (1) CONFIRMED, (2) REFUTED, (3) REFUTED. 6 and 54 SETTLED (representable).
- **Open question (integral Hasse principle).** By THM-4826 (7), n ≡ 5 (mod 8) is the only local obstruction. Is every
  rational n ≢ 5 (mod 8) a sum of four fourth powers of Eisenstein integers? The open cases below 200 are 60, 75, 78,
  124 and 140.
- **What would settle it.**
  - A representation settles each case.
  - Ruling one out needs a global argument, for instance an integral Brauer–Manin obstruction on the affine threefold
    x₁⁴ + x₂⁴ + x₃⁴ + x₄⁴ = n over Z[ω].
- **The heuristic.** It predicts infinitely many representations as B grows, at a rate c(n)·log B that is small for
  low-density n. This is not a prediction of where they appear.
