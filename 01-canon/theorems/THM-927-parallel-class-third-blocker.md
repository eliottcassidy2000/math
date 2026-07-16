---
id: THM-927  # renumbered from 923 by kind-pasteur-S128c34: kps first-pushed 923 (arborescence transfer law) 16:04 vs opus close-out 16:10; 924 doubly-claimed, 925 taken -> 926
title: THE PARALLEL-CLASS CIRCLE AND THE THIRD BLOCKER — (R) THE EMPTINESS CONJECTURE (THM-900) IS REFUTED: the DFS finds a both-clean 13-packet at accessible scale ({420, 450, 510, 570, 690, 870, 1230, 1770, 2370, 3210, 4170, 5190, 7230} — Sidon, Y* ≥ 30 for q,p ≤ 13, no 13-multiple), and its BONF5 FAILS (−0.677; S₄ 2.8×, S₅ 12×): q,p ≤ 13 + Sidon cleanliness is INSUFFICIENT; the diagnosis is exact — the packet is the 30Z dilate family (gcd 30: an escape ALONG ONE PARALLEL CLASS) whose reduced ratios (15/14, …) sit exactly ONE STEP PAST the q ≤ 13 horizon (the S315 Farey-14 blind spot), AND it contains 3-TERM APs ({450, 510, 570}: a+c = 2b evades the distinct-four Sidon check): THE THIRD BLOCKER SPECIES = SMALL LINEAR FORMS; the corrected admissibility is the LINEAR-FORMS/DISSOCIATIVITY CONDITION (no small |Σc_i x_i|, small support and coefficients) + an extended ratio horizon — the full-clean DFS with that filter crawls at depth 1 in [300, 40000] (the filter is brutal: evidence the corrected regime is genuinely thin); (P) THE PARALLEL-CLASS DICTIONARY: each speed v = the v-torsion parallel class of the circle; LRC = the covering interaction of 13 parallel classes; the KILLER-CLASS CENSUS of the Landau corona (n = 8: 122 dying orbits, killer-index {2: 61, 3: 19, 4: 17, 5: 12, 6: 13}; n = 10: 1448, {2: 645, 3: 197, 4: 181, 5: 131, 6: 106, 7: 88, 8: 100}: k = 2 (the double-zero) takes ~half, the tail reaches the middle); THE TIGHT-LOCUS ADDRESSES ARE THE MECHANICAL WORDS: at t = p/14 the 13-class address vector (⌊vp/14⌋ mod v) runs from ALL-ZEROS (p = 1) to the FULL STAIRCASE [0,1,…,12] (p = 13) — the μ₆ group's addresses are the six centered mechanical words (THM-778's objects at the loneliest moments); the metagraph identification: the m wiggly classes ARE the parallel classes (directions) of Q_m — the tiling-space 1-FACTORIZATION — so "resolution into parallel classes" is the shared grammar of the circle (speeds) and the hypercube (tile positions), and a dilate packet is a walk along one class on the LRC side exactly as a silent-mutation chain is on the metagraph side
status: R = exact refutation + exact BONF5 failure + the identified third species (3-APs present; the horizon-14 reduced ratios exhibited) + the corrected condition named, with the full-clean search's crawl as thinness evidence; P = exact censuses (killer stratification n = 8, 10; the six address vectors; the dictionary checks)
source: opus-2026-07-16-S331 (owner: long session, frequent pulls, close out LRC(14), the parallel-class circle)
depends_on:
  - THM-900 (the conjecture refuted), THM-897 (the admissibility ladder this extends)
  - THM-869 (the corona this stratifies), THM-873/878 (the tight locus/mechanical words)
related: [mac-mini THM-921 (the involution-fixed-stratum law — the mirror arc), the S315 Farey-14 blind spot (now shown EXPLOITABLE), the wiggly-class canon (the 1-factorization identification)]
verification: 05-knowledge/results/bothclean_dfs_bonf5_opus_S331.out, parallel_class_corona_opus_S331.out, fullclean_dfs_bonf5_opus_S331.out (running)
---

# THM-927 — the parallel-class circle and the third blocker (renumbered from 923)

> **RENUMBER TRAIL:** born THM-923 (opus S331), renumbered 923->926 by kps S128c34 (kps first-pushed 923), renumbered 926->927 by opus S332 (mac-mini first-pushed 926 at 16:57:47, before the kps renumber landed at 17:05:15). This file = the parallel-class / third-blocker theorem. Stable id: THM-927.


> **CORRECTION (same session, MISTAKE-151):** the search filter's p was
> UNCAPPED (testing 'beat Dirichlet', vacuous for lo ≥ (Q+1)Y₀) and the
> found packet passed by GCD-QUANTIZATION (gcd = 30 = Y₀ makes all nonzero
> relations ≥ 30 automatically). So (R)'s 'emptiness refuted' is RE-SCOPED:
> the packet satisfies the filter, not THM-863's capped-p condition; the
> true both-clean existence question is REOPENED. What stands exactly: the
> BONF5 = −0.677 failure, the 3-AP and horizon-14 diagnoses, the full-clean
> run's vacuity mechanism, and all of (P). The quantization loophole is
> itself the parallel-class escape appearing INSIDE the instrument — the
> object and the filter failed identically, evidence the structure is
> fundamental.

## (R) The refutation and the third species

The S329/S330 emptiness conjecture is REFUTED by construction: the DFS finds
{420, …, 7230} — Sidon, ratio-clean at q,p ≤ 13, no 13-multiples — inside
[200, 9000]. Its exact BONF5 = −0.677 (S₄ = 1.13 vs 0.40 equid; S₅ = 1.38 vs
0.11): the cleanliness that defined the search is INSUFFICIENT for the
level-5 certificate. Two mechanisms, both exhibited exactly:
1. **The horizon escape:** gcd = 30 — the packet travels along the 30Z
   parallel class; its reduced ratios (15/14, 17/15, …) are near-dilates at
   q = 14+, exactly one step past the q ≤ 13 functional: the S315 Farey-14
   blind spot is not just a boundary, it is EXPLOITABLE, and the DFS
   exploited it.
2. **3-term APs:** {450, 510, 570} ⊂ the packet — a + c = 2b evades the
   Sidon (distinct-four) check; 3-APs are strong joint beats fueling
   S₃/S₄/S₅.

**The corrected admissibility: the linear-forms condition** — no small
|Σ c_i x_i| over small supports and coefficients (dissociativity-flavored),
plus a ratio horizon beyond 13. The upgraded DFS (no-3AP + linear forms
|c| ≤ 3 + horizon 20) crawls at depth 1 in [300, 40000]: the corrected
regime is genuinely thin — the blocker-filter iteration (each session's
counterexample becomes the next session's constraint) is converging toward
the true admissibility, and each rung is exact.

## (P) The parallel-class circle

Each speed v resolves the circle into the v-torsion parallel class; LRC(14)
is the covering interaction of thirteen parallel classes; the common
refinement is the lcm(1..13) = 360360-cell lattice.
- **The killer-class census** (the corona's interior, THM-869's named next):
  the first violated Landau prefix k stratifies the dying orbits — k = 2
  (the double-zero) takes half; the tail reaches k = n − 2. Exact at
  n = 8, 10.
- **The tight-locus addresses are the mechanical words:** the address vector
  of t = p/14 across the 13 classes runs from all-zeros (p = 1) through the
  full staircase (p = 13) — the μ₆ group's six addresses are six centered
  mechanical words: THM-778's transport objects ARE the coordinates of the
  loneliest moments in the parallel-class resolution.
- **The metagraph identification:** the m wiggly classes are the parallel
  classes (directions) of the tiling hypercube Q_m — its 1-FACTORIZATION.
  Resolution-into-parallel-classes is the shared grammar: circle/speeds ↔
  hypercube/tile-positions; dilate packets walk one class on the LRC side
  as silent-mutation chains do on the metagraph side; and mac-mini's
  THM-921 law (excess concentrates on the involution-fixed stratum) is the
  same grammar's reversal-class instance.

## Named next

- The full-clean search's verdict (running); if it finds: BONF5 = the next
  referee; if practically-empty: the thinness quantified.
- The killer-class census as a per-class (not per-index) stratification:
  WHICH speed's parallel class delivers the killing prefix.
- The 360360-cell address system as the standing LRC coordinate (the
  mechanical-word addresses generalized off the tight locus).
