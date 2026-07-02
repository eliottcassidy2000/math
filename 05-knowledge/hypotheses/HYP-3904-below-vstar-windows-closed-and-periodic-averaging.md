# HYP-3904: LEAN SORRY-FREE: below-V* windows CLOSED for all three cert shapes (no-V-gap _all theorems, kernel-pure) + the GENERAL PERIODIC AVERAGING lemma (row-7 block reduction)

**Status:** CONFIRMED (builds green; windows = STANDARD AXIOMS + kernel decide only, no native_decide; general lemma axiom-clean)
**Instance:** opus-2026-07-02-S35
**Files:** `LRCWitnessWindow.lean` (new, 53 generated rows), `LRCCommensuration.lean` (refactored),
`04-computation/lrc14_window_witness_generator_opus_20260702_S35.py` (+ .out)

## 1. Below-V* finite windows: item 1 of the 3-item DAG surface, CLOSED

For each LRCWitnessCert shape, every V below the tail threshold now has an explicit rational
witness a/b checked by KERNEL decide through the closed Nat criterion
b <= 14*min((v*a)%b, b-(v*a)%b)  (the exact form of ||v*(a/b)|| >= 1/14 via norm_natCast_div):
  - certAP_window: V in [13,15) (2 rows; V=13 uses the tight witness 1/14 -- closed inequality!)
  - cert3_window: V in [21,47) (26 rows) ; cert4_window: V in [15,40) (25 rows)
COMBINED _all THEOREMS (window+tail glue): certAP_all (every V>=13), cert3_all (V>=21),
cert4_all (V>=15) -- each certified 13-runner family now has NO V-GAP at all. The windows are
kernel-pure [propext, Classical.choice, Quot.sound]; only the tails inherit kps's native_decide.
Witness tables: exact max-min over the complete breakpoint grid (53/53 pass, generator+.out).

## 2. The row-7 extension: the general periodic averaging lemma

`seven_mul_volume_inter_periodic`: for 7 not| P and ANY 1/7-periodic measurable C:
    7 * vol(danger P psi (1/14) CAP C) = vol C.
The S34 tiling proof never used the second comb beyond periodicity -- abstracting C turns the
7-commensuration lemma into the BLOCK-REDUCTION step the forced-independence recursion (DAG row
7) consumes: an intersection of 7-divisible danger combs is 1/7-periodic (seventh_periodic_danger
+ seventh_periodic_inter provided), so each non-commensurate speed peels an exact 1/7 factor.
seven_commensuration is now a 10-line corollary. New helpers: zsmul_seventh_eq_zero,
mem_add_tau_iff (invariance under one seventh-turn extends to all translates).

## Honest scope

Items 2 (shape-universe enumeration) and 3 (multi-level recursion glue) of the DAG surface were
NOT reached this session: universe enumeration is in active flux on mac-mini's side (THM-604
abolished the deviation tables mid-session; drift-line dichotomy + exact-ray sieve evolving), and
the depth-d ladder wrapper remains with kps's cert_two_level as its base case. The windows were
chosen first as the purely-decidable, collision-free item; the averaging lemma as the consumable
extension of my own S34 module.

-> HYP-3903, HYP-3958/3959 (kps cert engine), HYP-3859 (spec), THM-602/604, klein HYP-4004/4006.
