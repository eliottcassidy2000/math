---
id: THM-4517
title: "No root-uniform positive lower density of Collatz basins exists: the basins B(a_i) of the trunk-entry points a_i = (4^i - 1)/3 are pairwise disjoint, so no c > 0 bounds the lower density of every predecessor set B(a), a != 0 (mod 3), from below, and Krasikov-Lagarias X^0.84 is the only uniform type; a sheet-blind cycle-uniform lower bound is at most the smallest 3x-1 basin; the 3x-1 basins on [1,2^29] are 0.3269 / 0.3248 / 0.3484 and on [1,2^30] 93.8% of 3x+1 orbits enter the trunk through 5, with e_i = 0 exactly when 3 | i"
status: >
  PROVED (elementary, unconditional) for (1)-(2) + FINITE-EXACT for (3) + CITED
  (Krasikov-Lagarias 2003, Barina 2^68). NOT independently audited. Setting: T-form
  3x+1 on positive integers, B(a) = {n : T^k(n) = a for some k >= 0}, a_i = (4^i - 1)/3.
  (1) The sets B(a_i), i >= 2, are pairwise disjoint (after a_i the orbit is
  2^(2i-1), ..., 2, 1, 2, ... and contains no a_(i') >= 5); a_i is prime to 3 iff
  3 does not divide i, and for 3 | i the basin is the doubling ray of a_i (density 0).
  Hence there is no c > 0 with liminf |B(a) cap [1,X]|/X >= c for every root
  a != 0 (mod 3): M > 1/c disjoint basins of lower density >= c would give a union of
  lower density > 1. A root-uniform statement can only be sublinear, as X^0.84 is.
  (2) The pasted argument "basin density > 1/2 for every root would prove Collatz" has
  a void premise; even dens B({1,2}) > 1/2 alone shows only that every other cycle basin
  and divergent class has upper density < 1/2. A cycle-uniform lower bound c proved
  by a sheet-blind method (invariant under n -> -n, THM in the counterexample portrait)
  must hold for the three disjoint positive cycle basins of 3x-1, so c <= 1/3, and
  c <= min of their densities if these exist.
  (3) FINITE-EXACT. 3x-1 on [1, 2^29]: basins of {1}, {5,7,10}, {17..91} have counting
  densities 0.3268569, 0.3247598, 0.3483833 (top dyadic range 0.326787, 0.324876,
  0.348337); to 2^32 (2-bit sieve): 0.32676070, 0.32499135, 0.34824796, dyadic rows
  [2^28..2^32): {1} flat at 0.32675, {5,7,10} rising 5e-5 per doubling (0.3248756 ->
  0.3250474), {17..} falling 5e-5 per doubling (0.3483374 -> 0.3482025). 3x+1 on [1, 2^30], first trunk
  element 2^(2i-1) entered from a_i: e_2 = 0.93796, e_4 = 0.023647, e_5 = 0.037789,
  e_7 = 8.0e-5, e_8 = 4.85e-4, e_10 = 3.2e-5, e_11 = 2.1e-6, e_13 = 4e-7, e_14 = 1e-7,
  e_3, e_6, e_9, e_12 = 0 to seven decimals (their basins are the doubling rays of 21, 1365, 87381, 5592405: 26, 20, 14, 8 elements below 2^30) [corrected by the independent audit `collatz_necklace_20260929_audit.md`, opus S23, 2026-09-29]; so dens B(5) ~ 0.938, dens B(32) ~ 0.062, and 341
  owns a larger basin than 85. Ray-entry spectrum of the minus sheet in the note.
  NOT claimed: existence of any of these densities (HYP-9165); anything about
  divergent orbits; that the 3x-1 basins cover everything (their sum 0.99999+ is a
  count below 2^29 only).
source: collatz-necklace-20260929 session (mac-mini), 2026-09-29; owner seed: the pasted Krasikov-Lagarias / 1729 / 3x-1 argument
depends_on: []
related:
  - 05-knowledge/results/collatz_mod6_20260922_counterexample_portrait.md (sheet criterion T_+(-n) = -T_-(n); sheet-blind conditions)
  - 05-knowledge/results/collatz_procgen_20260924_inverse_tree_mod192.md (basins of 3x-1 seen by residue types)
  - 05-knowledge/results/mazur_positive_density_20260928.md (harmonic mass of predecessor trees; Krasikov-Lagarias row of the barrier atlas)
  - 05-knowledge/hypotheses/HYP-9165-basin-densities-exist-sheet-cap.md
note: 05-knowledge/results/collatz_necklace_20260929_fair_splits_power_clocks_basins.md
scripts:
  - 04-computation/experiments/collatz_necklace_20260929_basins_minus.c
  - 04-computation/experiments/collatz_necklace_20260929_entries_plus.c
outputs:
  - 05-knowledge/results/collatz_necklace_20260929_basins_minus_2e29.out
  - 05-knowledge/results/collatz_necklace_20260929_entries_plus_2e30.out
  - 05-knowledge/results/collatz_necklace_20260929_basins_minus_2e32.out (script 04-computation/experiments/collatz_necklace_20260929_basins_minus_2bit.c)
extra_sha256: 2d947096fd356fe3335523f5035d04f365a186fba633485de5250d2f83aab3ff (basins_minus_2bit.c), cfac7d1933da9bb2f63131cb6d74633a306a8bf32dbc2606950e666970f06304 (basins_minus_2e32.out)
script_sha256: f2d9f0326d206f969c41e85b9b4b5cb52810a2c0466a94555a86103f676ae1ad (basins_minus.c), 06094371c3c6d7c824fc37c06252d51a674688b37e30c3edd76703b6844d41e9 (entries_plus.c)
output_sha256: f86129b8e0e8c69b9a9d760138ab09cb0d876171f0c23dd1d4e9bc2a63a485fb (basins_minus_2e29), 744ba1caa9aaef711ebea41962e234c4f42a30f7cef3ac7970ae088cf7a0e019 (entries_plus_2e30)
hash_basis: raw LF bytes
audit: NOT independently audited; the two sieves memoise on the first orbit value below the start, use 128-bit orbits, and their small-range rows agree with the 10^6 test runs; the counts for 3 | i equal the ray lengths 26, 20, 14, 8, an exact check of (1) [corrected by the independent audit `collatz_necklace_20260929_audit.md`, opus S23, 2026-09-29].
---

# THM-4517 -- no root-uniform positive density of Collatz basins; the sheet cap; the trunk-entry spectrum

**PROVED (1)-(2), FINITE-EXACT (3); not independently audited.** Full note: [collatz_necklace_20260929_fair_splits_power_clocks_basins](../../05-knowledge/results/collatz_necklace_20260929_fair_splits_power_clocks_basins.md), section 3.

## 1. Disjoint trunk-entry basins

`T(n) = n/2` (`n` even), `(3n+1)/2` (`n` odd). `B(a) = {n >= 1 : T^k(n) = a, k >= 0}`. The trunk is the set of powers of two; an orbit that is not on the trunk enters it only through `a_i = (4^i - 1)/3 -> T(a_i) = 2^{2i-1}` (`i >= 2`; `(2^k - 1)/3` is an integer iff `k` is even). After `a_i` the orbit is `2^{2i-1}, 2^{2i-2}, ..., 2, 1, 2, 1, ...`, which contains no `a_{i'} >= 5`; so no orbit contains two distinct `a_i`, and the `B(a_i)` are pairwise disjoint. `a_i = 0 (mod 3)` iff `4^i = 1 (mod 9)` iff `3 | i`; a multiple of three has no odd preimage, so then `B(a_i)` is the doubling ray of `a_i`, of density zero. For `3` not dividing `i` the `a_i` are odd roots prime to `3`, the roots to which Krasikov–Lagarias's `|B(a) cap [1,X]| >= X^{0.84}` applies.

**Proposition.** There is no `c > 0` such that `liminf_X |B(a) cap [1,X]|/X >= c` for every root `a != 0 (mod 3)`.

*Proof.* Lower density is superadditive on disjoint sets. `M > 1/c` of the disjoint basins `B(a_i)`, `3` not dividing `i`, would have a union of lower density `>= Mc > 1`. ∎

So a root-uniform statement about predecessor sets can only be sublinear; `X^{0.84}` is of the only admissible type, and the premise "for any root `r` the basin has proportion `c > 1/2`" of the pasted argument is void before any sheet is invoked.

## 2. What density bounds can and cannot prove

* `dens B({1,2}) > 1/2` alone gives: every other cycle basin and every divergent class has upper density `< 1/2`. It does not exclude them.
* A bound `dens B(C) >= c` for *every* positive cycle `C` gives "at most `1/c` cycles with positive-density basins", nothing about density-zero cycles or divergent orbits.
* A sheet-blind proof of such a cycle-uniform bound (one invariant under `n -> -n`, which conjugates `3x+1` to `3x-1`: `T_+(-n) = -T_-(n)`) must also hold for the three disjoint positive cycle basins of `3x-1`, so `c <= 1/3`; and `c <=` the smallest of their densities when these exist — `0.3248` on the evidence of (3). This is the repo's SHEET control with a number attached; the pasted "cannot prove more than `1/3`" is right in spirit and off by the inequality between the three basins.

## 3. FINITE-EXACT densities

`3x-1` on `[1, 2^29]` (`T_-(n) = (3n-1)/2` for odd `n`; the three known cycles `{1}`, `{5,7,10}`, `{17,25,37,55,82,41,61,91,136,68,34}`): counting densities `0.3268569 / 0.3247598 / 0.3483833`; top dyadic range `[2^28, 2^29)`: `0.326787 / 0.324876 / 0.348337`; the `{1}` basin drifts down from `0.3277` at `2^20`, the `{5,7,10}` basin up from `0.3243`, the seven-odd-step cycle's basin is flat at `0.348` to three decimals from `2^20` (rows `0.34804 .. 0.34854`) [corrected by the independent audit `collatz_necklace_20260929_audit.md`, opus S23, 2026-09-29]. To `2^32` (2-bit sieve): cumulative `0.32676070 / 0.32499135 / 0.34824796`; dyadic rows `[2^28..2^32)`: `{1}` flat at `0.32675`, `{5,7,10}` rising `5e-5` per doubling (`0.3248756 -> 0.3250474`), `{17,...}` falling at the same rate (`0.3483374 -> 0.3482025`). The three basins are not equal, and the five-cycle basin is still gaining from the seven-odd-step basin at `2^32`. First ray element hit (with the odd number entering it): `7 x 2^2` (from 19) `0.1984`, `61 x 2^2` (163) `0.1914`, `2^6` (43) `0.1878`, `2^4` (11) `0.1336`, `7 x 2^8` (1195) `0.0526`, `17 x 2` (23) `0.0460`, `25 x 2^2` (67) `0.0433`, `37 x 2^4` (395) `0.0335`, `5 x 2^5` (107) `0.0282`, `5 x 2^7` (427) `0.0264`; the ray of `1` is never entered at `2^2, 2^8, 2^14` (`(2^{k+1}+1)/3` a multiple of three).

`3x+1` on `[1, 2^30]`, `e_i = ` counting density of `B(a_i)`: `e_2 = 0.93796` (`5`), `e_4 = 0.023647` (`85`), `e_5 = 0.037789` (`341`), `e_7 = 8.0e-5`, `e_8 = 4.85e-4`, `e_10 = 3.2e-5`, `e_11 = 2.1e-6`, `e_13 = 4e-7`, `e_14 = 1e-7`, and `e_3 = e_6 = e_9 = e_12 = 0` exactly. Top range `[2^29, 2^30)`: `0.93794 / 0.02366 / 0.03780`. Hence `dens B(5) ~ 0.938`, `dens B(32) ~ 0.062` (the two children of `16` partition `B(1)` minus the trunk; as "everything but the trunk" the claim would be the Collatz conjecture) [corrected by the independent audit `collatz_necklace_20260929_audit.md`, opus S23, 2026-09-29], and `341` owns a larger basin than `85`; the tail beyond `i = 11` is not converged (the basin of `a_i` is invisible below `~ 4^i`).

## 4. Dual-graph reading and boundary

In the dual of the planar unicyclic component the trees are loops at the outer face; these entry statistics are the loops, invisible to the `K`-bond that carries the necklace (THM-4515). Existence of the densities is HYP-9165. Nothing is claimed about divergent orbits or about the completeness of the `3x-1` cycle list.
