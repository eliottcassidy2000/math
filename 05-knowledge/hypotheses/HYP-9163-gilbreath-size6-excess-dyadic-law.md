---
id: HYP-9163
title: "Dyadic law of the size-6 extinction excess in Gilbreath's automaton: p_6(F) = F 2^(1-F) + 2^-F c_6(F), where c_6(F) = sum_z c_z(F) over the number z of zeros before the wall, c_0(F) = (2/3)(1 - 4^-(2^k-1)) and c_1(F) = c_2(F) = (2/15)(1 - 16^-(2^(k-1)-1)) for 2^k < F <= 2^(k+1), c_(z+2^m)(F) = c_z(F - 2^m) while fewer than 2^m hit rows matter, and c_6(F) ~ 1.4 log_2 F; hence p_6(F)/(F 2^(1-F)) -> 1 like 1 + 0.7 log_2 F / F"
status: >
  FINITE-EXACT: exact rationals for F <= 11 (absorbing Markov chain, Theorem A of the
  note), floating values to F = 17, exact per-z decomposition for F <= 11 and floating
  to F = 16; the closed forms for z = 0, 1, 2 match every computed F, the shift law
  c_(z+8) = c_z(F-8) holds for z = 0, 1, 2, 4 and F <= 16, the halving from odd F to
  F + 1 holds for all computed pairs. CONJECTURE beyond that range; mechanism identified
  (defect-column phases, Lucas hit toggling of the wall, released 4s every two or four
  rows, dyadic refresh of the left region) but not turned into a proof.
source: opus-2026-09-26 session gilbreath6-collatz-precision-20260926
related: [THM-4511 (wall theorem, size 4), S10 note gilbreath_fermat_platonic_20260926.md (finite-context table), S9 note (Lucas kernels)]
verification: 04-computation/experiments/gilbreath_size6_exact_20260926.py -> .out
---

# HYP-9163 -- the size-6 excess is a dyadic tower

See `05-knowledge/results/gilbreath_size6_extinction_20260926.md`, sections 2-4.

**Statement.** With `E(F) = p_6(F) - F 2^(1-F)` and `c_6(F) = 2^F E(F)`:

1. `E(F + 1) = E(F)/2` for every odd `F >= 3`.
2. `c_6(F) = sum_(z >= 0) c_z(F)` with `c_z(F) = 0` for `z > F - 3`, and
   `c_0(F) = (2/3)(1 - 4^-(2^k - 1))`, `c_1(F) = c_2(F) = (2/15)(1 - 16^-(2^(k-1) - 1))`
   for `2^k < F <= 2^(k+1)`, `k >= 2` (`c_0 = 1/2` for `F = 3, 4`).
3. `c_(z + 2^m)(F) = c_z(F - 2^m)` whenever `F - 2^m <= 2^m` (the hit
   pattern of `z + 2^m` coincides with that of `z` for `2^m` rows).
4. `c_6(F) = 1.4 log_2 F + O(1)`; consequently `p_6(F)/(F 2^(1-F)) = 1 +
   (0.7 + o(1)) log_2 F / F`.

**Evidence.** Exact: `c_6(F) = 1/2, 1/2, 29/32, 29/32, 45/32, 45/32,
12013/8192, 12013/8192, 16109/8192` (`F = 3..11`); per-`z` exact values in the
note; floating `c_6 = 1.9664, 2.3727, 2.3727, 2.8727, 2.8727, 2.8732` for
`F = 12..17`.

**What a proof needs.** The defect column's three phases have geometric
lengths; in the `6`-phase the wall column toggles at the Lucas hits
`z subset s`; in the `4`-phase the wall alternates `0/4` and releases a `4`
every two rows; each released `4` is destroyed iff its path through the
evolved sea is zero, a unit-triangular condition on the left word *at the
release time*; the dyadic truncation is the point where that path leaves
the part of the region determined by the initial word. Summing the
release probabilities over release times gives the geometric series; the
bookkeeping of which release times are compatible with which words is the
missing step.
