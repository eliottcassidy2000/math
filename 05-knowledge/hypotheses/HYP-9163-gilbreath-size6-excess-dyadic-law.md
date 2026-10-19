---
id: HYP-9163
title: "Dyadic law of the size-6 extinction excess in Gilbreath's automaton: p_6(F) = F 2^(1-F) + 2^-F c_6(F), where c_6(F) = sum_z c_z(F) over the number z of zeros before the wall, c_0(F) = (2/3)(1 - 4^-(2^k-1)) and c_1(F) = c_2(F) = (2/15)(1 - 16^-(2^(k-1)-1)) for 2^k < F <= 2^(k+1), c_(z+2^m)(F) = c_z(F - 2^m) while fewer than 2^m hit rows matter, and c_6(F) ~ 1.4 log_2 F; hence p_6(F)/(F 2^(1-F)) -> 1 like 1 + 0.7 log_2 F / F"
status: >
  PARTLY REFUTED, PARTLY OPEN; updated 2026-10-19. The logarithmic excess
  c6(F)=1.4 log2(F)+O(1) and the resulting ratio p6(F)/(F 2^(1-F))->1
  are REFUTED by an all-F lower bound. Dyadic rationality is now PROVED
  for every F by the alphabet-three shield and zero-run formula.
  The seven fixed-age counting formulas for ages 2,4,6,8,10,12,14 are
  PROVED by independently checked finite-state certificates. All-F
  paired halving, the historical per-z formulas, the complete all-age
  recurrence, and the precise asymptotic remain OPEN. Historical exact
  rationals through F=11 survive; later floating values are superseded
  by exact tables through F=23.
source: opus-2026-09-26 session gilbreath6-collatz-precision-20260926
related: [THM-4511 (wall theorem, size 4), S10 note gilbreath_fermat_platonic_20260926.md (finite-context table), S9 note (Lucas kernels)]
verification: 04-computation/experiments/gilbreath_size6_exact_20260926.py -> .out
---

# HYP-9163 -- the size-6 excess is a dyadic tower

> **Current-truth correction, 2026-10-19.** Item 4 below is false. For every F, the exact number of binary left words whose first possible failure occurs after exactly two appended zeros is `floor((F+1)/4)`. Thus `c6(F) >= (1/2)floor((F+1)/4)` and `liminf p6(F)/(F 2^(1-F)) >= 17/16 > 1`. Seven independently verified fixed-age formulas strengthen the latter bound to `285677/262144 = 1.089771270751953125`. This refutes both logarithmic growth of c6 and asymptotic exactness of the front-only probability.
>
> **What survives.** The finite absorbing-chain formulation and the historical exact rationals through F=11 survive. A new shield theorem proves dyadic rationality for every F and reduces the full risk to one finite-state zero-run threshold. The historical all-F paired-halving and per-z statements remain open; they are not proved or refuted by the new fixed-age formulas. The stronger unproved replacement is linear growth of c6, not logarithmic growth. The new age variable counts right-sea zeros; the historical z counts zeros before the left wall. These are different decompositions and must not be identified without a map.
>
> **Numerical correction.** In particular, `c6(13)=19437/8192=2.3726806640625`, replacing the old floating estimate 2.372697. Exact Markov-chain and zero-run computations agree through F=15; the independently cycle-certified zero-run computation extends through F=23. The historical floating program is absent from this checkout, so the cause of that numerical discrepancy is not established.
>
> See [the complete proof and exact certificates](../results/gilbreath_alphabet_three_size_six_20261019.md). The original conjecture below is retained as historical context.

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
