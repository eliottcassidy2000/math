---
id: THM-4793
title: "The Fermat line under Collatz: for n >= 2 the Syracuse orbit of F_n = 2^(2^n)+1 begins with 2^(n-1)-1 steps of valuation two and one of valuation three, with S^k(F_n) = 3^k 2^(2^n - 2k) + 1 and landing point (3^(2^(n-1))+1)/2 (the dual of THM-4556's valuation-one Mersenne line; the valuation-two run count is THM-4601's J = floor((m-1)/2) for 2^m+1); the orbits of 2^m+1 and 2^(m+1)+1, m = 2L+1, first merge after exactly L+2 odd steps when L is even and exactly L+3 when L = 3 mod 4, by the classical 4x+1 successor rule applied to the landing points 2·3^L+1 and 4·3^L+1 (the Fermat twin of THM-4556 (ii), an explicit THM-4581 pair-chain merge); and Pepin's test reads: F_n is prime iff the affine Collatz branch x -> 3x+1 on Z/F_n has exactly two cycles, the fixed point 2^(2^n - 1) and one cycle of length 2^(2^n). The line is prime-blind: F_4 and the composite F_5 follow the same law."
status: >
  PROVED (klein, 2026-10-25; elementary, proofs in the source note). FINITE-EXACT: the run law for all odd x < 20000 and
  L <= 19, the line for F_2..F_5 and 2^m+1 for m = 2..40, the exact merge depths for L = 2..19 (asserted in
  check_merge_points.py), the affine cycle types for n = 1..5 (order of 3 modulo F_5 is 11167360), and the Paley
  spectra at 5 and 17. The survey (2026-10-25) found no prior analysis of 2^m+1 or Fermat numbers as a Collatz family;
  the ingredients THM-4601 (valuation-two counts), THM-4556 (Mersenne line and its equal-time merge), THM-4581
  (pair-chain coalescence) and the 4x+1 rule are inherited. Not independently audited. Completeness of the five Fermat
  primes is OPEN and not addressed by this theorem.
session: klein-2026-10-25 (codex/session-tournament-twist-unify)
source: 05-knowledge/results/fermat_five_pepin_collatz_20261025/PEPIN_COLLATZ.md
scripts:
  - 05-knowledge/results/fermat_five_pepin_collatz_20261025/check_fermat_line.py (+ fermat-normal.json / fermat-optimized.json)
  - 05-knowledge/results/fermat_five_pepin_collatz_20261025/check_merge_points.py (+ merge-normal.json / merge-optimized.json)
related:
  - 01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md
  - 01-canon/theorems/THM-4601-two-anchor-reduction-of-the-residual-first-reset-two-branch.md
  - 01-canon/theorems/THM-4581-haar-coalescence.md
  - 01-canon/theorems/THM-4749-the-repunit-degree-law-phi-2-of-2-is-the-triangle-phi-3-of-2-is-the-gap-and-fermat-products-reverse-the-order.md
  - 01-canon/theorems/THM-4644-fermat-primes-rigid-circulant-tournaments-wieferich-orbit-explosion.md
  - 05-knowledge/results/fermat_five_pepin_collatz_20261025/FIVE_LADDER.md
---

# THM-4793 — The Fermat line under Collatz

## Statements

Let `S(x) = (3x+1)/2^v(x)` on odd `x`, `v(x) = v_2(3x+1)`.

1. **Run law.** The first `L` Syracuse steps of odd `x` all have valuation two iff `x ≡ 1 (mod 2^(2L+1))`,
   and then `S^k(x) = 3^k (x-1)/4^k + 1` for `k <= L`.
2. **Fermat line.** For `n >= 2`: `F_n` has `2^(n-1) - 1` valuation-two steps, then one valuation-three
   step, landing on the odd number `(3^(2^(n-1)) + 1)/2`. For general `2^m + 1` the run has
   `floor((m-1)/2)` steps.
3. **Exact merge.** For `m = 2L+1`, `L >= 2`, the orbits of `2^m+1` and `2^(m+1)+1` land on
   `a = 2·3^L + 1`, `b = 4·3^L + 1`; `S(a) = 3^(L+1)+2` and, for even `L`, `S(b) = (3^(L+1)+1)/4`, so
   `S(a) = 4 S(b) + 1` and `S^2(a) = S^2(b)`: first coincidence at odd index exactly `L + 2`. For
   `L ≡ 3 (mod 4)` the `4x+1` relation appears one step later (index `L + 3`); for `L ≡ 1 (mod 4)` not
   within `L + 3` steps.
4. **Pépin, affine form.** For `n >= 1`, `F_n` is prime iff `x -> 3x+1` on `Z/F_n` has cycle type
   `(1, 2^(2^n))`: it is conjugate to `x -> 3x` by `x -> x + 1/2`, and `3` generates `(Z/F_n)^x` iff
   `3^((F_n-1)/2) ≡ -1`.

## Proofs

In the source note, §§1–2a. The run law is the induction `x = 1 + 2^(2L+1)s ⟹ S(x) = 1 + 3·2^(2L-1)s`;
the landing valuation is three because `3^(2^(n-1)) ≡ 1 (mod 8)`; the merge is the identity
`3(4q+1)+1 = 4(3q+1)`; the affine form is the translation conjugacy plus Pépin.

## Mechanism and boundary

Valuations one and two bracket the Collatz slope `log_2 3`: the Mersenne residue `-1 (mod 2^p)` buys
the longest pure ascent, the Fermat residue `+1 (mod 2^(2^n))` the longest pure descent. The line
carries no primality information (hostile: `F_5`), which is the exact boundary of "Pépin for Collatz":
primality lives in the multiplicative tower `3^(2^k) mod F_n`, the line in the additive tower
`3^k 2^(2^n-2k) + 1`, and only the base `3` and the power-of-two certificate are shared.
