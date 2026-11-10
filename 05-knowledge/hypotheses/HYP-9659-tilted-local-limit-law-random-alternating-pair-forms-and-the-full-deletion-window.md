---
id: HYP-9659
title: "Three falsifiable steps of a local-limit proof of the {7, 21} spectrum. (A) Tilted local limit law: there are n_0 and c > 0 such that for n >= n_0, every tilt p in [1/2, p_n] and every odd h within one standard deviation of E_p[H_n], P_p(H_n = h) >= c / sd_p(H_n); this implies that every odd h > (n_0 - 1)!/2^(n_0 - 2) is a Hamiltonian path count. (B) Random-alternating law: for R uniform on m vertices the rank of beta_R = (mu_R / 2 mod 2) is distributed as that of a uniformly random alternating matrix over F_2 up to o(1), so the mod-4 discrepancy of H at order m+1 is O(2^(-m/2)). (C) Window law: for every odd h >= 237 that is a path count, every odd n in [3h, 5h) is H(T + z) for some tournament T with H(T) = h; this implies the Chan-Pak conjecture."
status: >
  OPEN.
  (A): no evidence beyond the bulk measurements of THM-5127 (order 10: at least 1,504 classes for every odd h within
  one class-standard-deviation of the class mean; near-Poisson local fluctuation; residues balanced to 1.2%). The
  implication is PROVED (source note, Proposition 5.1).
  (B): FINITE-EXACT agreement within 1.2 percentage points at orders 8 and 9 (all 6,880 and 191,536 classes).
  (C): FINITE-EXACT for 237 <= h <= 1023 with bases of order <= 10 and all out-sets. The implication is PROVED
  (source note, Proposition 6.1). (C) is the all-out-sets form of klein's HYP-9453 (2).
session: mac-mini-2026-11-10-ensemble
source: 05-knowledge/results/tilted_ensemble_local_limit_20261110/README.md
related:
  - THM-5127
  - HYP-9453
  - HYP-9474
  - HYP-9640
  - THM-4885
  - THM-4876
  - THM-5124
---

# HYP-9659 - a local-limit route to the {7, 21} spectrum

Definitions and proofs of the implications are in the
[source note](../results/tilted_ensemble_local_limit_20261110/README.md) and in
[THM-5127](../../01-canon/theorems/THM-5127-heat-flow-reading-of-the-w-polynomial-tilted-ensemble-mean-and-the-mod-four-quadratic-form-law-of-a-vertex-cube.md).

## (A) Tilted local limit law

**Statement.** Let `P_p` orient each pair `i < j` as `i -> j` with probability `p`, independently, and let `p_n` be
the tilt with `E_{p_n}[H_n] = (n-1)!/2^(n-2)`. There are `n_0` and `c > 0` such that for all `n >= n_0`, all
`p in [1/2, p_n]` and all odd `h` with `|h - E_p[H_n]| <= sd_p(H_n)`: `P_p(H_n = h) >= c / sd_p(H_n)`.

**Consequence (PROVED).** Every odd `h > (n_0 - 1)!/2^(n_0 - 2)` is a Hamiltonian path count.

**Why it is the right regime.** Over all tournaments the number of structures with count at most `X` is
super-polynomial in `X`, and by THM-5127 (B) every target is the mean of exactly one tilt at each admissible order.
The word families studied so far have bounded, logarithmic or polynomial multiplicity.

**Falsifier.** An order `n` and an odd `h` in the stated window with `P_p(H_n = h) = 0`, for `n` beyond any proposed
`n_0`; or, at finite range, class counts in the bulk at order 11 that are not smooth.

**Route.** Characteristic function of `H_n` under `P_p`: Gaussian regime from the Walsh levels (THM-071, THM-589);
rational frequencies from THM-5127 (C), (C4) and part (B) below; the remaining range needs a decoupling over several
exposed vertices and is untouched.

## (B) Random-alternating law for the pair form

**Statement.** For `R` uniform on `m` labelled vertices, the distribution of `rank(beta_R)` is within `o(1)` in total
variation of the rank distribution of a uniformly random alternating `m x m` matrix over `F_2`. In particular
`E_R[2^(-rank(beta_R)/2)] = O(2^(-m/2))`, and by THM-5127 (C) (4) the mod-4 discrepancy of `H` at order `m+1` is
`O(2^(-m/2))`.

**Evidence.** Orders 8 and 9, all classes: `0.4183 / 0.5516 / 0.0301` against `0.4205 / 0.5585 / 0.0210`, and
`0.8278 / 0.1699 / 0.0023` against `0.8394 / 0.1593 / 0.0014` (top rank / next / lower). Low ranks are
over-represented by structured classes.

**Falsifier.** A rank distribution at order 10 or 11 that moves away from the random-alternating one.

## (C) Window law

**Statement.** For every odd `h >= 237` that is a Hamiltonian path count and every odd `n` with `3h <= n < 5h` there
is a tournament `T` with `H(T) = h` and a one-vertex extension `T + z` with `H(T + z) = n`.

**Consequence (PROVED).** Every odd number other than 7 and 21 is a Hamiltonian path count.

**Evidence.** Exact for `237 <= h <= 1023` with bases of order at most 10; with the same bases the window `[3h, 9h)` is
complete for `365 <= h <= 1023`. Below 237 there are 48 exceptional `h` (largest 235).

**Falsifier.** An odd `h >= 237` and an odd `n` in `[3h, 5h)` that is not an extension of any realiser of `h`.

**Limits.** By THM-4876 the witnesses cannot come from one realiser; by THM-5124 they cannot come from additive
presence cubes. No structural control of the extension sets is known; the statement is stronger than the conjecture.
