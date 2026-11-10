---
id: HYP-9659
title: "Three falsifiable steps of a local-limit proof of the {7, 21} spectrum. (A) Tilted local limit law: there are n_0 >= 4 and c > 0 such that for n >= n_0, every tilt p in [1/2, p_n] and every odd h within one standard deviation of E_p[H_n], P_p(H_n = h) >= c / sd_p(H_n); this implies that every odd h > (n_0 - 1)!/2^(n_0 - 2) is a Hamiltonian path count. (B) Random-alternating law: for R uniform on m vertices the rank of beta_R = (mu_R / 2 mod 2) is distributed as that of a uniformly random alternating matrix over F_2 up to o(1); (B') quantitatively E_R[2^(-rank(beta_R)/2)] <= C 2^(-m/2), which gives a mod-4 discrepancy O(2^(-m/2)) at order m+1. (C) Window law: for every odd h >= 237 that is a path count, every odd n in [3h, 5h) is H(T + z) for some tournament T with H(T) = h; this implies the Chan-Pak conjecture."
status: >
  OPEN.
  (A): no evidence beyond the bulk measurements of THM-5127 (order 10: at least 1,504 classes for every odd h within
  one class-standard-deviation of the class mean; local variance about twice Poisson; residues balanced to 1.2%). The
  implication is PROVED (source note, Proposition 5.1).
  (B), (B'): FINITE-EXACT agreement within 1.2 percentage points at orders 8 and 9; pre-registered P1-P2 at order 10
  (9,733,056 classes) CONFIRMED, agreement within 0.08 points. (B) alone gives only that the mod-4 discrepancy tends
  to 0; (B') gives the exponential rate (audit M3). (B'') (deletion-pair form, added after the audit): labelled
  rho_m = 1.244, 1.156 at m = 8, 9. Order 10 pre-registered (P3-P6), run pending at this commit.
  (C): FINITE-EXACT for 237 <= h <= 1023 with bases of order <= 10 and all out-sets. The implication is PROVED
  (source note, Proposition 6.1). (C) is the all-out-sets form of klein's HYP-9453 (2).
  AUDIT: PASS WITH CORRECTIONS (the packet's audit/AUDIT.md; M3 and M4 concern this file; applied; MISTAKE-753).
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
the tilt with `E_{p_n}[H_n] = (n-1)!/2^(n-2)`. There are `n_0 >= 4` and `c > 0` such that for all `n >= n_0`, all
`p in [1/2, p_n]` and all odd `h` with `|h - E_p[H_n]| <= sd_p(H_n)`: `P_p(H_n = h) >= c / sd_p(H_n)`.

**Consequence (PROVED).** Every odd `h > (n_0 - 1)!/2^(n_0 - 2)` is a Hamiltonian path count.

**Why it is the right regime.** Over all tournaments the number of structures with count at most `X` is
super-polynomial in `X`, and by THM-5127 (B) every target is the mean of exactly one tilt at each admissible order.
The word families studied so far have bounded, logarithmic or polynomial multiplicity.

**Falsifier.** An order `n` and an odd `h` in the stated window with `P_p(H_n = h) = 0`, for `n` beyond any proposed
`n_0`; or, at finite range, class counts in the bulk at order 11 that are not smooth.

**Route.** Characteristic function of `H_n` under `P_p`: Gaussian regime at `p = 1/2` from the Walsh levels
(THM-4757 section 5, THM-4761; CONDITIONAL on Moran's CITED theorem);
rational frequencies from THM-5127 (C), (C4) and part (B) below; the remaining range needs a decoupling over several
exposed vertices and is untouched. No variance formula or limit law is known for `p != 1/2` (at `p = 1/2`: THM-4761,
conditional on Moran's theorem).

## (B) Random-alternating law for the pair form

**Statement.** For `R` uniform on `m` labelled vertices, the distribution of `rank(beta_R)` is within `o(1)` in total
variation of the rank distribution of a uniformly random alternating `m x m` matrix over `F_2`. By THM-5127 (C) (4) this gives that the mod-4 discrepancy of `H` at order `m+1` tends to 0.

**(B') Quantitative form.** There is an absolute constant `C` with `E_R[2^(-rank(beta_R)/2)] <= C 2^(-m/2)` (for a
uniformly random alternating matrix the left side is about `1.62 * 2^(-m/2)` for even `m` and `1.65 * 2^(-m/2)` for
odd `m`). By THM-5127 (C) (4), (B') gives a
mod-4 discrepancy `O(2^(-m/2))` at order `m+1`. (B) does not imply (B'): closeness `o(1)` in total variation only
gives `o(1)`, and the low-rank tail is exactly where the census deviates (rank-2 classes are over-represented by
factors 7.2, 66, 145 at orders 8, 9, 10).

**Evidence.** Orders 8 and 9, all classes: `0.4183 / 0.5516 / 0.0301` against `0.4205 / 0.5585 / 0.0210`, and
`0.8278 / 0.1699 / 0.0023` against `0.8394 / 0.1593 / 0.0013` (top rank / next / lower). Order 10, pre-registered
(P1: fractions in `[0.405, 0.430]`, `[0.550, 0.575]`, `[0.020, 0.032]`; P2: control): `0.4189 / 0.5594 / 0.0217`
against `0.4197 / 0.5590 / 0.0213`, CONFIRMED. For (B') with labelled weights (audit):
`E_R[2^(-rank/2)] = 0.10239, 0.07348, 0.05077` at `m = 8, 9, 10` against `0.10135, 0.07271, 0.05072`. The exact
labelled mod-4 discrepancy is much smaller than either bound: `0.0415, 0.00752, 0.00242` at orders 8, 9, 10.

**(B'') Deletion-pair form (added after the audit; not audited).** Let `T` be uniform on `m+2` labelled vertices and
`v != v'` two fixed vertices. Then the mod-4 classes of the two deletions decorrelate:
`E_T (-1)^((H(T-v) - H(T-v'))/2) -> 0`; quantitatively `rho_m = 2^m E_T (-1)^((H(T-v) - H(T-v'))/2)` is bounded (and
tends to 1). By THM-5127 (C) (5) the expectation equals `E_R[b_R^2]` with `b_R^2 in {0, 2^(-rank(beta_R))}`, and it is
`E_T (-1)^(N_1(T))`, `N_1` the number of odd cycles of `T` through exactly one of `v, v'`. Relations (PROVED): (B)
implies the qualitative form; bounded `rho_m` implies a mod-4 discrepancy `O(2^(-m/2))` at order `m+1`
(Cauchy-Schwarz); neither of (B') and bounded `rho_m` implies the other (zero-bias forms count in (B') only; (B')
bounds `E_R[b_R^2]` only by `C 2^(-m/2)`). The value `rho_m = 1` is the contribution of twin extensions and the exact
value for a uniformly random quadratic form. This is the smallest open statement on the mod-4 route: no bound below
the trivial 1 is proved here for `E_T (-1)^(N_1(T))` at large `m`.

**Evidence for (B'').** Labelled, exact: `rho_m = 1, 2, 2.875, 2.625, 1.583, 1.758, 1.244, 1.156` at `m = 2, ..., 9`. Order 10: pre-registered P3-P6 (`PREREGISTRATION.md`), run pending at this commit.

**Falsifier.** A rank distribution at order 11 that moves away from the random-alternating one; for (B'), a sequence
of orders along which `2^(m/2) E_R[2^(-rank(beta_R)/2)]` is unbounded; for (B''), a sequence along which `rho_m` is
unbounded.

## (C) Window law

**Statement.** For every odd `h >= 237` that is a Hamiltonian path count and every odd `n` with `3h <= n < 5h` there
is a tournament `T` with `H(T) = h` and a one-vertex extension `T + z` with `H(T + z) = n`.

**Consequence (PROVED).** Every odd number other than 7 and 21 is a Hamiltonian path count.

**Evidence.** Exact for `237 <= h <= 1023` with bases of order at most 10; with the same bases the window `[3h, 9h)` is
complete for `365 <= h <= 1023`. Below 237 there are 48 exceptional `h` (largest 235).

**Falsifier.** An odd `h >= 237` and an odd `n` in `[3h, 5h)` that is not an extension of any realiser of `h`.

**Limits.** In the census no single realiser covers more than 44% of `[3h, 5h)` for `h >= 237` (FINITE-EXACT, bases
of order `<= 10`, audit); by THM-5124 the witnesses cannot come from additive presence cubes. No structural control of the extension sets is known; the statement is stronger than the conjecture.
