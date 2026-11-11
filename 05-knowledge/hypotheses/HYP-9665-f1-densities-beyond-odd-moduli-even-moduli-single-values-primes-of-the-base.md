---
id: HYP-9665
title: "F1 local densities beyond THM-5134: (A) even moduli 2^j m carry the product of the flat 2-adic law and L_m; (B) the multiplicity of a single odd N follows 2 kappa_s (8/pi^2) sigma(n)/n ln N on average over short ranges; (C) at the primes of the base b the second digit of N_b is equidistributed; (D) the first-order constants A_(m,d) are computable from the phase clocks of the one-apex count"
status: >
  OPEN (all four parts). PRE-AUDIT DRAFT (2026-11-11). NUMERICAL evidence only where stated; nothing here is
  pre-registered. (B) would follow from a version of THM-5134 uniform in the modulus, which is HYP-9660 (B), OPEN.
session: opus-2026-11-11-S92
date: 2026-11-11
source: 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/README.md
related:
  - 01-canon/theorems/THM-5134-* (the proved law at odd moduli and in every base)
  - 01-canon/theorems/THM-5128-f1-by-value-local-densities-first-order-class-law.md ((E): powers of 2)
  - 05-knowledge/hypotheses/HYP-9660-* ((B) uniform gap)
  - 05-knowledge/results/odd_tables_signed_carries_20261111/prime_layers/PHASE_CLOCKS.md
---

# HYP-9665 - what THM-5134 leaves open about the local densities of F1

Notation of THM-5134.

**(A) Even moduli.** For `j >= 1` and odd `m >= 3`, the limiting class factor of an odd class `c` modulo `2^j m` is
`2 L_m(c - c_sigma)` on each sheet (flat on the odd classes modulo `2^j`, THM-5128 (E), times the odd law). What is
missing is the joint equidistribution on cells of `q` modulo `2^(j-1)` and of the state modulo `m`, i.e.
`rho(T_k) < 2` for the non-trivial characters of `Z/(2^(j-1) m)`. Test: exact census modulo 12, 36, 20, 72 with the
remainder of the product law; a remainder that does not decrease from `2^24` to `2^28` at 12 would refute it.

**(B) Single values.** Let `r_sigma(N)` be the number of triples of the sheet `sigma` with value `N`. Question:
for `X^eps <= K <= X`, is

    sum_{N in (X, X + K], N odd} | r_sigma(N) - 2 kappa_sigma S_sigma(N) ln N | = o(K ln X),      S_sigma(N) = (8/pi^2) sigma(n)/n ?

NUMERICAL (exploratory, series_vs_counts_25_24.out; odd `N` in `[2^24, 2^25)`): mean multiplicities `5.918` (sheet
0) and `7.482` (sheet 1); correlation of `r` with `sigma(n)/n`: `0.4016` and `0.4275` (with the squarefree-level
product of THM-5128 (D3): `0.3929`, `0.4181`); mean of `r` by quintile of `sigma(n)/n`: from `4.30` to `8.55` (sheet 0)
and from `5.44` to `10.77` (sheet 1). The ratios by `v_3(n)` and `v_5(n)` are above their limits (`1.455, 1.722` at
3 against `1.333, 1.444`), as the `1/ln X` term of THM-5134 (D) requires. The statement is the F1 analogue of
Ingham's asymptotic for a sum of two products; it is false if the remainder of THM-5134 (B) is not `o(1)`
uniformly for `m` up to a power of `ln X`.

**(C) The primes of the base.** `N_b = 1 + alpha_1 b + alpha_2 b^2`, so `N_b = 1 (mod b)` and the second base-`b`
digit of `N_b` is `alpha_1 = psi1 + psi2 + (2^s y1 y2 - [s = 0])/2 (mod b)`. Conjecture: it is equidistributed
modulo `b` by value on each sheet, for every `b >= 2`. NUMERICAL (every_base.out, sheet `s >= 1`, `N_b <= 2^24`):
deviation from flatness at most `0.003` for `b = 2, 4, 5, 10`, but `0.040` for `b = 3` and `0.052` for `b = 6`. A
deviation at 3 and 6 that does not decrease at `2^28` would refute it. For odd `b` the 2-adic law of `N_b` is not
known either.

**(D) First-order constants.** The constants `A_(m,d)^sigma` of THM-5134 (C) are convergent sums over the states of
type `d`. For `p | 2^t + 1` the first balanced state of type `p` is the fan `h = 2^(t-1)`, `t = ord_p(2)/2`
(bridges_exact.out S8: it is the first state of type `p` for `p = 3, 5, 11, 13, 17, 43, 113` and not for
`p = 19, 29, 37, ...`); for `p` with `ord_p(2)` odd there is no balanced state. Question (from the sweep, not
tested): are the positions of the states of type `p`, hence `A_(p,p)`, determined by the phase data `(i, a_i)` of
the odd-tables packet, which describe the one-apex count along the doubling orbits modulo `p`?
