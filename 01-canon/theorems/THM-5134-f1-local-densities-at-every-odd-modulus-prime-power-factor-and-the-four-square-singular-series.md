---
id: THM-5134
title: "Local densities of the one-arc join family F1 at every odd modulus. (A) For a fixed smaller state with d = gcd(x1, y1, m), the value M = x1 x2 + 2^s y1 y2 is uniform on the multiples of d modulo m, up to a power saving. (B) First-order law: #{N = c mod m} = sum over d | gcd(m, c - c_s) of (d/m) #{triples whose smaller state has gcd(g, m) = d} + O(X^(1 - theta_m/2)). (C) The shares of the types tend to d^-2 prod_{p | m/d} (1 - p^-2) at the rate 1/ln X. (D) Limit: the class factor m * share tends to L_m(c - c_s) = sum_{d | gcd(m, e)} d^-1 prod_{p | m/d} (1 - p^-2), multiplicative over prime powers, equal to (1 - p^-2) sigma(p^v)/p^v at a class of valuation v < k modulo p^k. (E) Singular series: along odd moduli divisible by ever higher powers of ever more primes the class factor of N tends to (8/pi^2) sigma(n)/n, n the odd part of N (sheet s >= 1) or of N + 1 (sheet s = 0). (F) This is the singular series of Jacobi's four-square theorem, r_4(n)/(pi^2 n) for odd n: at every odd prime the gluing form of F1 is the norm form of the quaternions. (G) A seam of rank r with uniform states has factor (1 - p^-r) sigma_(r-1)(p^v)/p^((r-1)v); r = 1, 2, 4 are the split composition algebras, and r = 4 gives the eight-square series (96/pi^4) sigma_3(n)/n^3."
status: >
  PRE-AUDIT DRAFT (2026-11-11). (A)-(E) PROVED, by the argument of THM-5128 (C), (D) carried out for an odd modulus
  with THM-5128 (B2), (B4)(iv) (the step that THM-5128 (D3) left "not written out" and typed NUMERICAL after its
  audit). (F): the identification is PROVED (equivalence of two unimodular quaternary forms over Z_p, p odd) with
  Jacobi's theorem CITED. (G): PROVED (an exact count for uniform states); the identities for r_8 are CITED.
  The exact factor for uniform states modulo p^k and the series (8/pi^2) sigma(n)/n were found by the audit of
  THM-5128 (auditors proof-5128-CF, citations, prior-art-wording; 2026-11-11), which also enumerated 9, 25, 27, 49;
  NEW here: the by-value proof for F1 at every odd modulus, the first-order law (B) with its types, and (F), (G).
  NUMERICAL / FINITE-EXACT support: exact remainder of (B) at X = 2^22, 2^24, 2^26, 2^28 for twelve moduli; four
  PRE-REGISTERED tests (window [2^27, 2^28) and the remainder at 2^28; registration commit 73c8437aac precedes the
  runs), all CONFIRMED, the fourth at its threshold (8 of 12). Not uniform in the modulus (HYP-9660 (B) is OPEN);
  even moduli beyond powers of 2 are not treated (remark at the end). Chan-Pak completeness remains OPEN: local
  densities are positive in every class and decide nothing about it.
session: opus-2026-11-11-S92
date: 2026-11-11
source: 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/README.md
scripts:
  - 04-computation/experiments/s92_prime_powers_20261111/odd_modulus_census.py (+ cum_22.out, cum_24.out, cum_26.out, cum_28.out, window_insample_*.out, window_prereg_28_27.out), PREREGISTRATION.md, prereg_window_28_27.json
  - 04-computation/experiments/s92_prime_powers_20261111/seam_rank_exact.py (+ .out)
related:
  - 01-canon/theorems/THM-5128-f1-by-value-local-densities-first-order-class-law.md (the prime and squarefree case; (B2), (B4)(iv), (D1), (D2) are used)
  - 05-knowledge/results/s91_f1_local_densities_20261110/audit/AUDIT_REPORT.md (findings proof-5128-CF#0, citations#2, prior-art-wording#12: the prime-power factor)
  - 05-knowledge/hypotheses/HYP-9660-* (uniform gap, OPEN), HYP-9665-* (this session)
---

# THM-5134 - local densities of F1 at every odd modulus, and the four-square singular series

## Setting

As in THM-5128. One-apex states `x(h) = 2 psi(h) + 1`, `y(h) = 2h + 1`, `g(h) = gcd(x(h), y(h))`. F1:
`N = M - [s = 0]`, `M = x1 x2 + 2^s y1 y2`, over ordered triples `(h1, h2, s)`; sheet 0 is `s = 0` (`c_0 = -1`), sheet
1 is `s >= 1` (`c_1 = 0`); `T_sigma(X)` is the set of triples of the sheet `sigma` with `N <= X`; `Y = X/2^s`.
`m >= 3` is an odd modulus, `rho_m` and `theta_m = ln(1/rho_m)/ln(2/rho_m^2)` are as in THM-5128 (B). For `d | m`:

- a state is of *type* `d` if `gcd(g(h), m) = d`; a triple with `h1 != h2` is of type `d` if its state of smaller
  word is; `D_sigma^(d)(X)` is the number of triples of type `d` in `T_sigma(X)`;
- `c_d = d^-2 prod_{p | m/d} (1 - p^-2)`; this is the share of the residue vectors `(alpha, beta)` modulo `m` with
  `gcd(alpha, beta, m) = d` (Jordan's totient `J_2(m/d)` over `m^2`), and `sum_{d | m} c_d = 1`;
- `chi_d(h) = [gcd(g(h), m) = d] - c_d`;
- `L_m(e) = sum_{d | gcd(m, e)} d c_d = sum_{d | gcd(m, e)} d^-1 prod_{p | m/d} (1 - p^-2)` for a class `e` modulo `m`;
- `sigma(n)` is the sum of the divisors of `n`, `sigma_j(n)` the sum of their `j`-th powers.

## Statements

**(A) Fixed smaller state.** Let `h1`, `s` be fixed, `d = gcd(g(h1), m)` and `G = #{h2 > h1 : M <= X}`. For every
class `c'` modulo `m`:

    #{h2 > h1 : M <= X, M = c' (mod m)} = [d | c'] (d/m) G + O_m((Y/y1)^(1 - theta_m)).

If `d = m` then `M = 0 (mod m)` for every `h2`.

**(B) First-order law.** For every class `c` modulo `m` and each sheet:

    #{t in T_sigma(X) : N = c (mod m)} = sum over d | gcd(m, c - c_sigma) of (d/m) D_sigma^(d)(X) + O_m(X^(1 - theta_m/2)).

So all the classes with the same `g0 = gcd(m, c - c_sigma)` have the same share up to the error, and the class
factor is `m * share = Phi(g0) + O(X^(-theta_m/2))`, `Phi(g0) = sum_{d | g0} d F_d`, `F_d = D^(d)/#T` the share of the
triples of type `d`. Every triple of type `d` has `d | N - c_sigma` exactly.

**(C) The types.** `D_sigma^(d)(X) = c_d #T_sigma(X) + 2 A_(m,d)^sigma X + O_m(X^(1 - theta_m/2))`, with
`A_(m,d)^sigma = lim_{H -> infinity} sum_{h < H} chi_d(h) w^sigma(h)` and the weights `w^0`, `w^1` of THM-5128 (D2)
(the limit exists as there). Hence `F_d -> c_d` at the rate `1/ln X`.

**(D) Limit.** For every class `c` modulo `m`:

    m #{N = c} / #T_sigma(X) = L_m(c - c_sigma) + (2/kappa_sigma) (sum_{d | gcd(m, c - c_sigma)} d A_(m,d)^sigma) / ln X + O(1/ln^2 X).

`L_m` is multiplicative: `L_m(e) = prod_{p^k || m} L_(p^k)(e)`, and with `v = v_p(e)`

    L_(p^k)(e) = (1 - p^-2) sigma(p^v)/p^v = (1 + 1/p)(1 - p^(-v-1))                  if v < k,
    L_(p^k)(e) = (1 - p^-2) sigma(p^(k-1))/p^(k-1) + p^-k = 1 + 1/p - p^(-k-1)        if v >= k.

For `k = 1` this is the factor `1 - 1/p^2 + [p | e]/p` of THM-5128 (D3); the mean of `L_m` over the classes is 1.

**(E) Singular series.** Let `N` be odd, `sigma` a sheet with `N + [sigma = 0] != 0`, and `n` the odd part of
`N + [sigma = 0]`. For every odd `m` in which every prime `p | m` occurs to a power higher than `v_p(n)`, the
limiting class factor of the class of `N` modulo `m` is `prod_{p | m} (1 - p^-2) sigma(p^(v_p(n)))/p^(v_p(n))`. Along
any sequence of such `m` whose prime factors exhaust the odd primes it tends to

    S_sigma(N) = prod_{p odd} (1 - p^-2) sigma(p^(v_p(n)))/p^(v_p(n)) = (8/pi^2) sigma(n)/n.

`S` has mean 1 over the odd `N` of either sheet and minimum `8/pi^2 = 0.8106` (at `n = 1`: the powers of 2 on sheet
1, and `N + 1` a power of 2 on sheet 0); it is unbounded, of maximal order `(4 e^gamma/pi^2) ln ln N` (Gronwall's
theorem for odd `n`; CITED). Over both sheets the density is `w_0 S_0 + w_1 S_1` (triple shares of THM-5128 (D1)),
not an Euler product.

**(F) The four-square series.** For odd `n`, `(8/pi^2) sigma(n)/n = r_4(n)/(pi^2 n)`, where `r_4(n) = 8 sigma(n)` is
the number of representations of `n` as a sum of four squares (Jacobi; CITED): `S_sigma(N)` is the singular series
of the four-square problem at `n`. The reason is local: for an odd prime `p` and a unit `u`, the quaternary forms
`a1 a2 + u b1 b2` and `a^2 + b^2 + c^2 + d^2` are equivalent over `Z_p`, so they take every class modulo every `p^k`
equally often; and the states of F1 are equidistributed modulo `p^k` (THM-5128 (B2)). `a1 a2 + b1 b2` is the
determinant of `[[a1, -b1], [b2, a2]]`, the norm form of the split quaternions `M_2`; the sum of four squares is the
norm form of Hamilton's quaternions; the two algebras are isomorphic over `Q_p` for every odd `p` and differ only at
2 and at infinity. There the two problems differ as well: the four-square problem has 2-adic density 1 at odd `n`
and real density `pi^2 n`, F1 has the flat law on the odd classes modulo `2^j` (THM-5128 (E)) and the mass
`2 kappa_sigma ln N` per odd `N` (THM-5128 (D1)). Equivalently `sigma(n)` is the number of sublattices of index `n`
in `Z^2`, i.e. of integer matrices of determinant `n` up to `SL_2(Z)` (Hermite normal forms; CITED).

**(G) Seams of rank `r`.** Let `u, v` be independent and uniform in `(Z/p^k)^r` (`p` any prime) and `B` a bilinear
form on `(Z/p^k)^r` with unit determinant. For a class `c` with `v = v_p(c) < k`:

    p^k P(B(u, v) = c) = (1 - p^-r) sigma_(r-1)(p^v)/p^((r-1)v) = (1 - p^-r)(1 + p^-(r-1) + ... + p^-(r-1)v).

So the series of a rank-`r` seam with equidistributed states is `prod_p (1 - p^-r) sigma_(r-1)(p^(v_p(n)))/p^((r-1)v_p(n))`,
that is `sigma_(r-1)(n)/(zeta(r) n^(r-1))` over all primes for `r >= 2`. The form `sum_{i <= r} a_i b_i` is the norm
form of a composition algebra exactly for `r = 1, 2, 4` (split complex numbers, split quaternions, split octonions
in Zorn's vector-matrix form `a b - <u, v>`; Hurwitz; CITED). For `r = 2` this is (F). For `r = 4` and odd `n`,
restricted to the odd primes, it is `(96/pi^4) sigma_3(n)/n^3 = r_8(n)/(pi^4 n^3/6)`, the singular series of the
eight-square problem (`r_8(n) = 16 sigma_3(n)` for odd `n`, Jacobi; CITED). For `r = 1` the product diverges to 0:
the count of `a b = n` is the divisor function and has no singular series.

**Numerical support.** Exact remainder of (B), `m max_c |#{N = c} - sum_{d | gcd(m, c - c_sigma)} (d/m) D^(d)| / #T`,
sheet 0 / sheet 1 (cum_*.out; the triple counts equal those of THM-5128 at `2^24`, `2^26`):

| `m` | `X = 2^22` | `2^24` | `2^26` | `2^28` | noise at `2^28` |
|---|---|---|---|---|---|
| 9 | `0.035 / 0.028` | `0.021 / 0.024` | `0.024 / 0.020` | `0.024 / 0.019` | `0.0001` |
| 27 | `0.046 / 0.044` | `0.022 / 0.029` | `0.029 / 0.022` | `0.029 / 0.027` | `0.0002` |
| 81 | `0.062 / 0.055` | `0.024 / 0.042` | `0.034 / 0.023` | `0.031 / 0.031` | `0.0003` |
| 25 | `0.004 / 0.008` | `0.004 / 0.005` | `0.004 / 0.005` | `0.004 / 0.004` | `0.0002` |
| 125 | `0.011 / 0.015` | `0.008 / 0.008` | `0.005 / 0.007` | `0.005 / 0.005` | `0.0003` |
| 49 | `0.182 / 0.140` | `0.163 / 0.141` | `0.130 / 0.125` | `0.093 / 0.103` | `0.0002` |
| 45 | `0.083 / 0.057` | `0.040 / 0.056` | `0.039 / 0.043` | `0.038 / 0.035` | `0.0002` |
| 75 | `0.097 / 0.065` | `0.050 / 0.056` | `0.034 / 0.031` | `0.027 / 0.029` | `0.0003` |

The remainder is far above the sampling noise and decreases slowly, with the oscillation at 3 and the slow rate at
7 already seen for primes (THM-5128 (C)); at `m = 9, 27, 81` it does not decrease from `2^26` to `2^28` on one sheet
each. Every triple of type `d` lies in a class with `d | c - c_sigma` (exact, all runs).
Class factors at `X = 2^28` (measured mean over the group; first-order value `Phi(g0)` with the measured types;
limit `L_m`), sheet 1: at 9, `g0 = 1, 3, 9`: see cum_28.out; the factors are still far from their limits (at 9 the
class `0` has about `1.45` against `35/27 = 1.296`), as (D) predicts with its `1/ln X` term.
Pre-registered (PREREGISTRATION.md; window of odd `N` in `[2^27, 2^28)`, a-priori predictions of `Phi(g0)` from the
smaller states alone): P1 group means within the registered tolerances (largest deviations `0.027` at 9, `0.028` at
27, `0.032` at 81, `0.004` at 25, `0.004` at 125, `0.091` at 49, `0.030` at 45, `0.021` at 75) - CONFIRMED; P2
flatness inside each group - CONFIRMED; P3 the measured factor is closer to the a-priori prediction than to the
limit for every group with `g0 > 1` at 9, 27, 25 (14 of 14) - CONFIRMED; P4 the remainder of (B) is smaller at
`2^28` than at `2^26` for at least 8 of 12 pairs - CONFIRMED with exactly 8 (it grew at 9 sheet 0, 27 sheet 1,
81 sheet 1, 125 sheet 0). The ordering of the factors along divisor chains fails at finite size and was not
registered. (G): seam_rank_exact.out enumerates `r = 1, 2, 3` at `9, 27, 25` and `r = 4` at 9.

## Proofs

**(A).** `M = a x2 + b y2` with `(a, b) = (x1, 2^s y1)`, and `gcd(a, b, m) = gcd(x1, y1, m) = d` because `2^s` is a
unit modulo `m`. The homomorphism `(Z/m)^2 -> Z/m`, `w -> a w_1 + b w_2`, has image `d Z/m`, of order `m/d`, so it
takes each multiple of `d` on `m d` residue vectors and no other value. By THM-5128 (B4)(ii) in the form (iv), each
residue vector modulo `m` is taken by `(x2, y2)` on `m^-2 #(B_n ∩ I ∩ Omega) + O_m(2^(n(1 - theta_m)))` states of a
block; sum over the `m d` vectors, and over the blocks `2^n <= Y/y1` with `I = (h1, infinity)`, as in the proof of
THM-5128 (C1). If `d = m` both coefficients vanish modulo `m`.

**(B).** As for THM-5128 (C2). `M` is symmetric in `(h1, h2)`; a triple with `h1 != h2` is counted through its
smaller state, and the triples with `h1 = h2` number `O(Y^(1/2))` for each `s`. `N = c` is `M = c - c_sigma` on the
sheet `sigma`. By (A) the triples of type `d` contribute `(d/m) D^(d)` to every class with `d | c - c_sigma`, and
nothing to the others, up to errors that add up to `sum_{y1 < Y^(1/2)} (Y/y1)^(1 - theta_m) = O(Y^(1 - theta_m/2))`;
the sum over `s` keeps the bound. The exact statement is `d | x1 x2 + 2^s y1 y2` when `d | g(h1)`.

**(C).** The condition `gcd(g(h), m) = d` is a condition on `(x(h), y(h))` modulo `m`, met by `c_d m^2` residue
vectors: the vectors with `gcd(alpha, beta, m) = d` are `d` times the vectors modulo `m/d` that are primitive at
every prime of `m/d`, and there are `J_2(m/d) = (m/d)^2 prod_{p | m/d} (1 - p^-2)` of those. So `chi_d` is a `chi_S`
of THM-5128 (B4)(iv), and the proof of THM-5128 (D2) applies word for word with `chi_d` in place of `chi`.

**(D).** From (B), (C) and THM-5128 (D1): the share of the class is `sum_{d | gcd(m, c - c_sigma)} (d/m) F_d` and
`F_d = c_d + 2 A_(m,d) X/#T + ...`, `#T = kappa_sigma X ln X + O(X)`. Multiplicativity: `d -> d c_d` is a product
over the primes of `m` of the local terms `p^-i (1 - p^-2)^[i < k]` (`p^i || d`, `p^k || m`), and the divisors of
`gcd(m, e)` are the products of the `p^i` with `i <= min(v_p(e), k)`. For `m = p^k`:
`L = sum_{i <= min(v, k)} p^-i (1 - p^-2)^[i < k]`, which is the stated value. The mean over the classes is
`sum_d (d/m)(m/d) c_d = sum_d c_d = 1`.

**(E).** By (D) with `v_p(N - c_sigma) = v_p(n) < k` for every `p^k || m`. The Euler product:
`prod_{p odd} (1 - p^-2) = 8/pi^2` and `prod_p sigma(p^v)/p^v = sigma(n)/n`. Mean: the Dirichlet series of
`sigma(n)/n` over odd `n` is `zeta_odd(s) zeta_odd(s + 1)`, so its mean over odd `n` is `zeta_odd(2) = pi^2/8`; on
sheet 0 the odd part of `N + 1` is distributed over the odd integers with the weights `2^-a`, which does not
change the mean. Minimum: `sigma(n)/n >= 1` with equality at `n = 1`.

**(F).** Over `Z_p`, `p` odd, a unimodular quadratic form is determined up to equivalence by its rank and by its
discriminant modulo squares of units. The sum of four squares has discriminant 1. The form `a1 a2 + u b1 b2` is the
orthogonal sum of the hyperbolic planes `a1 a2` and `u b1 b2`, with Gram determinant `(-1/4)(-u^2/4) = (u/4)^2`, a
square. Equivalent forms have the same number of solutions of `Q = c` modulo every `p^k`. The count for the first
form with uniform variables is (G) with `r = 2`; so it is also the `p`-adic density of the four-square form, and
Jacobi's `r_4(n) = 8 sigma(n)` (odd `n`) with real density `pi^2 n` and 2-adic density 1 gives the same product.

**(G).** A linear change of `v` by the matrix of `B` (a unit determinant) reduces to `B(u, v) = u . v`. With
probability `p^(-ri) (1 - p^-r)` the vector `u` is `p^i` times a primitive vector, `i < k`, and then `v -> u . v` is
`p^i` times a surjective homomorphism onto `Z/p^k`, uniform on the multiples of `p^i`. So for `v_p(c) = v < k`:
`P(c) = sum_{i <= v} p^(-ri) (1 - p^-r) p^(i-k)`. The product over the primes is
`(1/zeta(r)) prod_p sigma_(r-1)(p^v)/p^((r-1)v)`.

## Remarks

1. *What changed with respect to THM-5128.* Its (D3) product `prod_p (1 - 1/p^2 + [p | n]/p)` is the squarefree
   level: the mean of `L_(p^k)` over the classes divisible by `p`. The series at `N` is `(8/pi^2) sigma(n)/n`:
   smaller by `p^-3` at each prime dividing `n` exactly once, larger from `p^2 | n` on.
2. *Finite size.* At every computed size the factors are governed by the `1/ln X` term of (D), i.e. by the few
   smallest states of each type: at 25 the classes with `v = 1` and `v >= 2` have nearly the same factor (`1.25`),
   because the state `h = 2` (`x = y = 5`) is of type 5 and no state of type 25 is small.
3. *Even moduli.* For `2^j m` with `j >= 2` and `m > 1` odd the natural statement is the product of the flat law on
   the odd classes modulo `2^j` and of `L_m`. It needs the joint equidistribution of `q` modulo `2^(j-1)` and of the
   states modulo `m` on cells, that is `rho(T_k) < 2` for the characters of `Z/(2^(j-1) m)`, which is not proved
   here (HYP-9665).
4. *No uniformity.* All constants depend on `m`; nothing here is a level of distribution.
