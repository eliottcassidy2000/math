---
id: THM-5134
title: "Local densities of the one-arc join family F1 at every odd modulus. (A) For a fixed smaller state with d = gcd(x1, y1, m), the value M = x1 x2 + 2^s y1 y2 is uniform on the multiples of d modulo m, with a power saving. (B) First-order law: #{N = c mod m} = sum over d | gcd(m, c - c_s) of (d/m) #{triples whose smaller state has gcd(g, m) = d} + O(X^(1 - theta_m/2)). (C) The shares of the types tend to d^-2 prod_{p | m/d} (1 - p^-2) at the rate 1/ln X. (D) Limit: the class factor m * share tends to L_m(c - c_s) = sum_{d | gcd(m, e)} d^-1 prod_{p | m/d} (1 - p^-2), multiplicative over prime powers, equal to (1 - p^-2) sigma(p^v)/p^v at a class of valuation v < k modulo p^k. (E) Singular series (8/pi^2) sigma(n)/n, n the odd part of N (sheet s >= 1) or of N + 1 (sheet s = 0). (F) It is the series of Jacobi's four-square theorem, r_4(n)/(pi^2 n) for odd n, and Ingham's series for a sum of two products with the prime 2 removed. (H) Every base: the base-b packing counts of the same tournaments obey the same law at every odd modulus prime to b, with the class -b/2 in place of -1."
status: >
  PRE-AUDIT DRAFT (2026-11-11; independent audit pending, statuses below are the author's).
  PROVED: (A)-(E) and (H), by the argument of THM-5128 (C), (D) carried out for an odd modulus with THM-5128
  (B2), (B4)(iv); this is the step that THM-5128 (D3) left "not written out" and typed NUMERICAL after its audit.
  (F): the identification is PROVED (two unimodular quaternary forms of the same discriminant over Z_p, p odd);
  Jacobi's and Ingham's theorems are CITED.
  KNOWN, not claimed: the fixed-smaller-state law of (A) at every fixed odd modulus, composite and prime-power
  included, without a rate (Session206, weighted parameter capacity: uniform on the residues allowed by
  g = gcd(R, d, p), error o(X)); the exact factor (1 - p^-2) sigma(p^v)/p^v for uniform states and the series
  (8/pi^2) sigma(n)/n (audit of THM-5128, 2026-11-11, which enumerated 9, 25, 27, 49 and called the law "the
  classical determinant-form density"); valuation-resolved prime-power factors for another family of path counts
  (Session223: two root changes, 3-adic f(v), fixed length); the base-b packing polynomial of the joins
  (THM-5070 (2), kind-pasteur S18) and N_b = 1 (mod b) (THM-4940).
  NEW here: the power saving in (A) and the sum over the smaller state, i.e. the by-value law (B)-(E) at every
  odd modulus; the identification (F) inside the repository (Jacobi's r_4 was cited in THM-4826 and Session116,
  never for path counts); (H).
  NUMERICAL / FINITE-EXACT support: exact remainder of (B) at X = 2^22, 2^24, 2^26, 2^28 for twelve moduli; four
  PRE-REGISTERED tests (window [2^27, 2^28) and the remainder at 2^28; registration commit 73c8437aac precedes the
  runs), all CONFIRMED, the fourth at its threshold (8 of 12); base-b censuses (exploratory, not pre-registered).
  Not uniform in the modulus (HYP-9660 (B) is OPEN); even moduli beyond powers of 2 are not treated (HYP-9665).
  Chan-Pak completeness remains OPEN: the local densities are positive in every class, 7 and 21 included, and
  decide nothing about it.
session: opus-2026-11-11-S92
date: 2026-11-11
source: 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/README.md
scripts:
  - 04-computation/experiments/s92_prime_powers_20261111/odd_modulus_census.py (+ cum_22.out, cum_24.out, cum_26.out, cum_28.out, window_insample_*.out, window_prereg_28_27.out), PREREGISTRATION.md, prereg_window_28_27.json
  - 04-computation/experiments/s92_prime_powers_20261111/seam_rank_exact.py (+ .out)
  - 04-computation/experiments/s92_prime_powers_20261111/every_base.py (+ .out), bridges_exact.py (+ .out, sections S1, S4, S10), series_vs_counts.py (+ series_vs_counts_25_24.out)
related:
  - 01-canon/theorems/THM-5128-f1-by-value-local-densities-first-order-class-law.md (the prime and squarefree case; (B2), (B4)(iv), (D1), (D2) are used)
  - 05-knowledge/results/s91_f1_local_densities_20261110/audit/AUDIT_REPORT.md (findings proof-5128-CF#0, citations#2, prior-art-wording#12: the prime-power factor)
  - session206/weighted-parameter-capacity/WEIGHTED_CAPACITY.md (fixed-source law at a fixed odd modulus)
  - 05-knowledge/results/session223_root_changes_and_local_factors_20261107/README.md (valuation-resolved 3-adic factors)
  - 01-canon/theorems/THM-5070-port-chain-transfer-for-odd-cycle-packing-counts-and-ternary-coverage-to-a-million.md, THM-4940-odd-cycle-packing-counts-in-every-base.md (base b)
  - 01-canon/theorems/THM-4826-which-grid-carries-a-norm-shell-one-to-two-split-cubic-character-of-two-centred-units-complete-the-quartic-designs-four-eisenstein-fourth-powers-miss-five-mod-eight.md (Jacobi, Hurwitz, Lipschitz in canon)
  - 01-canon/theorems/THM-5135-* (class laws of bilinear seams; the algebra behind (F))
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
- `sigma(n)` is the sum of the divisors of `n`.

## Statements

**(A) Fixed smaller state.** Let `h1`, `s` be fixed, `d = gcd(g(h1), m)` and `G = #{h2 > h1 : M <= X}`. For every
class `c'` modulo `m`:

    #{h2 > h1 : M <= X, M = c' (mod m)} = [d | c'] (d/m) G + O_m((Y/y1)^(1 - theta_m)).

If `d = m` then `M = 0 (mod m)` for every `h2`. (Without the power saving, with error `o(X)`, this is Session206's
weighted-capacity theorem; KNOWN.)

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
`m^3 L_m(c)` is the number of `(a1, a2, b1, b2)` modulo `m` with `a1 a2 + b1 b2 = c`, the additive square of the
multiplication table modulo `m` (THM-5135 (1); checked for eleven odd `m <= 105`).

**(E) Singular series.** Let `N` be odd, `sigma` a sheet with `N + [sigma = 0] != 0`, and `n` the odd part of
`N + [sigma = 0]`. For every odd `m` in which every prime `p | m` occurs to a power higher than `v_p(n)`, the
limiting class factor of the class of `N` modulo `m` is `prod_{p | m} (1 - p^-2) sigma(p^(v_p(n)))/p^(v_p(n))`. Along
any sequence of such `m` whose prime factors exhaust the odd primes it tends to

    S_sigma(N) = prod_{p odd} (1 - p^-2) sigma(p^(v_p(n)))/p^(v_p(n)) = (8/pi^2) sigma(n)/n.

`S` has mean 1 over the odd `N` of either sheet and minimum `8/pi^2 = 0.8106` (at `n = 1`: on sheet 1 only `N = 1`;
on sheet 0 exactly when `N + 1` is a power of 2); it is unbounded, of maximal order `(4 e^gamma/pi^2) ln ln N`
(Gronwall's theorem for odd `n`; CITED). Over both sheets the density is `w_0 S_0 + w_1 S_1` (triple shares of
THM-5128 (D1)), not an Euler product: `0.876` at `N = 7` and `1.084` at `N = 21`.
This is a statement about limits of class factors at fixed moduli. It is not an asymptotic for the number of
representations of a single `N` (which would need uniformity in the modulus, HYP-9660 (B), OPEN). NUMERICAL
(exploratory, series_vs_counts_25_24.out): on the odd `N` of `[2^24, 2^25)` the correlation of the multiplicity
with `sigma(n)/n` is `0.40` (sheet 0) and `0.43` (sheet 1), and the mean multiplicity doubles from the lowest to
the highest quintile of `sigma(n)/n`.

**(F) The four-square series and Ingham's series.** For odd `n`, `(8/pi^2) sigma(n)/n = r_4(n)/(pi^2 n)`, where
`r_4(n) = 8 sigma(n)` is the number of representations of `n` as a sum of four squares (Jacobi; CITED) and `pi^2 n`
the real density of that problem: `S_sigma(N)` is the product of the densities of the four-square form at the odd
primes. The reason is local: for an odd prime `p` and a unit `u`, the quaternary forms `a1 a2 + u b1 b2` and
`a^2 + b^2 + c^2 + d^2` are equivalent over `Z_p`, so they take every class modulo every `p^k` equally often; and
the states of F1 are equidistributed modulo `p^k` (THM-5128 (B2)). `a1 a2 + b1 b2` is the determinant of
`[[a1, -b1], [b2, a2]]`, the norm form of the split quaternions `M_2`; the sum of four squares is the norm form of
Hamilton's quaternions, which split at every odd prime and ramify only at 2 and at the real place (MISTAKE-506
records the same fact). At those two places the problems differ: the four-square problem has 2-adic density 1 at
odd `n` and real density `pi^2 n`; F1 has the flat law on the odd classes modulo `2^j` (THM-5128 (E)) and the mass
`2 kappa_sigma ln N` per odd `N` (THM-5128 (D1)).
The classical theorem with the same form is Ingham's (1927; CITED): the number of solutions of
`a1 a2 + b1 b2 = n` in positive integers is `sum_{m < n} d(m) d(n - m) ~ (6/pi^2) sigma(n) ln^2 n`, with singular
series `sigma(n)/(zeta(2) n)` over all primes; `(8/pi^2) sigma(n)/n` is that series without its factor at 2
(`(6/pi^2)/(1 - 1/4) = 8/pi^2`). F1 is a thin analogue (one parameter per state, mass `ln N` instead of
`n ln^2 n`); Ingham's theorem does not apply to it and is not used.

**(H) Every base.** For an integer `b >= 2` let `N_b(h1, h2, s) = I(Omega(J(h1, h2, s)), b)` be the value at `b` of
the odd-cycle packing polynomial of the join (the independence polynomial of its conflict graph; `N_2 = N` by OCF;
THM-4940). Then (THM-5070 (2) in the present coordinates; checked on 314 joins)

    N_b = X1 X2 + b (2^s y1 y2 - [s = 0])/2,          X_i = 1 + b psi(h_i),

a polynomial `1 + alpha_1 b + alpha_2 b^2` with `alpha_1 = psi1 + psi2 + (2^s y1 y2 - [s = 0])/2`, `alpha_2 = psi1 psi2`.
For every odd `m` prime to `b`, (A)-(E) hold for `N_b` with `x` replaced by `X` (types by `gcd(X(h), y(h), m)`),
the same `c_d` and `L_m`, and the distinguished classes `c_1 = 0`, `c_0 = -b/2 (mod m)`: the limiting class factor
of `N_b = c` is `L_m(c)` on the sheet `s >= 1` and `L_m(2c + b)` on the sheet `s = 0`; the constants `kappa`, `A`
depend on `b`. At the primes of `b`: `N_b = 1 (mod b)` (THM-4940). The states with `X = y` (which are of type `p`
for every `p | y`) are `h = 0` and the powers of 2 for `b = 2` (THM-5126 (A)); for `b = 3` they are
`h = 0, 3, 45, 117, 309, 9045, ...` (12 below `2^24`; FINITE-EXACT); for `b >= 4` only `h = 0`, because
`psi(h) > h/2` for `h >= 1`. Tao's digit condition (the note, section 5): for `b > max(alpha_1, alpha_2)` the
base-`b` numeral of `N_b` is `(alpha_2)(alpha_1)(1)`.

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
The factors are still far from their limits: at 9, `X = 2^28`, the classes of valuation `0, 1, >= 2` have
`0.850 / 1.225 / 1.451` (sheet 0) and `0.847 / 1.229 / 1.462` (sheet 1) against the limits
`0.889 / 1.185 / 1.296`, as (D) predicts with its `1/ln X` term.
Pre-registered (PREREGISTRATION.md; window of odd `N` in `[2^27, 2^28)`, a-priori predictions of `Phi(g0)` from the
smaller states alone): P1 group means within the registered tolerances (largest deviations `0.027` at 9, `0.028` at
27, `0.032` at 81, `0.004` at 25, `0.004` at 125, `0.091` at 49, `0.030` at 45, `0.021` at 75) - CONFIRMED; P2
flatness inside each group - CONFIRMED; P3 the measured factor is closer to the a-priori prediction than to the
limit for every group with `g0 > 1` at 9, 27, 25 (14 of 14) - CONFIRMED; P4 the remainder of (B) is smaller at
`2^28` than at `2^26` for at least 8 of 12 pairs - CONFIRMED with exactly 8 (it grew at 9 sheet 0, 27 sheet 1,
81 sheet 1, 125 sheet 0). The ordering of the factors along divisor chains fails at finite size and was not
registered.
(H), exploratory: sheet `s >= 1`, `N_b <= 2^24`, `b = 2, 3, 4, 5, 6, 10`, primes `3..13` not dividing `b`: remainder of
the first-order law at most `0.13` (largest at 7, as for `b = 2`); sheet `s = 0`, `N_b <= 2^22`, `b = 2, 3, 5, 10`:
every triple whose smaller state is `p`-divisible has `2 N_b + b = 0 (mod p)` (exact), remainder at most `0.16`
(every_base.out, bridges_exact.out S10). In base 3 the prime 7 is boosted early by the balanced state `h = 3`
(`X = y = 7`), in base 5 the prime 13.

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
form with uniform variables is THM-5135 (1) with `r = 2`; so it is also the `p`-adic density of the four-square
form, and Jacobi's `r_4(n) = 8 sigma(n)` (odd `n`) with real density `pi^2 n` and 2-adic density 1 gives the same
product (class counts compared modulo 9, 27, 25, 49, 121; `r_4 = 8 sigma` for odd `n <= 400`).

**(H).** The formula: `I(Omega(R(h)), z) = 1 + psi(h) z` and
`I(Omega(J), z) = (1 + psi1 z)(1 + psi2 z) + z (2^s y1 y2 - [s = 0])/2` (every odd cycle of a one-apex tournament
passes through the apex, so cycles of the same side meet; the cycles through the reversed arc meet every other
cycle); it is THM-5070 (2) in the coordinates `a(2h) = psi(h)` and is checked by brute force on 314 joins of at
most 13 vertices. For `s >= 1`, `N_b = X1 X2 + b 2^(s-1) y1 y2`; for `s = 0`, `2 N_b + b = 2 X1 X2 + b y1 y2`. Since
`X = 1 + b (x - 1)/2`, in both cases the value is `a x2 + b' y2 + const` with `a > 0`, `b' > 0` depending on the
smaller state, a region of the kind of THM-5128 (B4) (which holds for all `a > 0`, `b' >= 0`), and the linear form
in `(X2, y2)` has coefficients `(X1, b 2^(s-1) y1)` resp. `(2 X1, b y1)`, of content `gcd(X1, y1, m)` modulo an odd
`m` prime to `b`. `(X, y) = (1 + b psi, 2h + 1)` is the image of `(q, h)` under an affine bijection of `(Z/m)^2`
(`b` and 2 are units), so THM-5128 (B2) gives its equidistribution on cells, and (B4)(ii)-(iv) hold for it. The
proofs of (A)-(E) then apply; the weights of (C) are again smooth functions of `y` and `rho` (the ratio `X/y` is
`b rho/2 + (1 - b/2)/y`). The class of `N_b` on sheet 0 is that of `2 N_b + b`. Balanced states: `q(h) < h/2` for
`h >= 1` (induction on the recursion: `q(1) = 0`, `q(2h) = q(h)`, `q(2h + 1) = q(h) + ceil(h/2)`), so
`psi(h) = h - q(h) > h/2` and `1 + b psi(h) > 2h + 1` for `b >= 4`, `h >= 1`.

## Remarks

1. *What changed with respect to THM-5128.* Its (D3) product `prod_p (1 - 1/p^2 + [p | n]/p)` is the squarefree
   level: the mean of `L_(p^k)` over the classes divisible by `p`. The series at `N` is `(8/pi^2) sigma(n)/n`:
   smaller by `p^-3` at each prime dividing `n` exactly once, larger from `p^2 | n` on.
2. *Finite size.* At every computed size the factors are governed by the `1/ln X` term of (D), i.e. by the few
   smallest states of each type: at 25 the classes with `v = 1` and `v >= 2` have nearly the same factor (`1.25`),
   because the state `h = 2` (`x = y = 5`) is of type 5 and no state of type 25 is small.
3. *Where the divisors come from.* `sigma(n)/n = sum_{d | n} 1/d`: the divisor `d` is the type of the smaller state,
   present with density `d^-2` and concentrated on a share `d` of the classes. `sigma(n)` is also the number of
   sublattices of index `n` in `Z^2` (Hermite normal forms of determinant `n`; classical).
4. *Even moduli.* For `2^j m` with `j >= 2` and `m > 1` odd the natural statement is the product of the flat law on
   the odd classes modulo `2^j` and of `L_m`. It needs the joint equidistribution of `q` modulo `2^(j-1)` and of the
   states modulo `m` on cells, that is `rho(T_k) < 2` for the characters of `Z/(2^(j-1) m)`, which is not proved
   here (HYP-9665).
5. *No uniformity.* All constants depend on `m`; nothing here is a level of distribution.
6. *Session265's "prime-power lifts".* They are reachability statements for the moving-mark source state modulo
   `p^k` (exceptions 5 and 761 of that proof); they carry no measure and are unrelated to the densities above. F1
   has no exceptional odd prime: 25 and 125 follow the same law as 9, 27, 81.
