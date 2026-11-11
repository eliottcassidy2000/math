---
id: THM-5134
title: "Local densities of the one-arc join family F1 at every odd modulus. (A) For a fixed smaller state with d = gcd(x1, y1, m), the value M = x1 x2 + 2^s y1 y2 is uniform on the multiples of d modulo m, with a power saving. (B) First-order law: #{N = c mod m} = sum over d | gcd(m, c - c_s) of (d/m) #{triples whose smaller state has gcd(g, m) = d} + O(X^(1 - theta_m/2)). (C) The shares of the types tend to d^-2 prod_{p | m/d} (1 - p^-2) at the rate 1/ln X. (D) Limit: the class factor m * share tends to L_m(e) = sum_{d | gcd(m, e)} d^-1 prod_{p | m/d} (1 - p^-2), e = c - c_s, multiplicative over prime powers, equal to (1 - p^-2) sigma_1(p^v)/p^v at a class of valuation v < k modulo p^k. (E) Singular series (8/pi^2) sigma_1(n)/n, n the odd part of N (sheet s >= 1) or of N + 1 (sheet s = 0). (F) It is the series of Jacobi's four-square theorem, r_4(n)/(pi^2 n) for odd n, and Ingham's series for a sum of two products with the prime 2 removed. (G) Every base: the base-b packing counts N_b of the same tournaments obey (A)-(D) at every odd modulus prime to b, with the class -b/2 in place of -1 on the sheet s = 0, and (E) with the product over the primes not dividing 2b; for odd b the second base-b digit of N_b is equidistributed by value on each sheet, with a power saving."
status: >
  AUDITED 2026-11-11: independent audit by eight auditors, PASS WITH CORRECTIONS (6 MAJOR, 64 MINOR, 41 NOTE; all
  applied; MISTAKE-763); report: 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/AUDIT_REPORT.md.
  PROVED: (A)-(E), by the argument of THM-5128 (C), (D) carried out for an odd modulus with THM-5128 (B2),
  (B4)(iv); this is the step that THM-5128 (D3) left "not written out" and typed NUMERICAL after its audit.
  PROVED: (G). Its class law in base b is a transfer of the proofs of (A)-(E): the region, the content of the
  linear form, the equidistribution of the base-b states, the weights with their regularity, the treatment of the
  additive constant and the mass constants kappa_sigma(b) are written out (the explicit region, the weights, their
  regularity, the treatment of the additive constant and the mass constants were supplied by the audit, auditor
  proof-5134-FH, which found no error and confirmed the constants against complete censuses); the
  block-by-block summations are those of THM-5128 (D1), (D2) with the base-b weights and are referred to, not
  repeated line by line. Its second-digit clause (odd b) is proved in full; statement and proof were found by the
  audit (auditor "hypotheses"; it was the conjecture HYP-9665 (C) of the draft). For b = 2 the second digit is
  THM-5128 (E) with j = 2; for even b > 2 it is HYP-9665 (C), OPEN (not written out).
  (F): the identification is PROVED (two unimodular quaternary forms of the same discriminant over Z_p, p odd).
  CITED: the classification of unimodular quadratic forms over Z_p, p odd (rank and discriminant modulo squares;
  Cassels, Rational Quadratic Forms, ch. 8; O'Meara, Introduction to Quadratic Forms, section 92); Jacobi's
  four-square theorem; Siegel's formula for a genus of one class; Ingham's theorem (1927) and Estermann's secondary
  terms (1930); Gronwall's theorem (1913), used in its form for odd n.
  KNOWN, not claimed: the fixed-smaller-state law of (A) at every fixed odd modulus, composite and prime-power
  included, without a rate (Session206, weighted parameter capacity: uniform on the residues allowed by
  g = gcd(R, d, p), error o(X); the dictionary is written in (A)); (A)-(D) at primes and squarefree odd moduli,
  with the power saving (THM-5128 (C1), (C2) with its product law, (D2), (D3)); the regions lemma with a power
  saving at every odd modulus (THM-5128 (B4)(iv)) and the summation over the smaller state (THM-5128 (C2), (D2));
  the mechanism at p^k (stated in THM-5128 (D3), from its audit, finding proof-5128-CF#0); the exact factor
  (1 - p^-2) sigma_1(p^v)/p^v for uniform states and the series (8/pi^2) sigma_1(n)/n (audit of THM-5128,
  2026-11-11, which enumerated 9, 25, 27, 49; its report calls the law "the determinant-form density", the
  docstring of its script audit/citations/check_singular_series_prime_powers.py "the classical determinant-form
  density"); valuation-resolved prime-power factors for another family of path counts (Session223: two root
  changes, 3-adic f(v), fixed length); the base-b packing polynomial of the joins (THM-5070 (2), kind-pasteur S18);
  its form for a fixed smaller state, the base-x tail law h = h0 + (x/2)[R0 N(I) - h_u Q(I)] with
  (h0, R0, h_u) = (x_b(h1) + (b/2)(2^s y1 - [s = 0]), 2 (x_b(h1) + 2^s y1), x_b(h1)) (THM-5085 (2), kind-pasteur
  S19); N_b = 1 (mod b) (THM-4940); the balanced states of base 2 (THM-5126 (A), THM-5128 (F)).
  NEW here: the written proof of the by-value law (A)-(D) at odd moduli that are not squarefree, which THM-5128
  (D3) announced as not written out: the fibre count for a smaller state of content d (image of order m/d), the
  types d | m with the shares c_d = J_2(m/d)/m^2, the closed form and the multiplicativity of L_m with its two
  prime-power cases (the value at v >= k included); hence (E), as a statement about limits along m; the censuses;
  the identification (F) inside the repository (Jacobi's r_4 was cited in THM-4826, in Session116 and in an old
  theta-series output, 05-knowledge/results/theta_lattice_tournament.out; never for path counts); in (G): the class
  law at odd moduli prime to b, the classes 0 and -b/2, the series over the primes not dividing 2b, the balanced
  states in bases 3 and >= 4 and, both found by the audit, the constants kappa_sigma(b) and the second-digit law
  for odd b.
  NUMERICAL / FINITE-EXACT support: exact remainder of (B) at X = 2^22, 2^24, 2^26, 2^28 for twelve moduli (and at
  2^30 in the audit's census); four PRE-REGISTERED tests (window [2^27, 2^28) and the remainder at 2^28; the
  registration commit 73c8437aac precedes the commit of the outputs), all CONFIRMED as registered, the fourth at
  its threshold (8 of 12). The registration was blind only in part (prime-level numbers of the same window were in
  the repository), P1 has no power at 49, the comparison of P3 fails at 49 (where it was not registered), and the
  count of P4 is not stable under the choice of sizes: see "Numerical support". Base-b censuses (exploratory, not
  pre-registered; complete since the audit: the two censuses of the draft were truncated).
  Not uniform in the modulus (the uniform gap of HYP-9660 (B) is OPEN, and a statement about single values needs a
  second-moment input beyond it); the by-value law at even moduli 2^j m, j >= 2, m > 1 odd, is not written out
  (HYP-9665 (A); every ingredient is available; for j = 1 it is (D), N being odd).
  Chan-Pak completeness remains OPEN: the local densities are positive in every class, 7 and 21 included, and
  decide nothing about it.
session: opus-2026-11-11-S92
date: 2026-11-11
source: 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/README.md
scripts:
  - 04-computation/experiments/s92_prime_powers_20261111/odd_modulus_census.py (+ cum_22.out, cum_24.out, cum_26.out, cum_28.out, window_insample_*.out, window_prereg_28_27.out), PREREGISTRATION.md, prereg_window_28_27.json
  - 04-computation/experiments/s92_prime_powers_20261111/seam_rank_exact.py (+ .out; its two parts are labelled THM-5135 (1) and THM-5134 (F))
  - 04-computation/experiments/s92_prime_powers_20261111/every_base.py (+ .out; the receipt of (G), a complete census since the audit), bridges_exact.py (+ .out, sections S1, S4, S10; S10 is a complete census since the audit), series_vs_counts.py (+ series_vs_counts_25_24.out)
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/ (AUDIT_REPORT.md; the audit receipts quoted in this file are in its folders proof-5134-AE, proof-5134-FH, proof-5135, numerics-prereg, hypotheses, citations, prior-art; the checks made after the corrections were applied are in its folder verify-THM-5134)
related:
  - 01-canon/theorems/THM-5128-f1-by-value-local-densities-first-order-class-law.md (the prime and squarefree case; (B2), (B3), (B4), (C2), (D1), (D2), (E) are used)
  - 05-knowledge/results/s91_f1_local_densities_20261110/audit/AUDIT_REPORT.md (findings proof-5128-CF#0, citations#2, prior-art-wording#12 on the prime-power factor)
  - 05-knowledge/results/s91_f1_local_densities_20261110/audit/citations/check_singular_series_prime_powers.py (the phrase quoted in the status)
  - session206/weighted-parameter-capacity/WEIGHTED_CAPACITY.md (fixed-source law at a fixed odd modulus)
  - 05-knowledge/results/session223_root_changes_and_local_factors_20261107/README.md (valuation-resolved 3-adic factors)
  - 01-canon/theorems/THM-5070-port-chain-transfer-for-odd-cycle-packing-counts-and-ternary-coverage-to-a-million.md (base b; the packing polynomial of the joins)
  - 01-canon/theorems/THM-5085-fixed-port-tails-in-every-base-and-the-ternary-bank-to-ten-million.md (base b; the tail law for a fixed smaller state)
  - 01-canon/theorems/THM-4940-odd-cycle-packing-counts-in-every-base.md (base b; N_b = 1 mod b)
  - 01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md (balanced states of base 2)
  - 01-canon/theorems/THM-4826-which-grid-carries-a-norm-shell-one-to-two-split-cubic-character-of-two-centred-units-complete-the-quartic-designs-four-eisenstein-fourth-powers-miss-five-mod-eight.md (Jacobi, Hurwitz, Lipschitz in canon)
  - 05-knowledge/results/theta_lattice_tournament.out (an earlier citation of Jacobi's r_4, not for path counts)
  - 01-canon/theorems/THM-5135-* (class laws of bilinear seams; the algebra behind (F))
  - 05-knowledge/hypotheses/HYP-9660-* (uniform gap, OPEN), HYP-9665-* (even moduli, single values, the second digit in even bases)
  - 05-knowledge/results/session193_triangulated_proof_20261105/modular_sync/PROOF.md (equidistribution on cells at every modulus, used in Remark 4)
---

# THM-5134 - local densities of F1 at every odd modulus, and the four-square singular series

## Setting

As in THM-5128. One-apex states `x(h) = 2 psi(h) + 1`, `y(h) = 2h + 1`, `q(h) = h - psi(h)`, `rho(h) = x(h)/y(h)`,
`g(h) = gcd(x(h), y(h))`; `x1, y1, psi1` and `x2, y2, psi2` are the values at `h1` and `h2`. F1:
`N = M - [s = 0]`, `M = x1 x2 + 2^s y1 y2`, over ordered triples `(h1, h2, s)`; sheet 0 (the head sheet) is `s = 0`
(`c_0 = -1`), sheet 1 (the tail sheet) is `s >= 1` (`c_1 = 0`); "sheet" alone always means one of these two;
`T_sigma(X)` is the set of triples of the sheet `sigma` with `N <= X`; `Y = X/2^s`.
`m >= 3` is an odd modulus, `rho_m` and `theta_m = ln(1/rho_m)/ln(2/rho_m^2)` are as in THM-5128 (B). For `d | m`:

- a state is of *type* `d` if `gcd(g(h), m) = d`; a triple with `h1 != h2` is of type `d` if its state of smaller
  word is; `D_sigma^(d)(X)` is the number of triples of type `d` in `T_sigma(X)`;
- `c_d = d^-2 prod_{p | m/d} (1 - p^-2)`; this is the share of the residue vectors `(alpha, beta)` modulo `m` with
  `gcd(alpha, beta, m) = d` (Jordan's totient `J_2(m/d)` over `m^2`), and `sum_{d | m} c_d = 1`;
- `chi_d(h) = [gcd(g(h), m) = d] - c_d`;
- `L_m(e) = sum_{d | gcd(m, e)} d c_d = sum_{d | gcd(m, e)} d^-1 prod_{p | m/d} (1 - p^-2)` for a class `e` modulo `m`;
- `sigma_1(n)` is the sum of the divisors of `n`; `sigma` alone always denotes a sheet (`c_sigma`, `T_sigma`,
  `S_sigma`, `kappa_sigma`), as in THM-5128;
- `X` always denotes the bound (`N <= X`). In (G) the base-`b` state is written `x_b(h) = 1 + b psi(h)` (so
  `x_2 = x`), with `g_b(h) = gcd(x_b(h), y(h))`.

## Statements

**(A) Fixed smaller state.** Let `h1`, `s` be fixed, `d = gcd(g(h1), m)` and `G = #{h2 > h1 : M <= X}`. For every
class `c'` modulo `m`:

    #{h2 > h1 : M <= X, M = c' (mod m)} = [d | c'] (d/m) G + O_m((Y/y1)^(1 - theta_m)),

uniformly in `h1`, `s`, `X` and `c'`. If `d = m` then `M = 0 (mod m)` for every `h2`. (Without the power saving,
with error `o(X)`, this is Session206's weighted-capacity theorem; KNOWN. Dictionary, checked by the audit
(prior-art/p2_session206_dictionary.out): with `Q = 2q`, the map `h2 -> N(h1, h2, s)` is the value
`h_S + R_S h2 - d_S Q(h2)` of the physical source `W O^s J O` of the O/Z/J grammar (for `h1 >= 1`, `W` is the root
O/Z word ending in `Z` that corresponds to `h1`), where Session206's source parameters `(h, R, e, d)` are
`(h_S, R_S, e_S, d_S) = (x1 + 2^s y1 - [s = 0], 2 (x1 + 2^s y1), 0, x1)`; and `gcd(R_S, d_S, m) = gcd(x1, y1, m)`.)

**(B) First-order law.** For every class `c` modulo `m` and each sheet:

    #{t in T_sigma(X) : N = c (mod m)} = sum over d | gcd(m, c - c_sigma) of (d/m) D_sigma^(d)(X) + O_m(X^(1 - theta_m/2)).

So all the classes with the same `g0 = gcd(m, c - c_sigma)` have the same share up to the error, and the class
factor is `m * share = Phi(g0) + O(X^(-theta_m/2))`, `Phi(g0) = sum_{d | g0} d F_d`, `F_d = D^(d)/#T` the share of the
triples of type `d`. Every triple of type `d` has `d | N - c_sigma` exactly.

**(C) The types.** `D_sigma^(d)(X) = c_d #T_sigma(X) + 2 A_(m,d)^sigma X + O_m(X^(1 - theta_m/2))`, with
`A_(m,d)^sigma = lim_{H -> infinity} sum_{h < H} chi_d(h) w^sigma(h)` and the weights `w^0`, `w^1` of THM-5128 (D2).
This is a sum over ALL states of the centred indicator `chi_d`, conditionally convergent: it is not absolutely
convergent, and the sum over the states of type `d` alone diverges like `c_d kappa_sigma ln H`; the limit exists as
in THM-5128 (D2). Hence `F_d -> c_d` at the rate `1/ln X`.

**(D) Limit.** For every class `c` modulo `m`:

    m #{N = c} / #T_sigma(X) = L_m(c - c_sigma) + (2/kappa_sigma) (sum_{d | gcd(m, c - c_sigma)} d A_(m,d)^sigma) / ln X + O(1/ln^2 X).

`L_m` is multiplicative: `L_m(e) = prod_{p^k || m} L_(p^k)(e)`, and with `v = v_p(e)`

    L_(p^k)(e) = (1 - p^-2) sigma_1(p^v)/p^v = (1 + 1/p)(1 - p^(-v-1))                  if v < k,
    L_(p^k)(e) = (1 - p^-2) sigma_1(p^(k-1))/p^(k-1) + p^-k = 1 + 1/p - p^(-k-1)        if v >= k.

For `k = 1` this is the factor `1 - 1/p^2 + [p | e]/p` of THM-5128 (D3); the mean of `L_m` over the classes is 1.
`m^3 L_m(c)` is the number of `(a1, a2, b1, b2)` modulo `m` with `a1 a2 + b1 b2 = c`, the additive square of the
multiplication table modulo `m` (THM-5135 (1), which gives both cases, the class 0 modulo `p^k` included; checked
for eleven odd `m <= 105`, bridges_exact.out S1). That identity holds for every `m`, even `m` included, with `L_m`
defined by the same formula (audit: every `m <= 210`, proof-5135/a1_rank_law.out [L]); the present theorem needs
odd `m` only because the states are odd, hence not equidistributed, modulo 2.

**(E) Singular series.** Let `N` be odd, `sigma` a sheet with `N + [sigma = 0] != 0`, and `n` the odd part of
`N + [sigma = 0]`. For every odd `m` in which every prime `p | m` occurs to a power higher than `v_p(n)`, the
limiting class factor of the class of `N` modulo `m` is `prod_{p | m} (1 - p^-2) sigma_1(p^(v_p(n)))/p^(v_p(n))`.
Along any sequence of such `m` in which every odd prime divides all but finitely many terms (for instance a chain
`m_1 | m_2 | ...` whose prime factors exhaust the odd primes) it tends to

    S_sigma(N) = prod_{p odd} (1 - p^-2) sigma_1(p^(v_p(n)))/p^(v_p(n)) = (8/pi^2) sigma_1(n)/n.

It is not enough that the prime factors of the terms exhaust the odd primes: along `m_j = p_j^k` (`p_j` the `j`-th
odd prime, `k > v_(p_j)(n)`) the class factor `(1 - p_j^-2) sigma_1(p_j^v)/p_j^v` tends to 1.
`S_sigma` has mean 1 over the odd `N`, for either sheet. As a function of odd `N`, `S_1 >= 8/pi^2 = 0.8106` with
equality only at `N = 1`, which is not a value of sheet 1 (there `N >= 3`): over the odd `N >= 3`, `8/pi^2` is an
infimum of `S_1`, approached along the primes `N`, not a minimum. `S_0 = 8/pi^2` exactly when `N + 1` is a power of
2. These are statements about the functions `S_sigma`: which odd `N` are values of a sheet is not addressed here.
`S_sigma` is unbounded, of maximal order `(4 e^gamma/pi^2) ln ln N`: over odd `n`,
`limsup sigma_1(n)/(n ln ln n) = e^gamma/2` (CITED: Gronwall's theorem, Trans. Amer. Math. Soc. 14 (1913), 113-122,
is the statement with `e^gamma` over all `n`; for odd `n` the same proof applies with the prime 2 removed from the
Mertens product; the reference is the audit's, not re-read, and the odd value `e^gamma/2 = 0.8905` is confirmed
numerically by the audit, citations/c3_ingham_gronwall_constants.out G1). Over both sheets the density is
`w_0 S_0 + w_1 S_1` (triple shares of THM-5128 (D1)), not an Euler product: `0.876` at `N = 7` and `1.084` at
`N = 21`.
This is a statement about limits of class factors at fixed moduli. It is not an asymptotic for the number of
representations of a single `N`: that would need at least a version of (B) uniform in the modulus - the uniform
gap of HYP-9660 (B) with a norm bound polynomial in the modulus, OPEN - and, beyond it, a second-moment (minor-arc)
input, since a law uniform in the modulus gives at most statements on average over `N`; and F1 has misses at every
scale tested (THM-5128 Remark 3). NUMERICAL (exploratory, series_vs_counts_25_24.out): on the odd `N` of
`[2^24, 2^25)` the correlation of the multiplicity `r` with `sigma_1(n)/n` is `0.40` (sheet 0) and `0.43` (sheet 1),
and the mean multiplicity doubles (ratio `1.99 / 1.98`) from the lowest to the highest quintile of `sigma_1(n)/n`,
whose own means are `1.001` and `1.609` (ratio `1.61`): the dependence is steeper than the series at this size.
The least-squares slope of `r/mean(r)` against `S_sigma(N)` is `1.44 / 1.43` instead of 1 (audit,
proof-5134-AE/series_vs_counts_audit_25_24.out), as the positive `1/ln X` terms of (D) at 3 and 5 predict; the
series explains 16 % and 18 % of the variance of `r`, and `8.07 % / 4.86 %` of the odd `N` of the window have
`r = 0` on the sheet although `S_sigma >= 0.81`: single values are dominated by fluctuations that the class laws
do not see.

**(F) The four-square series and Ingham's series.** For odd `n`, `(8/pi^2) sigma_1(n)/n = r_4(n)/(pi^2 n)`, where
`r_4(n) = 8 sigma_1(n)` is the number of representations of `n` as a sum of four squares (Jacobi; CITED) and
`pi^2 n` the real density of that problem: `S_sigma(N)` is the product of the densities of the four-square form at
the odd primes. The reason is local: for an odd prime `p` and a unit `u`, the quaternary forms `a1 a2 + u b1 b2` and
`a^2 + b^2 + c^2 + d^2` are equivalent over `Z_p`, so they take every class modulo every `p^k` equally often; and
the states of F1 are equidistributed modulo `p^k` (THM-5128 (B2)). `a1 a2 + b1 b2` is the determinant of
`[[a1, -b1], [b2, a2]]`, the norm form of the split quaternions `M_2`; the sum of four squares is the norm form of
Hamilton's quaternions, which split at every odd prime and ramify only at 2 and at the real place (MISTAKE-506
records the same fact). At those two places the problems differ: the four-square problem has 2-adic density 1 at
odd `n` and real density `pi^2 n`, and `r_4(n) = pi^2 n * 1 * (8/pi^2) sigma_1(n)/n` is the product of its real
density and of its local densities (Siegel's formula for a genus of one class; CITED); F1 has the flat law on the
odd classes modulo `2^j` (THM-5128 (E)) and the mass `2 kappa_sigma ln N` per odd `N` (THM-5128 (D1)).
The classical theorem with the same form is Ingham's (CITED): the number of solutions of `a1 a2 + b1 b2 = n` in
positive integers is

    sum_{k < n} d(k) d(n - k) = (6/pi^2) sigma_1(n) ln^2 n + O(sigma_1(n) ln n ln ln n)

(A. E. Ingham, Some asymptotic formulae in the theory of numbers, J. London Math. Soc. 2 (1927), 202-208; the
secondary terms are due to T. Estermann, On the representations of a number as the sum of two products, Proc.
London Math. Soc. (2) 31 (1930), 123-133; both references through secondary sources, the originals were not read).
Ingham's proof is elementary; the Euler-product reading of his constant,
`sigma_1(n)/(zeta(2) n) = prod_p (1 - p^-2) sigma_1(p^(v_p(n)))/p^(v_p(n))` over all primes, is THM-5135 (1) with
`r = 2` at every prime, and for odd `n` the series `(8/pi^2) sigma_1(n)/n` is that product without its factor `3/4`
at 2 (`(6/pi^2)/(1 - 1/4) = 8/pi^2`). The convergence is slow, the relative error being of order `ln ln n / ln n`:
the ratio of the count to the main term is `0.70` at `n = 720720` and `1.17` at `n = 999983` (audit,
citations/c3_ingham_gronwall_constants.out I1). F1 is a thin analogue (one parameter per state, mass `ln N` instead
of `n ln^2 n`); Ingham's theorem does not apply to it and is not used.

**(G) Every base.** The parts are: the tournaments, the count, the class law, the series, the primes of `b`, the
second digit, balanced states, digits.

*The tournaments.* `R(h)` is the one-apex tournament of the word `h`: a transitive tournament on `c_1, ..., c_k`
(`c_i -> c_j` for `i < j`; `k` the bit length of `2h`) with an apex `u`, where `u -> c_i` if and only if the `i`-th
binary digit of `2h`, most significant first, is 1 (`R(0)` is the single vertex `u`). `J(h1, h2, s)` is `R(h1)`,
then a transitive tournament on `s` vertices, then the converse of `R(h2)`: every vertex beats every vertex of a
later block, except that the apex of the last block beats the apex of the first. `J(h1, h2, s)` has
`N(h1, h2, s)` Hamiltonian paths.

*The count (KNOWN).* For an integer `b >= 2` let `N_b(h1, h2, s) = I(Omega(J(h1, h2, s)), b)` be the value at `b` of
the odd-cycle packing polynomial of the join (the independence polynomial of its conflict graph; `N_2 = N` by OCF;
THM-4940). Then, with `x_b(h) = 1 + b psi(h)`,

    N_b = x_b(h1) x_b(h2) + b (2^s y1 y2 - [s = 0])/2,

a polynomial `1 + alpha_1 b + alpha_2 b^2` with `alpha_1 = psi1 + psi2 + (2^s y1 y2 - [s = 0])/2`, `alpha_2 = psi1 psi2`.
This is THM-5070 (2): `J(h1, h2, s)` is the apex bridge `T (.) B` of its one-apex sides with the words
`W = 2^(s1) (2 h1 + 1) - 1` and `W' = 2^(s2) (2 h2 + 1) - 1`, `s1 + s2 = s` (the padding enters as `s` trailing ones
shared between the two words); `a(W) = psi(h1)`, `a(W') = psi(h2)` (in particular `a(2h) = psi(h)`), and the bracket
of THM-5070 (2) equals `(2^s y1 y2 - [s = 0])/2`. For a fixed smaller state it is the base-`x` tail law of
THM-5085 (2), `h = h0 + (x/2) [R0 N(I) + Delta0 [1 in I] - h_u Q(I)]` at `x = b`, with
`(h0, R0, Delta0, h_u) = (x_b(h1) + (b/2)(2^s y1 - [s = 0]), 2 (x_b(h1) + 2^s y1), 0, x_b(h1))`, `N(I) = h2` and
`Q(I) = 2 q(h2)`. (Checked by brute force on 314 joins of at most 13 vertices, every_base.out; by the audit on
2,816 joins, with the isomorphism to `T (.) B` in 212 cases, proof-5134-FH/h1_packing_bruteforce.out and
h1b_session_construction.out, and on 134,400 cases of the THM-5085 dictionary, prior-art/p5_thm5085_dictionary.out.)

*The class law.* Let `m` be odd and prime to `b`. The *value* of a triple is `V = N_b` on the sheet `s >= 1` and
`V = 2 N_b + b` on the sheet `s = 0`; a state is of type `d` if `gcd(g_b(h), m) = d`; `T^b_sigma(X)` is the set of
triples of the sheet `sigma` with `N_b <= X`. Then (A)-(D) hold for the base-`b` family with `V` in place of `M`,
`N_b <= X` in place of `M <= X`, these types, the same `c_d` and `L_m`, the mass
`#T^b_sigma(X) = kappa_sigma(b) X ln X + O_b(X)` and the weights

    w^0_b(h) = (2/(b x_b)) ln(2 (x_b + y)/(x_b + 2y)),        w^1_b(h) = (2/(b x_b)) ln(1 + x_b/(2y)),
    kappa_0(b) = (4/b^2)(l_1 - 2 l_2 + l_3),        kappa_1(b) = (4/b^2)(l_2 - l_3),        l_j = -Li_2(-b/2^j)

(`x_b = x_b(h)`, `y = y(h)`; `kappa_0, kappa_1 = 0.091805, 0.132470` for `b = 3` and `0.013661, 0.028354` for
`b = 10`; for `b = 2` these are the weights and the constants of THM-5128 (D)). All error terms and the constants
`A_(m,d)` depend on `b` and on `m`. In particular every triple of type `d` has `d | V` exactly, the distinguished
classes are `c_1 = 0` and `c_0 = -b/2 (mod m)`, and the limiting class factor of `N_b = c` is `L_m(c)` on the sheet
`s >= 1` and `L_m(2c + b)` on the sheet `s = 0`.

*The series.* (E) holds with the product over the primes not dividing `2b`. Let `n` be the odd part of `N_b`
(sheet `s >= 1`), resp. of `2 N_b + b` (sheet `s = 0`); it is prime to `b`. For every odd `m` prime to `b` in which
every prime `p | m` occurs to a power higher than `v_p(n)`, the limiting class factor of the class of `N_b` modulo
`m` is `prod_{p | m} (1 - p^-2) sigma_1(p^(v_p(n)))/p^(v_p(n))`, and along any sequence of such `m` in which every
prime not dividing `2b` divides all but finitely many terms it tends to

    prod_{p not dividing 2b} (1 - p^-2) * sigma_1(n)/n = (8/pi^2) prod_{p | b, p odd} (1 - p^-2)^-1 * sigma_1(n)/n,

that is `(9/pi^2) sigma_1(n)/n` for `b = 3, 6` and `(25/(3 pi^2)) sigma_1(n)/n` for `b = 5, 10`. The other sentences
of (E) (the mean, the infimum, the maximal order, the sheet weights and the values `0.876`, `1.084`) are statements
about `b = 2` and are not claimed for other bases.

*The primes of `b`.* `N_b = 1 (mod b)` (THM-4940), hence `2 N_b + b = 2 (mod b)`.

*Second digit.* The second base-`b` digit of `N_b` is `(N_b - 1)/b = alpha_1 (mod b)`. For odd `b`, on each sheet
and for every class `c` modulo `b`,

    #{t in T^b_sigma(X) : alpha_1 = c (mod b)} = #T^b_sigma(X)/b + O_b(X^(1 - theta_b/2)),

`theta_b` being the exponent of the Setting for the modulus `b`: the second digit is equidistributed by value on
each sheet, with a power saving (PROVED; a corollary of THM-5128 (B4)(iv); statement and proof found by the
audit, auditor "hypotheses"; in the draft it was the conjecture HYP-9665 (C)). For `b = 2` the second digit is
`N mod 4` and the flat law is THM-5128 (E) with `j = 2`. For even `b > 2` it is HYP-9665 (C), OPEN.

*Balanced states.* A state is balanced in base `b` if `x_b(h) = y(h)`; then `g_b = y`, and the state is
`p`-divisible for every prime `p | y`. The balanced states are `h = 0` and the powers of 2 for `b = 2` (THM-5126 (A),
THM-5128 (F)); for `b = 3` they are `h = 0, 3, 45, 117, 309, 9045, ...` (12 below `2^24`, FINITE-EXACT; the audit
finds 12 below `2^28` and 13 below `2^30`, the thirteenth being `h = 659,903,829`; whether the list is finite is
OPEN); for `b >= 4` only `h = 0`, because `psi(h) > h/2` for `h >= 1`.

*Digits.* Tao's digit condition (the note, section 5; from the puzzle "The magic money plant" set by Terence Tao
for Alex Bellos's column, The Guardian, 26 October 2026; the attribution was confirmed by the audit on syndicated
copies, citations/sources.txt item 13): for `b > max(alpha_1, alpha_2)` the base-`b` numeral of `N_b` is
`(alpha_2)(alpha_1)(1)`.

**Numerical support.** Exact remainder of (B), `m max_c |#{N = c} - sum_{d | gcd(m, c - c_sigma)} (d/m) D^(d)| / #T`,
sheet 0 / sheet 1 (cum_*.out; the triple counts equal those of THM-5128 at `2^24`, `2^26`). The noise column is
`sqrt(m/#T)` with `#T` of sheet 1; the value of sheet 0 is 13 % larger (ratio `1.126`, from the triple counts
854,415,478 and 1,083,261,684 of cum_28.out) and rounds to the same entry except at 125, where it is `0.0004`:

| `m` | `X = 2^22` | `2^24` | `2^26` | `2^28` | noise at `2^28` (sheet 1) |
|---|---|---|---|---|---|
| 9 | `0.035 / 0.028` | `0.021 / 0.024` | `0.024 / 0.020` | `0.024 / 0.019` | `0.0001` |
| 27 | `0.046 / 0.044` | `0.022 / 0.029` | `0.029 / 0.022` | `0.029 / 0.027` | `0.0002` |
| 81 | `0.062 / 0.055` | `0.024 / 0.042` | `0.034 / 0.023` | `0.031 / 0.031` | `0.0003` |
| 25 | `0.004 / 0.008` | `0.004 / 0.005` | `0.004 / 0.005` | `0.004 / 0.004` | `0.0002` |
| 125 | `0.011 / 0.015` | `0.008 / 0.008` | `0.005 / 0.007` | `0.005 / 0.005` | `0.0003` |
| 49 | `0.182 / 0.140` | `0.163 / 0.141` | `0.130 / 0.125` | `0.093 / 0.103` | `0.0002` |
| 45 | `0.083 / 0.057` | `0.040 / 0.056` | `0.039 / 0.043` | `0.038 / 0.035` | `0.0002` |
| 75 | `0.097 / 0.065` | `0.050 / 0.056` | `0.034 / 0.031` | `0.027 / 0.029` | `0.0003` |

*The remainder.* It is far above the sampling noise. It decreases at 5 and 7 (49: `0.163 / 0.141` at `2^24`,
`0.093 / 0.103` at `2^28`); at `m = 9, 27, 81` it oscillates between `0.018` and `0.042` over the octaves from `2^24`
to `2^28`, without a net decrease on sheet 0 (`0.021, 0.022, 0.024` at `2^24` against `0.024, 0.029, 0.031` at
`2^28`; on sheet 1 it is `0.024, 0.029, 0.042` against `0.019, 0.027, 0.031`; audit, numerics-prereg/octaves.out).
This is the oscillation at 3, and the slow rate at 7, already seen for primes (THM-5128 (C)).
The proved exponents are small: `theta_m/2` is below `0.074` at 9, 27, 81, below `0.093` at 25, 125,
below `0.043` at the multiples of 7 and below `0.062` at 45, 75 (audit, proof-5134-AE/rates_and_states.out: rates
`0.866025`, `0.815263`, `0.932039`, `0.892191`, those of THM-5128, which give the suprema `0.0733`, `0.0927`,
`0.0422`, `0.0619` of the admissible `theta_m/2`). At `X = 2^30` (census of the audit,
proof-5134-AE/summary_tables.out; 3,658,114,303 and 4,649,378,323 triples) the remainder is `0.012 / 0.021` at 9,
`0.016 / 0.023` at 27, `0.017 / 0.026` at 81, `0.003 / 0.003` at 25, `0.003 / 0.004` at 125, `0.058 / 0.079` at 49,
`0.021 / 0.029` at 45 and `0.016 / 0.021` at 75; from `2^28` to `2^30` it decreases in 23 of the 24 pairs (modulus,
sheet) of the twelve moduli, all but 9 on sheet 1.
Every triple of type `d` lies in a class with `d | c - c_sigma` (exact, in every cumulative run; on the registered
window the session's script computes this flag without printing it, and the audit confirmed it there,
numerics-prereg/window_exactness.out). The pass criterion of the cumulative runs (every remainder below `0.25`) is
not a test of the law.
*The factors.* They are still far from their limits: at 9, `X = 2^28`, the classes of valuation `0, 1, >= 2` have
`0.850 / 1.225 / 1.451` (sheet 0) and `0.847 / 1.229 / 1.462` (sheet 1) against the limits
`0.889 / 1.185 / 1.296`, as (D) predicts with its `1/ln X` term. (NUMERICAL; the session evaluated no constant
`A_(m,d)`, the receipt is the audit's, proof-5134-AE/a_constants.out: with the partial sums of `A_(9,d)` over
`h < 2^27`, (D) gives `0.853 / 1.220 / 1.441` and `0.850 / 1.222 / 1.454` at `X = 2^28`; at 25 it gives
`0.9412 / 1.2339 / 1.2399` against the measured `0.9415 / 1.2328 / 1.2392` on sheet 0 and `0.9394 / 1.2414 / 1.2467`
against `0.9389 / 1.2432 / 1.2480` on sheet 1; at 49 the constants are not determined, the partial sums of
`A_(49,7)^0` being `-0.0151, +0.0017, +0.0080, +0.0054` at `H = 2^8, 2^14, 2^20, 2^27`, as for `A_7` in THM-5128.)
*Pre-registered tests* (PREREGISTRATION.md; window of odd `N` in `[2^27, 2^28)`, a-priori predictions of `Phi(g0)`
from the smaller states alone).
P1: group means within the registered tolerances (largest deviations `0.027` at 9, `0.028` at 27, `0.032` at 81,
`0.004` at 25, `0.004` at 125, `0.091` at 49, `0.030` at 45, `0.021` at 75) - CONFIRMED. The tolerances were set from
the in-sample maxima of the same statistic, so P1 tests the persistence of the in-sample size of the remainder,
not an error model. At 9, 27, 81, 45 they separate the first-order values from the flat law, from the limit and
from a law without the valuation split by factors `3.5` to `15.7`. At 25, 125, 75 a law without the valuation split
would also have passed (at 25 the factors at `v = 1` and `v >= 2` differ by `0.008` on sheet 0 and `0.006` on sheet
1), and at 49 the tolerance `0.115` is of the size of the predicted excess of the favoured classes (`0.097` to
`0.119`), so the test has no power there. The valuation-resolved factor is tested at 3 only (audit,
numerics-prereg/eval_prereg.out). On the window the class-level remainder of (B), which the session's window mode
does not print, is `0.0234 / 0.0268` (9), `0.0276 / 0.0369` (27), `0.0316 / 0.0431` (81), `0.0028 / 0.0035` (25),
`0.0038 / 0.0041` (125), `0.0761 / 0.0908` (49), `0.0376 / 0.0459` (45), `0.0283 / 0.0347` (75) (audit,
numerics-prereg/window_exactness.out).
P2: flatness inside each group - CONFIRMED.
P3: the measured factor is closer to the a-priori prediction than to the limit for every group with `g0 > 1` at
9, 27, 25 (14 of 14) - CONFIRMED. Of the 14 comparisons 12 are distinct: the two at (27, `g0 = 3`) are the same
sets of `N` as the two at (9, `g0 = 3`). P3 was registered for 9, 27, 25 only. At 49 the comparison fails in every
group with `g0 > 1`: in the four in-sample windows, on the registered window and cumulatively at `2^22, ..., 2^28`,
0 of 4 each time (on the window, sheet 0, `g0 = 7`: `1.1691` measured, `1.1098` first-order, `1.1195` limit), and
still at `2^30` in the census of the audit; at `X = 2^28` it holds in 6 of 10 groups at 63 and in 4 of 6 at 21
(audit, proof-5134-AE/summary_tables.out T3). So at 49 the remainder of (B) is still larger than the first-order
correction at every computed size: in every group with `g0 > 1` the measured factor is farther from the
first-order value than the first-order value is from the limit. At 63 and 21, at `X = 2^28`, the same holds in the
groups where the comparison fails (`g0 = 7`, and `g0 = 21` at 63; cum_28.out). This does not refute (B): for the
fixed smaller state `h1 = 0`, `s = 0`, the relative deviation from uniformity modulo 49 is `0.298, 0.222, 0.169` at
`X = 2^20, 2^24, 2^28`, a factor `0.93` per binary digit against the spectral rate `0.932` (audit,
proof-5134-AE/fixed_state_A.out).
P4: the remainder of (B) is smaller at `2^28` than at `2^26` for at least 8 of the 12 pairs - CONFIRMED with exactly
8 (it grew at 9 sheet 0, 27 sheet 1, 81 sheet 1, 125 sheet 0; the decreases at 27 sheet 0 and at 25 sheet 0 are of
`0.00013`, below the noise scale, and the increase at 125 sheet 0 is of `0.00007`). The verdict is read off
cum_26.out and cum_28.out; no committed program prints it (recomputed by the audit, numerics-prereg/eval_prereg.out).
The count is not stable under the choice of sizes: 4 of 12 from `2^26` to `2^27`, 10 of 12 from `2^27` to `2^28`,
11 of 12 from `2^28` to `2^30`; 8 or more of 12 has probability `0.19` under independent fair coin flips.
The ordering of the factors along divisor chains fails at finite size and was not registered. The registration
text says that the tolerance at 49 is below two of the four in-sample values; it is below three of them (`0.151`,
`0.129`, `0.118`; the fourth is `0.107`). The registered files are unchanged.
*Blindness.* The registration (commit 73c8437aac, 11:29) precedes the commit of the outputs of the registered
window and of `2^28` (commit 97d6a9f6ca, 11:43); the in-sample censuses belong to the registration commit. It was
not blind at the prime level: the audit of THM-5128 had enumerated the same window at prime moduli (its
recompute-5128/a6_window_27_28.out; its report prints the excess at 5), and THM-5128 (C) gives the remainders at
`2^26` and `2^28` for 3, 5, 7. These numbers fix the `g0 = 1` group means of P1 at 9, 27, 81, 25, 125, 49 (12 of the
68 registered group values) and the mean over the favoured classes, and the prime-level direction agrees with the
prime-power one in 9 of the 12 pairs of P4 (copying it would have predicted 9 decreases; registered threshold 8,
outcome 8). The blind content is the split of the favoured class by valuation (decisive only at 3), the flatness
inside the groups, and the composite moduli 45, 75 (audit, numerics-prereg/prior_exposure.out). The verdicts stand
as registered.
*(G), exploratory.* The two base-`b` censuses of the draft were incomplete (sheet 0: words below `2^12` only; sheet
`s >= 1`: words below `2^21` and an unsound skip of smaller states). The scripts were corrected and the outputs
regenerated after the audit; the numbers below are those of the regenerated outputs and agree with the complete
counts of the auditors.
Sheet `s >= 1`, `N_b <= 2^24`, `b = 2, 3, 4, 5, 6, 10` (57,818,155 / 35,564,915 / 24,947,862 / 18,838,778 /
14,919,260 / 7,634,213 triples), primes `3..13` not dividing `b`: the remainder of the first-order law is at most
`0.133` (at 7: `0.133, 0.067, 0.125, 0.042, 0.105, 0.060` in these bases, the first being the value of THM-5128; at
the other primes the largest printed value is `0.051`), and the second digit is flat to
`0.000, 0.040, 0.001, 0.002, 0.052, 0.003` (every_base.out).
Sheet `s = 0`, `N_b <= 2^22`, `b = 2, 3, 5, 10` (10,532,579 / 6,030,660 / 2,896,626 / 1,042,874 triples): every triple
whose smaller state is `p`-divisible has `2 N_b + b = 0 (mod p)` (exact, 0 violations); the remainder at 7 is
`0.167, 0.101, 0.078, 0.118`, and at the other primes the largest printed value is `0.057` (bridges_exact.out S10).
Second digit, odd `b` (NUMERICAL, audit, complete censuses, hypotheses/second_digit_b3.out): `max_c |b share - 1|`,
sheet 0 / sheet 1, for `b = 3` is `0.131 / 0.127`, `0.071 / 0.073`, `0.041 / 0.040`, `0.024 / 0.023` at
`X = 2^16, 2^20, 2^24, 2^28`: a factor `0.867` per binary digit, the rate `sqrt(3)/2 = 0.8660` of the digit chain at 3.
Constants (NUMERICAL, audit, proof-5134-FH/h5_kappa_base_b.out and h6_first_order_constants.out): the slopes of
`#T^b_sigma(X)/X` against `ln X` between `2^24` and `2^26` agree with `kappa_sigma(b)` to within `0.00001` for
`b = 2, 3, 4, 5, 10` on both sheets; the census values `(D^(p) - #T/p^2)/(2X)` at `X = 2^26` are close to the partial
sums over `h < 2^12` of the series with the weights `w_b`, for instance `-0.01328` against `-0.01338` (`b = 3`,
`p = 5`, sheet 0), `-0.01889` against `-0.01949` (sheet 1), `+0.03548` against `+0.03517` (`b = 3`, `p = 7`, sheet 0),
`+0.00508` against `+0.00496` (`b = 10`, `p = 7`, sheet 0). At 7 the partial sums keep drifting beyond `2^12`
(`b = 3`, sheet 0: `0.0276` over `h < 2^24`, same output), as for `A_7` in THM-5128, so these partial sums do not
determine the constants at 7. The limits of (D), (E) in base `b` are not visible at these sizes, only the
first-order structure and its constants.
Early boosts. The state `h = 3` (`psi = 2`, `y = 7`) is 7-divisible for every `b = 3 (mod 7)`: `x_b = y = 7` in base 3
(balanced) and `x_b = 21 = 3y` in base 10 (not balanced), with `F p^2 = 3.37` and `3.28` on sheet 1 at `2^24` (`F`
the share of the triples whose smaller state is `p`-divisible, as printed in every_base.out); the state `h = 6`
(`x_b = 26 = 2y`, `y = 13`; not balanced) boosts 13 in base 5 (`F p^2 = 3.38`). An early boost needs a small
divisible state, not a balanced one.

## Proofs

**(A).** `M = a x2 + b y2` with `(a, b) = (x1, 2^s y1)` (the coefficients of THM-5128 (B4); in (G) the letter `b` is
the base and the second coefficient is written `b'`), and `gcd(a, b, m) = gcd(x1, y1, m) = d` because `2^s` is a
unit modulo `m`. The homomorphism `(Z/m)^2 -> Z/m`, `w -> a w_1 + b w_2`, has image `d Z/m`, of order `m/d`, so it
takes each multiple of `d` on `m d` residue vectors and no other value. Let `S_n` be the set of the `h2` of the
block `B_n` that lie in the interval `I = (h1, infinity)` and in the region `Omega = {a x2 + b y2 <= X}`. By
THM-5128 (B4)(ii) in the form (iv), each residue vector modulo `m` is taken by `(x2, y2)` on
`m^-2 #S_n + O_m(2^(n(1 - theta_m)))` states of `S_n`; sum over the `m d` vectors and over the blocks
`2^n <= Y/y1`, as in the proof of THM-5128 (C1). The constant of the error term is uniform in `h1`, `s`, `X` and
`c'`, because THM-5128 (B4) is uniform in `a`, `b`, `X`, `I`, and here `b = 2^s y1 >= a = x1 > 0`. If `d = m` both
coefficients vanish modulo `m`.

**(B).** As for THM-5128 (C2). `M` is symmetric in `(h1, h2)`; a triple with `h1 != h2` is counted through its
smaller state, and the triples with `h1 = h2` number `O(Y^(1/2))` for each `s`. `N = c` is `M = c - c_sigma` on the
sheet `sigma`. On sheet 0, `N <= X` is `M <= X + 1`, and (A) is applied with `X + 1`. By (A) the triples of type
`d` contribute `(d/m) D^(d)` to every class with `d | c - c_sigma`, and nothing to the others, up to errors that
add up to `sum_{y1 < Y^(1/2)} (Y/y1)^(1 - theta_m) = O(Y^(1 - theta_m/2))` (here the uniformity of (A) in `h1`, `s`
and `X` is used); the sum over `s` keeps the bound. The exact statement is `d | x1 x2 + 2^s y1 y2` when `d | g(h1)`.

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

**(E).** By (D) with `v_p(N - c_sigma) = v_p(n) < k` for every `p^k || m`. The limit: let `P` be at least the
largest prime factor of `n`. If every odd prime divides all but finitely many terms `m_j`, then for large `j` the
class factor of `m_j` is the product of the factors of all odd primes `p <= P` and of factors `1 - p^-2` at primes
`p > P`; the latter product lies between `prod_{p > P} (1 - p^-2)` and 1, and `prod_{p > P} (1 - p^-2) -> 1` as
`P -> infinity`. The Euler product: `prod_{p odd} (1 - p^-2) = 8/pi^2` and `prod_p sigma_1(p^v)/p^v = sigma_1(n)/n`.
Mean: the Dirichlet series of `sigma_1(n)/n` over odd `n` is `zeta_odd(s) zeta_odd(s + 1)`, so its mean over odd `n`
is `zeta_odd(2) = pi^2/8`; on sheet 0 the odd part of `N + 1` is distributed over the odd integers with the weights
`2^-a`, which does not change the mean. Lower bound: `sigma_1(n)/n >= 1` with equality only at `n = 1`, that is at
`N = 1` for sheet 1 (not a value of that sheet, where `N >= 1 + 2 = 3`) and at `N + 1` a power of 2 for sheet 0;
along the primes, `S_1(p) = (8/pi^2)(1 + 1/p) -> 8/pi^2`. Maximal order: `prod_{2 < p <= x} (1 - 1/p)^-1` is
asymptotic to `(e^gamma/2) ln x` (Mertens' theorem without the prime 2), which in the proof of Gronwall's theorem
gives `limsup sigma_1(n)/(n ln ln n) = e^gamma/2` over odd `n`; times `8/pi^2` this is `4 e^gamma/pi^2 = 0.7218`.

**(F).** Over `Z_p`, `p` odd, a unimodular quadratic form is determined up to equivalence by its rank and by its
discriminant modulo squares of units (CITED: Cassels, Rational Quadratic Forms, ch. 8; O'Meara, Introduction to
Quadratic Forms, section 92). The sum of four squares has discriminant 1. The form `a1 a2 + u b1 b2` is the
orthogonal sum of the hyperbolic planes `a1 a2` and `u b1 b2`, with Gram determinant `(-1/4)(-u^2/4) = (u/4)^2`, the
square of a unit. Equivalent forms have the same number of solutions of `Q = c` modulo every `p^k`. The count for
the first form with uniform variables is THM-5135 (1) with `r = 2`; so it is also the `p`-adic density of the
four-square form, and Jacobi's `r_4(n) = 8 sigma_1(n)` (odd `n`) is the product of the real density `pi^2 n`, of the
2-adic density 1 and of these densities at the odd primes, as Siegel's formula says for a genus of one class
(CITED). Receipts: the class counts of the two forms agree modulo 9, 27, 25, 49, 121 for the unit `u = 2`
(seam_rank_exact.out), and for every unit at the fourteen moduli 3, 5, 7, 9, 11, 13, 15, 21, 25, 27, 33, 35, 45, 49,
with sampled units at 63, 75, 81, 105, 121, 125 (audit, proof-5134-FH/f1_local_counts.out); `r_4 = 8 sigma_1` for
odd `n <= 400` (bridges_exact.out S4). The real density `pi^2 n` and the 2-adic density 1 at the odd classes are
computed by the audit (citations/c1_classical_counts.out C7, c2_local_class_counts.out L5); the two lines of
bridges_exact.out S4 about them record identities, not tests.

**(G).** *The formula.* Every odd cycle of a one-apex tournament passes through its apex, so its odd cycles meet
pairwise and `I(Omega(R(h)), z) = 1 + psi(h) z` (there are `psi(h)` of them, since `x = 1 + 2 psi` by OCF). In
`J(h1, h2, s)` each of the `s` middle vertices is beaten by the first apex and beats the second, so it can be
counted as a chain vertex of either side: `J` is the apex bridge `T (.) B` of THM-5070 (2) for the words
`W = 2^(s1) y1 - 1` (the digits of `2 h1`, most significant first, followed by `s1` ones) and `W' = 2^(s2) y2 - 1`
(read least significant first: `s2` ones, then the digits of `2 h2`), for any `s1 + s2 = s`. By the recursion of
THM-5070 (2), `a(2W + 1) = a(W)` and `a(2W) = a(W) + ceil(W/2)`: trailing ones do not change `a`, so
`a(W) = a(2 h1)`, and `h -> a(2h)` satisfies the recursion of `psi` (`a(4h) = a(2h) + h`,
`a(4h + 2) = a(2h + 1) + h + 1 = a(2h) + floor(h/2) + 1`), so `a(2h) = psi(h)`. For `s1 = 0`, `ceil(W/2) = h1` and
`1 + floor(W/2) = h1 + 1`; for `s1 >= 1` both equal `2^(s1 - 1) y1`. So the bracket
`ceil(W/2)(1 + floor(W'/2)) + (1 + floor(W/2)) ceil(W'/2)` of THM-5070 (2) is `h1 (h2 + 1) + (h1 + 1) h2 = (y1 y2 - 1)/2`
for `s = 0` and `2^(s-1) y1 y2` for `s >= 1`, which is the formula for `N_b`. The form of THM-5085 (2) is this
formula with `x_b(h2) = 1 + b (h2 - q(h2))` and `y2 = 2 h2 + 1` expanded.
*The region and the value.* Fix the smaller state `(h1, s)`. Since `x_b(h) = (1 - b/2) + (b/2) x(h)`,

    N_b = a x2 + b' y2 + k,     a = b x_b(h1)/2,     b' = b 2^(s-1) y1  (s >= 1),   b' = b y1/2  (s = 0),
    k = x_b(h1) (1 - b/2) - (b/2) [s = 0],

so `{h2 : N_b <= X}` is the region `{a x2 + b' y2 <= X - k}` of THM-5128 (B4), with `a > 0`, `b' > 0` and the
threshold `X - k > 0`; the lemma is uniform in the threshold and holds for all `a > 0`, `b' >= 0`. The value is a
linear form in the base-`b` state of the larger word: `V = x_b(h1) x_b(h2) + b 2^(s-1) y1 y2` for `s >= 1` and
`V = 2 x_b(h1) x_b(h2) + b y1 y2` for `s = 0`, with the coefficients `(x_b(h1), b 2^(s-1) y1)`, resp.
`(2 x_b(h1), b y1)`, of content `gcd(x_b(h1), y1, m) = gcd(g_b(h1), m)` modulo an odd `m` prime to `b` (`b`, 2 and
`2^(s-1)` are units). `(x_b, y) = (1 + b (h - q), 2h + 1)` is the image of `(q, h)` under an affine bijection of
`(Z/m)^2` (determinant `-2b`, a unit), as `(x, y)` is; so a set of residue vectors of `(x_b, y)` modulo `m` is a set
of as many residue vectors of `(x, y)`, and THM-5128 (B4)(ii)-(iv) apply to it as they stand.
*(A), (B) in base `b`.* The proof of (A) applies to `V` with the region above; the blocks are bounded by
`2^n <= (2X + b)/(b 2^s y1)`, because `N_b >= b 2^(s-1) y1 y2 - b/2`. The proof of (B) applies: `N_b` is symmetric
in `(h1, h2)`, the diagonal has `O_b(Y^(1/2))` triples for each `s`, and `G = #{h2 > h1 : N_b <= X}` vanishes unless
`b 2^(s-1) y1^2 < X + b/2`. The class of `N_b` is read from `V`: `V = N_b` on the sheet `s >= 1`, and on sheet 0
`V = c'` is `N_b = (c' - b)/2 (mod m)`, so the distinguished class `V = 0` is `N_b = -b/2`.
*Mass and types, (C), (D) in base `b`* (details supplied by the audit, auditor proof-5134-FH). As in the proof of
THM-5128 (D1), `G = Gamma(h1) + O((Y/y1)^(1/2) + 1)` with `Gamma(h1) = int_0^1 (Z(v) - max(h1, 1))^+ dv`,
`Z(v) = (X - k)/(2 (a (1 - v/2) + b'))`, and `int_0^1 Z = (X - k) w_s^b(h1)`, where

    w_s^b(h) = (2/(b x_b)) ln((x_b + 2^s y)/(x_b/2 + 2^s y)),      w_0^b = w^0_b,      sum_{s >= 1} w_s^b = w^1_b

(the ratios are `2 (x_b + 2^s y)/(x_b + 2^(s+1) y)`, and their product over `1 <= s <= S` tends to `1 + x_b/(2y)`).
The additive constant is harmless: `t -> (t - c)^+` is 1-Lipschitz, so replacing `X - k` by `X` changes `Gamma(h1)`
by at most `|k| w_s^b(h1) <= ln 2`. (For `s >= 1` it is below `(1 - 2/b) ln 2`; for `s = 0` it is
`(1 - 2/b + 1/x_b) ln(2 (x_b + y)/(x_b + 2y))`, which is below `2 ln(4/3)` at `h = 0` and below `ln 2` for `h >= 1`,
where `x_b >= 1 + b`.) Over the smaller states this adds `O(Y^(1/2))`. So
`G = X w_s^b(h1) - Xi_X(h1) + O((Y/y1)^(1/2) + 1)` with `0 <= Xi_X(h1) <= max(h1, 1)`, as for `b = 2`.
Regularity of the weights: `w_s^b(h) = F(r)/y` with `r = x_b/y = (b/2) rho + (1 - b/2)/y` and
`F(r) = (2/(b r)) ln((r + 2^s)/(r/2 + 2^s)) <= 2^(-s)/b`. For `h >= 1`, `r` lies in `((b + 4)/12, b/2]` (because
`1/2 < rho <= 1` and `y >= 3`), where `F` is smooth with derivative `O_b(2^-s)`. On a cell of `B_n` with `j` top and
`j` bottom digits prescribed, `y` varies by the factor `1 + O(2^-j)` and `rho` by `O(2^-j)` (proof of THM-5128
(B4)), so `w_s^b` varies by `O_b(W_0 2^-j)` with `W_0 = 2^(-s-n)`: this is the hypothesis of THM-5128 (B4)(iii).
The same holds for `Xi_X = int_0^1 min(Z, max(h1, 1))` with `W_0 = 2^n`: `Z` and `max(h1, 1)` vary by factors
`1 + O_b(2^-j)`, hence so does their minimum, which is less than `2^n`.
With these two facts the proofs of THM-5128 (D1), (D2) apply with `w_s^b` in place of `w_s` and
`chi_d^b(h) = [gcd(g_b(h), m) = d] - c_d` in place of `chi`; `chi_d^b` is a `chi_S` of THM-5128 (B4)(iv), with
`#S = c_d m^2`, by the affine bijection. The block sums are
`sum over B_n of w_s^b = ln 2 * int_{1/2}^{1} F((b/2) rho) d rho + O_b(2^(-s) 2^(-n/2))` (the term `(1 - b/2)/y` of
`r` changes a block sum by `O_b(2^(-s) 2^(-n))`), and

    int_{1/2}^{1} F((b/2) rho) d rho = (4/b^2) [Li_2(-b rho/2^(s+2)) - Li_2(-b rho/2^(s+1))]_{1/2}^{1}
                                     = (4/b^2) (l_(s+1) - 2 l_(s+2) + l_(s+3)),        l_j = -Li_2(-b/2^j),

which is `kappa_0(b)` at `s = 0` and telescopes to `kappa_1(b)` over `s >= 1`; the smaller states run over
`(1/2) log_2 Y + O_b(1)` blocks. This gives the mass law, (C) and (D) in base `b`. The summations over the blocks
and over the smaller states are those of THM-5128 (D1), (D2); they are referred to and not repeated line by line
here. (The audit confirmed `kappa_sigma(b)` and the first-order constants against complete censuses: Numerical
support.)
*(E) in base `b`.* `N_b = 1 (mod b)` and `2 N_b + b = 2 (mod b)`, so the odd part `n` of the value is prime to `2b`
and `v_p(V) = v_p(n)` for every prime `p` not dividing `2b`. By (D) in base `b` the limiting class factor at an
admissible `m` is the product of the factors `(1 - p^-2) sigma_1(p^(v_p(n)))/p^(v_p(n))` over `p | m`, and the limit
along `m` is taken as in the proof of (E), over the primes not dividing `2b`; `prod_{p not dividing 2b} (1 - p^-2)`
is `8/pi^2` divided by the factors `1 - p^-2` of the odd primes of `b`.
*Second digit, odd `b`* (proof found by the audit, auditor "hypotheses"). `N_b = 1 + alpha_1 b + alpha_2 b^2`, so
the second digit is `alpha_1 (mod b)`. Fix the smaller state `(h1, s)`. Since `psi2 = (x2 - 1)/2`, with `2^(-1)` the
inverse of 2 modulo `b`:

    alpha_1 = 2^(-1) x2 + (2^(s-1) y1) y2 + (psi1 - 2^(-1))   (mod b)      for s >= 1,
    alpha_1 = 2^(-1) x2 + (2^(-1) y1) y2 + (psi1 - 1)         (mod b)      for s = 0.

This is an affine form in the residue vector `(x2, y2)` modulo `b` of the larger state, and its coefficient at `x2`
is a unit. So for every class `c` it takes the value `c` on exactly `b` of the `b^2` residue vectors: one for each
residue of `y2`. The `h2` with `N_b <= X` form the region `{a x2 + b' y2 <= X - k}` above. By THM-5128 (B4)(ii) in
the form (iv), for the odd modulus `b`, each residue vector modulo `b` is taken by `(x2, y2)` on
`b^-2 #S_n + O_b(2^(n (1 - theta_b)))` states of `S_n`, the set of the `h2 > h1` of the block `B_n` in the region.
Summing over the `b` vectors and over the blocks `2^n <= (2X + b)/(b 2^s y1)`:

    #{h2 > h1 : N_b <= X, alpha_1 = c (mod b)} = G/b + O_b((Y/y1)^(1 - theta_b)),

uniformly in `h1`, `s`, `X`, `c`. `alpha_1` and `N_b` are symmetric in `(h1, h2)`, so the summation over the smaller
states is that of THM-5128 (C2): the diagonal has `O_b(Y^(1/2))` triples for each `s`, `G = 0` unless
`b 2^(s-1) y1^2 < X + b/2`, the errors add up to `O_b(Y^(1 - theta_b/2))`, and the sum over `s` is geometric. This
is the statement, on each sheet. For `b = 2` the second digit of `N` is `N mod 4`, and the flat law is THM-5128 (E)
with `j = 2`. For even `b > 2` this proof does not apply as it stands (2 is not a unit modulo `b`, and by value the
states are not equidistributed modulo powers of 2: THM-5128 (E)); the audit indicates how the by-value argument of
THM-5128 (E) would cover `b = 2^a` and mixed `b`, which is not written out (HYP-9665 (C)).
*Balanced states.* `q(h) < h/2` for `h >= 1` (induction on the recursion: `q(1) = 0`, `q(2h) = q(h)`,
`q(2h + 1) = q(h) + ceil(h/2)`), so `psi(h) = h - q(h) > h/2` and `1 + b psi(h) > 2h + 1` for `b >= 4`, `h >= 1`.

## Remarks

1. *What changed with respect to THM-5128.* Its (D3) product `prod_p (1 - 1/p^2 + [p | n]/p)` is the squarefree
   level: the mean of `L_(p^k)` over the classes divisible by `p`. The series at `N` is `(8/pi^2) sigma_1(n)/n`:
   smaller by `p^-3` at each prime dividing `n` exactly once, larger from `p^2 | n` on.
2. *Finite size.* At every computed size the factors are governed by the `1/ln X` term of (D), i.e. by the few
   smallest states of each type: at 25 the classes with `v = 1` and `v >= 2` have nearly the same factor (`1.233`
   and `1.239` on sheet 0, `1.243` and `1.248` on sheet 1 at `X = 2^28`, cum_28.out), because the state `h = 2`
   (`x = y = 5`) is of type 5 and no state of type 25 is small (the first is `h = 512`, `x = y = 1025`; audit,
   proof-5134-AE/fixed_state_A.out).
3. *Where the divisors come from.* `sigma_1(n)/n = sum_{d | n} 1/d`: the divisor `d` is the content of the smaller
   state. The states with `d | g` have density `d^-2` (THM-5128 (B2)), and they put their triples on the share `1/d`
   of the classes, those with `d | N - c_sigma`, where their density is therefore multiplied by `d`:
   `d^-2 * d = 1/d`. (At a fixed `m` the exact type `d` has density `c_d = d^-2 prod_{p | m/d} (1 - p^-2)`; its formal
   limit over all moduli, `(8/pi^2) d^-2`, is not proved: by THM-5128 (F) the density of `g = 1` is not proved to be
   `8/pi^2`.) `sigma_1(n)` is also the number of sublattices of index `n` in `Z^2` (Hermite normal forms of
   determinant `n`; classical).
4. *Even moduli.* For `2^j m` with `j >= 2` and `m > 1` odd the natural statement is the product of the flat law on
   the odd classes modulo `2^j` and of `L_m` (for `j = 1` the odd classes modulo `2m` are the classes modulo `m`
   and the statement is (D), `N` being odd). It needs the joint equidistribution of `q` modulo `2^(j-1)` and of the
   states modulo `m` on cells, that is `rho(T_k) < 2` for the non-trivial characters of `Z/(2^(j-1) m)`, which
   follows from THM-5128 (B1) for the characters of odd order and from the unitarity computation of THM-5128 (E),
   carried out for an even modulus, for those of even order (audit, hypotheses/proof_sketches_and_tables.txt
   section 1 and mixed_gap.out; for the count alone, equidistribution on cells at every modulus is KNOWN:
   Session193 modular_sync, Theorems 1, 3, 4). The by-value statement is not written out (HYP-9665 (A)). NUMERICAL
   (audit, exact censuses, in HYP-9665 (A)): inside each class modulo `m` the odd classes modulo `2^j` are flat at
   the level of the sampling noise, while the remainder of the first-order product law at `2^j m` lies between
   twice the remainder of (B) at `m` and that value plus twice the flatness defect times the largest class factor
   modulo `m` (HYP-9665 (A), with the proof of the two bounds), so it oscillates with the remainder at `m`; it is
   therefore not a test of the even law.
5. *No uniformity.* All constants depend on `m` (and on `b` in (G)); nothing here is a level of distribution.
6. *Session265's "prime-power lifts".* They are reachability statements for the moving-mark source state modulo
   `p^k` (exceptions 5 and 761 of that proof); they carry no measure and are unrelated to the densities above. F1
   has no exceptional odd prime: 25 and 125 follow the same law as 9, 27, 81.
