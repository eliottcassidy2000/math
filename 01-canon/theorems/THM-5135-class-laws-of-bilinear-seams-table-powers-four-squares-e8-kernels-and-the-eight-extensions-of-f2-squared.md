---
id: THM-5135
title: "Class laws of bilinear seams with uniform states. (1) A unimodular seam of rank r takes a class of valuation v < k modulo p^k with factor (1 - p^-r) sigma_(r-1)(p^v)/p^((r-1)v), and the class 0 with factor (1 - p^-r) sigma_(r-1)(p^(k-1))/p^((r-1)(k-1)) + p^(-(r-1)k): it is the r-fold additive power of the multiplication table modulo p^k; for every modulus m the additive square of the table is m^3 L_m. (2) At an odd prime the sum of 2r squares has the same law iff r is even or p = 1 (mod 4). (3) At 2 the four-square form differs from the rank-2 seam (already modulo 2), while the E8 form has the law of the rank-4 seam modulo every 2^k: over all primes sigma_3(n)/(zeta(4) n^3) times (8 pi^4/3) n^3 is 240 sigma_3(n). (4) A kernel with Gram matrix G has the law (1) at every p not dividing det G, whatever its signature, and at every prime its law is the additive convolution of the laws of the scaled tables p^(e_i) a_i b_i, p^(e_i) the elementary divisors of G (Smith form); two rank-4 forms of determinant -1 occur in the join calculus; the one with a proved law (the K-form on the root-change words of THM-4916, letter-count model) does not follow the uniform law, its states being confined; the corner of the fusion law is unimodular at the odd primes only, in the path-mass coordinates of its states, and its law there is open (HYP-9666 (A)). (5) The eight central extensions of F_2^2 by F_2 are the sets of units of square -1: 0 = R^4, 1 = split quaternions (the dihedral group of order 8, D_4; transitive triangle), 2 = bicomplex numbers (C_4 x C_2, signed triangle), 3 = Hamilton (Q_8, cyclic triangle); commutative iff even; the norm counts on the Lipschitz order Z[i, j, k] convolve r_2 additively, those on Z[i] x Z[i], the maximal order of the bicomplex numbers, multiplicatively. (6) Brahmagupta identity for self-joins."
status: >
  AUDITED 2026-11-11: independent audit by eight auditors, PASS WITH CORRECTIONS (6 MAJOR, 64 MINOR, 41 NOTE; all
  applied; MISTAKE-765); report: 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/AUDIT_REPORT.md.
  PROVED (elementary; proofs in this file): (1), for every prime p and every k and r, the class 0 included, with
  the table identity for every modulus m; (2), for every odd prime p and every k and r, from the CITED
  classification; (3)(b) at the prime 2, modulo every 2^k, by the lifting lemma (supplied by the audit) from the
  enumerated counts modulo 2; in (3)(d) the step from the zero counts to an even unimodular lattice (supplied by
  the audit; its input at 2 is the audit's enumeration of the 1024 quaternary forms over F_2); (4)(a), the law (1)
  at every p not dividing det G, modulo every p^k; (4)(b), the Smith-form law at every prime, 2 included, modulo
  every p^k (proved by the audit; it was the conjecture HYP-9666 (C)); in (4)(c) the two determinants -1 and the
  parities of the fusion states (elementary; the factor of the K-form is KNOWN, THM-4916, letter-count model;
  nothing is claimed for the fusion corner at the odd primes, HYP-9666 (A), OPEN); in (5) the rule of squares and
  the counts on the orders Z[i] x Z[i] and Z[1, i, j, ij] (the count on Z[i, j, k] is Jacobi's theorem, CITED);
  (6).
  FINITE-EXACT, with the ranges (S1-S6 are the sections of bridges_exact.out): (3)(a) four squares against
  a1 a2 + b1 b2 modulo 4 and 8 (S3; audit: every class modulo 2^k, k <= 6); (3)(b) the E8 form against the rank-4
  seam modulo 2, 4, 8, all residue vectors (S3; audit: modulo 16, all 16^8 vectors, and modulo 32, 64, 128 by an
  exact dynamic programme on the Dynkin tree); (3)(c) eight squares against the rank-4 seam modulo 4 and 8 (S3;
  audit: every class modulo 2^k, k <= 6, with 128 zeros against 136 modulo 2); the shells of E8 of norm 2n for
  n <= 3 (S4; audit: n <= 180);
  r_4(n) = 8 sigma_1(n) for odd n <= 400, the Hurwitz count 24 sigma_1(n) for odd n < 60, r_8(n) = 16 sigma_3(n)
  for odd n <= 39 and r_8(2) = 112 (S4; audit: the general formulas for n <= 3000 for four squares, n <= 700 for
  the Hurwitz order, n <= 1500 for eight squares); the census of (5) over the 2^9 normalized cochains: 16
  cocycles, 8 classes, groups, sign tables, explicit zero divisors in seven classes (S6; that Hamilton's algebra
  has none is classical, not a search; the four algebras of the table were identified by the audit, by trace
  forms, centres and explicit models, not by the session's script). Checks of the proved statements:
  (1) at p^k = 9, 27, 25, 49 for r <= 4 and at 8, 16 for r <= 3 (seam_rank_exact.out; audit: 17 prime powers up
  to 125, r <= 6, class 0 included); the table identity for eleven odd m <= 105 (S1; audit: every m <= 210);
  (2) for r <= 4 modulo 3, 5, 7, 11, 13, 9, 25 (S2; audit: r <= 6 for the odd p <= 41, k <= 3, p^k <= 1400, and
  r <= 4 modulo 81); (4) for three kernels modulo 9, 25, 49 (S5), for [[1, 2], [3, 5]] modulo 27 and Session261's
  kernel modulo 5, 7, 47 (seam_rank_exact.out), and by the audit for all 2 x 2 matrices modulo 9 and modulo 8;
  (6) on one-apex states (S9 of bridges_exact.out; audit: 375,200 cases with s <= 6, no failure).
  CITED (classical, not proved here): Jacobi (1829, 1834; general forms r_4(n) = 8 sum_{d | n, 4 not | d} d and
  r_8(n) = 16 sum_{d | n} (-1)^(n+d) d^3, hence 8 sigma_1(n) and 16 sigma_3(n) for odd n); Hurwitz, Vorlesungen
  ueber die Zahlentheorie der Quaternionen (1919) (24 times the sum of the odd divisors of n, hence 24 sigma_1(n)
  for odd n); the composition algebras in the dimensions 1, 2, 4, 8 (Hurwitz, Nachr. Goettingen 1898, for sums of
  squares; Albert 1942 and Jacobson 1958 for arbitrary non-degenerate forms); A. E. Ingham, Some asymptotic
  formulae in the theory of numbers, J. London Math. Soc. 2 (1927), 202-208 (the sum of two products), and
  T. Estermann, On the representations of a number as the sum of two products, Proc. London Math. Soc. (2) 31
  (1930), 123-133 (its secondary terms); the theta series of E8 (OEIS A004009: 1 + 240 sum sigma_3(n) q^(2n),
  vectors of norm 2n) and H. S. M. Coxeter, Integral Cayley numbers, Duke Math. J. 13 (1946), 561-578 (his order
  is one of seven isomorphic maximal orders); Siegel's formula for a genus of one class (Z^4, E8); the
  classification of unimodular quadratic forms over Z_p, p odd (rank and discriminant modulo squares; Cassels,
  Rational Quadratic Forms, ch. 8; O'Meara, Introduction to Quadratic Forms, sec. 92), and for p = 2 the even
  unimodular case (rank and Arf invariant; O'Meara, sec. 93; Conway-Sloane, Sphere Packings, Lattices and Groups,
  ch. 15 sec. 7), which (3)(b) no longer needs; the rank of a positive definite even unimodular lattice is
  divisible by 8 (Serre, A Course in Arithmetic, ch. V), so none has rank 4; H^2(F_2^2; F_2) as the space of
  quadratic forms on F_2^2; the multiplicativity of the sum of four squares. The books and the papers (Jacobi;
  1898, 1919, 1927, 1930, 1942, 1946, 1958) are cited at the level of the audit's sources.txt (bibliographic
  records and secondary sources; no original was re-read).
  KNOWN in the repository, credited, not claimed: the factor of (1) for r = 2 (audit of THM-5128) and, for every r
  modulo a prime, for any non-degenerate bilinear form, hence for every kernel with p not dividing det G and every
  signature (THM-5128 (A); quoted in HYP-9660 (D)); the table multiplicities N_m(t) (modular centres packet,
  PROOFS.md P4, eq. (4)) and the box version of the additive square of the product table
  (overnight_research_20261029/TWO_PORTS.md P5: D_L(t) = sum_k f_L(k) f_L(k - t) counts the 2 x 2 matrices of
  determinant t with entries in a box); the counts 16, 8, 4, 2 of (5), the class of Q_8 and the three classes of
  the dihedral group of order 8 (D_4; written D8 in THM-2887), and the rule of squares, which is the polarization
  identity of THM-2887 eq. (4), stated there for the determinant commutator (the commutative classes are implied
  there only by 16 - 8); the line rule "cyclic = Hamilton, transitive = split quaternions, the unit of square -1
  is the middle vertex" (THM-5133 (D1), kind-pasteur S23, audited 2026-11-11, the same day, PASS WITH CORRECTIONS,
  MISTAKE-762; its (A), (B), (D1) are unchanged by that audit; found here independently by enumeration and
  credited to it); the norm-composing sign tables (Session112); Jacobi's and Hurwitz's counts in THM-4826 and
  Session116; the K-form (Session199 physical-adjoint, recorded in THM-4916 (1d)) and the local factors proved in
  THM-4916 (letter-count model); the kernel x x' + 2 y y' (Session260, QUASIPERFECT_NORMS.md eq. (1)) and the
  kernel of determinant 329 (Session261); the norm reading of self-substitutions (Session260 eq. (2); Session261
  BILINEAR_KERNEL.md eq. (4); THM-5133 (A)); Session264's fusion law (3) with its path masses; MISTAKE-484 (Z^8 is
  not E8), MISTAKE-506 (Hamilton's algebra ramifies only at 2 and at the real place).
  NEW in the repository: (1) at prime powers for every r with its table form; (2); (3) (the E8 reading of the
  rank-4 law at every prime); (4) at prime powers, its Smith-form law (audit) and its checks on the kernels of
  THM-5133 (B), Session260 and Session261; the commutative half of (5) with the bicomplex numbers (which no file
  names) and the additive / multiplicative convolution statement; (6) for F1 states (the identity is the classical
  two-square identity; no earlier statement for F1 states was found in the repository). None of this is new
  mathematics outside the repository.
  No statement here is about Chan-Pak, which remains OPEN.
session: opus-2026-11-11-S92
date: 2026-11-11
source: 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/README.md
scripts:
  - 04-computation/experiments/s92_prime_powers_20261111/bridges_exact.py (+ bridges_exact.out; sections S1-S6, S9)
  - 04-computation/experiments/s92_prime_powers_20261111/seam_rank_exact.py (+ seam_rank_exact.out; its two headings are THM-5135 (1), which is statement (1) here, and THM-5134 (F); the kernel checks of (4) follow them)
  - 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/audit/ (AUDIT_REPORT.md and the auditors' receipts in proof-5135/, citations/, hypotheses/, bridges-note/, numerics-prereg/; the re-derivations made when the corrections were verified are in verify-THM-5135/)
related:
  - 01-canon/theorems/THM-5134-* (the by-value law of F1; this theorem is its algebra)
  - 01-canon/theorems/THM-5128-f1-by-value-local-densities-first-order-class-law.md (parent; its (A) is the rank-r law modulo a prime for any non-degenerate bilinear form)
  - 01-canon/theorems/THM-5133-flip-and-spin-kernel-signature-sign-tournaments-of-twisted-algebras-and-the-pg32-extremal-theorems.md (audited 2026-11-11, PASS WITH CORRECTIONS, MISTAKE-762)
  - 01-canon/theorems/THM-2887-quaternionic-arf-lift-of-the-semantic-v4-and-global-carry-no-go.md
  - 01-canon/theorems/THM-4916-root-change-local-factor-identities-confinement-multipliers-product-family.md, 05-knowledge/results/root_change_local_factors_20261107/PROOF.md, session199/physical-adjoint/PROOF.md (the K-form)
  - 01-canon/theorems/THM-4826-which-grid-carries-a-norm-shell-one-to-two-split-cubic-character-of-two-centred-units-complete-the-quartic-designs-four-eisenstein-fourth-powers-miss-five-mod-eight.md
  - 01-canon/theorems/THM-4756-forward-step-polynomial-mod-3-is-h-times-eulerian.md (module states at 3, Remark 3)
  - 05-knowledge/results/session264_constructive_20261111/fusion/PROOF.md, 05-knowledge/results/session261_creative_chanpak_20261111/constructors/BILINEAR_KERNEL.md, 05-knowledge/results/session260_quasiperfect_chanpak_20261111/arithmetic/QUASIPERFECT_NORMS.md
  - 05-knowledge/results/modular_centers_collatz_pell_20261111/PROOFS.md, 05-knowledge/results/overnight_research_20261029/TWO_PORTS.md, 05-knowledge/results/session112_counting_interfaces_20261030/norms/NORM_SIGNS.md
  - 05-knowledge/hypotheses/HYP-9666-* (rank-4 seams, (A) OPEN; its (C) is proved as (4)(b) here), 05-knowledge/hypotheses/HYP-9660-* ((D) quotes THM-5128 (A)), 05-knowledge/hypotheses/HYP-9506-* (the by-value limit of THM-4916, OPEN)
---

# THM-5135 - class laws of bilinear seams, and the algebras behind them

A *seam of rank `r`* is a bilinear form `B(u, w)` in two state vectors of length `r`; its value distribution for
independent uniform states is its *class law*. THM-5134 proves that the one-arc join family F1 has, by value, the
class law of the rank-2 seam at every odd modulus. This theorem records the algebra: which laws belong to which
forms, where the classical definite forms agree with them, and which algebras of dimension 4 are involved.

Notation. The state vectors are `u, w`; the letter `v` always denotes a `p`-adic valuation. `sigma_j(n)` is the sum
of the `j`-th powers of the divisors of `n`: `sigma_1(n)` is the sum of the divisors (as in THM-5128; `sigma` alone,
which is the sheet index in THM-5128 and THM-5134, is not used here) and `d(n) = sigma_0(n)` is their number.
`r_j(n)` is the number of representations of `n` as a sum of `j` squares of integers.
`L_m(e) = sum_{d | gcd(m, e)} d^-1 prod_{p | m/d} (1 - p^-2)` is the function of THM-5134, taken here for every
modulus `m >= 1` by the same formula. The *factor* of a class modulo `m` is `m` times its probability.

## Statements

**(1) The rank-`r` law.** Let `p` be any prime, `u, w` independent and uniform in `(Z/p^k)^r`, and `B` a bilinear
form with unit determinant. For a class `c` modulo `p^k` with `v = v_p(c) < k`:

    p^k P(B(u, w) = c) = (1 - p^-r) sigma_(r-1)(p^v)/p^((r-1)v) = (1 - p^-r)(1 + p^-(r-1) + ... + p^-(r-1)v),

and for the class 0 (`v >= k`):

    p^k P(B(u, w) = 0) = (1 - p^-r) sigma_(r-1)(p^(k-1))/p^((r-1)(k-1)) + p^(-(r-1)k).

For `r = 1` this is the multiplication table: the number of `(a, b)` modulo `m` with `a b = t` is
`N_m(t) = sum_{d | gcd(t, m)} d phi(m/d)` (KNOWN), factor `(1 - 1/p)(v + 1)` at a class of valuation `v < k` modulo
`p^k`. The law of rank `r` is the `r`-fold additive convolution power of `N_m`. For every modulus `m`,

    (N_m * N_m)(c) = m^3 L_m(c)       for every class c modulo m,

also with one table multiplied by a unit. (THM-5134 uses odd `m` only, because the states of F1 are odd at 2; at
`2^k` the identity gives the factor `3/4` at the odd classes, the factor that THM-5134 (F) removes from Ingham's
series.) For `r >= 2` the product over all primes of the factors of an integer `n >= 1` (with `v = v_p(n) < k` at
each prime) is `sigma_(r-1)(n)/(zeta(r) n^(r-1))`; for `r = 1` it diverges to 0 (the count of `a b = n` is the
divisor function).

**(2) Sums of squares at odd primes.** For an odd prime `p` the sum of `2r` squares has the class counts of the
hyperbolic form `sum_{i <= r} a_i b_i` modulo every `p^k` if and only if `(-1)^r` is a square modulo `p`, that is
`r` even or `p = 1 (mod 4)`. So four and eight squares follow the laws of rank 2 and 4 at every odd prime; two and
six squares do so only at `p = 1 (mod 4)` (the character `chi_4`).

**(3) The prime 2 and the definite models.**
(a) Modulo 2, 4 and 8 the four-square form and `a1 a2 + b1 b2` have different class counts (already modulo 2:
`8, 8` against `10, 6`; FINITE-EXACT: modulo 4 and 8 in bridges_exact.out S3, modulo 2 and in every class modulo
`2^k`, `k <= 6`, by the audit). Jacobi (CITED; checked for odd `n <= 400`):
`r_4(n) = 8 sigma_1(n) = pi^2 n * (8/pi^2) sigma_1(n)/n` for odd `n`: real density, times the rank-2 factors at the
odd primes, times a 2-adic density 1. That the count is this product is Siegel's formula for a genus of one class
(CITED).
(b) The form `x^T C x / 2` of the lattice E8 (`C` the Cartan matrix) has the class counts of the rank-4 seam at
every odd prime (by the classification used in (2): unimodular of rank 8 and discriminant 1) and modulo every
`2^k`. The latter is PROVED by the lifting lemma of the Proofs section (audit): for an integer-valued quadratic
form whose polar form is unimodular at `p`, the class counts modulo `p^k` are determined by those modulo `p`; and
modulo 2 the E8 form and the rank-4 seam both have 136 zeros among the 256 vectors (Arf invariant 0; enumerated).
Equivalently (CITED): an even unimodular `Z_2`-lattice is a sum of planes `x y` and `x^2 + x y + y^2`, two of the
second kind being equivalent to two of the first, so E8 tensor `Z_2` is the sum of four hyperbolic planes.
FINITE-EXACT: modulo 2, 4, 8 (bridges_exact.out S3); modulo 16, 32, 64, 128 by the audit. Hence, by Siegel's
formula for the genus of E8, which has one class (CITED), over ALL primes the rank-4 series times the real
density is

    (8 pi^4/3) n^3 * sigma_3(n)/(zeta(4) n^3) = 240 sigma_3(n),

the number of vectors of norm `2n` in E8, i.e. of integral octonions of norm `n` in Coxeter's order, one of seven
isomorphic maximal orders (CITED; `240, 2160, 6720` recounted for `n = 1, 2, 3`).
(c) `Z^8` is not that model: `r_8(n) = 16 sigma_3(n)` for odd `n` only, `r_8(2) = 112` (MISTAKE-484; checked for odd
`n <= 39` and at `n = 2`). "Only" follows from Jacobi's general formula `r_8(n) = 16 sum_{d | n} (-1)^(n+d) d^3`
(CITED), which for even `n` equals `16 sigma_3(n) - 32 sum_{d | n, d odd} d^3 < 16 sigma_3(n)`. The eight-square
form differs from the rank-4 seam modulo 4 and 8 (S3) and already modulo 2 (128 zeros against 136; audit).
(d) For rank 2 the series over all primes is `sigma_1(n)/(zeta(2) n) = prod_p (1 - p^-2) sigma_1(p^v)/p^v`, that is
(1) with `r = 2` at every prime. It is the Euler-product reading of the main term of Ingham's theorem (1927;
CITED) on the number of solutions of `a1 a2 + b1 b2 = n` in positive integers:

    sum_{m < n} d(m) d(n - m) = (6/pi^2) sigma_1(n) ln^2 n + O(sigma_1(n) ln n ln ln n).

Ingham's proof is elementary and uses no singular series; the secondary terms are Estermann's (1930; CITED); the
convergence is slow (the ratio of the count to the main term is `0.70` at `n = 720720` and `1.17` at `n = 999983`;
audit). An integer-valued quaternary quadratic form with the seam's number of zeros modulo every prime is an even
unimodular lattice (Proofs; step supplied by the audit), and no positive definite even unimodular lattice of rank
4 exists (CITED: the rank is divisible by 8). So a definite quaternary model must differ from the seam at some
prime. Hamilton's form (the norm on the Lipschitz order `Z[i, j, k]`) and the norm form of the Hurwitz order do so
at 2 only: they have `8 sigma_1(n)` and `24 sigma_1(n)` elements of norm `n` at odd `n` (Jacobi, Hurwitz; CITED;
checked for odd `n <= 400` and odd `n < 60`), with the 2-adic densities 1 and `3/2` at the odd classes against
`3/4` for the seam (audit). For a positive definite even unimodular form in `2r` variables the real density is
`(2 pi)^r n^(r-1)/(r-1)!`, and its product with `sigma_(r-1)(n)/(zeta(r) n^(r-1))` is `(2r/|B_r|) sigma_(r-1)(n)`
for even `r` (`B_r` the Bernoulli number): `240 sigma_3(n)` at `r = 4`. At `r = 2`, where no such form exists, the
same expression is `24 sigma_1(n)`, which is Hurwitz's count at odd `n`: the Hurwitz order has the real density
`2 pi^2 n`, half of `4 pi^2 n`, and the 2-adic density `3/2`, twice `3/4` (audit).

**(4) Kernels.** Let `G` be an integer `r x r` matrix and `B(u, w) = u^T G w`, with `u, w` independent and uniform
modulo `p^k`.
(a) For every prime `p` not dividing `det G` the class law of `B` modulo every `p^k` is (1). It depends neither on
the signature of `G` nor on the sign of `det G`. Modulo `p` this is THM-5128 (A), which holds for any
non-degenerate bilinear form (KNOWN); the prime powers are a corollary of (1).
(b) Smith form (PROVED; the proof was supplied by the audit; it was the conjecture HYP-9666 (C)). For every prime
`p`, 2 included, let `p^(e_1), ..., p^(e_r)` be the elementary divisors of `G` over `Z_p` (with `e_i = infinity` for
the divisors equal to 0, which occur when `G` is singular). The class law of `B` modulo `p^k` is the law of
`sum_i p^(e_i) a_i b_i` with `a_1, b_1, ..., a_r, b_r` independent and uniform modulo `p^k`: it is the additive
convolution of the laws of the scaled tables `p^(e_i) a_i b_i`, and it depends on `G` only through the
`min(e_i, k)`. (a) is the case `e_1 = ... = e_r = 0`.
Checks. The elliptic and hyperbolic two-port kernels of THM-5133 (B) have the same law at the primes of (a): the
Gaussian kernel of order 4 (state `(1, 1)` in the notation `(r, M)` of THM-5133 (B), where `r` is not a rank;
visit matrix `[[5, 2], [2, 1]]`, determinant 1), the first hyperbolic kernel `[[11, 6], [6, 3]]` (order 5,
determinant `-3`) and Session260's kernel `diag(1, 2)` (the form `x x' + 2 y y'`) have the rank-2 law modulo 9,
25, 49, except the second modulo 9: there `3 | det G`, the elementary divisors are `1, 3` and the law is that of
`a_1 b_1 + 3 a_2 b_2` (numbers of solutions for a class of valuation `0, 1, 2` modulo 9: `486, 1134, 1377` against
`648, 864, 945` for the rank-2 law; audit). A non-symmetric kernel of determinant `-1`, `[[1, 2], [3, 5]]`, has the
rank-2 law modulo 27. Session261's kernel `[[149, 72], [72, 37]]`, of determinant `329 = 7 * 47`, has the
elementary divisors `1, p` at `p = 7` and `p = 47`: modulo `p` it has the rank-1 law (`p^2` times the
multiplication table), modulo 49 the law of `a_1 b_1 + 7 a_2 b_2` (factors `6/7, 90/49, 97/49` at the valuations
`0, 1, 2`; audit), and it has the rank-2 law modulo 3, 5, 11, 13.
FINITE-EXACT (audit): among all 6561 matrices modulo 9, all 4096 modulo 8, the matrices with entries in `[-4, 4]`
modulo 27 and in `[-3, 3]` modulo 25, each Smith class has one law, the convolution of its scaled tables.
(c) The law is a statement about uniform states. Two forms of rank 4 and determinant `-1` occur in the join
calculus: the K-form (Session199; recorded in THM-4916 (1d)), `U U' + E E' - d d' + a d' + a' d` (determinant `-1`
in `(U, E, d, a)`), and the corner of Session264's fusion law,
`H(A * B) = h h' + (H - h) h' + h (H' - h') + (U - h)(V' - h') + (d - h)(e' - h')` (determinant `-1` in the
coordinates `(h, H - h, U - h, d - h)` of the left state of `A` and `(h', H' - h', V' - h', e' - h')` of the right
state of `B`). For uniform states both have the law (1) with `r = 4`. What is known for the states of
tournaments:
- *K-form.* THM-4916 computes the class law of the values `K(t.alpha, t.beta) = F(t.beta alpha^dagger)` over words
  with at most two root changes (letter-count model, with the weighting `M_K` of THM-4916 at `K = 2`, which is
  not the algebra `M_2` of (5); one word, hence one free state, per value; the by-value limit is OPEN,
  HYP-9506 (A), HYP-9660 (A)): factor `1 + (3p - 2)(p - 1)/(2 p^3)` at the class 0 modulo an odd prime `p`. This
  is neither the value `1 + p^-3 - p^-4` of two independent uniform states nor the value 1 of one uniform state
  under a non-zero linear functional (at `p = 3`: `34/27` against `83/81` and 1): the states reached with few root
  changes are confined (THM-4916, cascade).
- *Fusion corner.* The states are not uniform at 2: `h = H(A - z)` (`z` the root of `A`) and `H = H(A)` are path
  counts, hence odd (Redei), so the first coordinate is odd and the second even; and `U = d (mod 2)`,
  `V' = e' (mod 2)`, so the last two coordinates have the same parity. In the path-mass coordinates
  `a_o = (U + d)/2 - h`, `a_e = (U - d)/2` (and `b_o`, `b_e` for `B`; defined in the Proofs) the last two terms of
  the corner are `2 (a_o b_o + a_e b_e)`: the form is unimodular in the coordinates `(U - h, d - h)`, whose last two
  the states do not fill (`U = d (mod 2)` is a sublattice of index 2); in the coordinates `(h, H - h, a_o, a_e)` the
  corner has determinant `-4`, so it is unimodular at the odd primes only. At 2 the uniform law fails in every family: `H(A * B)` is odd, while the uniform rank-4 law gives
  the class 0 modulo 2 the share `17/32`. At the odd primes nothing is claimed: whether the left and right states
  of some family are equidistributed is HYP-9666 (A), OPEN (the audit found no affine relation modulo 3, 5 or 7
  satisfied by the left states, or by the right states, of all labelled pointed tournaments of one order, at each
  of the orders 4 to 7).

So a unimodular Gram matrix of rank 4 does not by itself give the law of (3)(b); see HYP-9666.

**(5) The eight extensions of `F_2^2` by `F_2`.** A normalized 2-cocycle `f` on `V = F_2^2` gives the group on
`V x F_2` with `(a, z)(b, z') = (a + b, z + z' + f(a, b))` and the real algebra with basis `e_a`,
`e_a e_b = (-1)^f(a,b) e_(a+b)`. There are 16 cocycles and 8 classes modulo the one non-trivial coboundary, which
changes the sign of the six products `e_a e_b`, `a != b` non-zero (three unordered pairs, both orders). The
squaring map `a -> e_a^2` is a quadratic form on `F_2^2` whose polar form is the commutator, and `H^2(F_2^2; F_2)`
is the space of these forms (classical; CITED): the class is the support `S` of the form, i.e. the set of the
non-zero `a` with `e_a^2 = -1`, and every subset of the three points occurs:

| `#S` | classes | group | algebra | commutative | sign table on the three points | zero divisors |
|---|---|---|---|---|---|---|
| 0 | 1 | `C_2^3` | `R^4` | yes | signed triangle | yes |
| 1 | 3 | `D_4`, the dihedral group of order 8 | `M_2(R)`, split quaternions | no | transitive tournament, middle vertex `S` | yes |
| 2 | 3 | `C_4 x C_2` | `C + C`, bicomplex numbers | yes | signed triangle | yes |
| 3 | 1 | `Q_8` | Hamilton's quaternions | no | cyclic tournament | no |

The product of the three squares is `+1` in the commutative classes and `-1` in the others, so the algebra is
commutative iff `#S` is even: the polar form is non-zero iff `#S` is odd, and then the Arf invariant is 0 for
`D_4` and 1 for `Q_8` (the language of THM-2887, where `D_4` is written D8). In a non-commutative class the two
cocycles are a tournament on the three points (`a -> b` iff `e_a e_b = + e_(a+b)`) and its converse; in a
commutative class they are two signed triangles and there is no orientation. Quaternions and bicomplex numbers
share two units of square `-1` and differ in the third (`-1`, anticommuting; `+1`, commuting).
Norms. Hamilton: the sum of four squares (the reduced norm). `M_2`: the determinant, which is the seam of F1; it
is the norm on the maximal order `M_2(Z)` (on the order spanned by the units of the table the norm form is
`x0^2 + x1^2 - x2^2 - x3^2`, which has the seam law at the odd primes and not at 2). `C + C`: the multiplicative
norm `|z_1|^2 |z_2|^2` (the algebra norm, of degree 4). On the Lipschitz order `Z[i, j, k]` (the Z-span of the units
of the table; 8 units, the group `Q_8`) and on the maximal order `Z[i] x Z[i]` of `C + C` (idempotent coordinates
`(z_1, z_2)`; 16 units, `C_4 x C_4`):

    #{elements of Z[i, j, k] of norm n} = sum_{a + b = n} r_2(a) r_2(b) = 8 sigma_1(n)   (n odd),
    #{elements of Z[i] x Z[i] of norm n} = sum_{a b = n} r_2(a) r_2(b).

The Z-span `Z[1, i, j, ij]` of the units of the table of `C + C` is the suborder `z_1 = z_2 (mod 2)` of
`Z[i] x Z[i]`, of index 4, with 8 units (the group `C_4 x C_2` of the table). On it the count is half of the
second count at every odd `n`, and it is not a convolution of `r_2` at even `n`: it is 0 at every `n = 2 (mod 4)`,
and 16 against 48 at `n = 4`. Both counts on `C + C` vanish for `n = 3 (mod 4)`.
So, on these orders, the quaternion norm is the additive convolution of two Gaussian norms and the norm of
`C + C` on its maximal order their multiplicative convolution. In the join calculus the ordered join multiplies
path counts (rank 1; commutative) and the one-arc join adds two products (rank 2). For uniform factors the first
has the divisor-function law of (1); this is a comparison (the note, section 8): no equidistribution of the path
counts of ordered joins is claimed (path counts are odd, and 7 and 21 are not path counts). The second has the
law of rank 2, with the series `sigma_1(n)/n`, which F1 has at every odd modulus (THM-5134).

**(6) Brahmagupta identity.** With `M_s(A, B) = x_A x_B + 2^s y_A y_B` on one-apex states (for `s >= 1` the F1 value,
on the sheet 1; for `s = 0` the F1 value plus one, on the sheet 0):

    M_s(A, A) M_s(B, B) = M_s(A, B)^2 + 2^s (x_A y_B - x_B y_A)^2.

`M_s(A, B)` is the real part of `z_A conj(z_B)`, `z = x + sqrt(-2^s) y`, and minus the cross-determinant
`x_A y_B - x_B y_A` of the two states is the coefficient of `sqrt(-2^s)` (the imaginary part divided by `2^(s/2)`).
The identity is the classical two-square identity; the norm reading of self-substitutions is in Session260
(eq. (2)), Session261 (eq. (4)) and THM-5133 (A); the statement for F1 states was not found in the repository.

## Proofs

**(1).** A linear change of `w` by the matrix of `B` (a unit determinant) reduces to `B(u, w) = u . w`. With
probability `p^(-ri) (1 - p^-r)` the vector `u` is `p^i` times a primitive vector, `i < k`, and then `w -> u . w` is
`p^i` times a surjective homomorphism onto `Z/p^k`, uniform on the multiples of `p^i`; the vector `u = 0`, of
probability `p^(-rk)`, contributes to the class 0 only. So for `v_p(c) = v < k`:
`P(c) = sum_{i <= v} p^(-ri) (1 - p^-r) p^(i-k)`, and `P(0) = sum_{i < k} p^(-ri) (1 - p^-r) p^(i-k) + p^(-rk)`, which
are the two displayed values. For `r = 1` and a modulus `m`: `a` with `gcd(a, m) = d` (there are `phi(m/d)` of them)
takes each multiple of `d` exactly `d` times; `N_m(t)` depends only on `gcd(t, m)`, so the table multiplied by a
unit is the same function. The seam `sum a_i b_i` is a sum of `r` independent table entries. The identity with
`L_m`: by the Chinese remainder theorem `(N_m * N_m)(c)/m^4` is the product over the prime powers `p^k || m` of the
probabilities modulo `p^k`; and `d^-1 prod_{p | m/d} (1 - p^-2)` is the product over the primes of `m` of the local
terms `p^-i (1 - p^-2)^[i < k]` (`p^i || d`, `p^k || m`), so `L_m` is the product of the
`L_(p^k)(e) = sum_{i <= min(v, k)} p^-i (1 - p^-2)^[i < k]`, `v = v_p(e)` (the computation of THM-5134 (D), which
does not use that `m` is odd). For `v < k` this is `(1 - p^-2) sigma_1(p^v)/p^v`, and for `v >= k` it is
`(1 - p^-2) sigma_1(p^(k-1))/p^(k-1) + p^-k = 1 + 1/p - p^(-k-1)`: the two displayed values with `r = 2`. So the
probability modulo `p^k` is `L_(p^k)(c)/p^k`, the product of the probabilities is `L_m(c)/m`, and
`(N_m * N_m)(c) = m^4 L_m(c)/m`.
Checks: the law at `p^k = 9, 27, 25, 49` for `r <= 4` and at `8, 16` for `r <= 3` (seam_rank_exact.out); the
identity with `L_m` for eleven odd `m <= 105`, with the units 1, 2, 4 (bridges_exact.out S1); by the audit, the law
with the class 0 at the 17 prime powers `2, 4, 8, 16, 32, 3, 9, 27, 81, 5, 25, 125, 7, 49, 11, 121, 13` for
`r <= 6`, and the identity for every `m <= 210`, with the units 2 to 7. The Euler product is
`(1/zeta(r)) prod_p sigma_(r-1)(p^v)/p^((r-1)v)`.

**(2).** Over `Z_p`, `p` odd, a unimodular quadratic form is determined by its rank and its discriminant modulo
squares (CITED; the statement is false at `p = 2`: `Z^8` and E8 have the same rank and discriminant and different
class counts, (3)(c)). The sum of `2r` squares has discriminant 1, the hyperbolic form of rank `r` has `(-1/4)^r`,
in the class of `(-1)^r`. If the classes differ, the numbers of solutions of `Q = 0` modulo `p` differ: they are
`p^(2r-1) + eta (p^r - p^(r-1))`, `eta = +1` for the hyperbolic form and `-1` for the other class. (Checked for
`r <= 4` modulo 3, 5, 7, 11, 13 and modulo 9, 25: bridges_exact.out S2, in which the cases `r = 3, 4` modulo 25
were skipped before the audit and have been added since. The audit: `r <= 6` for the odd primes `p <= 41` and
`k <= 3` with `p^k <= 1400`, 162 cases, with the zero counts modulo `p`; `r <= 4` modulo 81.)

**(3).** (a), (c) and the shells are enumerations (bridges_exact.out S3, S4; the comparison modulo 2 and modulo
`2^k` for `k <= 6` is the audit's). The proofs of (b) at the prime 2 and of the step of (d) follow.

*Lifting lemma (audit).* In this lemma and its proof `n` is the number of variables, `B` is a polar form (not a
seam) and `N_k`, indexed by the exponent `k`, is not the table multiplicity `N_m` of (1). Let `Q` be an
integer-valued quadratic form in `n` variables, `B(x, y) = Q(x + y) - Q(x) - Q(y)` its polar form, and `p` a prime
not dividing the determinant of `B` (any prime; for `p = 2` this means that the form is even unimodular at 2, and
`n` is even). For `k >= 0` and a class `c` modulo `p^k` let `N_k(c)` be the number of `x` modulo `p^k` with
`Q(x) = c (mod p^k)`; so `N_0 = 1`. Then for `k >= 2`

    N_k(c) = p^((n-1)(k-1)) (N_1(c mod p) - [p | c]) + [p^2 | c] p^n N_(k-2)(c/p^2),

so that every `N_k` is determined by `N_1`, the class counts modulo `p`.

Proof. For an integer `t` and vectors `x, y`: `Q(x + t y) = Q(x) + t B(x, y) + t^2 Q(y)`. With `t = p^k` this shows
that `Q(x)` modulo `p^k` depends only on `x` modulo `p^k`, so `N_k` is well defined. Write `N_k = P_k + I_k`, where
`P_k` counts the `x` not divisible by `p` (primitive) and `I_k` the `x` divisible by `p`.
(i) Primitive vectors. Let `k >= 1` and let `x` modulo `p^k` be primitive with `Q(x) = c (mod p^k)`. Its lifts
modulo `p^(k+1)` are the `p^n` vectors `x + p^k y`, `y` modulo `p`, and, since `2k >= k + 1`,
`Q(x + p^k y) = Q(x) + p^k B(x, y) (mod p^(k+1))`. The Gram matrix of `B` is invertible modulo `p` and `x` is not
zero modulo `p`, so `y -> B(x, y)` is a non-zero linear form on `(Z/p)^n`: it takes each value of `Z/p` exactly
`p^(n-1)` times. Hence each of the `p` classes modulo `p^(k+1)` above `c` contains exactly `p^(n-1)` lifts of `x`.
Every primitive vector modulo `p^(k+1)` is a lift of a primitive vector modulo `p^k`, so
`P_(k+1)(c') = p^(n-1) P_k(c' mod p^k)`, and by induction `P_k(c) = p^((n-1)(k-1)) P_1(c mod p)`, where
`P_1(c) = N_1(c) - [c = 0]` (modulo `p` the only vector divisible by `p` is 0).
(ii) Vectors divisible by `p`. Let `k >= 2`. They are the `p x'`, `x'` modulo `p^(k-1)`, and `Q(p x') = p^2 Q(x')`.
So `I_k(c) = 0` unless `p^2 | c`; and if `c = p^2 c''`, then `I_k(c)` is the number of `x'` modulo `p^(k-1)` with
`Q(x') = c'' (mod p^(k-2))`, which is `p^n N_(k-2)(c'')`, because each vector modulo `p^(k-2)` has `p^n` lifts
modulo `p^(k-1)`.
Adding (i) and (ii) gives the formula, and an induction on `k` (in steps of 2, from `N_0` and `N_1`) gives the
last assertion. QED.

(b). At odd `p`, E8 is unimodular of rank 8 and determinant 1, hence hyperbolic by the classification used in
(2). At 2: the polar form of `x^T C x / 2` is `x^T C y`, of determinant `det C = 1`; the polar form of the rank-4
seam `a_1 b_1 + a_2 b_2 + a_3 b_3 + a_4 b_4`, a quadratic form in the 8 variables `(a, b)`, has the Gram matrix
`[[0, I], [I, 0]]`, of determinant 1. Modulo 2 both forms have the class counts `136, 120` (enumeration of the 256
vectors; bridges_exact.out S3). By the lemma with `p = 2`, `n = 8` they have the same class counts modulo every
`2^k`; those of the seam are given by (1) with `r = 4`. Receipts: the counts of the two forms agree modulo 2, 4, 8
by enumeration of all `2^8`, `4^8`, `8^8` residue vectors (S3), modulo 16 by enumeration of all `16^8` vectors and
modulo 32, 64, 128 by an exact dynamic programme on the Dynkin tree (audit, a3_e8.out); the recursion of the
lemma, started from `136, 120`, reproduces the enumerated counts modulo 2, 4, 8 and the law (1) with `r = 4`
modulo `2^k` for `k <= 14` (audit, c8_hensel_recursion_2adic.out). The number 136 is the Arf invariant 0: the
other even unimodular class of rank 8 at 2, three hyperbolic planes and one plane `x^2 + x y + y^2`, has 120 zeros
modulo 2, and `Z^8` has 128 (audit). The lemma holds at every prime: at an odd `p` it reduces the comparison of
two unimodular forms modulo `p^k` to their class counts modulo `p` (audit: four squares against the rank-2 law
modulo `3^k`, `5^k`, `7^k`, `11^k` for `k <= 8, 6, 5, 4`).

The displayed identity of (b) is `zeta(4) = pi^4/90`; the two lines of bridges_exact.out S4 about the densities
are identities, not tests. The real densities `pi^2 n` and `(8 pi^4/3) n^3` (the derivatives in `n` of the volumes
of `{x . x <= n}` in `R^4` and of `{x . x <= 2n}` in `R^8`, the covolumes being 1) and the 2-adic density 1 of the
four-square form at the odd classes (exactly 1 modulo 2, 4, 8, 16, 32) were computed by the audit
(c1_classical_counts.out C7, c2_local_class_counts.out L5).

(d). Let `Q` be an integer-valued quaternary quadratic form with `p^3 + p^2 - p` zeros modulo `p`, the number of
zeros of the seam, for every prime `p`. Modulo an odd `p`, a form of rank 0, 1, 2, 3 has `p^4`; `p^3`;
`2 p^3 - p^2` or `p^2`; `p^3` zeros, and a non-degenerate form has `p^3 + eta (p^2 - p)`, `eta` as in (2). None of
the first values is `p^3 + p^2 - p` (the equality `2 p^3 - p^2 = p^3 + p^2 - p` would mean `p (p - 1)^2 = 0`), so
`Q` is non-degenerate modulo `p`, of square discriminant. Modulo 2: of the 1024 quaternary quadratic forms over
`F_2` exactly 280 have 10 zeros, and all of them have a non-degenerate polar form (enumeration by the audit,
a2_squares.out [Q]; they are the hyperbolic ones). So the determinant of the polar form of `Q` is prime to every
prime, i.e. it is `+-1`, and the polar form has the even diagonal `B(x, x) = 2 Q(x)`: `Q` is an even unimodular
lattice. A positive definite even unimodular lattice has rank divisible by 8 (CITED), which excludes rank 4.
(Audit: over `F_3` the 21,060 quaternary forms with the seam's number of zeros are all non-degenerate; for
`p = 3, 5, 7, 11, 13` the diagonal forms with `p^3 + p^2 - p` zeros are exactly the non-degenerate ones of square
discriminant; a2_squares.out [Q], a7_lemmas.out [L].) Hamilton's form differs from the seam at 2 by (a) and agrees
with it at every odd prime by (2) with `r = 2`. The Lipschitz order has index 2 in the Hurwitz order, so the two
norm forms are equivalent over `Z_p` for every odd `p`; at 2 the factors at the odd classes are 1 for the
Lipschitz order, `3/2` for the Hurwitz order (audit, c6_followups.out F5, modulo 16) and `3/4` for the seam, by
(1). The value `2r/|B_r|` is `(2 pi)^r/((r-1)! zeta(r))`, by `zeta(r) = (2 pi)^r |B_r|/(2 r!)` for even `r` (audit,
c6_followups.out F6: 24, 240, 504, 480 for `r = 2, 4, 6, 8`). Ingham's asymptotic was compared with exact counts
by the audit (c3_ingham_gronwall_constants.out I1: the count is the number of positive solutions for `n <= 60`;
ratios to the main term up to `n = 10^6`).

**(4).** (a) If `p` does not divide `det G`, `w -> G w` is a bijection of `(Z/p^k)^r` and `u^T G w = u . (G w)`; so
the law is that of `u . w`, which is (1).

(b) (Audit.) `Z_p` is a principal ideal domain, so `G` has a Smith normal form over it: `G = P D Q` with `P, Q` in
`GL_r(Z_p)` and `D = diag(p^(e_1), ..., p^(e_r))`, `0 <= e_1 <= ... <= e_r <= infinity` (`p^infinity = 0`; the unit
parts of the diagonal entries have been moved into `Q`). Equivalently, with the Smith normal form
`G = P_0 D_0 Q_0` over `Z` (`P_0, Q_0` in `GL_r(Z)`, `D_0 = diag(d_1, ..., d_r)`), `e_i = v_p(d_i)`. The entries of
`P`, `Q` and of their inverses are `p`-adic integers, so modulo `p^k` the matrices `P^T` and `Q` are invertible:
`u -> a = P^T u` and `w -> b = Q w` are bijections of `(Z/p^k)^r`. Hence, for `u, w` independent and uniform, `a`
and `b` are again independent and uniform, i.e. the `2r` coordinates `a_1, b_1, ..., a_r, b_r` are independent and
uniform modulo `p^k`, and

    u^T G w = u^T P D Q w = (P^T u)^T D (Q w) = sum_i p^(e_i) a_i b_i      (mod p^k).

The `r` summands are independent, so the law of the sum is the additive convolution of their laws, and
`p^(e_i) a_i b_i` modulo `p^k` depends only on `min(e_i, k)` (it is 0 if `e_i >= k`). Nothing in the argument uses
that `p` is odd, or that `G` is symmetric or non-singular. If all `e_i = 0` the sum is `a . b`: this is (a) again.
Examples: `[[11, 6], [6, 3]]` has determinant `-3` and the unit entry 11, hence the elementary divisors `1, 3` at
3; `[[149, 72], [72, 37]]` has determinant `7 * 47` and is not zero modulo 7 or 47 (it is `2 [[1, 1], [1, 1]]`
modulo 7), hence the elementary divisors `1, p` at these primes, and modulo `p` its law is that of `a_1 b_1`, with
`p^2` free pairs `(a_2, b_2)`. Receipts: bridges_exact.out S5 and seam_rank_exact.out (session); a4_kernels.out
[K] and smith_r2.out (audit).

(c) The two determinants are computed from the displayed forms. K-form: in `(U, E, d, a)` the Gram matrix is
`[[1,0,0,0],[0,1,0,0],[0,0,-1,1],[0,0,1,0]]`, of determinant `-1`. Fusion corner: in the coordinates
`(h, H - h, U - h, d - h)` and `(h', H' - h', V' - h', e' - h')` the Gram matrix is
`[[1,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]]`, of determinant `-1`. Parities of the fusion states. Notation of
Session264 (PROOF.md, section 2): `A` is a pointed tournament with root `z`, `H = H(A)`, `h = H(A - z)`; the *mass*
of a path `R` of `A` is `H(A - V(R))`, the number of Hamiltonian paths of the tournament left when the vertices of
`R` are removed (1 if nothing is left); `U` is the sum of the masses of the paths that start at `z`, and `d` their
alternating sum, with the sign `+` for a path with an odd number of vertices and `-` for an even number; `V` and
`e` are the same sums over the paths that end at `z`; primes denote the quantities of `B`. The one-vertex path
has mass `h`. Let `a_o`, `a_e` be the total masses of the non-trivial root-start paths of `A` with an odd, even
number of vertices, and `b_o`, `b_e` those of the root-end paths of `B` (Session264, proof of the fusion law).
Then `U - h = a_o + a_e`, `d - h = a_o - a_e`, `V' - h' = b_o + b_e`, `e' - h' = b_o - b_e`; so `U = d` and
`V' = e' (mod 2)`, and `(U - h)(V' - h') + (d - h)(e' - h') = 2 (a_o b_o + a_e b_e)`. In the coordinates
`(h, H - h, a_o, a_e)` and `(h', H' - h', b_o, b_e)` the Gram matrix is
`[[1,1,0,0],[1,0,0,0],[0,0,2,0],[0,0,0,2]]`, of determinant `-4`. The share `17/32` is (1) with `p^k = 2`, `r = 4`:
the factor of the class 0 is `15/16 + 1/8 = 17/16`. The values at `p = 3` are
`1 + (3p - 2)(p - 1)/(2 p^3) = 34/27` and `1 + p^-3 - p^-4 = 83/81`. Receipts (audit): the two determinants, in
these and in the raw coordinates, and the law (1) of both forms for uniform states modulo 2, 3, 4, 5, 9
(a4_kernels.out [D]); the factor of THM-4916 against the uniform values at `p = 3, 5, 7, 11, 13` (a4_kernels.out
[F]); the fusion law on 5625 pairs of pointed tournaments of order at most 4, the parities for all labelled
pointed tournaments of order at most 7, and the absence of affine relations modulo 3, 5, 7 at the orders 4 to 7
(fusion_states.out).

**(5).** The census is an enumeration (bridges_exact.out S6): cocycle condition on the `2^9` normalized cochains;
groups by commutativity and the number of elements of order 4 (`0, 4, 2, 6` for `C_2^3`, `C_4 x C_2`, `D_4`, `Q_8`;
the exponent divides 4, because the square of `(a, z)` is `(0, f(a, a))`, of order at most 2, so the group is not
cyclic of order 8); explicit zero divisors with coordinates in `{-1, 0, 1}` in the seven classes other than
Hamilton's. Hamilton's algebra has none because the sum of four squares is multiplicative and definite (classical;
CITED), not because a search in a box found none (the audit checked `det L_x = (x_0^2 + x_1^2 + x_2^2 + x_3^2)^2`
for the left multiplication `L_x`). Rule of squares, with `s_a = e_a^2` and `c = a + b`:
`e_c^2 = (e_a e_b)(e_a e_b) = eps e_a^2 e_b^2` with `eps = +1` if `e_a, e_b` commute and `-1` if they anticommute
(associativity is the cocycle condition), so `s_a s_b s_c = eps`; in additive notation this is the polarization
identity `q(a + b) + q(a) + q(b) = commutator` (THM-2887 eq. (4), stated there for the determinant commutator).
The algebras: for `#S = 0`, `R[V]`; for `#S = 2`, `C` tensor `R[C_2] = C + C`; for `#S = 1`, `(1, -1)_R = M_2(R)`;
for `#S = 3`, `(-1, -1)_R`. These identifications are not checked by bridges_exact.py (it checks commutativity and
zero divisors); the audit checked them (trace forms of signatures `(4, 0)`, `(3, 1)`, `(2, 2)`, `(1, 3)` and
centres of dimensions 4, 1, 4, 1 for `#S = 0, 1, 2, 3`; explicit idempotents and matrix models;
a5_extensions.out, c4_extensions_census.out).
The counting formulas. Write an element of either algebra as `g + g' j` with `g, g'` in `C = R[i]`: the quaternion
norm is `|g|^2 + |g'|^2` and the norm of `C + C` is `|g + i g'|^2 |g - i g'|^2`, where
`(z_1, z_2) = (g + i g', g - i g')` are the idempotent coordinates (`j` acts as `i` and as `-i` on the two
factors). For Gaussian integers `g, g'` the first formula is the definition of the additive convolution, and it is
`8 sigma_1(n)` at odd `n` by Jacobi's theorem. On `Z[i] x Z[i]` the second formula is the definition of the count;
`r_2(n) = 0` for `n = 3 (mod 4)` and one of `a, b` is `3 (mod 4)` when `a b` is. The Z-span of `1, i, j, ij` is the
set of the `g + g' j` with `g, g'` Gaussian integers, i.e., as `z_1 - z_2 = 2 i g'` and `z_1 + z_2 = 2 g`, the
suborder `z_1 = z_2 (mod 2)`, of index 4. Modulo 2 a Gaussian integer of odd norm is in the class of 1 or of `i`,
and one of norm `2 (mod 4)` in the class of `1 + i`. For odd `n`, multiplication of `z_1` by `i` exchanges the
elements of norm `n` with `z_1 = z_2 (mod 2)` and the others, so the suborder contains half of them; for
`n = 2 (mod 4)` one of `z_1, z_2` has odd norm and the other is in the class of `1 + i`, so the suborder contains
none. (Audit: half at every odd `n <= 400`, zero at `n = 2 (mod 4)`, and the values at even `n`; a2_squares.out
[B], c6_followups.out F3, a3_algebras_fano_e8.out part b. The form `x0^2 + x1^2 - x2^2 - x3^2` equals the seam
modulo 3, 9, 5, 7 and differs from it modulo 2, 4, 8, 16; a7_lemmas.out [O].)

**(6).** `(x^2 + t y^2)(x'^2 + t y'^2) = (x x' + t y y')^2 + t (x y' - x' y)^2`.

## Remarks

1. *What (5) says about tournaments.* Only the non-commutative classes are tournaments; that is the dictionary of
   THM-5133 (D), restricted to one line. The bicomplex numbers sit in the commutative half, where the sign table
   is a signed graph. On one line the orientation and the set of the units of square `-1` determine each other,
   and the norm form of a non-commutative class does not see the orientation at odd primes: the determinant
   (transitive triangle, `H = 1`) and the sum of four squares (cyclic triangle, `H = 3`) have the same law there,
   by (2); they differ at the prime 2 and at the real place. This does not extend beyond one line: a norm form
   sees only which units have square `-1`. The 128 tournaments on the Fano plane with all lines cyclic (16
   alternative ones with `H = 189`, 112 others with `H = 171`) have the same eight-square form, which is
   multiplicative only for the 16 (audit, a3_algebras_fano_e8.out part d): their orientation is invisible to the
   form at every place. And with Coxeter's order the definite form E8 agrees with the rank-4 seam at 2 as well
   ((3)(b)), so that there only the real place separates the definite form from the split one. No statement about
   path counts follows.
2. *Composition algebras.* The hyperbolic form of rank `r` is the norm form of a composition algebra exactly for
   `r = 1, 2, 4`, i.e. in the dimensions 2, 4, 8 (split complex numbers, split quaternions, split octonions;
   Hurwitz 1898 for sums of squares, Albert 1942 and Jacobson 1958 for arbitrary non-degenerate forms; CITED). The
   laws of (1) hold for every `r`, and rank 3 occurs in the join calculus: kernels with three visits (THM-5133
   (A)) have `3 x 3` visit matrices, for instance `[[5, 3, 0], [3, 7, 2], [0, 2, 2]]` at order 4 (determinant 32)
   and `[[15, 13, 2], [13, 14, 4], [2, 4, 2]]` at order 5 (determinant `-6`); these were found by the audit
   (a6_identity_rank3.out [V]), THM-5133 (A) exhibits none. The ranks 1, 2, 4 are not singled out by the gluing
   calculus.
3. *Scope.* (1)-(4) are statements about uniform independent states. THM-5134 supplies the equidistribution for
   F1 at every odd modulus. For module states `(H, Pi_2)` it fails at 3 (THM-4756: `Pi_2 = (-1)^n H (mod 3)` for
   every tournament on `n` vertices), and so it does for every kernel of THM-5133 (A), whatever its number of
   visits: modulo 3, `H(K[A, B]) = c H(A) H(B)` with a constant `c` that depends only on the kernel and on the
   parities of `|A|` and `|B|`, because `P_2 = ((-1)^n + 1) H` and `P_m = 0 (mod 3)` for `m >= 3` (THM-4756; audit:
   all kernels of orders 3, 4, 5 with every ordered pair of ports, kernel_mod3.out). Here `P_m`, the number of
   ordered `m`-tuples of disjoint non-empty paths covering the tournament (THM-5133; `P_1 = H`,
   `P_2 = Pi_2 - 2 H`, `Pi_2` counting the ordered pairs with empty paths allowed), is `m!` times the number of
   unordered ones, so `P_2` and `P_3` are even and the prime 2 is excluded as well. So for module states the laws
   of (1) and (4) can concern only the primes `p >= 5` not dividing the determinant of the visit matrix.
