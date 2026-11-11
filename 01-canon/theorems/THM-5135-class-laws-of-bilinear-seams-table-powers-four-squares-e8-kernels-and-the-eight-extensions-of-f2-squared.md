---
id: THM-5135
title: "Class laws of bilinear seams with uniform states. (1) A unimodular seam of rank r takes a class of valuation v < k modulo p^k with factor (1 - p^-r) sigma_(r-1)(p^v)/p^((r-1)v): it is the r-fold additive power of the multiplication table modulo p^k; for odd m the additive square of the table is m^3 L_m. (2) At an odd prime the sum of 2r squares has the same law iff r is even or p = 1 (mod 4). (3) At 2 the four-square form differs from the rank-2 seam, while the E8 form has the law of the rank-4 seam: over all primes sigma_3(n)/(zeta(4) n^3) times (8 pi^4/3) n^3 is 240 sigma_3(n). (4) A kernel with Gram matrix G has the law (1) at every p not dividing det G, whatever its signature; two unimodular rank-4 forms of the join calculus do not have it, their states being confined. (5) The eight central extensions of F_2^2 by F_2 are the sets of units of square -1: 0 = R^4, 1 = split quaternions (D_4, transitive triangle), 2 = bicomplex numbers (C_4 x C_2, signed triangle), 3 = Hamilton (Q_8, cyclic triangle); commutative iff even; quaternion norms convolve r_2 additively, bicomplex norms multiplicatively. (6) Brahmagupta identity for self-joins."
status: >
  PRE-AUDIT DRAFT (2026-11-11; independent audit pending, statuses below are the author's).
  PROVED (elementary): (1), (2), (4) first part, (6), and the rule of squares in (5).
  FINITE-EXACT: (3) (class counts modulo 2, 4, 8; shells of E8 for n <= 3; r_4, r_8, Hurwitz counts on the stated
  ranges) and the census of (5) (16 cocycles, 8 classes, groups, zero divisors, sign tables).
  CITED (classical, not proved here): Jacobi (r_4(n) = 8 sigma(n), r_8(n) = 16 sigma_3(n) for odd n), Hurwitz
  (24 sigma(n); the dimensions 1, 2, 4, 8), Ingham 1927 (sum of two products), the theta series of E8 and
  Coxeter's integral octonions (240 sigma_3(n)), the classification of unimodular forms over Z_p, the absence of
  positive definite even unimodular lattices of rank 4.
  KNOWN in the repository, credited, not claimed: the factor of (1) for r = 2 (audit of THM-5128) and for every r
  modulo a prime (HYP-9660 (D)); the table multiplicities N_m(t) (modular centres packet, PROOFS.md P4, eq. (4));
  the counts 16, 8, 4, 2 and the Q_8 / D_8 classes of (5) (THM-2887); the line rule "cyclic = Hamilton,
  transitive = split quaternions, the unit of square -1 is the middle vertex" (THM-5133 (D1), kind-pasteur S23,
  PRE-AUDIT, same day; found here independently by enumeration and credited to it); the norm-composing sign
  tables (Session112); Jacobi's and Hurwitz's counts in THM-4826 and Session116; the unimodular K-form of
  THM-4916 (1d) and its proved local factors; Session264's fusion law (3); MISTAKE-484 (Z^8 is not E8),
  MISTAKE-506 (Hamilton's algebra ramifies only at 2 and at the real place).
  NEW in the repository: (1) at prime powers for every r with its table form, (2), (3) (the E8 reading of the
  rank-4 law at every prime), (4), the commutative half of (5) with the bicomplex numbers (which no file names)
  and the additive / multiplicative convolution statement. None of this is new mathematics outside the repository.
  No statement here is about Chan-Pak, which remains OPEN.
session: opus-2026-11-11-S92
date: 2026-11-11
source: 05-knowledge/results/s92_prime_powers_quaternions_bases_20261111/README.md
scripts:
  - 04-computation/experiments/s92_prime_powers_20261111/bridges_exact.py (+ bridges_exact.out; S1-S6, S9)
  - 04-computation/experiments/s92_prime_powers_20261111/seam_rank_exact.py (+ .out)
related:
  - 01-canon/theorems/THM-5134-* (the by-value law of F1; this theorem is its algebra)
  - 01-canon/theorems/THM-5133-flip-and-spin-kernel-signature-sign-tournaments-of-twisted-algebras-and-the-pg32-extremal-theorems.md (PRE-AUDIT)
  - 01-canon/theorems/THM-2887-quaternionic-arf-lift-of-the-semantic-v4-and-global-carry-no-go.md
  - 01-canon/theorems/THM-4916-root-change-local-factor-identities-confinement-multipliers-product-family.md
  - 01-canon/theorems/THM-4826-which-grid-carries-a-norm-shell-one-to-two-split-cubic-character-of-two-centred-units-complete-the-quartic-designs-four-eisenstein-fourth-powers-miss-five-mod-eight.md
  - 05-knowledge/results/session264_constructive_20261111/fusion/PROOF.md, 05-knowledge/results/session261_creative_chanpak_20261111/constructors/BILINEAR_KERNEL.md, 05-knowledge/results/session260_quasiperfect_chanpak_20261111/arithmetic/QUASIPERFECT_NORMS.md
  - 05-knowledge/results/modular_centers_collatz_pell_20261111/PROOFS.md, 05-knowledge/results/session112_counting_interfaces_20261030/norms/NORM_SIGNS.md
  - 05-knowledge/hypotheses/HYP-9666-* (rank-4 seams, OPEN)
---

# THM-5135 - class laws of bilinear seams, and the algebras behind them

A *seam of rank `r`* is a bilinear form `B(u, v)` in two state vectors of length `r`; its value distribution for
independent uniform states is its *class law*. THM-5134 proves that the one-arc join family F1 has, by value, the
class law of the rank-2 seam at every odd modulus. This theorem records the algebra: which laws belong to which
forms, where the classical definite forms agree with them, and which algebras of dimension 4 are involved.
`L_m` is as in THM-5134.

## Statements

**(1) The rank-`r` law.** Let `p` be any prime, `u, v` independent and uniform in `(Z/p^k)^r`, and `B` a bilinear
form with unit determinant. For a class `c` with `v = v_p(c) < k`:

    p^k P(B(u, v) = c) = (1 - p^-r) sigma_(r-1)(p^v)/p^((r-1)v) = (1 - p^-r)(1 + p^-(r-1) + ... + p^-(r-1)v).

For `r = 1` this is the multiplication table: the number of `(a, b)` modulo `m` with `a b = t` is
`N_m(t) = sum_{d | gcd(t, m)} d phi(m/d)` (KNOWN), factor `(1 - 1/p)(v + 1)`. The law of rank `r` is the `r`-fold
additive convolution power of `N_m`; for odd `m`, `(N_m * N_m)(c) = m^3 L_m(c)`, also with one table multiplied by
a unit. For `r >= 2` the product of the factors over all primes is `sigma_(r-1)(n)/(zeta(r) n^(r-1))`; for `r = 1`
it diverges to 0 (the count of `a b = n` is the divisor function).

**(2) Sums of squares at odd primes.** For an odd prime `p` the sum of `2r` squares has the class counts of the
hyperbolic form `sum_{i <= r} a_i b_i` modulo every `p^k` if and only if `(-1)^r` is a square modulo `p`, that is
`r` even or `p = 1 (mod 4)`. So four and eight squares follow the laws of rank 2 and 4 at every odd prime; two and
six squares do so only at `p = 1 (mod 4)` (the character `chi_4`).

**(3) The prime 2 and the definite models.**
(a) Modulo 4 and 8 the four-square form and `a1 a2 + b1 b2` have different class counts. Jacobi (CITED):
`r_4(n) = 8 sigma(n) = pi^2 n * (8/pi^2) sigma(n)/n` for odd `n`: real density, times the rank-2 factors at the odd
primes, times a 2-adic density 1.
(b) The form `x^T C x / 2` of the lattice E8 (`C` the Cartan matrix) has the class counts of the rank-4 seam
modulo 2, 4, 8 (FINITE-EXACT) and at every odd prime (by (2): unimodular of discriminant 1); CITED: the same
holds modulo every `2^k`. Hence over ALL primes the rank-4 series times the real density is

    (8 pi^4/3) n^3 * sigma_3(n)/(zeta(4) n^3) = 240 sigma_3(n),

the number of vectors of norm `2n` in E8, i.e. of integral octonions of norm `n` in Coxeter's order (CITED;
`240, 2160, 6720` recounted for `n = 1, 2, 3`).
(c) `Z^8` is not that model: `r_8(n) = 16 sigma_3(n)` for odd `n` only, `r_8(2) = 112` (MISTAKE-484).
(d) For rank 2 the series over all primes is `sigma(n)/(zeta(2) n)`, Ingham's series for `a1 a2 + b1 b2 = n`
(CITED). No positive definite even unimodular lattice of rank 4 exists (CITED), so a definite quaternary model
must differ from the seam at some prime; Hamilton's does so at 2 only.

**(4) Kernels.** Let `G` be an integer `r x r` matrix and `B(u, v) = u^T G v`. For every prime `p` not dividing
`det G` the class law of `B` modulo every `p^k` is (1). It depends neither on the signature of `G` nor on the sign
of `det G`: the elliptic and hyperbolic two-port kernels of THM-5133 (B) have the same law at such `p` (checked for
`[[5, 2], [2, 1]]`, `[[11, 6], [6, 3]]`, `diag(1, 2)` modulo 9, 25, 49; the second differs at 9, where `3 | det G`;
Session261's kernel of determinant `329 = 7 * 47` has the rank-1 law modulo 7 and 47).
The law is a statement about uniform states. Two unimodular forms of rank 4 occur in the join calculus: the
K-form of THM-4916 (1d), `U U' + E E' - d d' + a d' + a' d` (determinant `-1` in `(U, E, d, a)`), and the corner of
Session264's fusion law, `H(A * B) = h h' + (H - h) h' + h (H' - h') + (U - h)(V' - h') + (d - h)(e' - h')`
(determinant `-1`). For the first, THM-4916 proves the factor `1 + (3p - 2)(p - 1)/(2 p^3)` at the class 0 (two root
changes), not the uniform rank-4 value `1 + p^-3 - p^-4`: its states are confined. So a seam of rank 4 does not by
itself carry the law of (3)(b); see HYP-9666.

**(5) The eight extensions of `F_2^2` by `F_2`.** A normalized 2-cocycle `f` on `V = F_2^2` gives the group on
`V x F_2` with `(a, z)(b, w) = (a + b, z + w + f(a, b))` and the real algebra with basis `e_a`,
`e_a e_b = (-1)^f(a,b) e_(a+b)`. There are 16 cocycles and 8 classes modulo the one non-trivial coboundary, which
changes the sign of the three products `e_a e_b`, `a != b` non-zero. The class is determined by the set `S` of the
non-zero `a` with `e_a^2 = -1`, and every subset of the three points occurs:

| `#S` | classes | group | algebra | commutative | sign table on the three points | zero divisors |
|---|---|---|---|---|---|---|
| 0 | 1 | `C_2^3` | `R^4` | yes | signed triangle | yes |
| 1 | 3 | `D_4` (order 8) | `M_2(R)`, split quaternions | no | transitive tournament, middle vertex `S` | yes |
| 2 | 3 | `C_4 x C_2` | `C + C`, bicomplex numbers | yes | signed triangle | yes |
| 3 | 1 | `Q_8` | Hamilton's quaternions | no | cyclic tournament | no |

The product of the three squares is `+1` in the commutative classes and `-1` in the others, so the algebra is
commutative iff `#S` is even. In a non-commutative class the two cocycles are a tournament on the three points
(`a -> b` iff `e_a e_b = + e_(a+b)`) and its converse; in a commutative class they are two signed triangles and
there is no orientation. Quaternions and bicomplex numbers share two units of square `-1` and differ in the third
(`-1`, anticommuting; `+1`, commuting).
Norms. Hamilton: the sum of four squares. `M_2`: the determinant, which is the seam of F1. `C + C`: the
multiplicative norm is `|z_1|^2 |z_2|^2`. On the orders `Z[i, j, k]` and `Z[i]^2`:

    #{quaternions of norm n} = sum_{a + b = n} r_2(a) r_2(b) = 8 sigma(n)   (n odd),
    #{bicomplex integers of norm n} = sum_{a b = n} r_2(a) r_2(b),          which is 0 for n = 3 (mod 4).

So the quaternion norm is the additive convolution of two Gaussian norms and the bicomplex norm their
multiplicative convolution. In the join calculus the ordered join multiplies path counts (rank 1; commutative) and
the one-arc join adds two products (rank 2): the first has the divisor-function law of (1), the second the law
`sigma(n)/n`.

**(6) Brahmagupta identity.** With `M_s(A, B) = x_A x_B + 2^s y_A y_B` on one-apex states (the F1 value on the sheet
`s >= 1`, and plus one on the sheet 0):

    M_s(A, A) M_s(B, B) = M_s(A, B)^2 + 2^s (x_A y_B - x_B y_A)^2.

The join value is the real part of `z_A conj(z_B)`, `z = x + sqrt(-2^s) y`, and the cross-determinant of the two
states its imaginary part.

## Proofs

**(1).** A linear change of `v` by the matrix of `B` (a unit determinant) reduces to `B(u, v) = u . v`. With
probability `p^(-ri) (1 - p^-r)` the vector `u` is `p^i` times a primitive vector, `i < k`, and then `v -> u . v` is
`p^i` times a surjective homomorphism onto `Z/p^k`, uniform on the multiples of `p^i`. So for `v_p(c) = v < k`:
`P(c) = sum_{i <= v} p^(-ri) (1 - p^-r) p^(i-k)`. For `r = 1` and a modulus `m`: `a` with `gcd(a, m) = d` (there are
`phi(m/d)` of them) takes each multiple of `d` exactly `d` times. The seam `sum a_i b_i` is a sum of `r` independent
table entries; the identity with `L_m` is the case `r = 2` and the Chinese remainder theorem (checked for eleven
odd `m <= 105`, with the units 1, 2, 4). The Euler product is `(1/zeta(r)) prod_p sigma_(r-1)(p^v)/p^((r-1)v)`.

**(2).** Over `Z_p`, `p` odd, a unimodular quadratic form is determined by its rank and its discriminant modulo
squares. The sum of `2r` squares has discriminant 1, the hyperbolic form of rank `r` has `(-1/4)^r`, in the class of
`(-1)^r`. If the classes differ, the numbers of solutions of `Q = 0` modulo `p` differ: they are
`p^(2r-1) + eta (p^r - p^(r-1))`, `eta = +1` for the hyperbolic form and `-1` for the other class. (Checked for
`r <= 4` modulo 3, 5, 7, 11, 13, and modulo 9, 25.)

**(3).** (a), (c) and the shells are enumerations (bridges_exact.out S3, S4). (b): at odd `p`, E8 is unimodular of
rank 8 and determinant 1, hence hyperbolic by the classification used in (2); at 2 the class counts were
enumerated over all `2^8`, `4^8`, `8^8` residue vectors. The displayed identity is `zeta(4) = pi^4/90`.

**(4).** If `p` does not divide `det G`, `v -> G v` is a bijection of `(Z/p^k)^r` and `u^T G v = u . (G v)`. The two
determinants are computed from the displayed forms: in the coordinates `(h, H - h, U - h, d - h)` and
`(h', H' - h', V' - h', e' - h')` the fusion corner has Gram matrix `[[1,1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]]`.

**(5).** The census is an enumeration (bridges_exact.out S6): cocycle condition on the `2^9` normalized cochains;
groups by their numbers of elements of order 4 (`0, 4, 2, 6` for `C_2^3`, `C_4 x C_2`, `D_4`, `Q_8`); zero divisors
searched with coordinates in `{-1, 0, 1}`. Rule of squares: `e_c^2 = (e_a e_b)(e_a e_b) = eps e_a^2 e_b^2` with
`eps = +1` if `e_a, e_b` commute and `-1` if they anticommute (associativity is the cocycle condition), so
`s_a s_b s_c = eps`. The algebras: for `#S = 0`, `R[V]`; for `#S = 2`, `C` tensor `R[C_2] = C + C`; for `#S = 1`,
`(1, -1)_R = M_2(R)`; for `#S = 3`, `(-1, -1)_R`. The two counting formulas are the definitions (a quaternion is a
pair of Gaussian integers with norm the sum of the norms; an element of `Z[i]^2` has norm the product), with
Jacobi's theorem; `r_2(n) = 0` for `n = 3 (mod 4)` and one of `a, b` is `3 (mod 4)` when `a b` is.

**(6).** `(x^2 + t y^2)(x'^2 + t y'^2) = (x x' + t y y')^2 + t (x y' - x' y)^2`.

## Remarks

1. *What (5) says about tournaments.* Only the non-commutative classes are tournaments; that is the dictionary of
   THM-5133 (D), restricted to one line. The bicomplex numbers sit in the commutative half, where the sign table
   is a signed graph. The norm form of a non-commutative class does not see the orientation at odd primes: the
   determinant (transitive triangle, `H = 1`) and the sum of four squares (cyclic triangle, `H = 3`) have the same
   law there, by (2). Orientation is a matter of the prime 2 and of the real place.
2. *Composition algebras.* The hyperbolic form of rank `r` is the norm form of a composition algebra exactly for
   `r = 1, 2, 4` (split complex numbers, split quaternions, split octonions; Hurwitz, CITED). The laws of (1) hold
   for every `r`; rank 3 occurs in the join calculus (kernels with three visits, THM-5133 (A)). The ranks 1, 2, 4
   are not singled out by the gluing calculus.
3. *Scope.* (1)-(4) are statements about uniform independent states. THM-5134 supplies the equidistribution for
   F1; for module states `(H, Pi_2)` it fails at 3 (THM-4756: `Pi_2 = (-1)^n H (mod 3)`).
