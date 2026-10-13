---
id: THM-4628
title: "Residue colourings of the unit-distance graph of L = Q(zeta_55) at a conjugation-fixed prime with residue field F_(q^2) are proper colourings of the finite Euclidean graph Cay(F_(q^2), N^1) on the norm-one torus, whose eigenvalues are minus the Kloosterman sums over F_q, so they need at least 1 + (q+1)/(2 sqrt q) colours (Hoffman and Weil): every fixed prime has even residue degree and q >= 109, so no residue colouring at any fixed prime uses six or fewer colours; exactly chi >= 7 at 109 and chi >= 12, 15, 16, 20, 21, 23, 24 at 439, 659, 769, 1319, 1429, 1759, 1979; at least 788, 2266, 27737 at the degree-ten residue fields of 19, 29, 79"
status: >
  CORRECTED 2026-10-13 (the session's audit and, independently, the concurrent codex note
  incoming_cylinder_prime_repair_20261013.md, which enumerates all forty unit classes modulo 55): the conjugation-fixed
  unramified primes have residue degree two (class 54) or ten (classes 19, 24, 29, 39) only; the earlier degree-twenty
  example at 2 was false (2^10 = 34 mod 55); the ramified primes 5 and 11 are not fixed either.
  PROVED: (a) the parity of the residue degree, the list of fixed residues, and the image of the unit vectors (Hilbert 90
  for L/L^+, Frobenius = conjugation on the residue field); (b) the torus-to-Kloosterman identity, by the four-line count
  below, valid for every odd q (Medrano, Myers, Stark, Terras, J. Comput. Appl. Math. 68 (1996), CITED, treat the form
  x^2 + y^2, i.e. q = 3 mod 4), verified to 1e-13 for every prime q = -1 mod 55 with q <= 1979 and, by the audit, for every
  character a at 109 and 439 and by full adjacency spectra at q = 7, 11, 23; (c) given Hoffman's bound (CITED) and Weil's
  bound for Kloosterman sums over any finite field (CITED). NUMERICAL (float64 character sums, 1e-13 agreement with the independent Kloosterman computation; the integer
  conclusions beyond the Weil floor rest on margins of at least 0.13): the eight least eigenvalues and Hoffman bounds; the
  Weil floors themselves are PROVED.
  Extends S24's single computation at 109 (HYP-9247) to every conjugation-fixed prime, closing its "unexamined
  residue fields of 19, 29, 79". Independently audited 2026-10-13 (04-computation/experiments/audit_S25_20261013.md):
  CORRECT WITH FIXES (the claimed fixed primes of order twenty do not exist, since l^10 = 1 mod 11; p^+ is inert
  because unramified; the citation scoped; MISTAKE-597).
session: opus-2026-10-13-S25
source: 05-knowledge/results/five_manuscripts_chowla_collatz_residue_bsd_20261013.md (section 4)
scripts:
  - 04-computation/experiments/hn_residue_kloosterman_20261013.py (+ .out)
related:
  - HYP-9247 (every finite unit graph in L three-colourable; residue colourings need at least 7 at 109)
  - THM-4614 (doubling extension; three-colourable exteriors cannot change the six-colour port relation)
  - 05-knowledge/results/six_colour_core_doubling_20261012.md (S24, section 4: the 109 computation)
  - Hoffman (1970) eigenvalue bound; Weil (1948) for Kloosterman sums; Medrano-Myers-Stark-Terras (1996) finite Euclidean graphs
  - the Deligne-Drinfeld manuscript (OpenAI, September 2026; PREPRINT CLAIM): the degeneration whose colouring analogue this theorem refutes
---

# THM-4628 — residue colourings are Kloosterman-forced

## Setting

`L = Q(zeta_55)`, `L^+` its maximal real subfield, `p` a prime of `O_L` fixed by complex conjugation. The decomposition group of `p` is `<l> mod 55` (`l` the rational prime below `p`), and conjugation is `-1 mod 55`; so `p` is fixed iff `-1 in <l>`, iff the residue degree `f` is even and `l^(f/2) = -1 mod 55`. Write `f = 2f'`, `q = l^(f')`, so the residue field is `F_(q^2)` and conjugation induces `x -> x^q`. A *residue colouring* colours the `p`-integral points of `L` through `O_(L,p) -> F_(q^2)` followed by a colouring of `F_(q^2)`; it is proper on the unit-distance graph of `L` iff the colouring of `F_(q^2)` is proper on `G_q = Cay(F_(q^2), S)`, `S` the set of residues of `p`-integral unit vectors.

## Statement

(a) `S = N^1 = {z : z^(q+1) = 1}`, the norm-one torus, of order `q+1`. The conjugation-fixed primes of `L` are exactly the primes above `l = 54 mod 55` (`f' = 1`, `q = l >= 109`) and above `l = 19, 24, 29, 39 mod 55` (`f' = 5`, `q = l^5 >= 19^5`); there are none of degree four or twenty, and `q >= 109` always.
(b) The eigenvalue of `G_q` at the additive character `psi_a`, `a != 0`, depends only on `c = N(a) in F_q^x` and equals `-K_q(c)`, `K_q(c) = sum_(x in F_q^x) e((x + c/x)/q)` for `q` prime (and the corresponding trace form in general).
(c) `chi(G_q) >= 1 + (q+1)/|lambda_min| >= 1 + (q+1)/(2 sqrt q)`; hence no residue colouring of `L` at a conjugation-fixed prime uses six or fewer colours.
(d) Values (`l`, `lambda_min`, Hoffman bound, `chi >=`): `109, -19.6447, 6.599, 7`; `439, -40.0880, 11.976, 12`; `659, -50.2584, 14.132, 15`; `769, -53.5492, 15.379, 16`; `1319, -71.7060, 19.409, 20`; `1429, -74.3776, 20.226, 21`; `1759, -83.1833, 22.158, 23`; `1979, -88.3613, 23.408, 24`; Weil floors for `q = 19^5, 29^5, 79^5`: `787.8, 2265.5, 27736.6`.

## Proof

(a) A unit vector `u in L` (`u ubar = 1`) is `w/wbar` for some `w in L^x` by Hilbert 90 for the quadratic extension `L/L^+`. Since `p` is fixed by conjugation, the prime `p^+ = p cap O_(L^+)` is inert in `L`, and `v_p(w) = v_p(wbar)`; scaling `w` by a power of a uniformizer of `p^+` (an element of `L^+`, so invariant under conjugation) makes `w` a `p`-unit. Conjugation acts on `F_(q^2) = O_L/p` as the Frobenius over `F_q = O_(L^+)/p^+`, so `u mod p = w^(1-q)`. As `w` ranges over `p`-units its residue ranges over `F_(q^2)^x`, and the `(q-1)`-th powers form the kernel of the norm to `F_q`, of order `q+1`. Conversely every element of `N^1` is such a residue. The decomposition-group statements: `(Z/55)^x = Z/4 x Z/10`; `-1 = (-1 mod 5, -1 mod 11)`; `<l>` contains `-1` iff the order of `l` is even and `l^(ord/2) = -1 mod 55`. Order `2`: `l = -1 mod 55`. Order `4`: `l^2 = -1 mod 11` would be needed, impossible since `-1` is a non-residue modulo `11`. Order `20`: `l^10 = -1 mod 11` would be needed, impossible since `l^10 = 1 mod 11` for every `l` prime to `11` (so the prime `2`, of order twenty with `2^10 = 34 mod 55`, is not fixed). Order `10`: `l^5 = -1 mod 55`, i.e. `l = -1 mod 5` and `l^5 = -1 mod 11`, which holds exactly for the residues `19, 24, 29, 39`; the least such prime is `19`. Hence `q = l >= 109 = 2 55 - 1` or `q = l^5 >= 19^5`. Also `p^+` is inert in `L` because it is unramified (`l != 5, 11`) and `p` is fixed by the nontrivial automorphism.
(b) For `z in N^1`, `z^q = z^(-1)`, so `Tr(az) = az + a^q/z`; the sum `lambda_a = sum_(z in N^1) psi(Tr(az))` is invariant under `a -> au`, `u in N^1`, hence depends only on `c = N(a)`, and `lambda_a = sum_(w : N(w) = c) psi(Tr w)`. Count the `w in F_(q^2)` with `Tr w = t` and `N w = c`: they are the roots of `X^2 - tX + c`, so their number is `1 - eta(t^2 - 4c)` (`eta` the quadratic character of `F_q`: two roots in `F_(q^2) \ F_q` when the discriminant is a non-square, one root `t/2 in F_q` when it is zero, none in `F_(q^2) \ F_q` and two in `F_q`, which have trace `2w` and norm `w^2` with zero discriminant, when it is a nonzero square). Hence `lambda_a = sum_t psi(t) - sum_t psi(t) eta(t^2 - 4c) = -sum_t psi(t) eta(t^2 - 4c)`. On the other side, `K_q(c) = sum_(x != 0) psi(x + c/x) = sum_t psi(t) #{x : x + c/x = t} = sum_t psi(t) (1 + eta(t^2 - 4c)) = sum_t psi(t) eta(t^2 - 4c)`. So `lambda_a = -K_q(N(a))` for every odd `q` (the count uses `2` invertible). Verified by direct evaluation of both sides for the eight primes in (d) (maximum discrepancy `4.7e-13`) and by the audit at every character for `109`, `439` and by full spectra at `q = 7, 11, 23`.
(c) Hoffman: for a `d`-regular graph with least eigenvalue `lambda_min`, `chi >= 1 - d/lambda_min`; here `d = q+1`. Weil: `|K_q(c)| <= 2 sqrt q`. The floor `1 + (q+1)/(2 sqrt q)` exceeds `6` for `q >= 101`, and `q >= 109` by (a).
(d) `hn_residue_kloosterman_20261013.py`. ∎

## Scope and boundary

- The theorem concerns colourings that factor through a residue field at a conjugation-fixed prime; colourings periodic under other sublattices (powers of `p`, composite moduli, primes not fixed by conjugation, where the unit circle does not reduce) are not covered.
- Hoffman's bound is a lower bound; the exact chromatic numbers of `G_q` are not computed (S24 noted the same for `109`).
- For `q` a prime power the identity is proved by the same count; the eigenvalues were not evaluated numerically there (only the Weil floor is used).

## Reading

Every periodic colouring of the plane's unit graph that one can build inside `L` by reduction modulo a fixed prime needs far more than seven colours, and the need grows like `sqrt(q)/2`: the finite Euclidean plane is a worse colouring problem than the real one. This is the precise reason the Deligne-Drinfeld degeneration (linear equations, saturated lattices, an injective leading projection) has no colouring analogue: reduction modulo a prime is not a faithful degeneration of a colouring problem, and the plane's number is invisible from every finite residue field of `L`.
