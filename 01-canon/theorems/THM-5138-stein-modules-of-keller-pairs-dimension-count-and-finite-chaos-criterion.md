---
id: THM-5138
title: "Stein modules of Keller pairs. For a Keller pair (f, g) on C^2 with dual fields d_u, d_v, a centre z = (a, b), U = f - a, V = g - b: the Stein kernel M_z = Im(d_u - V) + Im(d_v - U) is the image of d - d(UV) ^ on 1-forms, and the Stein space S_z = C[x,y]/M_z is the fibre at z of C[x,y] as a module over the polynomial ring C[E_a, E_b], E_a = f - d_v, E_b = g - d_u. (A) Identities: [U^i V^j] = i! [i = j] [1]; [h V^m] = [d_u^m h]; [d_u^m d_v^n h] = [h H_(n,m)(U, V)]; for an automorphism S_z is one-dimensional and GMC(2) says that M_z is a Mathieu-Zhao subspace. (B) Dimension count: dim S_z = b_1 of the generic fibre of (f - a)(g - b) = c + sum_i delta_i (deg_u S_i + deg_v S_i) = c + ((deg f + deg g)/gcd) sum_i delta_i rho_i >= 2d - 1 over the components S_i of the non-properness curve; it is 1 iff the pair is an automorphism. (C) Finite chaos at one centre (all classes [d_u^m d_v^n h] vanish for m + n large) holds iff S_z = 0 or the map is an isomorphism; for a Keller pair, iff it is an automorphism. (D) The module is torsion-free: a nonzero polynomial with coefficients in k_0 does not lie in M_z when a, b are algebraically independent over k_0 (injectivity of the Stein transform). (E) Detection: at such a centre a Keller pair is an automorphism iff some iterates d_u^r x and d_u^s y lie in M_z, iff SN, iff CV, iff M_z is Mathieu-Zhao; hence JC(2) holds iff the Stein kernel of every Keller pair is a Mathieu-Zhao subspace (the instance q = fg of the statement of the Image Conjecture of Zhao, which is false in general), and GMC(2) is the pair (x, y). (F) At every centre: isomorphism => K-GMC => CV => SN => first hits. (G) Exact examples on three non-proper planar etale maps: Stein spaces of dimension 2 (double cover on C^* x C), 11 (slice of the C^3 Keller map of THM-1300; linear potentials 2, 4, 6) and 3 (inclusion of C^2 minus a hyperbola); finite chaos fails wherever S_z != 0, SN fails at every centre in general position, and at special centres of the third map SN and CV hold without an isomorphism."
status: >
  First audit (audit/AUDIT.md in the source packet): PASS WITH CORRECTIONS, every number reproduced by an
  independent engine. Parts (D), (E), item (A4), the present form of (C) and the special-centre example were written
  after it; second audit (audit/AUDIT-2.md, of the commit c440c4a4b7): PASS WITH CORRECTIONS, no mathematical error
  found in them (the valuation lemma tested on exact series on five maps; no torsion found by brute force on five
  maps), one major correction (the Image Conjecture was already refuted in canon, THM-2801), sixteen minor, five
  strengthenings. All corrections of both audits are applied (MISTAKE-764).
  (A) PROVED (elementary). The dual fields, their Weyl-algebra module and the untwisted statements are KNOWN
  (THM-4620, THM-4626; Bavula). The automorphism case is a transport of THM-2022 (GMC(2), PROVED in canon). GMC as
  a statement about images of first-order operators is Zhao's; "S_z is finite dimensional" is Zhao 2010,
  Corollary 2.11. The description of the family of Stein spaces as one module is the first audit's.
  (B) PROVED from a CITED comparison theorem (Dimca-Saito 1993, Theorem 0.3, for C^n; Sabbah 1996, Theorem 1.1 and
  Corollary 1.2, for smooth affine varieties; already cited in the repo for the potential W = f, the case of
  Gurjar-Masuda-Miyanishi 2013), not re-proved here, plus an Euler characteristic count, the connectedness of the
  generic fibre of a Keller polynomial (THM-4620 (5)) and THM-4122 (4) (Nguyen Van Chau) for the second formula.
  (C) PROVED. The sharp statement and the trace argument are the first audit's; the Weyl-module map is the
  author's. For Keller pairs "S_z != 0" comes from (B) (CITED) at an arbitrary centre and from (D) at a centre in
  general position.
  (D) PROVED (a valuation at infinity of the target plane; elementary). The reduction of membership at a centre in
  general position to a rational certificate is Session269's; the reduction of one-direction detection to
  torsion-freeness was sketched by the first audit, which left torsion-freeness unproved in general. With the CITED
  holonomicity of the twisted module the exceptional set of each h is a proper algebraic subset (second audit).
  (E) PROVED from (D), (F), THM-2022 and the classical criterion THM-4620 (5) (a Keller pair with a locally nilpotent
  dual field is an automorphism: KNOWN). It does not use (B). CONCURRENT: Session269
  (05-knowledge/results/session269_centres_ports_switches_20261111/jacobian/proof.md, pushed 2026-11-11
  13:03-13:08, built on the pre-audit checkpoint of this theorem) proved first, by rational certificates, that first
  hits in both pure directions on an open set of centres give mixed finite chaos at very general centres, and
  derived the collapse of the ladder conditionally on (B) and (C); first-pusher credit for moving-centre detection
  in two directions and for the collapse is theirs. What is added here: one direction, two test polynomials, no
  use of (B), and the surface case; testing at one centre algebraically independent over the coefficients is the
  formulation of Theorem 2a of Session269. K-GMC is the instance q = (f - a)(g - b) of the statement of the Image
  Conjecture of Zhao (J. Algebra 324 (2010), Conjecture 1.3); that conjecture is REFUTED in general in
  characteristic zero (THM-2801: SIC(n) false for n >= 2), so the instance inherits no support from it. The nearest
  published equivalence is van den Essen-Wright-Zhao 2011, Theorem 4.3; in the repo the vacuum-Mathieu programme
  HYP-8350. Two literature sweeps (the audits) did not find the product potential, the count, torsion-freeness of
  this module, or these equivalences; no priority is claimed (torsion-freeness may be folklore).
  (F) PROVED.
  (G) FINITE-EXACT (exact linear algebra modulo the prime 2147483629; dimensions stable in the truncation and
  reproduced by the audit with other primes and filtrations; a nonzero normal form certifies non-membership
  relative to the truncation only). PROVED by hand: the two-dimensional space and the failure of K-GMC and SN on
  C^* x C for centres with a != 0, and S_z = 0 at a = 0; for the inclusion of C^2 minus {1 + xy = 0}: dimension 2,
  class of 1 nonzero, SN and CV at the centre (0, 0), one-sided SN on the lines ab = 0 (second audit). PROVED from
  (C), (D): on the slice finite chaos fails where S_z != 0, and SN, CV, K-GMC fail at every centre algebraically
  independent over Q (all classes of the iterates of x, y, w nonzero: second audit). NOT proved: failure of SN, CV,
  K-GMC at the rational centre (2, 3) of the slice (no vanishing observed); K-GMC at the special centre (0, 0) of
  the inclusion (not decided).
  Pre-registered predictions S2-S6 (commits d701fad1cd, 90e0a6855d): all CONFIRMED.
  JC(2) remains OPEN. Nothing here excludes a Keller pair; (E) is a reformulation, and at centres in general
  position the Gaussian conditions are equivalent to the classical criterion. Detection at a prescribed centre is
  OPEN on C^2 and false on surfaces (HYP-9668).
session: kind-pasteur-2026-11-11-S24
source: 05-knowledge/results/kps_s24_stein_modules_keller_20261111/README.md
scripts:
  - 05-knowledge/results/kps_s24_stein_modules_keller_20261111/scripts/stein.py, stein_runs.py, stein_more.py, stein_centres.py, stein_inject.py, slice_pred.py, c3check.py
related:
  - 01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md
  - 01-canon/theorems/THM-4620-keller-dual-field-module-kernel-ladder-and-de-rham-cokernel.md
  - 01-canon/theorems/THM-4626-surjective-dual-fields-are-locally-nilpotent-and-hyp-9255-is-the-planar-jacobian-conjecture.md
  - 01-canon/theorems/THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction.md
  - 05-knowledge/hypotheses/HYP-9668-gaussian-ladder-for-keller-pairs-between-gmc2-and-jc2.md
  - 05-knowledge/results/session269_centres_ports_switches_20261111/jacobian/proof.md
---

Notation. `(f, g)` a Keller pair on `C^2` with `{f, g} = 1`; `d_u = {., g}`, `d_v = {f, .}` the dual fields
(`d_u f = d_v g = 1`, `d_u g = d_v f = 0`, commuting; THM-4620). Centre `z = (a, b)`, `U = f - a`, `V = g - b`.
`T_u = d_u - V`, `T_v = d_v - U`, `M_z = T_u C[x,y] + T_v C[x,y]`, `S_z = C[x,y] / M_z`, classes `[h]_z`.
`d = [C(x,y) : C(f,g)]`. The same definitions apply to an etale map `(u, v) : X -> C^2` from a connected smooth
affine surface, `B = O(X)`. Proofs: sections 2-6 of the source note.

## (A) The Stein module

1. `(d - d(UV) ^)(alpha df + beta dg) = (T_u beta - T_v alpha) df ^ dg`; so `S_z = H^2(C^2, d - d(UV) ^)`.
2. `[V k] = [d_u k]`, `[U k] = [d_v k]`; `[U^i V^j] = i! [i = j] [1]`; `[h V^m] = [d_u^m h]`, `[h U^m] = [d_v^m h]`;
   `[d_u^m d_v^n h] = [h H_(n,m)(U, V)]` with `H_(n,m)(U, V) = sum_j (-1)^j j! C(n,j) C(m,j) U^(n-j) V^(m-j)`;
   `(U d_u - V d_v) C[x,y]` lies in `M_z`.
3. If `(f, g)` is an automorphism, `S_z = C [1]`, `M_z` is the kernel of `U^i V^j -> i! [i = j]`, and `M_z` is a
   Mathieu-Zhao subspace (this is GMC(2), THM-2022).
4. `E_a = f - d_v` and `E_b = g - d_u` commute, `M_z = (E_a - a) C[x,y] + (E_b - b) C[x,y]`, and
   `[d_u, E_a] = [d_v, E_b] = 1`, `[d_u, E_b] = [d_v, E_a] = 0`: the spaces `S_z` are the fibres of `C[x,y]` over
   `C[E_a, E_b]`, with `d_u, d_v` acting as `d/da, d/db`. For the standard pair `e^(d_x d_y)` trivializes the module
   and `[h]_z = (e^(d_x d_y) h)(a, b) [1]_z`.

## (B) Dimension count

Let `S` be the non-properness set of `F = (f, g)`, with components `S_i`, `delta_i = d - #F^(-1)(p)` for generic
`p` on `S_i`, `G = gcd(deg f, deg g)` and `rho_i >= 1` the integers of THM-4122 (4). For every centre:

    dim S_z = b_1(Y) = c + sum_i delta_i (deg_u S_i + deg_v S_i) = c + ((deg f + deg g)/G) sum_i delta_i rho_i >= 2d - 1,

`Y` a generic fibre of `(f - a)(g - b)` and `c >= 1` the number of its components (`c = 1` for generic `z`).
`dim S_z = 1` iff `F` is an automorphism; `S_z != 0` at every centre. For an etale map from a smooth affine surface
`X`: `dim S_z = b_1(Y) - rank(H^1(X) -> H^1(Y)) + b_2(X)`, for centres off the lines through the boundary points
`(a, infinity)`, `(infinity, b)` of the closure of `S`.

## (C) Finite chaos

Let `F : X -> C^2` be etale, `X` connected. If for every `h` there is `N` with `[d_u^m d_v^n h]_z = 0` for
`m + n >= N`, then `S_z = C [1]`, and if `[1]_z != 0` then `F` is an isomorphism onto `C^2`. So finite chaos at `z`
holds iff `S_z = 0` or `F` is an isomorphism. A Keller pair on `C^2` has finite chaos at one centre iff it is an
automorphism.

## (D) Injectivity of the Stein transform

1. Every extension `w` to `Frac(B)` of the valuation `-deg_(u,v)` of `C(u, v)` satisfies
   `w(d_u phi) >= w(phi) + 1`, `w(d_v phi) >= w(phi) + 1`.
2. For `P(a, b) != 0` of degree `n` and `phi != 0`: `w(P(E_a, E_b) phi) = w(phi) - n`. So `B` is a torsion-free
   module over `C[E_a, E_b]`.
3. If `X`, `F` and `h != 0` are defined over a subfield `k_0` of `C` and `a, b` are algebraically independent over
   `k_0`, then `h` is not in `M_(a,b)`. For every `h != 0` the set of centres with `[h]_z = 0` lies in a countable
   union of proper algebraic subsets (in a proper algebraic subset, by the CITED holonomicity of the module).

## (E) Detection

For a Keller pair with coefficients in `k_0` and a centre `z` with coordinates algebraically independent over
`k_0`, the following are equivalent: (1) automorphism; (2) `M_z` is a Mathieu-Zhao subspace; (3) CV at `z`;
(4) `d_u^m x` and `d_u^m y` lie in `M_z` for `m` large; (5) some `d_u^r x` and some `d_u^s y`, `r, s >= 1`, lie in
`M_z`; (6) `d_u` is locally nilpotent. The same with `d_v`, and with "at `z`" replaced by "at every centre of a set
not contained in a countable union of proper algebraic subsets". At such a centre a polynomial `h` over `k_0` has
some iterate `d_u^r h` in `M_z` iff `h` lies in `k_0[f, g]` (THM-4620 (1)): item (5) says that `x` and `y` are
polynomials in `f` and `g`.

For an etale map of a connected smooth affine surface defined over `k_0`, such a centre and `h` defined over `k_0`:
`d_u^m h` lies in `M_z` iff `d_u^m h = 0`; SN in both directions, or K-GMC, holds at `z` iff `F` is an isomorphism
onto `C^2`.

Hence the following are equivalent: JC(2); for every Keller pair `(f, g)`,
`Im(d/dx - (fg)_x) + Im(d/dy - (fg)_y)` is a Mathieu-Zhao subspace of `C[x, y]`; CV, or SN, or SN in the direction
`d_u`, for all Keller pairs and centres; and (using (B) and holonomicity, CITED) `C[x, y]` has rank one over
`C[f - d_v, g - d_u]` for every Keller pair.

## (F) Ladder at one centre

Isomorphism => (K-GMC) `M_z` is a Mathieu-Zhao subspace => (CV) for every `h`,
`[h U^k V^(m+k)] = [h U^(m+k) V^k] = 0` for `m >= N(h)`, `0 <= k <= m` => (SN) for every `h`, `d_u^m h` and
`d_v^m h` lie in `M_z` for `m` large => first hits. Finite chaos => K-GMC.

## (G) Examples (exact)

- `C^* x C`, `(u, v) = (x^2, y/x)`: for `a != 0`, `dim S_z = 2`, functionals
  `lambda_+-(x^k y^l) = l! Res_(x = +-sqrt(a)) [x^(k+l+1) e^(b(x^2-a)) / (x^2-a)^(l+1)]` (one half of the local
  Gaussian functionals at the two preimages of the centre), `lambda_+-(x^(-2) V^m) = (-1)^m m! / (2 a^(m+1))` while
  `[V^m] = 0`: K-GMC and SN in the direction `d_u` fail; `d_v` is locally nilpotent. For `a = 0`: `S_z = 0`.
- Slice `F_1 = 1` of the map of THM-1300: `X = C^2 \ {xy = -1}`, `u = y + 3xw`, `v = xw + xw^2 - x^3 w^3`,
  `w = 1/(1+xy)`, Jacobian `-2w^3`, generic degree 3, non-properness curve
  `27 v^2 - 18 uv + 16 + u^3 v - u^2 = 0` of bidegree (3, 2) with two sheets missing. `dim S_z = 11` at the centres
  (2, 3), (5, -7), (1, 1), (-3, 4), (7, 2), (-1, -2) and 9 at (0, 0) (centre on an excluded line); for the linear
  potentials `u`, `v`, `u + v`: 2, 4, 6. At (2, 3) the normal forms of `w V^m`, `w U^m` (`m <= 5`; audit `m <= 9`,
  `m <= 12`), `x V^m`, `x U^m` (`m <= 4`) and of all computed Hermite classes of `w` and `x` are nonzero. With
  `r = xw` the slice map is the restriction of the finite cubic map `(r, u) -> (u, 2r^3 - u r^2 + 2r)` to the
  complement of its critical curve `w = 1 - ru + 3r^2 = 0`; the dual fields raise the order of pole along that curve
  by exactly 2, so no iterate of `x`, `y`, `w` vanishes, and by (E) all their classes are nonzero at every centre
  algebraically independent over `Q`: SN fails there in both directions.
- Inclusion of `C^2 \ {1 + xy = 0}` in `C^2`, centre (0, 0): the Stein kernel is graded by the charge
  `deg x - deg y`, contains every element of nonzero charge and has codimension 2 in charge 0; `[1] != 0`; SN and CV
  hold, finite chaos fails, the map is not an isomorphism. On the lines `ab = 0` one-sided SN holds. At (2, 3) the
  dimension is 3 and nothing vanishes.
