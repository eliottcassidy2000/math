---
id: THM-5138
title: "Stein modules of Keller pairs. For a Keller pair (f, g) on C^2 with dual fields d_u, d_v, a centre (a, b), U = f - a, V = g - b: the Stein kernel M_F = Im(d_u - V) + Im(d_v - U) is the image of d - d(UV) ^ on 1-forms and the Stein space S_F = C[x,y]/M_F is the space of classes of the Gaussian integration-by-parts identities in the coordinates (f, g). (A) Identities: [U^i V^j] = i! [i = j] [1]; [h V^m] = [d_u^m h]; [d_u^m d_v^n h] = [h H_(n,m)(U, V)] (complex Hermite); the rotation field U d_u - V d_v maps into M_F. For an automorphism S_F is one-dimensional and GMC(2) says that M_F is a Mathieu-Zhao subspace. (B) Dimension count: for a generic centre, dim S_F = b_1 of the generic fibre of UV = c + sum_i delta_i (deg_u S_i + deg_v S_i) over the components S_i of the non-properness curve (delta_i missing sheets), which is at least 2d - 1 for generic degree d; so JC(2) is the statement that the Stein functional of every Keller pair is unique. (C) Finite chaos: if for every h the classes of d_u^m d_v^n h vanish for m + n large, then dim S_F <= d; hence a Keller pair has finite chaos iff it is an automorphism, and JC(2) <=> every Keller pair has finite Hermite expansions modulo its Stein kernel. (D) Ladder: automorphism <=> finite chaos => the Stein kernel is Mathieu-Zhao (K-GMC) => cone vanishing => one-sided Stein nilpotence; GMC(2) is K-GMC for automorphisms. (E) On the two planar non-proper etale maps in canon (a double cover on C^* x C; a slice of the C^3 Keller map of THM-1300 on C^2 minus a hyperbola) the Stein spaces have dimensions 2 and 11 (linear potentials on the slice: 2, 4, 6), as the count predicts, and one-sided Stein nilpotence fails."
status: >
  PRE-AUDIT DRAFT (independent audit pending; statuses below are the author's).
  (A) PROVED (elementary). The dual fields, their Weyl-algebra module and the untwisted statements are KNOWN
  (THM-4620, THM-4626; Bavula). The automorphism case is a transport of THM-2022 (GMC(2), PROVED in canon).
  (B) PROVED from a CITED comparison theorem (Dimca-Saito 1993 for C^n; Sabbah 1996 and Hien 2009 for smooth affine
  varieties), which was not re-proved here and is checked on ten examples, plus an Euler characteristic count and
  the connectedness of the generic fibre of a Keller polynomial (THM-4620 (5)). 'JC(2) iff the generic fibre of f
  is a line' is KNOWN; the count for the product potential and its reading as uniqueness of the Stein functional
  are specific here.
  (C) PROVED (module argument, elementary) given (B).
  (D) PROVED. The three weaker rungs for all Keller pairs are conjectures (HYP-9668).
  (E) FINITE-EXACT (exact linear algebra modulo the prime 2147483629, stable in the truncation); the two-dimensional
  space and the failure of the Mathieu property on C^* x C are also PROVED by hand (residue functionals).
  Pre-registered predictions S2-S6 (commits d701fad1cd, 90e0a6855d, pushed before the runs): all CONFIRMED.
  CONCURRENT: Session267 (same owner prompt, 2026-11-11) has bridges of a different kind (NC2 for families, a
  rotation operator in the row chart, a shear-word compiler, a wreath law); no Stein functional of a Keller pair, no
  twisted cohomology, no dimension count. JC(2) remains OPEN; no rung of the ladder is proved for a Keller pair
  that is not an automorphism, and none is known to imply JC(2).
session: kind-pasteur-2026-11-11-S24
source: 05-knowledge/results/kps_s24_stein_modules_keller_20261111/README.md
scripts:
  - 05-knowledge/results/kps_s24_stein_modules_keller_20261111/scripts/stein.py, stein_runs.py, stein_more.py, slice_pred.py, c3check.py
related:
  - 01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md
  - 01-canon/theorems/THM-4620-keller-dual-field-module-kernel-ladder-and-de-rham-cokernel.md
  - 01-canon/theorems/THM-4626-surjective-dual-fields-are-locally-nilpotent-and-hyp-9255-is-the-planar-jacobian-conjecture.md
  - 05-knowledge/hypotheses/HYP-9668-gaussian-ladder-for-keller-pairs-between-gmc2-and-jc2.md
---

Notation. `(f, g)` a Keller pair on `C^2` with `{f, g} = 1`; `d_u = {., g}`, `d_v = {f, .}` the dual fields
(`d_u f = d_v g = 1`, `d_u g = d_v f = 0`, commuting; THM-4620). Centre `(a, b)`, `U = f - a`, `V = g - b`.
`T_u = d_u - V`, `T_v = d_v - U`, `M_F = T_u C[x,y] + T_v C[x,y]`, `S_F = C[x,y] / M_F`, classes `[h]`.
`d = [C(x,y) : C(f,g)]`. The same definitions apply to an etale map `(u, v) : X -> C^2` from a smooth affine
surface. Proofs: sections 2-6 of the source note.

## (A) The Stein module

1. `(d - d(UV) ^)(alpha df + beta dg) = (T_u beta - T_v alpha) df ^ dg`; so `S_F = H^2(C^2, d - d(UV) ^)`.
2. `[V k] = [d_u k]`, `[U k] = [d_v k]`; `[U^i V^j] = i! [i = j] [1]`; `[h V^m] = [d_u^m h]`, `[h U^m] = [d_v^m h]`;
   `[d_u^m d_v^n h] = [h H_(n,m)(U, V)]` with `H_(n,m)(U, V) = sum_j (-1)^j j! C(n,j) C(m,j) U^(n-j) V^(m-j)`;
   `(U d_u - V d_v) C[x,y]` lies in `M_F`.
3. If `(f, g)` is an automorphism, `S_F = C [1]`, `M_F` is the kernel of `U^i V^j -> i! [i = j]`, and `M_F` is a
   Mathieu-Zhao subspace (this is GMC(2), THM-2022).

## (B) Dimension count

Let `S` be the non-properness set of `F = (f, g)`, with components `S_i`, and `delta_i = d - #F^(-1)(p)` for
generic `p` on `S_i`. For every centre `(a, b)` such that the closure of `S` in `P^1 x P^1` contains neither
`(a, infinity)` nor `(infinity, b)`:

    dim S_F = b_1(Y) = c + sum_i delta_i (deg_u S_i + deg_v S_i) >= 2d - 1,

`Y` a generic fibre of `(f - a)(g - b)` and `c` the number of its components. `dim S_F = 1` iff `F` is an
automorphism. For an etale map from a smooth affine surface `X`:
`dim S_F = b_1(Y) - rank(H^1(X) -> H^1(Y)) + b_2(X)`.

## (C) Finite chaos

If for every `h` there is `N` with `[d_u^m d_v^n h] = 0` for `m + n >= N`, then `dim S_F <= d` (for any etale map
from a smooth affine surface with `S_F` finite dimensional). A Keller pair on `C^2` has this property at one centre
allowed in (B) iff it is an automorphism. Hence JC(2) <=> every Keller pair has finite chaos <=> every Keller pair
has a one-dimensional Stein space.

## (D) Ladder

For a Keller pair and a centre: automorphism <=> finite chaos => (K-GMC) `M_F` is a Mathieu-Zhao subspace =>
(CV) for every `h`, `[h U^k V^(m+k)] = [h U^(m+k) V^k] = 0` for `m >= N(h)`, `0 <= k <= m` => (SN) for every `h`,
`d_u^m h` and `d_v^m h` lie in `M_F` for `m` large.

## (E) Examples (exact)

- `C^* x C`, `(u, v) = (x^2, y/x)`: `dim S_F = 2`, functionals
  `lambda_+-(x^k y^l) = l! Res_(x = +-sqrt(a)) [x^(k+l+1) e^(b(x^2-a)) / (x^2-a)^(l+1)]`,
  `lambda_+-(x^(-2) V^m) = (-1)^m m! / (2 a^(m+1))` while `[V^m] = 0`: K-GMC and SN fail.
- Slice `F_1 = 1` of the map of THM-1300: `X = C^2 \ {xy = -1}`, `u = y + 3xw`, `v = xw + xw^2 - x^3 w^3`,
  `w = 1/(1+xy)`, Jacobian `-2w^3`, generic degree 3, non-properness curve
  `27 v^2 - 18 uv + 16 + u^3 v - u^2 = 0` of bidegree (3, 2) with two sheets missing. `dim S_F = 11` at the centres
  (2, 3), (5, -7), (1, 1), (-3, 4) and 9 at (0, 0) (an excluded centre); for the linear potentials `u`, `v`,
  `u + v`: 2, 4, 6. At (2, 3) the classes of `w V^m`, `w U^m` (`m <= 5`), `x V^m`, `x U^m` (`m <= 4`) and all computed
  Hermite classes of `w` and `x` (`n, m <= 4`) are nonzero.
