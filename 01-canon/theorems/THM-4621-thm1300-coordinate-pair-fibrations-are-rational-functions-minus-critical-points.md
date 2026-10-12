---
id: THM-4621
title: "THM-1300's Keller map on its three coordinate-pair fibrations: each generic fibre is a punctured line on which the remaining coordinate is a degree-three function whose critical points are exactly the punctures — a cubic polynomial in y minus its two critical points for (F_1, F_2), 2aX + 2/X - c/X^2 on P^1 minus five points for (F_1, F_3), s(cs^2 - 2s + b)/2 on P^1 minus three points for (F_2, F_3); the fibre cubic is Gao's tangency polynomial w^3 - 2w^2 + bcw - 2ac^2 under w = 3ac/(b - y); the Jelonek set is {27a^2c^2 - 18abc + 16a + b^3c - b^2 = 0}, the image misses exactly the curve (4/(27t^2), 4/(3t), t), fibre sizes are 3, 1, 0; the Galois group of K(x,y,z)/K(F) is S_3; every ring of constants R_j is K[F_(i != j)], the locally nilpotent part of every dual field is K[F], x, y, z lie outside it, and no dual field is surjective (x^2 not in V_3 K[x], x not in V_2 K[x], x^2 not in V_1 K[x]); the Keller module has A_3-length at least three"
status: >
  PROVED: the parametrizations, restrictions, derivative identities, the discriminant, the image
  complement and the tangency-polynomial identity are exact symbolic identities (sympy, rational
  function arithmetic; 04-computation/experiments/logic4_20261012/thm1300_fibrations.py); the
  consequences (fibre topology, Galois group, rings of constants, nilpotent parts, non-surjectivity
  witnesses by residues, the length bound) are elementary arguments on top of THM-4620 and standard
  facts (minimal extensions of local systems, j_+ of a rank-one local system). FINITE-EXACT hostile
  controls: dim Ker V_3 on polynomials of degree <= 12 equals 4 = #{1, F_2, F_2^2, F_1}; iterate
  degrees deg V_j^m(x_i), m <= 6, grow linearly with slopes 7, 8, 10. CITED: Gao, arXiv:2608.00222,
  Theorem 3.4 (fibre stratification of F, the polynomial c_3, the cusp curve, the tangency polynomial;
  derived here independently from the (F_1, F_2)-fibre cubic and found to agree), THM-1310 (fibre
  sizes in {0, 1, 3} modulo odd primes), THM-4562 (no dual field is locally nilpotent), THM-4567.
  Not independently audited yet (audit requested this session).
session: mac-mini-2026-10-12-logic4
source: 05-knowledge/results/logic_lenses_jc2_collatz_hn_20261012.md
scripts:
  - 04-computation/experiments/logic4_20261012/thm1300_fibrations.py (+ .out)
  - 04-computation/experiments/logic4_20261012/thm1300_module_checks.py (+ .out)
related:
  - THM-1300 (the map; det JF = -2 re-verified here), THM-1310, THM-4562, THM-4567
  - THM-4620 (the general dictionary this instantiates)
  - HYP-9257 (the exact A_3-length)
  - THM-4625 (the same anatomy for every tangent-sweep Keller map; THM-1300 is the d = 2 member)
  - S. Gao, Counterexamples to the Jacobian conjecture in dimensions greater than two, arXiv:2608.00222 (tangent sweep; Theorem 3.4)
---

# THM-4621 — the three coordinate-pair fibrations of THM-1300's map

## The map

With `u = 1 + xy`,

    F_1 = u^3 z + y^2 u (4 + 3xy),   F_2 = y + 3x u^2 z + 3x y^2 (4 + 3xy),   F_3 = 2x - 3x^2 y - x^3 z,

`det JF = -2`, degrees `7, 6, 4`. The identity `(F_2 - y)(1 + xy) = 3x F_1` holds in `K[x, y, z]`.
Write `G_(12) = (F_1, F_2)`, `G_(13) = (F_1, F_3)`, `G_(23) = (F_2, F_3)`, with target coordinates
`(a, b)`, `(a, c)`, `(b, c)`. The dual fields `V_1, V_2, V_3` (THM-4620) have polynomial coefficients
and `V_j F_i = δ_(ij)`.

## Statement

1. **The `(F_1, F_2)`-fibration.** For `a ≠ 0`, put `Q(y) = y^2 - by + 3a`. The fibre `G_(12)^(-1)(a, b)`
   is the image of the bijective parametrization

       y ↦ ( x, y, z ) = ( (b - y)/Q,  y,  -Q · (2y^4 - by^3 + (30a - b^2) y^2 + 6ab y - 9a^2) / (27a^2) ),

   on which `u = 3a/Q`; no point of the fibre has `Q(y) = 0`. So the fibre is `A^1_y ∖ {Q = 0}`:
   `P^1` minus three points when `b^2 ≠ 12a`, minus two when `b^2 = 12a`. On it

       F_3 = -(y - b)(2y^2 - by + 18a - b^2)/(27a^2) = κ (y^3/3 - b y^2/2 + 3ay) + c_0,
       κ = -2/(9a^2),   c_0 = 2b/(3a) - b^3/(27a^2),       dF_3/dy = κ Q(y).

   The fibre of `F` over `(a, b, c)` is the set of roots of the cubic `F_3(y) = c` minus the roots of
   `Q`: `F` restricted to the fibre is a cubic polynomial with its critical points deleted.
2. **Fibres of `F`, Jelonek set, image (re-derived).** The discriminant of the cubic is
   `-4 c_3(a, b, c)/(729 a^6)` with

       c_3 = 27a^2 c^2 - 18abc + 16a + b^3 c - b^2,

   Gao's `c_3(v)`. Hence (for `a ≠ 0`): three preimages when `c_3 ≠ 0`; exactly one when `c_3 = 0` and
   `b^2 ≠ 12a` (a double root sits at a puncture); none on the curve `C = {(4/(27t^2), 4/(3t), t)}`
   (`Q = (y - 2/(3t))^2`, `F_3(2/(3t)) = t`, the cubic has a triple root at the puncture). No fibre has
   exactly two points. The Jelonek non-properness set is `{c_3 = 0}` and the image of `F` is
   `A^3 ∖ C` (agreeing with Gao, Theorem 3.4, and with THM-1310's `N ∈ {0, 1, 3}`).
3. **Tangency polynomial.** Under `w = 3ac/(b - y)` the equation `F_3(y) = c` becomes
   `(9c/2)(w^3 - 2w^2 + bc·w - 2ac^2) = 0`: the fibre cubic is Gao's tangency polynomial
   `w^3 - 2w^2 + Xw - 2Y` at the sweep target `(X, Y) = (bc, ac^2) = (F_2 F_3, F_1 F_3^2)`.
4. **Galois group.** `K(x, y, z) = K(F)(y)`, the cubic is the minimal polynomial of `y` over `K(F)`,
   `deg F = 3`, and since `c_3` is not a square in `K(a, b, c)` the Galois group of the Galois closure
   of `K(x, y, z)/K(F)` is `S_3`. Equivalently the monodromy of the finite étale triple cover over
   `A^3 ∖ {c_3 = 0}` is the full symmetric group, and the local monodromy around a generic point of the
   Jelonek surface is a transposition (`c_3` divides the discriminant exactly once).
5. **The `(F_1, F_3)`-fibration.** For `ac ≠ 0`, with `X = x/(1 + xy)` and `v = X/(aX^3 - X + c)`,
   the fibre `G_(13)^(-1)(a, c)` is parametrized bijectively by

       x = vX,   y = (v - 1)/x,   z = (2x - 3x^2 y - c)/x^3,

   so it is `P^1_X` minus `{0, ∞}` and the three roots of `aX^3 - X + c` (five points generically).
   On it `F_2 = 2aX + 2/X - c/X^2` and `dF_2/dX = 2(aX^3 - X + c)/X^3`.
6. **The `(F_2, F_3)`-fibration.** For `c ≠ 0`, with `w = xy`, the fibre `G_(23)^(-1)(b, c)` is the plane
   conic `b x^2 - (6 + 4w) x + 3c(1 + w)^2 = 0` (minus its point `(0, -1)`), parametrized by the lines
   `w = -1 + sx` through that point:

       x = 2/(3cs^2 - 4s + b),   w = -1 + sx,   z = (2x - 3x^2 y - c)/x^3,

   so it is `P^1_s` minus `{∞}` and the two roots of `3cs^2 - 4s + b` (three points). On it
   `F_1 = s(cs^2 - 2s + b)/2`, `dF_1/ds = (3cs^2 - 4s + b)/2`, and `u = 2s/(3cs^2 - 4s + b)`.
7. **Rings of constants and nilpotent parts.** Each generic fibre is geometrically connected, so
   (THM-4620(1)) `R_1 = K[F_2, F_3]`, `R_2 = K[F_1, F_3]`, `R_3 = K[F_1, F_2]`, the locally nilpotent
   part of every `V_j` is `N_j = K[F]`, and `Ker V_j^m = (+)_(i<m) F_j^i K[F_(i≠j)]`. The coordinates
   `x, y, z` lie outside `K[F]`, so `V_j^m(x), V_j^m(y), V_j^m(z) ≠ 0` for all `j, m` (computed through
   `m = 6`: degrees `1 + 7m`, `1 + 8m`, `1 + 10m` for `x` and `z`, two less for `y`).
8. **No dual field is surjective.** The generic ranks of the cokernels (THM-4620(3)) are
   `μ_3 = 2`, `μ_2 = 4`, `μ_1 = 2`. Explicit witnesses, by residues on the fibre curves:
   * `x^2 ∉ V_3 K[x, y, z]`: on the `(F_1, F_2)`-fibre with `Q = (y - r_1)(y - r_2)`, the residues of
     `x^2 dF_3` at `r_1, r_2` are `-2/(r_1^2 (r_1 - r_2))`, `2/(r_2^2 (r_1 - r_2))` (also `x^2 y`, `x^3`);
   * `x ∉ V_2 K[x, y, z]`: the residue of `x dF_2` at `X = 0` is `2` (also `y`, `z`, `y^2`, `x^2`);
   * `x^2 ∉ V_1 K[x, y, z]`: the residue of `x^2 dF_1` at a root `σ` of `3cs^2 - 4s + b` is `1/(3cσ - 2)`
     (also `x^3`, `x^2 y`).
   So the first-order sentence `∀v ∃w (θw = v)` fails for `θ = V_1, V_2, V_3` and holds for the
   standard module (THM-4620, reading).
9. **Length.** The twisted `A_3`-module `V_F = K[x, y, z]` has length at least three: the simple
   submodule `K[F]`, the minimal extension of the irreducible rank-two local system (the standard
   representation of the `S_3`-monodromy) on `A^3 ∖ {c_3 = 0}`, and a module supported on the Jelonek
   surface (at a generic point of it one sheet stays and the two escaping sheets are swapped by the
   local monodromy, so `V_F` is locally `O ⊕ O[1/s] ⊕ j_+ L_(-1)`, and `O[1/s]/O` is a delta module along
   the surface). The exact length is HYP-9257.

## Proofs

**(1).** With `x = (b - y)/Q` and `u = 1 + xy = (Q + y(b - y))/Q = 3a/Q`, the identity
`(F_2 - y) u = 3x F_1` shows that `F_1 = a` and `F_2 = b` are equivalent on `{u ≠ 0}` to `(b - y) u = 3ax`
together with `z = (a - y^2 u (4 + 3xy))/u^3`; the first equation is `x = (b - y)/Q` when `Q(y) ≠ 0`,
and at a root `r` of `Q` it reads `b - r = 0`, impossible since `Q(b) = 3a ≠ 0`. Points with `u = 0`
have `F_1 = 0 ≠ a`. The displayed `z` and `F_3` are the substitutions (exact, script). `F_3` is a
polynomial in `y` because the fibre is a smooth affine curve on which the polynomial vector field
`V_3` restricts to a derivation `c(y) d/dy` with `c = V_3(y)` regular; `V_3 F_3 = 1` forces
`F_3'(y) = 1/c(y)` to have poles only where `F_3` has poles, and the computation shows there are none.
`dF_3/dy = κQ` vanishes exactly at the punctures, as it must for an unramified `F`.

**(2).** The fibre over `(a, b, c)` is `{y : F_3(y) = c, Q(y) ≠ 0}`. A root of `F_3 - c` is multiple iff
it is a root of `F_3' = κQ`, i.e. a puncture. If `Q` has distinct roots `r_1 ≠ r_2` with
`F_3(r_1) = F_3(r_2) = c`, the cubic `F_3 - c` would have two double roots: impossible. So over
`{c_3 = 0}` with `b^2 ≠ 12a` one double root is deleted and one simple root remains; over `C`
(`b^2 = 12a`, `c = F_3(r)`) the cubic is `κ(y - r)^3/3` and the fibre is empty. The discriminant and the
identity `disc = -4c_3/(729a^6)` are exact (script); `c_3` is irreducible and appears to the first
power. The locus `a = 0` is excluded from the parametrization but carries no fibre of size two either
(Gao, Theorem 3.4, whose stratification is complete).

**(3)** is the substitution `y = b - 3ac/w` (exact; the factor `9c/2`). Gao's sweep coordinates are
`X = F_2 F_3`, `Y = F_1 F_3^2` (from his twisted form `F = (Y/C^2, X/C, C)`).

**(4).** `x` and `z` are rational in `y` over `K(a, b)` on the fibre, so `K(x, y, z) = K(F)(y)`; the
cubic has degree `3 = deg F` (three generic preimages) hence is irreducible; its Galois group is `S_3`
iff the discriminant is not a square in `K(F)`. The local-monodromy statement is the standard fact
that a simple zero of the discriminant of a family of simple polynomials gives a transposition.

**(5), (6).** Exact substitutions (script): the displayed points satisfy `F_1 = a, F_3 = c`
(resp. `F_2 = b, F_3 = c`) identically, the parameter is recovered as `X = x/(1 + xy)`
(resp. `s = (w + 1)/x`), and the remaining coordinate restricts as displayed with the displayed
derivative. For (5): on `F_3 = c ≠ 0` one has `x ≠ 0` and `z = (2x - 3x^2 y - c)/x^3`; then
`a = F_1 = [(1 + w)^3 (x(2 - 3w) - c) + x w^2 (1 + w)(4 + 3w)]/x^3` and the bracket collapses to
`(1 + w)[(1 + w)(2 + w) x - (1 + w)^2 c] + ...`, giving `a x^3 - (1 + w)(2 + w) x + c(1 + w)^3 = 0`; with
`X = x/(1 + w)` this is `v(aX^3 - X + c) = X`. The punctures are where `v → ∞` (roots of
`aX^3 - X + c`), `X = 0` (`v = 0`, i.e. `u = 0`, impossible for `a ≠ 0`) and `X = ∞`. For (6): on
`F_3 = c` the equation `F_2 = b` becomes `b x^2 = x(6 + 4w) - 3c(1 + w)^2` (the bracket
`w + 3(1 + w)^2 (2 - 3w) + 3w^2 (4 + 3w) = 6 + 4w`, exact); the conic passes through `(x, w) = (0, -1)`,
which is not on the fibre (`x = 0` gives `F_3 = 0`), and the pencil of lines through it gives the
parametrization; `x = ∞` at the two roots of `3cs^2 - 4s + b` and `(x, w) → (0, -1)` as `s → ∞`. In
both cases the derivative of the remaining coordinate vanishes exactly at the finite punctures.

**(7).** Connected generic fibres give `L_j = K(G_j)` (THM-4620(1)), so `R_j ⊆ K(G_j) ∩ K[x]`. If
`P(G_j)/Q(G_j) ∈ K[x]` with `P, Q` coprime in `K[X, Y]` and `q` an irreducible factor of `Q`, then
`q(G_j)` is a non-constant polynomial (the components of `G_j` are algebraically independent) and
vanishes somewhere; the image of `G_j` is open (`G_j` is smooth) and meets the curve `{q = 0}` in a
dense subset, on which `P` must vanish; so `q | P`, a contradiction unless `Q` is constant. Hence
`R_j = K[G_j]` and `N_j = R_j[F_j] = K[F]`. On the generic `(F_1, F_2)`-fibre, `x` has poles at the
punctures while every element of `K[F] = K[a, b][F_3]` restricts to a polynomial in `y`; `y` has
degree one, not a multiple of three; and `z` (degree six) is not a quadratic polynomial in
`F_3|_fibre` (the linear system for its coefficients has no solution, script). So `x, y, z ∉ K[F]`.
The kernel-dimension control at degree `≤ 12` (four-dimensional, spanned by `1, F_2, F_2^2, F_1`)
is a FINITE-EXACT consistency check.

**(8).** `V_j` restricts to each fibre of `G_j` as the derivation `d/dF_j` of its coordinate ring
(the fibre is a smooth affine curve `C`, `V_j` is tangent to it and `V_j F_j = 1`). If `q = V_j h`
with `h ∈ K[x, y, z]`, then on `C` the differential `q dF_j = dh|_C` is exact. On `C = P^1` minus a
finite set, a differential regular on `C` is exact iff its residues at the finite punctures vanish
(the residue at `∞` is then forced, and a rational differential on `P^1` with all residues zero is
exact). The listed residues are exact rational functions of the symbolic roots (script); a nonzero
value shows `q ∉ V_j K[x, y, z]`. The generic ranks are `2γ + r - 1 = 2, 4, 2` for `r = 3, 5, 3`.

**(9).** Over `U = A^3 ∖ {c_3 = 0}` the map is finite étale of degree three with monodromy group
`S_3` (4), so `V_F|_U = O_U ⊕ T` with `T` the rank-two bundle with flat connection attached to the
standard representation, irreducible; its minimal extension is a simple holonomic module distinct
from `K[F] ≅ O`. Near a generic point of the Jelonek surface `S = {c_3 = 0}` with local equation `s`,
the preimage of a small ball consists of one sheet mapping isomorphically and a connected double
cover `W` of the ball minus `S` (the two roots that coalesce at the puncture are swapped by the local
monodromy, (4)); the direct image of `O_W` is `j_+(O ⊕ L_(-1)) = O[1/s] ⊕ j_+ L_(-1)`, with `L_(-1)` the
sign local system. So `V_F` has a subquotient `O[1/s]/O`, a nonzero module supported on `S`. Three
pairwise non-isomorphic constituents give length `≥ 3`.  ∎

## Reading

* **The whole family.** THM-4625 proves statements (1), (5), (6) for every member of Gao's tangent-sweep family with `deg p = d`: `P^1` minus `d+1` and `d+3` points carrying a degree-`(d+1)` function minus its critical points, and a hyperelliptic `(F_1, F_2)`-fibre of genus `⌊(d-2)/2⌋` (here `d = 2`: the conic of (6) is the `t^2 = b^2 - 8aP(w)/w^2` model with `P/w^2 = 2 - w`), with Galois group `S_(d+1)`.
* **A uniform mechanism.** On every coordinate-pair fibre the picture is the same: a rational
  function of degree three on `P^1`, with its critical points deleted. Deleting the critical points
  removes the ramification and leaves an étale map which is not proper; over the critical values the
  fibre loses two of its three points. This is the fibrewise shadow of Gao's tangent sweep: the
  tangency polynomial is the fibre cubic after a Möbius change of variable (3).
* **Why not the plane.** For a planar Keller pair the analogous object is a single curve `{f = c}`
  with `g` a degree-`d` function on it whose critical points are deleted — but the curve must be
  a *plane* affine curve, all of whose fibres are smooth, and THM-4620(5) says the obstruction is
  exactly `b_1` of the generic fibre of `f`. THM-1300's fibres have `b_1 = 2, 4, 2`.
* **Provenance note.** Gao's paper attributes the July-2026 dimension-three counterexample to
  Alpöge, with an infinite family by Gallagher and a geometric explanation by Speyer. The repo's
  provenance ruling for THM-1300 (owner correction of 2026-07-20, MISTAKE-205) is unchanged by this
  file; the external attribution is reported, not adopted.
