---
id: THM-4718
title: "The generic Apollonius Galois group for labelled input circles with radii in the base field is the parity cube (Z/2)^3"
status: >
  PROVED over K=Q(x1,y1,x2,y2,x3,y3,r1,r2,r3), with algebraically independent coordinates and labelled input radii in K:
  the generic splitting field of the eight solution circles has Galois group (Z/2)^3 in its parity-cube action.
  Direct elimination proves the discriminant formula; valuations at six distinct irreducible defect divisors give squareclass rank three.
  The generic rank gap in the original 2026-10-19 account is closed by the 2026-10-20 addendum and independently audited result below.
  Exact uniform-group fixed-point weights are 1/8,6/8,1/8 for 8,4,0 points. Finite-field census values remain finite evidence;
  their specialization-frequency limits and error terms are not proved by that census. Original text is retained as historical provenance.
session: mac-mini-2026-10-19-apollonius
source: 05-knowledge/results/apollonius_enumerative_20261019_galois_parity_cube_ford_descartes_windows.md
scripts:
  - 04-computation/experiments/apollonius_20261019/apollonius_finite_field.py (+ .out, _large.out), apollonius_rational_test.py (+ .out), apollonius_lie_gram.py (+ .out), apollonius_galois.py, apollonius_normalized.py
related:
  - THM-4653 (Legendre bilinear forms over F_q), THM-4701 (finite unit-circle graphs), THM-4714 (Ford circles via Minkowski ?), HYP-9338
  - Erdos-1057 Carmichael paper (Howe-Huffman-Kohlmeyer-Kruer 2026): the same Chebotarev/character style of counting
---

## Addendum, 2026-10-20: generic squareclass rank is proved

For labelled input circles over `K=Q(x1,y1,x2,y2,x3,y3,r1,r2,r3)` with all nine variables algebraically independent, the generic splitting field has group exactly `(Z/2)^3`, acting by even swaps within four pairs. [The direct elimination, distinct-divisor squareclass proof, graph-code extension and independent audit](../../05-knowledge/results/apollonius_squareclass_code_20261019/APOLLONIUS_SQUARECLASS_CODE.md) close the original generic-rank gap without a census or an unproved Lie-coordinate normalization. The six defect polynomials are pairwise nonassociate absolutely irreducible quadrics; their valuation coordinates are independent, and the four discriminant incidence rows have rank three with their total sum as the only relation.

The uniform finite-group weights for 0,4,8 fixed circles are exactly `1/8,6/8,1/8`, with mean four. Applying those weights to finite-field specializations still requires a separate distribution theorem and explicit degeneracy control. In particular, the historical census paragraph below does not prove convergence or an `O(1/q)` error. Giving only squared radii, forgetting input labels, or specializing the parameters can change the field and group. No historical priority is claimed.

### Status at first publication, retained for provenance

> PROVED: the parity relation (the Gram determinant of three points of the Lie quadric has zero diagonal, so it is twice the
> product of the three pairwise Lie products, each a function of e_i e_j only) and its consequence for the Galois group, given
> the CITED Lie-sphere description of the problem (the two solutions of a sign class are the two points in which the Lie
> quadric meets the line Lie-orthogonal to the three oriented data circles; Knight, J. Geom. 2005; Cecil, Lie Sphere Geometry;
> Pedoe) and the elementary lemma disc(Q|W) disc(Q|W^perp) = disc(Q) mod squares. FINITE-EXACT: Delta_e Gram_e is a rational
> square for every class at 300 random integer triples, and the product of the four discriminants is a square at 400 more
> (apollonius_lie_gram.out, apollonius_rational_test.out); the finite-field census for q = 11..503 (apollonius_finite_field*.out).
> The identification of the full group (Z/2)^3 rather than a proper subgroup is read off the census (all three fixed-point
> types occur with the predicted weights) and is not separately proved. Literature: McKean, "Circles of Apollonius two ways"
> (arXiv:2210.13288) treats the problem as three quadrics with quadratically enriched counts and records Pedoe's real-count
> theorem; the Sottile-Yahl survey (arXiv:2108.07905) was not checked for this group; no priority claim. Pre-audit.

## Original 2026-10-19 text below: the addendum controls the generic proof and distribution scope


# THM-4718 — the Galois group of the problem of Apollonius is the parity cube

## Setting

* Three circles `C_i` with centres `c_i = (a_i, b_i)`, radii `r_i > 0`, `d_ij = |c_i − c_j|`, in general position (centres not collinear, no two tangent, no circle through the solutions' degeneracies). A circle `(x, y, r)` is tangent to `C_i` with sign `e_i ∈ {±1}` iff `|(x, y) − c_i|² = (r + e_i r_i)²`.
* The global flip `(e, r) ↦ (−e, −r)` preserves the circle, so the sign vectors fall into four classes `e ∈ {(+++), (++−), (+−+), (−++)}`. For a class, subtracting the equations gives two linear equations in `(x, y, r)`; substituting into one equation gives a quadratic in `r` with discriminant `Δ_e`. Generically each class has two solutions, eight in all (Bézout: three quadrics in `P³`).
* **Lie coordinates.** An oriented circle `(c, ρ)` (`ρ = e r` signed) has the Lie product `⟨(c, ρ), (c', ρ')⟩ = (|c − c'|² − (ρ − ρ')²)/2`; it vanishes iff the oriented circles are tangent, and `⟨v, v⟩ = 0`. The oriented circles are the points of the Lie quadric `Q ⊂ P⁴` (signature (3,2)); a circle oriented-tangent to the three data circles `C_i^(e_i)` is a point of `Q` on the projective line `ℓ_e = P(W_e^⊥)`, `W_e = span(C_i^(e_i))`.

## Statements

1. **Discriminant = oriented tangency defects (PROVED, modulo the cited Lie description).** For each class,

       Δ_e ≡ G_e := det( ⟨C_i^(e_i), C_j^(e_j)⟩ )_(i,j) = 2 ∏_(i<j) ( d_ij² − (e_i r_i − e_j r_j)² ) / 8     (mod squares),

   i.e. the two circles of class `e` are rational over a field `k` iff the product of the three pairwise oriented tangency defects `d_ij² − (e_i r_i − e_j r_j)²` is a square in `k` (up to the fixed constant). Numerically `Δ_e · G_e` is an exact rational square at every one of 300 random integer triples.

2. **The parity relation (PROVED).** Each factor `d_ij² − (e_i r_i − e_j r_j)² = d_ij² − r_i² − r_j² + 2 e_i e_j r_i r_j` depends only on `e_i e_j`. Over the four classes, each product `e_i e_j` takes each sign exactly twice, so

       Δ_1 Δ_2 Δ_3 Δ_4 ≡ ∏_(i<j) ( defect_ij^+ · defect_ij^- )²  ≡ 1     (mod squares).

   Hence the Galois (= monodromy, Harris) group of the eight solutions, which acts by swapping within the four conjugate pairs, lies in the even subgroup `{(s_1, …, s_4) ∈ {±1}^4 : ∏ s_i = 1} ≅ (Z/2)³` of the four independent swaps.

3. **Chebotarev census (FINITE-EXACT; the equality with (Z/2)³ read from it).** Over `F_q` a generic Apollonius problem has `2·#{e : Δ_e ∈ F_q^×²}` rational solution circles. If the group is the full even subgroup, the Frobenius is the identity with probability 1/8 (8 solutions), swaps exactly two pairs with probability 6/8 (4 solutions) or all four with probability 1/8 (0 solutions); the counts 2 and 6 never occur generically and the mean is 4. Random generic triples (3000 per `q`):

   | `q` | 0 | 2 | 4 | 6 | 8 | odd (degenerate) | mean |
   |---|---|---|---|---|---|---|---|
   | 31 | 283 | 138 | 1908 | 106 | 202 | 363 | 3.84 |
   | 101 | 368 | 39 | 2089 | 51 | 336 | 117 | 3.96 |
   | 211 | 347 | 18 | 2233 | 24 | 325 | 53 | 3.97 |
   | 503 | 387 | 8 | 2236 | 11 | 323 | 35 | 3.92 |

   The weights converge to `1/8, 3/4, 1/8` and the forbidden counts vanish like `O(1/q)` (degenerate linear or double-root classes), so no further relation among the `Δ_e` exists: the group is the parity cube.

4. **Real circles and degenerations (PROVED from 1; Pedoe's count CITED).** Over `R` the number of solution circles of a generic problem is `2·#{e : Δ_e > 0}`, i.e. twice the number of classes whose three oriented tangency defects have positive product; the odd counts of Pedoe's theorem ("anything but seven") arise only from degenerate classes. If the three circles are mutually externally tangent, the `(+++)` defects vanish: the class degenerates (Soddy circles).

## Proofs

**1.** In Lie coordinates the three conditions are `⟨v, C_i^(e_i)⟩ = 0`, `v ∈ Q`, so the oriented solutions of class `e` are `Q ∩ P(W_e^⊥)`, two points of a line, with discriminant `disc(Q|_(W_e^⊥))`. For a nondegenerate quadratic form, `disc(Q|_W) · disc(Q|_(W^⊥)) ≡ disc(Q)` mod squares (choose a basis adapted to `W ⊕ W^⊥`), and `disc(Q|_W) = det(⟨C_i^(e_i), C_j^(e_j)⟩) = G_e`. The diagonal of the Gram matrix vanishes (`⟨v, v⟩ = 0` on `Q`), and a symmetric 3×3 matrix with zero diagonal has determinant `2 m_12 m_13 m_23`. The elimination discriminant `Δ_e` of the setting is the same quantity in affine coordinates, up to the square of the elimination denominator; the exact numerical identity `Δ_e G_e ∈ Q²` confirms the constant. ∎

**2.** Immediate from `defect_ij(e) = d_ij² − r_i² − r_j² + 2 e_i e_j r_i r_j` and the class table: `e_1 e_2` is `+, +, −, −`, `e_1 e_3` is `+, −, +, −`, `e_2 e_3` is `+, −, −, +` over the classes `(+++), (++−), (+−+), (−++)`. The Galois group of the splitting field of the eight solutions acts on each class by the Galois group of `k(√Δ_e)/k`, and the relation `√(Δ_1 Δ_2 Δ_3 Δ_4) ∈ k` forces the product of the four signs to be `+1`. ∎

**3.** `apollonius_finite_field.py`: for each sign class solve the two linear equations for `x, y` as affine functions of `r`, substitute, and count the roots of the quadratic in `F_q` (via the square table), identifying `r` and `−r`; configurations with a singular linear system are skipped as degenerate. ∎

**4.** Over `R`, a class contributes two circles iff `Δ_e > 0`. The oriented tangency defect `d² − (r_i − r_j)²` is positive for disjoint or crossing circles and negative when one circle lies inside the other, so the real count is read off the three pairwise positions. ∎

## Remarks

* **Enumerative reading.** The classical count 8 = 2³ is the Bézout number; the "arithmetic" count over a finite field is 4 on average (Chebotarev), not 8 and not 1: the Galois group is neither the symmetric group `S_8` (average 1) nor the regular `(Z/2)³` (average 1), but the parity cube acting with four orbits of size 2. This is the discriminant-level refinement of McKean's quadratically enriched count.
* **Ford circles.** Three mutually tangent "circles" with one of them a line (Ford circles at Farey neighbours `a/b, c/d` and the real line, curvatures `2b², 2d², 0`) are the Soddy degeneration of this problem: Descartes' rule gives the two solutions `2(b ± d)²`, the Ford circles at the mediant `(a+c)/(b+d)` and the parent `(a−c)/(b−d)` (checked for all Farey neighbours with `q < 40`); the Stern–Brocot moves are the Apollonian reflections, and THM-4714's Minkowski `?` sends the tangency points to dyadic rationals.
