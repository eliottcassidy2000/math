---
id: THM-4910
title: "The seven-point sphere holds with the anticomplementary vertex. For an n-simplex V0...Vn (n >= 3) in THM-4742's regime (O, V0, G, O' not coplanar), the seven points A, B, C*, D, E, F, H' are cospherical when (C) is read as the anticomplementary vertex C* = V0 + (n+1)(G - V0) = (n+1)G - nV0, not as V0. They lie on THM-4742's sphere S2, together with an eighth point X, the inverse of the facet centroid G' in the facet circumsphere. S2 is the sphere whose power at each corner P of the skew quadrilateral O, G, G', O' is c_P (R^2 - |OP|^2), with c = n+1, n, n-1, n, and Vieta on the four sides yields the pairs {D, E}, {C*, A}, {H', X}, {F, B}. The literal V0 misses S2 by its power along the median, V0A*V0C* = sum_i |V0Vi|^2 (THM-4742's defect). Reduced mod 2, C* is the facet point {1..n} for odd n, which is the PG(3,2) point 'complementary' to V0 when n = 3"
status: >
  PROVED: statements 1-4 for every n >= 3 in THM-4742's regime. The proof is elementary (Vieta on the sides of the
  skew quadrilateral of centres). Statement 5 is STRUCTURAL: an exact congruence plus an interpretation.
  FINITE-EXACT:
  - symbolic check in THM-4742's Gram basis, all n;
  - 150 exact rational tetrahedra, and 40, 20, 10 exact rational simplices for n = 4, 5, 6: all eight labels on S2,
    pairwise distinct, corner powers, defect.
  NEGATIVE, generic: statement 6 at n = 4, 5, 6, by an exact rank test on random rational simplices (audit J1). The
  flag system is consistent, with a real sphere, on a real hypersurface of nondegenerate 4-simplices (audit J1,
  j1_flag_locus.py), so this is not a universal negative.
  KNOWN, credited:
  - the false verdict for the literal reading C = V0, the sphere S1 and the exceptional locus (session23b,
    05-knowledge/results/tetrahedral_inversion_spheres_20261021/);
  - the six-label sphere S2 through A, B, D, E, F, H', the pencil, the defect Sigma, the all-dimension form and the
    regime (THM-4742, opus-2026-10-21-S49);
  - power of a point, Vieta, the anticomplementary simplex and the Euler line (classical).
  New here:
  - C* in S2 (the pasted seven-point claim is true with the anticomplementary reading);
  - the eighth point X;
  - the corner powers at G, G', O' and the side pairs on GG' and G'O' (the power (n+1)R^2 at O, the iota-pairs {D,E},
    {F,B} and the circle {D,E,F,B} are THM-4742's);
  - the defect read as V0A*V0C*;
  - the mod-2 remark;
  - the negative flag test.
  AUDIT: audit J1 (05-knowledge/results/sevensphere_audit_J1_20261107.md), PASS WITH CORRECTIONS, all applied
  (MISTAKE-732). Statements 1-4 were re-derived symbolically from the corner powers alone, with exact checks for
  n = 3..7.
session: mac-mini-2026-11-07-sevensphere (owner directive 2026-11-07, re-sent from 2026-10-21: tournament structure vs open questions; the Fano plane vs the pasted seven-point sphere; merge diverse ideas)
source: 05-knowledge/results/sevensphere_anticomplementary_signature_simplex_20261107.md
depends_on: [THM-4742, session23b tetrahedral_inversion_spheres_20261021]
verification:
  - 04-computation/experiments/sevensphere_anticomplementary_20261107/eight_point_sphere_exact.py (+ .out): (P1) symbolic all n; (P2) n = 3 exact; (P3) n = 4, 5, 6 exact
  - 04-computation/experiments/sevensphere_anticomplementary_20261107/flag_negative.py (+ .out): statement 6
---

# THM-4910 — the seven-point sphere holds with the anticomplementary vertex

## Setting

This follows THM-4742's notation.
- An `n`-simplex `V0V1…Vn` (`n ≥ 3`) has circumsphere `S` (centre `O`, radius `R`) and centroid `G`.
- Its facet `V1…Vn` has circumcentre `O'`, centroid `G'` and circumradius `r'`, with `r'² = R² − |OO'|²`.
- `G` lies on the median `V0G'`: `G − V0 = (n/(n+1))(G' − V0)`.
- **Regime** (THM-4742): `O, V0, G, O'` are not coplanar, and the lengths used below are nonzero. Equivalently, the
  four centres `O, G, G', O'` are not coplanar; they form a *skew quadrilateral* `O → G → G' → O' → O`.

**The labels** (the pasted list, with (C) corrected):

| label | definition | type |
|---|---|---|
| `A` | second point of the line `V0G` on `S` | chord |
| `B` | pole of the facet hyperplane, `B = O + R²(O' − O)/|O' − O|²` | inversion of `O'` in `S` |
| `C*` | `V0 + (n+1)(G − V0) = (n+1)G − nV0`, the vertex of the anticomplementary simplex | `×(n+1)` from `V0` |
| `D` | `O + (n+1)(G − O)` | `×(n+1)` from `O` |
| `E` | inverse of `G` in `S` | inversion of `G` in `S` |
| `F` | `O + (n+1)(O' − O)` | `×(n+1)` from `O` |
| `H'` | `O' + n(G' − O')`, the facet's own `D`-point (for `n = 3`, the face orthocentre) | `×n` from `O'` |
| `X` | `O' + r'²(G' − O')/|G' − O'|²`, the inverse of `G'` in the facet circumsphere | the facet's own `E`-point |

- The anticomplementary simplex has vertices `Wi = (n+1)G − nVi`. Its facet opposite `Wi` has centroid `Vi`. It is
  the image of the simplex under the homothety of centre `G` and ratio `−n`.
- For `n = 3`, "four times" describes three labels alike: `D` from `O` through `G`, `F` from `O` through `O'`, and
  `C*` from `V0` through `G`.

## Statement

1. **Corner powers.** Let `S2` be the sphere in the 3-flat `O + span(V0 − O, G − O, O' − O)` whose power at the four
   corners is

       Pow_S2(P) = c_P · (R² − |OP|²),    c_O = n+1,  c_G = n,  c_G' = n−1,  c_O' = n.

   It exists and is unique, and it is THM-4742's sphere `S2`.
   - Note `R² − |OP|² = −Pow_S(P)`.
   - Writing the corners as `(α, β) ∈ {0,1}²` (`α` = centroid instead of circumcentre, `β` = facet instead of
     simplex), `c = n + 1 − α − β`.
2. **Eight points.** `S2` meets the four side lines of the quadrilateral exactly in:
   - `OG`: `D` and `E`, the simplex's `D`- and `E`-points;
   - `OO'`: `F` and `B`;
   - `GG'` (the median through `V0`): `C*` and `A`;
   - `O'G'` (the facet's Euler line): `H'` and `X`, the facet's own `D`- and `E`-points.
3. **The pasted seven-point sphere holds with `C = C*`.** `A, B, C*, D, E, F, H'` lie on one sphere, for every simplex
   in the regime and every `n ≥ 3`.
   - In the exact checks all eight labels are pairwise distinct. This is generic only. Labels on different sides never
     coincide in the regime, but each side pair can merge into a tangency point inside the regime (audit J1):
     - `D = E` and `F = B` on session23b's exceptional family;
     - `C* = A` for `V0 = (2, 2/3, 4/3)` over the facet `(0,0,0), (4,0,0), (0,2,0)`;
     - `H' = X` for an algebraic example.
   - Relative to `O`, `C* = s − (n−1)v` mirrors THM-4742's `H' = s − (n−1)w`, swapping the vertex `v = V0 − O` and
     the facet circumcentre `w = O' − O`.
4. **The literal reading.** `V0` lies on the median outside `S2`, and

       Pow_S2(V0) = V0A · V0C* = Σ_(i≥1) |V0Vi|².

   This is THM-4742's defect `Σ`, now read as the product of the distances from `V0` to the median's two points on
   `S2`. So `C = V0` fails because `V0` is the wrong point of the median; the other point `C*` works.
5. **Mod 2.**
   - In barycentric coordinates, `C* = −(n−1)V0 + Σ_(i≥1) Vi`.
   - For odd `n` this is congruent mod 2 to the facet point `Σ_(i≥1) Vi`, whose support is the complement of `{0}`.
   - For `n = 3` it is the PG(3,2) point `{1,2,3}`, where PG(3,2) is the nonempty subsets of the vertex set.
   - So over GF(2), `C*` is the face's own point `{1,2,3}`, the complement of `V0`'s point `{0}`.
   - In the same subset model, the point complementary to `{1,2,3}` is `{0} = V0`. So the GF(2) reduction does *not*
     select `C*` as "the vertex complementary to V1V2V3". Reading (C) as the anticomplementary vertex remains an
     interpretation, justified by statement 3, not by this congruence. (STRUCTURAL.)
6. **No flag extension, generically** (NEGATIVE, generic). The power rule does not continue down a full flag of faces
   of a generic simplex.
   - Take faces `F_n ⊃ F_(n−1) ⊃ … ⊃ F_2` with circumcentres `O_k` and centroids `G_k`.
   - For a generic simplex there is no sphere with `Pow(O_k) = (k+1)(R² − |OO_k|²)` and `Pow(G_k) = k(R² − |OG_k|²)` for
     all `k`: the `2(n−1)` conditions on `n+1` unknowns are inconsistent (exact rank test at `n = 4, 5, 6`, audit J1;
     floating-point residuals 0.16, 0.19, 0.12 here).
   - The exceptions form a proper algebraic subset. For `n = 4` it is a real hypersurface of nondegenerate simplices
     that carry a real flag sphere (audit J1). Whether `n = 5, 6` have exceptions is open.

## Proof

**Powers on a line.** Restrict `Pow_Σ` of any sphere `Σ` to a unit-speed line `ℓ(t)`. It is a monic quadratic
`t² + βt + γ`. So the two points `t1, t2` of `Σ ∩ ℓ` satisfy

    t1·t2 = Pow_Σ(ℓ(0)),    (t1 − L)(t2 − L) = Pow_Σ(ℓ(L)).

The pair is therefore determined by the powers at any two points of the line (Vieta).

**1. Existence and uniqueness.** A power function is `|x|² − 2u·x + κ`. In the regime the four corners are affinely
independent in the 3-flat, so the four prescribed values determine `(u, κ)` uniquely. The sphere is real because it
contains `D` (below). Statement 2 shows it contains `A, B, D, E, F, H'`. By THM-4742(1) (when `D ≠ E`), or in general
by comparing coefficients in the Gram basis ((P1) of `eight_point_sphere_exact.py`), the two coincide.

**2. The four sides.** Put `g = |OG|`, `c = |OO'|`, `γ = |O'G'|` and `m = |V0G'|`.
- **Side `OG`**, measured from `O`:
  - `t1 t2 = (n+1)R²` and `(t1 − g)(t2 − g) = n(R² − g²)`;
  - so `t1 + t2 = R²/g + (n+1)g`, giving `{t1, t2} = {(n+1)g, R²/g}`, i.e. `{D, E}`.
- **Side `OO'`**, measured from `O`:
  - `t1 t2 = (n+1)R²` and `(t1 − c)(t2 − c) = n(R² − c²)`;
  - so `{(n+1)c, R²/c}`, i.e. `{F, B}`.
- **Side `G'O'`**, measured from `O'` towards `G'`:
  - `t1 t2 = n r'²` and `(t1 − γ)(t2 − γ) = (n−1)(R² − |OG'|²) = (n−1)(r'² − γ²)`. This uses `|OG'|² = c² + γ²`,
    since `OO'` is perpendicular to the facet;
  - so `t1 + t2 = (r'² + nγ²)/γ`, giving `{nγ, r'²/γ}`, i.e. `{H', X}`.
- **Side `GG'`**, the median, measured from `V0`:
  - `G` is at `nm/(n+1)` and `G'` at `m`. `A` is at some `a > m`, and `C*` is at `nm`, since
    `C* − V0 = (n+1)(G − V0) = n(G' − V0)`.
  - On the chord `V0A` of `S`, `(nm/(n+1))(a − nm/(n+1)) = R² − g²` and `m(a − m) = R² − |OG'|²`.
  - Hence `(nm − nm/(n+1))(a − nm/(n+1)) = n(R² − g²) = Pow_S2(G)`, and
    `(nm − m)(a − m) = (n−1)(R² − |OG'|²) = Pow_S2(G')`.
  - So `{C*, A}` is the Vieta pair of this side.

All eight points lie on `S2`.

**3.** This is the special case of 2.

**4.** Both points of `S2` on the median lie beyond `G` as seen from `V0`, so `Pow_S2(V0) = a · nm`. On the other side,
`Σ_(i≥1)|V0Vi|² = n m² + Σ_i |Vi − G'|²` (parallel axes). Also `Σ_i |Vi − G'|² = n(r'² − γ²) = n(R² − |OG'|²)`, and
the chord gives `R² − |OG'|² = m(a − m)`. So `Σ = n(m² + m(a − m)) = n m a`.

**5.** Direct, since `(n+1)G = Σ_i Vi`.

**6.** The script `flag_negative.py` fits `Pow − |x|²` by an affine function on the `2(n−1)` flag corners (200 random
simplices per `n`). For `n ≥ 4` the system is overdetermined, and it is inconsistent for generic simplices. Audit J1's
exact rank test (`audit_J1/j1_flag_exact.py`) confirms this. Its locus search (`j1_flag_locus.py`) finds the real
exceptional hypersurface at `n = 4`. ∎

## Remarks

- **What changes against the October verdict.** Session23b and THM-4742 read (C) literally as `V0` and proved the
  pasted claim false. Over ordered fields the seventh incidence then needs `Σ = 0`, which never happens. With the
  anticomplementary reading there is no obstruction: statement 2 is an identity, valid over any field in THM-4742's
  regime (with `|G' − O'|² ≠ 0` for `X`).
- **The Fano plane, complementary against anticomplementary** (ANALOGY).
  - The real picture of the Fano plane is a triangle with its centroid and the vertices of its *complementary*
    (medial) triangle, the homothety of ratio `−1/2` about `G`. Its seventh line, through the three midpoints, is
    drawn as the circumcircle of the complementary triangle, i.e. the nine-point circle. It is straight only in
    characteristic 2 (THM-4742's Fano comparison).
  - The seven-point sphere's seventh point is a vertex of the *anticomplementary* simplex, the homothety of ratio `−n`
    about `G`.
  - Both incidences use the homothety about the centroid. In characteristic 2, where Fano's line is straight, the
    anticomplementary vertex reduces to the complementary point for odd `n` (statement 5). For the triangle itself
    (`n = 2`) it reduces to the centroid point `V0 + V1 + V2`.
- **The quadrilateral as a Boolean square** (STRUCTURAL).
  - The corners are `(α, β) ∈ {0,1}²`, and each side pairs a homothety point (`D, F, C*, H'`) with an inversion or
    chord point (`E, B, A, X`).
  - The four planes of consecutive sides cut `S2` in four circles: `{D,E,F,B}` at `O`, `{D,E,C*,A}` at `G`,
    `{C*,A,H',X}` at `G'` and `{H',X,F,B}` at `O'`. Each point lies on two of them.
  - The 'top' quadruple `{D,F,C*,H'}` is never coplanar: its determinant in THM-4742's basis is `n(n−1) ≠ 0`. So there
    is no Miquel cube. The 'bottom' quadruple `{E,B,A,X}` is generically non-coplanar.
- **The corner inversions.** Inversion in `S2`'s orthogonal sphere about each corner `P` (radius `√Pow_S2(P)`)
  preserves `S2` and swaps the pairs on the two sides through `P`. For `O` this is THM-4742's `ι`.
