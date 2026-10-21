---
id: THM-4742
title: "The seven-point sphere claim of a tetrahedron splits into two inversion-invariant spheres through a common circle: the six labels A, B, D, E, F, H' always lie on one sphere S2 (new; the concurrent session23b packet has the other sphere S1 through B, C, D, E, F), S1 and S2 meet in the circle through B, D, E, F in the plane O G O', and the vertex C misses S2 by exactly Sigma = |V0V1|^2 + |V0V2|^2 + |V0V3|^2 = -8 (V0 - O).(G - V0). Hence (generic regime) the seven labels are cospherical iff Sigma = 0 iff the median V0G is tangent to the circumsphere iff A = C: never over an ordered field, possible over non-real fields (Q(i), F_p) but only with A = C. This is the precise form of the Fano analogy: six incidences are identities and the seventh is controlled by one defect that is positive over R (Fano's seventh line: determinant -2)"
status: >
  PROVED: statements 1-4 (exact identities in the Gram basis v, s, w; symbolic check in seven_point_pencil.py (P1)).
  FINITE-EXACT: (P2) 150 random rational tetrahedra in actual coordinates; (P3) every valid random configuration over F_p,
  5 <= p <= 97 (no mismatch); (P4) an explicit tetrahedron over Q(i); (P5) the analogues in dimensions 2, 4, 5.
  CLASSICAL: Fano's seventh line (determinant -2, characteristic 2 only).
  CREDIT: the false-claim verdict, the four-point obstruction A,C,D,E, the real exceptional locus (D = E, B = F, five
  distinct points), the sphere S1 through B,C,D,E,F with its inversion invariance and the values Q1(H') = 6 - 2 v.s,
  Q1(A) = mu(16 - T) are the concurrent session23b packet 05-knowledge/results/tetrahedral_inversion_spheres_20261021/
  (commit 7798cbbecb). New here: S2, the pencil, the defect identity, the tangency/coincidence reading, the field dichotomy,
  the dimension check, the Fano comparison. Session opus-2026-10-21-S49; pre-audit.
session: opus-2026-10-21-S49 (owner: tournament structure vs open questions; the Fano plane vs a pasted 'seven-point sphere'; Gluck's paper)
source: 05-knowledge/results/seven_point_pencil_fano_and_converse_gluck_20261021.md
scripts:
  - 04-computation/experiments/seven_point_spheres_20261021/seven_point_pencil.py (+ .out)
related:
  - 05-knowledge/results/tetrahedral_inversion_spheres_20261021/README.md (concurrent: S1, exceptional locus, counterexample)
  - 05-knowledge/results/fano_tournament_move_quotients_20261021/session23b-fano/README.md (concurrent: tetrahedral Fano incidence, real non-representability of the Fano plane)
  - THM-4610 / MISTAKE-649 (the 3x - 1 barrier: an algebraic structure that cannot see the real sign)
---

# THM-4742 — two spheres through a common circle, and the seventh incidence

## Setting

* A tetrahedron `V0V1V2V3` has circumcentre `O` and circumradius `R`. Its centroid is `G`. The face `V1V2V3` has circumcentre `O'`, centroid `G'` and orthocentre `H'`.
* The seven labels of the pasted claim are:
  * `A`, the second point of the line `V0G` on the circumsphere;
  * `B`, the pole of the face plane, i.e. the common point of the tangent planes at `V1, V2, V3`;
  * `C = V0`;
  * `D = O + 4(G − O)`;
  * `E`, the inverse of `G` in the circumsphere;
  * `F = O + 4(O' − O)`;
  * `H'`.
* Normalise `O = 0`, `R = 1`. Put `v = V0`, `s = V1 + V2 + V3`, `w = O'`, `k = |w|²`, `t = v + s = 4G` and `T = |t|²`.
  * Every vertex lies on the circumsphere and the face lies in the plane `w·x = k`, so `s·w = 3k`.
  * In these terms `C = v`, `D = t`, `E = 4t/T`, `F = 4w`, `B = w/k`, `H' = s − 2w` and `A = v + μ(s − 3v)`, where `μ = −2 v·(s − 3v)/|s − 3v|²`.
* **Generic regime:**
  * `δ = v·(s × w) ≠ 0`, i.e. `O, V0, G, O'` are not coplanar;
  * `G ≠ O` and `O' ≠ O`;
  * for uniqueness, `D ≠ E` (i.e. `|OG| ≠ R/2`).
* `ι` denotes the inversion `x ↦ 4x/|x|²`, the inversion in the sphere of radius `2R` about `O`.

## Statements

1. **Two spheres.** Define spheres `|x|² − 2u·x + 4 = 0` by the values of `(2u·v, 2u·s, 2u·w)`:

       S1: (5, T − 1, 4k + 1),        S2: (2 v·s − 1, |s|² + 6, 4k + 1).

   Since `(v, s, w)` is a basis, each sphere is well defined.
   * `S1` contains `B, C, D, E, F` (session23b).
   * `S2` contains `A, B, D, E, F, H'`.
   * Both have power `4R²` at `O`, so both are `ι`-invariant. `ι` swaps `D ↔ E` and `B ↔ F`, and sends `C ↦ 4C`, `A ↦ 4A`, `H' ↦ ι(H')`. Hence `S1` also contains `4C`, and `S2` also contains `4A` and `ι(H')`.
   * With `D ≠ E`, any sphere through `D` and `E` has power `OD·OE = 4R²` at `O`. So `S1` and `S2` are the unique spheres through `{B, C, D, E, F}` and `{A, B, D, E, F, H'}` respectively.
2. **The pencil.**
   * `2(u1 − u2)` has components `(Σ, −Σ, 0)` on `(v, s, w)`, where `Σ = 6 − 2 v·s`.
   * Hence `u1 − u2` is orthogonal to `t` and to `w`, and the radical plane of `S1` and `S2` is the plane `O G O'`.
   * `S1 ∩ S2` is the circle `Γ` through `B, D, E, F`, which are always coplanar.
   * `S1` and `S2` are the members of the pencil of `ι`-invariant spheres through `Γ` that pass through `C` and through `H'` respectively.
3. **The defect.**
   * In actual units,

         Pow_S2(C) = Pow_S1(H') = Σ = |V0V1|² + |V0V2|² + |V0V3|² = −8 (V0 − O)·(G − V0).

   * Also `Pow_S1(A) = μ(16 − T)` (session23b), with `μ = Σ/|s − 3v|²`.
   * So `Σ = 0` if and only if the median line `V0G` is tangent to the circumsphere at `V0`, if and only if `A = C`.
4. **The seventh incidence over fields.** Work over a field of characteristic `≠ 2, 3`, in the generic regime with `D ≠ E`. The seven labels lie on one sphere iff `S1 = S2`, iff `Σ = 0`, iff `A = C`.
   * Over any ordered field, `Σ` is a sum of squares of edge coordinates and is positive for a nondegenerate tetrahedron, so the seventh incidence never holds.
   * Over `Q(i)` it does hold, for `V0 = 0`, `V1 = e1`, `V2 = e2`, `V3 = (0, i, i)`. Over `F_p`, `5 ≤ p ≤ 97`, it holds in every sampled configuration with `Σ = 0`, and only there.
   * Whenever it holds, `A = C`, so seven *distinct* cospherical labels never occur in the generic regime.
   * The real exceptional locus of session23b has `D = E` and `B = F`, and lies outside this regime.
5. **Three dimensions only (FINITE-EXACT, random exact trials).**
   * For simplices of dimension 4 and 5, the analogues (`D = O + (n+1)(G − O)`, `F = O + (n+1)(O' − O)`, `H' = O' + n(G' − O')`) have no cospherical set of more than `n + 1` of the seven labels.
   * In the plane, `A, B, C, E` and `B, D, E, F` are concyclic.
   * `D` and `E` do not lie on the circumsphere of the tangential tetrahedron (30 rational trials, scratch; recorded only as a negative).

## Proofs

**1–2.**
* Write `Q(x) = |x|² − 2u·x + 4`. Every point `X = αv + βs + γw` has `|X|²` given by the Gram matrix

      [[1, v·s, v·w], [v·s, |s|², 3k], [v·w, 3k, k]],

  and `2u·X = α(2u·v) + β(2u·s) + γ(2u·w)`.
* **`S2`.**
  * `Q(D) = T − (2v·s − 1) − (|s|² + 6) + 4 = 0`, using `T = 1 + 2v·s + |s|²`.
  * `Q(F) = 16k − 4(4k + 1) + 4 = 0`.
  * `Q(H') = |s|² − 8k − (|s|² + 6) + 2(4k + 1) + 4 = 0`, using `|s − 2w|² = |s|² − 12k + 4k`.
  * A constant term of 4 gives `Q(4x/|x|²) = 4Q(x)/|x|²`. So `E = ι(D)` and `B = ι(F)` lie on `S2`.
  * For `A = v + μ(s − 3v)`:

        Q(A) = 6 − 2 v·s − μ(T + 8 − 8 v·s) = 6 − 2 v·s − μ|s − 3v|² = 0,

    since `−2 v·(s − 3v) = 6 − 2 v·s` and `T + 8 − 8 v·s = |s|² − 6 v·s + 9 = |s − 3v|²`.
* **`S1`** is session23b's sphere. The same bookkeeping gives `Q1(B…F) = 0`.
* **Statement 2.** Subtract the two coefficient triples: `(5 − (2v·s − 1), (2v·s + |s|²) − (|s|² + 6), 0) = (Σ, −Σ, 0)`. On `t = v + s` this gives `Σ − Σ = 0`, and on `w` it gives 0.
* All identities are re-derived symbolically in (P1) and checked exactly in (P2). ∎

**3.**
* Since `|Vi| = 1`, `|V0 − Vi|² = 2 − 2 v·Vi`. Summing over `i = 1, 2, 3` gives `Σ = 6 − 2 v·s`.
* `Q2(C) = 1 − (2v·s − 1) + 4 = Σ`, and `Q1(H') = Σ` likewise.
* `(V0 − O)·(G − V0) = v·(t/4 − v) = (1 + v·s)/4 − 1 = −Σ/8`.
* Scaling back multiplies powers by `R²`, so in actual units the defect is the sum of squared edge lengths at `V0`.
* `μ = −2 v·(s − 3v)/|s − 3v|² = Σ/|s − 3v|²`, so `μ = 0` (`A = C`) iff `Σ = 0`. ∎

**4.**
* If the seven labels lie on one sphere `S` with `D ≠ E`, then `S` has power `4R²` at `O`.
  * `S ⊇ {B, C, D, E, F}` forces `S = S1`.
  * `S ⊇ {A, B, D, E, F, H'}` forces `S = S2`.
* By 2, `S1 = S2` iff `Σ = 0`. The identities are polynomial, so they hold over any field of characteristic `≠ 2, 3`, where the construction (division by 4 and 3) makes sense.
* Over an ordered field, `Σ` is a sum of squares of coordinates, zero only if `V1 = V2 = V3 = V0`.
* The `Q(i)` and `F_p` cases are exact computations ((P3), (P4)). ∎

## The Fano comparison (what is and is not claimed)

* **The two configurations.**
  * The real picture of the Fano plane (vertices, edge midpoints and centroid of a triangle) realises six of its seven lines as collinearities. The seventh line, through the three midpoints, has determinant `−2` in `0/1` coordinates. It exists exactly where `2 = 0`.
  * The pasted seven-label sphere realises its "lines" as two spheres. One incidence, `C ∈ S2`, has defect `Σ`, a sum of squares. It vanishes for a nondegenerate tetrahedron exactly over non-real fields, and then only through the coincidence `A = C`.
* **What is common.** In both, one positive quantity (`2` or `Σ`) separates real geometry from an algebraic incidence.
* **What is not claimed.**
  * The seven labels are not identified with the seven Fano points.
  * Session23b's tetrahedral Fano incidence (six edges plus a formal point) is a different, combinatorial bridge.
  * No Fano incidence structure on the seven metric labels is asserted.
