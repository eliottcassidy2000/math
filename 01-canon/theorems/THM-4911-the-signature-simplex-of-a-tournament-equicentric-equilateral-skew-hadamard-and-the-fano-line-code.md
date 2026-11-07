---
id: THM-4911
title: "The signature simplex of a tournament. The columns p_i = (I+M)e_i of a tournament's signed matrix are n vertices of the cube {+-1}^n on the sphere of radius sqrt(n), with |p_i - p_j|^2 = 4(1 + delta_ij), where delta_ij counts the vertices beating exactly one of i, j. THM-1950's inverse response s(T) = u^T(I+M^TM)^(-1)u is 1/dist(0, aff{p_i})^2, so the circumradius is r_T^2 = n - 1/s(T). Centroid = circumcentre (equicentric) iff M^2 u is parallel to u: for odd n iff T is regular (Pfaffian kernel argument); a non-strong T is equicentric iff T = R_a => R_b with both parts regular; an almost-regular T (scores n/2-1, n/2) is equicentric iff its two score classes are regular halves joined biregularly, which exists iff n = 2 (mod 4); and for even n, -M^2 u = lambda u forces lambda = n-1 (mod 4). The simplex is equilateral iff T is doubly regular (c = -1) or I+M is skew-Hadamard (c = 0), so the skew-Hadamard conjecture says exactly that every n = 0 (mod 4) carries a tournament whose n signature vectors are pairwise orthogonal. Equivalently, T is equicentric iff every vertex is an endpoint of the same number of directed 2-paths. For the Paley tournament QR7 the signature vectors are the +-1 incidence vectors of the seven Fano lines. THM-4910's eight-point sphere transfers to every tournament and vertex v in THM-4742's regime (never for an equicentric T, and not for some non-equicentric T such as source => C3 => sink), with Pow(O_T) = n(n - 1/s(T)) and Pow(p_v) = 4(n - 1 + #2-paths with an end at v)"
status: >
  PROVED: statements 1-6.
  FINITE-EXACT, all isomorphism classes (nauty gentourng):
  - (1), (4) for n <= 8;
  - (2) in floating point for n <= 8 (audit J1: exact for n <= 7, floating point for n <= 9);
  - (3a) for n <= 9;
  - (3b)-(3d) for n <= 10, with equicentric counts 1, 2, 1, 5, 3, 14, 15, 165 for n = 3..10;
  - (5) exact;
  - (6) exact on 24 random non-equicentric tournaments (n = 5..8).
  KNOWN in substance, credited:
  - doubly regular tournaments correspond to skew-Hadamard matrices (Reid and Brown, 1972);
  - (I+M)^T(I+M) = I + M^TM and the inverse response s with 0 < s <= n, s = n iff regular, and the velocity-addition
    law (THM-1950, klein-2026-07-21-S400);
  - Paley QR7's out-neighbourhoods are the Fano lines (THM-211; session23b fano packet; the arithmetic_braids Fano
    code note);
  - the skew-Sylvester double D (THM-447), which is the skew-Hadamard doubling;
  - even-order tournament Pfaffians are odd (THM-1475(1), klein-2026-07-20-S339; the same count in THM-3729), and the
    kernel vector of an odd-order tournament's skew matrix has all-odd entries
    (07-reflections/border-remembers-doubling-forgets-the-gleason-ring-from-two-paley-tournaments-cbx0611s3.md,
    2026-06-11). These are the ingredients of (3a);
  - the cases c = -1 (SS^T = nI - J) and c = 0 (SS^T = (n-1)I, skew conference, equivalently skew-Hadamard) are
    THM-472's equality cases;
  - the main-angle framework is Nozaki-Suda, "A characterization of skew Hadamard matrices and doubly regular
    tournaments" (arXiv:1202.5374, 2012). In their language, equicentric means at most two main eigenvalues +-sqrt(lambda)
    of sqrt(-1) M. Their Theorem 2.7 (order 4k+2, main eigenvalues +-1, i.e. lambda = 1, plus a spectral condition;
    equivalent to DRT(4k+3)) is the closest prior result to (3c). A brief search found no classification of (3a)-(3d)
    and no statement of the c in {-1, 0} dichotomy, which is likely folklore;
  - the 2-path characterisation in (3) and the Sonneborn-Berger form were supplied by audit J1.
  New here:
  - the signature simplex as a geometric home for s(T);
  - the equicentric classification (3a)-(3d);
  - the c in {-1, 0} dichotomy for equilateral simplices and the skew-Hadamard reformulation;
  - the transfer (6) of the eight-point sphere.
  AUDIT: audit J1 (05-knowledge/results/sevensphere_audit_J1_20261107.md), PASS WITH CORRECTIONS, all applied
  (MISTAKE-732). Every proof was re-checked, with an independent census of all classes n <= 10 and an exact check of
  (6) on 160 vertex choices.
session: mac-mini-2026-11-07-sevensphere (owner directive 2026-11-07: tournament structure vs open questions; Fano vs the seven-point sphere)
source: 05-knowledge/results/sevensphere_anticomplementary_signature_simplex_20261107.md
depends_on: [THM-4910, THM-4742, THM-1950, THM-447, THM-211]
verification:
  - 04-computation/experiments/sevensphere_anticomplementary_20261107/signature_simplex.py (+ .out): statements 1-6 (about 2 minutes, n <= 10)
---

# THM-4911 — the signature simplex of a tournament

## Setting

- `T` is a tournament on `n` vertices. `M = A − Aᵀ` with `A_ij = 1` iff `i → j`, so `M` is skew with `±1` off the
  diagonal. `u` is the all-ones vector, and `x = Mu` is the score-excess vector, `x_i = 2d⁺_i − (n−1)`.
- **Signature vectors.** `p_i = (I + M)e_i`. Its `k`-th entry is `+1` for `k = i` and for `k → i`, and `−1` for
  `i → k`.
- The **signature simplex** of `T` is the simplex with vertices `p_1, …, p_n`. Its Gram matrix is
  `(I+M)ᵀ(I+M) = I + MᵀM`, which is positive definite, so the simplex is nondegenerate (dimension `n − 1`).
- **Notation.**
  - `δ_ij` is the number of vertices `k ∉ {i, j}` that beat exactly one of `i, j`.
  - `s(T) = uᵀ(I + MᵀM)⁻¹u` is THM-1950's inverse response, which is also `F_T(1)` of THM-4900.
  - `G_T` is the centroid of the signature simplex and `O_T` its circumcentre.
  - `T` is *equicentric* if `G_T = O_T`, and *equilateral* if all edges are equal.

## Statement

1. **Cube and sphere.** Every `p_i` is a vertex of `{±1}^n` with `|p_i|² = n`. The Hamming distance between `p_i` and
   `p_j` is `1 + δ_ij`, so `|p_i − p_j|² = 4(1 + δ_ij)`.
2. **The inverse response is a distance.**
   - `dist(0, aff{p_i})² = 1/s(T)`.
   - The circumcentre is `O_T = (I+M)(I+MᵀM)⁻¹u / s(T)`, and the circumradius is `r_T² = n − 1/s(T)`.
   - Since `G_T ∈ aff{p_i}` and `|G_T|² = (n + |x|²)/n²`:

         n²/(n + |x|²) ≤ s(T) ≤ n,

     with equality on the left iff `T` is equicentric. The right bound is THM-1950's, with equality iff `T` is regular.
3. **Equicentric tournaments.** `T` is equicentric iff `M²u = −λu` for some `λ`, and then `λ = |x|²/n`.
   - **The 2-path form** (audit J1). `(M²u)_v = 2Δ_v − (n−1)²`, where `Δ_v` is the number of directed 2-paths with an
     end at `v`, as in 6.
   - So `T` is equicentric iff every vertex is an endpoint of the same number `Δ` of directed 2-paths, and then
     `λ = (n−1)² − 2Δ`.
   - Equivalently, `2SB(v) − (n−2)d⁺(v)` is constant, where `SB(v)` is the Sonneborn–Berger score (the sum of the
     scores of the vertices `v` beats).
   - By 6 this also says that all vertex powers `Pow(p_v)` are equal.
   - **(a) Odd n.** Equicentric iff regular.
   - **(b) Non-strong.** `T` is equicentric iff it has exactly two strong components and both are regular,
     `T = R_a ⇒ R_b`. Then `a, b` are odd, `n` is even and `λ = ab`.
   - **(c) Almost regular** (all scores in `{n/2 − 1, n/2}`). Equicentric iff:
     - the two score classes `S` (score `n/2`) and `S̄` (score `n/2 − 1`) each have `n/2` vertices;
     - both induce regular tournaments, so `n/2` is odd;
     - every vertex of `S` beats `(n/2 + 1)/2` vertices of `S̄`, and every vertex of `S̄` beats `(n/2 − 1)/2` of `S`.

     Such tournaments exist exactly when `n ≡ 2 (mod 4)`: the single arc for `n = 2`, and strong ones for every
     `n ≥ 6`. They have `λ = 1`.
   - **(d) Even n.** `λ ≡ n − 1 (mod 4)`.
   - **Counts** (FINITE-EXACT, isomorphism classes).

     | `n` | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
     |---|---|---|---|---|---|---|---|---|
     | equicentric | 1 | 2 | 1 | 5 | 3 | 14 | 15 | 165 |
     | non-strong (3b) | 0 | 2 | 0 | 3 | 0 | 8 | 0 | 37 |
     | almost regular (3c) | 0 | 0 | 0 | 2 | 0 | 0 | 0 | 88 |
     | values of `λ` | 0 | 3 | 0 | 1, 5, 9 | 0 | 3, 7, 11, 15 | 0 | 1, 5, …, 25 |
4. **Equilateral tournaments.** The simplex is equilateral iff `I + MᵀM = (n − c)I + cJ`, and then `c ∈ {−1, 0}`:
   - `c = −1` iff `T` is doubly regular (`n ≡ 3 mod 4`, with `δ_ij = (n−1)/2`);
   - `c = 0` iff `I + M` is a skew-Hadamard matrix (`n ≡ 0 mod 4`, with `δ_ij = (n−2)/2`). Then `T` is equicentric
     but not regular, and the signature vectors are pairwise orthogonal.

   Consequently the **skew-Hadamard conjecture** (OPEN) is equivalent to either of the following. The first form is the
   definition of a skew-Hadamard matrix; the second is Reid–Brown.
   - for every `n ≡ 0 (mod 4)` some tournament on `n` vertices has pairwise orthogonal signature vectors;
   - for every `n ≡ 3 (mod 4)` some tournament on `n` vertices has a regular signature simplex.
5. **The Fano line code.** For the Paley tournament QR7 (`i → j` iff `j − i ∈ {1, 2, 4}` mod 7), `p_i` is `−1`
   exactly on the out-neighbourhood `i + {1, 2, 4}`, which is a line of the cyclic Fano plane.
   - So the seven signature vectors are the `±1` incidence vectors of the seven Fano lines. Two lines meet in one
     point, so every Hamming distance is 4.
   - The signature simplex of QR7 is a regular 6-simplex in the 7-cube. This is case `c = −1` of 4.
   - More generally, the signature simplex of a doubly regular tournament is the incidence code of its out-neighbourhood
     design, a symmetric 2-`(n, (n−1)/2, (n−3)/4)` design (classical).
6. **The eight-point sphere transfers at vertices in the regime.** Let `n ≥ 4` and `v` be a vertex such that the
   signature simplex is in THM-4742's regime at `v` (`O_T, p_v, G_T, O'_v` not coplanar, where `O'_v` is the
   circumcentre of the facet opposite `p_v`). Then THM-4910 gives an eight-point sphere `S2(T, v)` with
   - `Pow(O_T) = n·r_T² = n(n − 1/s(T))`;
   - `Pow(p_v) = Σ_(j≠v)|p_v − p_j|² = 4(n − 1 + Δ_v)`, where `Δ_v = Σ_j δ_vj` is the number of directed 2-paths
     `j → k → v` plus the number of 2-paths `v → k → j`.

   The regime fails whenever `T` is equicentric. So every regular tournament (including the Fano tournament QR7), every
   `R_a ⇒ R_b` and every skew-Hadamard tournament has no seven-point sphere.
   - It also fails at every vertex of some non-equicentric tournaments, for example source ⇒ C3 ⇒ sink: 2 classes at
     `n = 5` and 6 at `n = 7` (audit J1, `j1_regime_failures.py`).

## Proofs

**1.**
- The coordinates `k ∉ {i, j}` where `p_i` and `p_j` differ are exactly those where `M_ki ≠ M_kj`, i.e. where `k`
  beats exactly one of `i, j`. That gives `δ_ij` coordinates.
- At the coordinates `i` and `j` themselves, `(p_i)_i = 1`, `(p_j)_i = M_ij`, `(p_i)_j = M_ji` and `(p_j)_j = 1`.
  Since `M_ij = −M_ji`, exactly one of these two coordinates differs.

**2.**
- All `p_i` lie on the sphere of radius `√n` about `0`, so the circumcentre is the orthogonal projection of `0` onto
  the affine hull.
- Minimising `|Pλ|²` subject to `uᵀλ = 1`, with `PᵀP = I + MᵀM`, gives `λ = (PᵀP)⁻¹u/s` and minimum `1/s`.
- Then `r_T² = n − |O_T|²`.
- For the bounds, `G_T` lies in the hull and `|G_T|² = |(I+M)u|²/n² = (n + |x|²)/n²`. Equality holds iff the projection
  is `G_T`, i.e. `λ = u/n`, i.e. `(I + MᵀM)u ∝ u`.

**3.** `MᵀM = −M²`, so `G_T = O_T` iff `M²u ∝ u`. Apply `uᵀ` to get `λn = −uᵀM²u = |x|²`.
- **The 2-path form.** Write `(M²u)_v = Σ_(v→k) x_k − Σ_(k→v) x_k` with `x_k = 2d⁺_k − (n−1)` and `d⁺_k = n − 1 − d⁻_k`.
  This gives

      (M²u)_v = 2Σ_(v→k) d⁺_k + 2Σ_(k→v) d⁻_k − (n−1)(d⁺_v + d⁻_v) = 2Δ_v − (n−1)².

- **Sonneborn–Berger.** With `SB(v) = Σ_(v→k) d⁺_k` and `Σ_k d⁺_k = C(n,2)`,

      Δ_v = 2SB(v) − (n−2)d⁺_v + (n−1)² − C(n,2).
- **(a)** Let `n` be odd.
  - `M` is skew of odd order, so `Mz = 0` for the vector of signed principal sub-Pfaffians,
    `z_i = (−1)^(i+1) Pf(M_(îî))`. This holds because `adj M = zzᵀ` and `M·adj M = 0`.
  - Each `M_(îî)` is skew of even order with `±1` off the diagonal. Its Pfaffian is a signed sum of `(n−2)!!` terms
    `±1`, hence odd.
  - So every `z_i` is odd, `z ≠ 0`, and `Σ z_i` is a sum of `n` odd numbers, hence odd and nonzero.
  - Now `zᵀM²u = −(Mz)ᵀ(Mu) = 0`, while `zᵀ(−λu) = −λ Σ z_i`. So `λ = 0`, hence `|x|² = 0` and `T` is regular. The
    converse is immediate.
- **(b)** Write `T = A ⇒ B` with `|A| = a` and `|B| = b`, so `M = [[M_A, J], [−Jᵀ, M_B]]`. Put `x_A = M_A u_A` and
  `x_B = M_B u_B`; note `u_Aᵀx_A = u_Bᵀx_B = 0`.
  - The two blocks of `M²u = −λu` read

        M_A x_A + b x_A = (ab − λ) u_A,    M_B x_B − a x_B = (ab − λ) u_B.

  - Dot the first with `x_A`: `x_Aᵀ M_A x_A = 0`, so `b|x_A|² = 0` and `x_A = 0`. Likewise `x_B = 0`.
  - Both parts are regular, hence strong, so they are the two strong components. Conversely, `R_a ⇒ R_b` gives
    `x = (b u_A, −a u_B)` and `M x = −ab·u`.
- **(c)** Let `S = {i : x_i = +1}` and `S̄ = {i : x_i = −1}`. They have equal size `k = n/2`, since `Σ x_i = 0`, and
  `λ = |x|²/n = 1`.
  - For every vertex `i`, let `e_S(i) = Σ_(j∈S) M_ij` and `e_S̄(i) = Σ_(j∈S̄) M_ij` be its excess of wins over losses
    against each class. Then `x_i = e_S(i) + e_S̄(i)` and `(Mx)_i = e_S(i) − e_S̄(i)`.
  - The condition `(Mx)_i = −1`, together with `x_i`, gives:
    - for `i ∈ S`: `e_S(i) = 0` and `e_S̄(i) = 1`;
    - for `i ∈ S̄`: `e_S(i) = −1` and `e_S̄(i) = 0`.
  - So `T[S]` and `T[S̄]` are regular. Since `e_S(i)` for `i ∈ S` is a sum of `k − 1` signs, this needs `k` odd, i.e.
    `n ≡ 2 (mod 4)`. The cross counts are `(k+1)/2` and `(k−1)/2`, as stated.
  - Conversely these conditions give `Mx = −u`.
  - For existence (`k ≥ 3` odd), take any regular tournaments on `S` and on `S̄` and a circulant cross pattern: `s_a`
    beats `t_b` iff `b − a mod k ∈ {0, 1, …, (k−1)/2}`. Every `t_b` then beats `(k−1)/2 ≥ 1` vertices of `S`, so the
    tournament is strong.
- **(d)** For even `n` all `x_j` are odd. Then

      (Mx)_i = 2Σ_j M_ij d_j − (n−1)x_i,    Σ_j M_ij d_j ≡ C(n,2) − d_i (mod 2).

  So `−λ ≡ 2C(n,2) + (n−1)² (mod 4)`, i.e. `λ ≡ 3` if `4 | n` and `λ ≡ 1` if `n ≡ 2 (mod 4)`. In both cases
  `λ ≡ n − 1 (mod 4)`.

**4.** Suppose `I + MᵀM = (n−c)I + cJ`.
- Then `MᵀM u = (n−1)(1+c)u`, and `MᵀM y = (n−c−1)y` for every `y ⊥ u`.
- The vector `x = Mu` is orthogonal to `u`, and `MᵀM` commutes with `M`. So
  `MᵀM x = M(MᵀM u) = (n−1)(1+c)x`.
- If `x ≠ 0`, comparing with `(n−c−1)x` gives `nc = 0`, i.e. `c = 0`. If `x = 0`, then `(n−1)(1+c) = 0`, i.e.
  `c = −1`.
- Case `c = −1`: `T` is regular with constant `δ`, equivalently constant common out-degree, i.e. doubly regular.
- Case `c = 0`: `(I+M)ᵀ(I+M) = nI`, so `I + M` is skew-Hadamard. Conversely every skew-Hadamard matrix
  `H = I + S` gives a tournament with `M = S`.
- The equivalences with the skew-Hadamard conjecture are Reid–Brown's bordering of a doubly regular tournament
  (KNOWN).

**5.** By the definition of `p_i`, the out-neighbourhood `i + {1, 2, 4}` is a translate of the difference set
`{1, 2, 4}`, i.e. a cyclic Fano line.

**6.** Apply THM-4910 to the `(n−1)`-simplex: `m = n − 1`, so `m + 1 = n`.
- `Pow(O) = (m+1)R² = n r_T²`.
- The defect is `Σ_(j≠v)|p_v − p_j|² = 4Σ_(j≠v)(1 + δ_vj)`.
- Count `Δ_v` by the separating vertex `k`. If `k → v`, the `j` separated from `v` by `k` are the in-neighbours of `k`.
  If `v → k`, they are the out-neighbours of `k`. Summing gives the 2-path count.
- If `G_T = O_T`, the basis vector `G − O` vanishes and the regime fails. ∎

## Remarks

- **Why this is a bridge** (STRUCTURAL).
  - The seven-point sphere is about how a simplex's centres separate: `O ≠ G`, `O' ≠ G'`, non-coplanar.
  - The signature simplex turns tournament regularity into exactly these centre conditions:
    - regular (odd `n`) means the centres coincide;
    - doubly regular or skew-Hadamard means the simplex is equilateral.
  - Precisely:
    - for odd `n`, regular is the same as equicentric;
    - for all `n`, equilateral is the same as doubly regular or skew-Hadamard;
    - for even `n`, equicentric is a much larger class with no classical name.
  - Where the seven-point sphere exists, its data (`s(T)`, the 2-path counts `Δ_v`, the power `n(n − 1/s)`) are
    classical second-order score quantities.
- **The Fano plane sits where the sphere degenerates** (STRUCTURAL).
  - QR7's signature simplex is the Fano line code. It is equilateral, so it is equicentric, so the seven-point sphere
    does not exist for it.
  - The Fano tournament is one of many tournaments for which the sphere's construction collapses (`G_T = O_T`). The
    others include every equicentric tournament, for example every regular one.
- **The open question.** By 4, the skew-Hadamard conjecture asks which orders carry a tournament with pairwise
  orthogonal signature vectors. The weaker HYP-9500(a), a strong equicentric tournament at every even order, turned out
  to be elementary (audit J1). Its part (b) is the open content.
