---
id: THM-4742
title: "The seven-point sphere claim splits into a pencil of two inversion-invariant spheres through a common circle. For an n-simplex (n >= 3) with circumcentre O, centroid G and facet circumcentre O', when O, V0, G, O' span a 3-flat (not coplanar), all seven labels lie in that flat and: the six labels A, B, D, E, F, H' lie on one sphere S2 (new); session23b's sphere S1 carries B, C, D, E, F; both have power (n+1)R^2 at O, so both are invariant under inversion of radius sqrt(n+1) R; they meet in a circle through B, D, E, F in the plane O G O'; and the vertex C misses S2 by exactly Sigma = sum_i |V0Vi|^2 = -2(n+1)(V0 - O).(G - V0). So the seven labels share a sphere iff Sigma = 0 iff the median V0G is tangent to the circumsphere iff A = C. This never happens over an ordered field; over fields where -1 is a short sum of squares (Q(i), F_p) it happens, but only with A = C. In the tetrahedron, D = O + 4(G - O), F = O + 4(O' - O) and H' is the face orthocentre; in general H' = O' + n(G' - O') is the facet's own D-point. The plane (n = 2) is degenerate"
status: >
  PROVED: statements 1-4 for every n >= 3, as exact identities in the Gram basis v, s, w, homogeneous in R^2
  (seven_point_pencil.py (P1), (P1n)). FINITE-EXACT: (P2) 150 rational tetrahedra in actual coordinates; (P3) every
  valid random configuration over F_p, 5 <= p <= 97 (50,025 configurations, no mismatch); (P4) an explicit tetrahedron
  over Q(i); (P5) n = 3..6 in exact Gram coordinates, and n = 2 degenerate. CLASSICAL: Fano's seventh line (determinant
  -2). AUDITED 2026-10-21 by the session's independent audit; statement 5 was wrong in the first version and is
  corrected, along with credits, the field condition and the regime (MISTAKE-650).
  CREDIT: concurrent session23b packet 05-knowledge/results/tetrahedral_inversion_spheres_20261021/ (commit 7798cbbecb):
  - the false-claim verdict and the real exceptional locus;
  - the four-point obstruction A,C,D,E;
  - the sphere S1 with its inversion invariance;
  - Q1(H') = 6 - 2 v.s > 0 (i.e. Sigma > 0 over R) and A != C for a genuine tetrahedron;
  - that D != E forces the constant term.
  New here:
  - S2;
  - the pencil and its radical plane;
  - Pow_S2(C) = Sigma;
  - Sigma = 0 iff A = C over general fields;
  - the all-dimension form.
  Session opus-2026-10-21-S49.
session: opus-2026-10-21-S49 (owner: tournament structure vs open questions; the Fano plane vs a pasted 'seven-point sphere'; Gluck's paper)
source: 05-knowledge/results/seven_point_pencil_fano_and_converse_gluck_20261021.md
scripts:
  - 04-computation/experiments/seven_point_spheres_20261021/seven_point_pencil.py (+ .out)
  - 04-computation/experiments/seven_point_spheres_20261021/audit_S49/ (independent audit scripts and outputs: geo_*.py)
related:
  - 05-knowledge/results/tetrahedral_inversion_spheres_20261021/README.md (concurrent: S1, exceptional locus, counterexample)
  - 05-knowledge/results/fano_tournament_move_quotients_20261021/session23b-fano/README.md (concurrent: tetrahedral Fano incidence, real non-representability of the Fano plane)
---

# THM-4742 — two spheres through a common circle, and the seventh incidence

## Setting

* **The simplex.** An `n`-simplex `V0V1…Vn` (`n ≥ 3`) has circumcentre `O`, circumradius `R` and centroid `G`. Its facet `V1…Vn` has circumcentre `O'` and centroid `G'`.
* **The labels.**
  * `A` is the second point of the line `V0G` on the circumsphere.
  * `B` is the pole of the facet hyperplane: `(B − O)·(Vi − O) = R²` for `i ≥ 1`.
  * `C = V0`.
  * `D = O + (n+1)(G − O)`.
  * `E` is the inverse of `G` in the circumsphere.
  * `F = O + (n+1)(O' − O)`.
  * `H' = O' + n(G' − O')`, the facet's own analogue of `D`. For `n = 3` this is the orthocentre of the face triangle.
* **The basis.** Put `v = V0 − O`, `s = Σ_(i≥1)(Vi − O)` and `w = O' − O`.
  * Each `Vi − O'` (`i ≥ 1`) is orthogonal to `w`, so `s·w = n|w|² =: nk`.
  * Also `t = v + s = (n+1)(G − O)` and `T = |t|²`.
* **Labels in the basis.** Relative to `O`:
  * `C = v`, `D = t`, `E = (n+1)R² t/T`, `F = (n+1)w`, `B = R² w/k`, `H' = s − (n−1)w`;
  * `A = v + μ(s − nv)`, where `μ = −2 v·(s − nv)/|s − nv|²`.

  All seven labels lie in the 3-flat `O + span(v, s, w)`.
* **Regime.**
  * `v, s, w` are independent (`O, V0, G, O'` not coplanar).
  * `R²`, `|G − O|²`, `|O' − O|²` and `|G − V0|²` are nonzero, and the facet hyperplane is non-isotropic. Over `R` this reduces to `G ≠ O` and `O' ≠ O`.
  * For uniqueness, `D ≠ E`.
* **The inversion.** `ι(x) = O + (n+1)R²(x − O)/|x − O|²` is the inversion of radius `√(n+1) R` about `O`.

## Statements

1. **Two spheres.** Within the 3-flat, define spheres `|x − O|² − 2u·(x − O) + (n+1)R² = 0` by the values of `(2u·v, 2u·s, 2u·w)`:

       S1: ((n+2)R², T − R², (n+1)k + R²),      S2: (2 v·s + (2 − n)R², |s|² + 2nR², (n+1)k + R²).

   * `S1` contains `B, C, D, E, F` (session23b, for `n = 3`). `S2` contains `A, B, D, E, F, H'`.
   * Both have power `(n+1)R²` at `O`, so both are `ι`-invariant. `ι` swaps `D ↔ E` and `B ↔ F`. So `S1` also contains `ι(C)`, and `S2` also contains `ι(A)` and `ι(H')`.
   * With `D ≠ E`, any sphere through `D` and `E` (on one ray from `O`) has power `OD·OE = (n+1)R²`. So `S1` and `S2` are the unique such spheres.
2. **The pencil.**
   * The coefficient differences are `2(u1 − u2) = (Σ, −Σ, 0)` on `(v, s, w)`, with `Σ = 2nR² − 2 v·s`.
   * So `u1 − u2` is orthogonal to `t` and `w`, and the radical plane of `S1` and `S2` is the plane `O G O'`.
   * `S1 ∩ S2` is the circle `Γ` through `B, D, E, F`.
3. **The defect.**
   * `Pow_S2(C) = Pow_S1(H') = Σ = Σ_(i≥1) |V0Vi|² = −2(n+1)(V0 − O)·(G − V0)`.
   * `μ = Σ/|s − nv|²`, so `Σ = 0` iff the median `V0G` is tangent to the circumsphere at `V0`, iff `A = C`.
4. **The seventh incidence over fields** (characteristic `≠ 2` and not dividing `n(n+1)`, in the regime).
   * The seven labels lie on one sphere iff `S1 = S2`, iff `Σ = 0`, iff `A = C`.
   * **Ordered fields.** `Σ` is a sum of `n²` squares of coordinates, so a nondegenerate simplex has `Σ > 0` and the seventh incidence never holds.
   * **Fields where it can hold.** A nondegenerate solution needs `−1` to be a sum of at most `n² − 1` squares (level `≤ n² − 1`; for the tetrahedron, `≤ 8`).
     * It does hold over `Q(i)`, for `V0 = 0`, `V1 = e1`, `V2 = e2`, `V3 = (0, i, i)`.
     * It holds over `F_p` (`5 ≤ p ≤ 97`) in every sampled configuration with `Σ = 0`, and only there.
     * Every time it holds, `A = C`, so seven distinct cospherical labels never occur in the regime.
   * The real exceptional locus of session23b has `D = E` and `B = F`, so it lies outside the regime.
5. **Dimensions.**
   * Statements 1–4 hold for every `n ≥ 3` ((P1n) symbolic; (P5) exact for `n = 3…6`).
   * In the plane (`n = 2`), `v, s, w` are always dependent, so the statement is void.

## Proofs

**1–3.**
* Write `Q(x) = |x|² − 2u·x + (n+1)R²` relative to `O`. For `X = αv + βs + γw`, `|X|²` comes from the Gram matrix

      [[R², v·s, v·w], [v·s, |s|², nk], [v·w, nk, k]],

  and `2u·X = α(2u·v) + β(2u·s) + γ(2u·w)`.
* **The coefficients of `S2`.** They are forced by `D`, `F` and `H'`:
  * `Q(D) = 0` gives `2u·t = T + (n+1)R²`.
  * `Q(F) = 0` gives `2u·w = (n+1)k + R²`.
  * `Q(H') = 0` gives `2u·s = |s|² + 2nR²`, using `|s − (n−1)w|² = |s|² − (n² − 1)k`.
  * Subtracting gives `2u·v = 2 v·s + (2 − n)R²`, since `T = R² + 2 v·s + |s|²`.
* **The remaining labels of `S2`.**
  * A constant term of `(n+1)R²` gives `Q(ι(x)) = (n+1)R² Q(x)/|x|²`, so `E = ι(D)` and `B = ι(F)` lie on `S2`.
  * For `A`, using `|A|² = R²`:

        Q(A) = 2nR² − 2 v·s − μ(2u·t − (n+1)2u·v) = Σ − μ |s − nv|² = 0,

    since `−2 v·(s − nv) = Σ` and `T + (n+1)R² − (n+1)(2 v·s + (2 − n)R²) = |s − nv|²`.
* **`S1`.** Replace `H'` by `C`: `Q(C) = 0` gives `2u·v = (n+2)R²` and `2u·s = T − R²`.
* **The defects.**
  * `Q2(C) = R² − 2 v·s − (2 − n)R² + (n+1)R² = 2nR² − 2 v·s`, and `Q1(H') = |s|² − T + (2n+1)R²` is the same number.
  * Since `|Vi − O| = R`, `Σ_i |V0 − Vi|² = 2nR² − 2 v·s`.
  * `(V0 − O)·(G − V0) = v·(t/(n+1) − v) = (R² + v·s)/(n+1) − R² = −Σ/(2(n+1))`.
* **The pencil.** Subtract the two coefficient triples.
* The identities are checked symbolically in `n` and `R²` in (P1n), at `n = 3` in (P1), and exactly in actual coordinates in (P2) and (P5). ∎

**4.**
* If one sphere `S` contains all seven labels with `D ≠ E`, then `S` has constant term `(n+1)R²`.
  * `S ⊇ {B, C, D, E, F}` forces `S = S1`, since the coefficients are determined (`v, s, w` independent).
  * `S ⊇ {A, B, D, E, F, H'}` forces `S = S2`.
* By 2, `S1 = S2` iff `Σ = 0`. The identities are polynomial and homogeneous, so they hold over any field of the stated characteristic.
* **Ordered fields.** `Σ` is a sum of squares of the `n²` coordinates of `V0 − Vi`, so it is zero only if every `Vi = V0`.
* **Other fields.** A nontrivial zero of a sum of `n²` squares exists iff `−1` is a sum of `n² − 1` squares.
* The `Q(i)` and `F_p` cases are exact computations ((P3), (P4)). ∎

**5.** For `n = 2`, `v, s, w` are three vectors in the plane, so they are dependent. ∎

## The Fano comparison (what is and is not claimed)

* **The two configurations.**
  * The real picture of the Fano plane (vertices, midpoints and centroid of a triangle) realises six of its seven lines. The seventh line, through the three midpoints, has determinant `−2` in `0/1` coordinates, so it exists exactly when the characteristic is 2. It fails over `Q(i)`.
  * The seven-label sphere realises its six-label sphere `S2` always. Its seventh incidence, `C ∈ S2`, has defect `Σ`, a sum of squares. It holds exactly when that sum of squares vanishes nontrivially, which never happens over an ordered field. Even then `A = C`.
* **What is common.** In both cases a single quantity, nonzero over `R`, decides the last incidence. For Fano it is the characteristic (`2 ≠ 0`). For the sphere it is positivity (`Σ > 0`).
* **What is not claimed.**
  * The seven labels are not identified with the seven Fano points.
  * Session23b's tetrahedral Fano incidence (six edges plus a formal point) is a different, combinatorial bridge.

## Audit

The independent audit of 2026-10-21 (`audit_S49/`) found:

1. **The first statement 5 was false.** It claimed the result failed for `n = 4, 5`. The test was vacuous: every label lies in a 3-flat, so the sphere systems in `R^n` were singular and silently skipped. The audit showed the theorem holds for every `n ≥ 3`, and (P5) and (P1n) now prove and check it.
2. **The regime was understated.** Non-coplanarity of `O, V0, G, O'` is needed. In session23b's counterexample `H', B, F` are collinear.
3. **The field condition was overstated.** It said "non-real fields"; the right condition is level `≤ n² − 1`.
4. **The Fano comparison was mis-typed.** Fano's obstruction is the characteristic, not positivity.
5. **The credits were incomplete.** Session23b already had `Q1(H') > 0` and `A ≠ C` over `R`.

The geometry was re-derived independently (`geo_exact.py`, `geo_symbolic.py`, `geo_qi.py`, `geo_fields.py`; p = 5 exhaustive: 368,640 in-regime configurations, no mismatch). See MISTAKE-650.
