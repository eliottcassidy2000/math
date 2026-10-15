---
id: THM-4650
title: "Guy D19 (a point at rational distance from all four corners of the unit square): no such point lies on the four side lines, the two midlines, the two diagonals, or the circumcircle. The midline case is Fermat's right-triangle theorem: the midline curve is 2-isogenous to the congruent-number curve y^2 = x^3 - x (conductor 32, rank 0, torsion Z/4), so a midline point would make 1 a congruent number; the side-line case is the conductor-40 curve y^2 = x^3 - 2x + 1 (rank 0, torsion Z/4). On the midline of a 1 x a rectangle the problem is exactly a rational triangle with base 2a and height 1 (Guy's base/height-ratio remark in D19), so the square's sides and midlines are the ratios 1 and 2"
status: >
  PROVED: statements 1-4. KNOWN in substance: the non-existence on the diagonals, midlines and side lines is Yang Ji
  (arXiv:2105.05250, 2021, Fermat descent), with further special cases in Ieno (arXiv:2111.02888); the rectangle family of
  statement 5 overlaps Bremner-Ulas (arXiv:1502.07312), and the base/height-ratio problem is MacLeod (SWJPAM 2 (2003) 9-17).
  New here (not found in those abstracts): the circumcircle via Ptolemy, the explicit elliptic-curve identifications
  (conductor 40; conductor 32 with the Velu 2-isogeny from y^2 = x^3 - x), the j-invariant of the rectangle family and its rank
  table. The two rank-0 facts are unconditional: L(E,1) = 0.6555... and 0.7422... are nonzero (audit B), so rank 0 follows from
  Kolyvagin (Coates-Wiles for the CM midline curve); PARI's ellrank 2-descent bound [0,0] agrees. The isogeny in 3 is the Velu
  2-isogeny with kernel <(-1,0)> (ellisogeny; the original isomorphism loop of d19_isogeny.gp aborted and was fixed after
  audit B, MISTAKE-610). FINITE-EXACT: the rank table of statement 5 (ellrank with equal bounds) and the brute-force
  cross-check of statements 2-3 (no solutions among the 5,472,915 reduced fractions p/q with |p| <= 3000, q <= 1500).
  CITED: Berry 1990; Corvaja-Turchet-Zannier 2025. The general D19 question remains OPEN. Audited (audit B, MISTAKE-610).
session: mac-mini-2026-10-14-wilson
source: 05-knowledge/results/d19_parity_shift_synthesis_20261014.md
scripts:
  - 04-computation/experiments/d19_square_distances_20261014/d19_special_curves.gp (+ .out)
  - 04-computation/experiments/d19_square_distances_20261014/d19_isogeny.gp (+ .out)
  - 04-computation/experiments/d19_square_distances_20261014/d19_rectangle_family.gp (+ .out)
  - 04-computation/experiments/d19_square_distances_20261014/d19_lines_check.py (+ .out)
related:
  - continuum_finite_witness_boundary_20261014f.md section 4 (codex: the integer form of D19 and its mod-24 local guard; this theorem is a global obstruction on special curves)
  - Guy, Unsolved Problems in Number Theory, D19 (rational distances from the corners of a square) and D27 (congruent numbers)
  - T. G. Berry, Points at rational distance from the corners of a unit square, Ann. Sc. Norm. Sup. Pisa (4) 17 (1990) 505-529
  - THM-4558, THM-4628 and the Hadwiger-Nelson notes (unit- and rational-distance graphs in field planes)
---

# THM-4650 — Guy D19: the special loci of the square carry no rational-distance point

## Setting

Let the square have corners `(0,0), (1,0), (0,1), (1,1)`. A **D19 point** is a point `P = (x, y)` whose distances `a, b, c, d` to the four corners (in that order) are all rational. Guy's D19 asks whether one exists; it is open.

* **Rationality of `P` (classical).** `x = (a² − b² + 1)/2` and `y = (a² − c² + 1)/2`, so a D19 point has rational coordinates. D19 is therefore a question about rational points on the surface `S` given by `a² = x² + y²`, `b² = (1−x)² + y²`, `c² = x² + (1−y)²`, `d² = (1−x)² + (1−y)²`.

## Statements

1. **Diagonals and circumcircle (PROVED, elementary).** No D19 point lies on the diagonal `y = x`, the anti-diagonal `x + y = 1`, or the circumcircle `(x − ½)² + (y − ½)² = ½`.
2. **Side lines (PROVED).**
   * No D19 point lies on the line `y = 0`. By symmetry the same holds for the lines `x = 0`, `x = 1` and `y = 1`.
   * The condition is that `x² + 1` and `(1 − x)² + 1` are both rational squares. With `x = (s² − 1)/(2s)` this becomes `w² = s⁴ − 4s³ + 6s² + 4s + 1`.
   * That quartic is a genus-1 curve with Jacobian `y² = x³ − 2x + 1`: conductor 40, rank 0, torsion `Z/4`.
   * Its four rational points are `s = 0` (two points) and `s = ∞` (two points). All of them are degenerate (`x = ∞`).
3. **Midlines (PROVED).**
   * No D19 point lies on `x = ½`, nor by symmetry on `y = ½`.
   * The condition is that `4y² + 1` and `4(1 − y)² + 1` are both squares. With `y = (t² − 1)/(4t)` this becomes `w² = t⁴ − 8t³ + 18t² + 8t + 1`.
   * Its Jacobian `E_mid : y² = x³ − 11x + 14` has conductor 32, `j = 287496 = j(2i)`, rank 0 and torsion `Z/4`. Its four rational points are again `t = 0, ∞`.
   * **Fermat's right triangle.** `E_mid` is 2-isogenous to `y² = x³ − 16x`, which is isomorphic over Q to `y² = x³ − x` via `(x, y) ↦ (4x, 8y)`. That is the congruent-number curve for `n = 1`.
   * Isogenous curves have equal rank. So a D19 point on a midline exists iff `rank E_1(Q) > 0`, iff 1 is a congruent number. Fermat's descent (no right triangle with rational sides and area 1, equivalently `x⁴ − y⁴ = z²` has no nontrivial solution) excludes it.
   * This links Guy's D19 to his D27.
4. **The parametrizations are complete (PROVED).**
   * Every rational `(y, u)` with `u² − 4y² = 1` has `t = u + 2y ≠ 0`, `u − 2y = 1/t`, so `y = (t² − 1)/(4t)`. The same holds for the side-line conic.
   * Each quartic has a rational point, so it is isomorphic to its Jacobian and has exactly `#E(Q) = 4` rational points.
5. **Rectangles and the base/height ratio (PROVED reduction; FINITE-EXACT table).**
   * On the midline `x = ½` of the `1 × a` rectangle, the D19 conditions are `1 + u² = □`, `1 + v² = □` and `u + v = 2a`, with `u = 2y` and `v = 2(a − y)`.
   * Equivalently: two rational right triangles with a common leg 1 glued along it. That is a rational triangle with base `2a` and height 1, i.e. base/height ratio `2a`, which is the ratio problem in Guy's D19 remarks (MacLeod).
   * The square's side lines are ratio 1, and its midlines are ratio 2.
   * The Jacobians `E_a` have `j = 16(16a⁴ + 16a² + 1)³/(a²(a² + 1))`, a non-constant `j`. So this is not a twist family, and only the square (`a = 1`) has CM by `Z[2i]`.
   * **Rank table.** Among the 46 reduced ratios `n/d` with `n ≤ 12`, `d ≤ 6`:
     * rank 1 for `1/3, 1/5, 3/4, 3/5, 4/3, 4/5, 5, 5/2, 6, 7/2, 7/4, 7/5, 7/6, 8, 9, 9/5, 10/3, 11/3, 12/5`;
     * rank 0 for the rest, including 1, 2, 3, 4 and 1/2.
   * Torsion `Z/2 × Z/4` (order 8; `elltors` gives `[8, [4, 2]]`) appears at `3/2, 5/6, 8/3`. These are the isosceles cases, where `P` is the centre of a rectangle with rational diagonal. (The pre-audit text said `Z/8`; the script printed only the order.)

## Proofs

**1.**
* On `y = x`, `a = √2·|x|` is rational only for `x = 0`. Then `d = √2`. The anti-diagonal is symmetric.
* On the circumcircle, let `P` lie on the arc `AB` of the square `ABCD`. Ptolemy for the cyclic quadrilateral `A, P, B, C` gives `AB·PC = AP·BC + PB·CA`, i.e. `PC = AP + √2·PB`. With all three distances rational this forces `PB = 0`, so `P = B`, whose distance to `D` is `√2`. ∎

**2–3.**
* Each condition `k² + z² = □` (with `k` fixed) is a conic with the rational point `(0, k)`. Its rational points are given by the stated parametrization (statement 4).
* The second condition then becomes the displayed quartic.
* `ellfromeqn` gives its Jacobian, and `ellminimalmodel`, `ellglobalred`, `elltors` and `ellrank` give conductor, torsion and the rank bounds `[0, 0]`. The analytic rank is also 0.
* The isogeny class of `E_mid` (`ellisomat`) contains `[a4, a6] = [−16, 0]` at degree 2, and `ellisisom` identifies it with `y² = x³ − x` by `(u, r, s, t) = (2, 0, 0, 0)`. Directly: Vélu's 2-isogeny from `y² = x³ − x` with kernel `⟨(−1, 0)⟩` has image exactly `y² = x³ − 11x + 14` (`ellisogeny`). Also `a_p(E_mid) = a_p(y² = x³ − x)` for every prime `3 ≤ p < 1000`. ∎

**4.**
* The conic map is the inverse of `t = u + 2y`, and similarly for the side line.
* A genus-1 curve with a rational point is an elliptic curve isomorphic to its Jacobian. ∎

**5.**
* Substitute `u = 2y` and `v = 2a − 2y`.
* For the triangle: the altitude of length 1 splits the base `u + v` into the legs of two right triangles with hypotenuses `√(1+u²)` and `√(1+v²)`.
* The table is `d19_rectangle_family.out`. ∎

## Reading

* **Where the general problem sits.**
  * `S` is a `(Z/2)⁴`-cover of the plane branched along the eight isotropic lines through the corners. Four of these lines pass through each circular point.
  * Dropping one corner gives Berry's K3 surface, on which three-distance points are dense.
  * The four-corner surface is of general type: Corvaja–Turchet–Zannier (Geom. Dedicata 219 (2025), arXiv:2403.02030) prove that `r >= 4` points give a simply connected surface of general type (here `K² = 8`, `p_g = 5`). So the Bombieri–Lang conjecture would put all but finitely many D19 points on finitely many curves of genus 0 or 1, the same mechanism as the conditional Erdős–Ulam results (Shaffaf; Tao). Bombieri–Lang is ineffective, so this gives no finite search. Corvaja–Turchet–Zannier also note that Berry's density hypotheses fail when the three points are rational.
  * Statements 2–3 eliminate the natural genus-1 curves.
* **Local versus global.**
  * codex's mod-24 guard (continuum note §4) is local: `12 | L`, and one coordinate divisible by 4 and one by 3.
  * Statements 2–3 are global Mordell–Weil obstructions on curves, and they come out as Fermat-type descents.
* **Arithmetic inside geometry.** The midline of D19 is Fermat's first descent; the sides are a sibling curve of conductor 40. A geometric rational-distance question about the unit square, a Hadwiger–Nelson-type object, is governed on its symmetry lines by congruent-number arithmetic.
