# Placing the 17-vertex Exoo-Ismailescu graph in Q(zeta55): excluded to coordinate height 1e63, the published field certificate refuted, exact decision OPEN

2026-10-12, opus-2026-10-12-S24 (addendum). **VERIFIED (numerical, with margins
stated):** no placement of EI17 in `L = Q(zeta55)` with coordinate height below
about `1e63`; the degree-20 polynomial published as the minimal polynomial of a
coordinate does not vanish at the published realization. **PROVED:** the
parallelogram structure (13 directions, 6 relations) and the exact facts about
the published polynomial (irreducible, Galois group `S20`). **NUMERICAL (lower
bound):** the census of at least eleven non-congruent real faithful realizations.
**OPEN:** the exact coordinate field of EI17 and hence the exact placement
decision. **CITED:** Exoo-Ismailescu, Geombinatorics 26 (2016) 49-64 (the graph);
House of Graphs 51375 (adjacency, API read 2026-10-12);
arXiv:2607.19995v3 (Costa Silva Filho, July 2026; the certificate refuted here).
Self-audited; no independent peer audit.

## Why this test, and what it decides

[HYP-9247](../hypotheses/HYP-9247-every-finite-unit-distance-graph-in-q-zeta55-is-three-colourable.md)
conjectures that every finite unit-distance graph in `L` is three-colourable;
[THM-4614](../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md)
makes that conjecture decisive for six-colour forcing at the eight ports. The
cheapest refutation would be a known triangle-free four-chromatic unit-distance
graph with all coordinates in `L`. The smallest known is EI17: 17 vertices, 31
edges (`2*17-3`, minimally rigid), girth 4, chromatic number 4
(Exoo-Ismailescu 2016; the smallest such graph per Soifer 2024).

A placement of a rigid graph in `L` is a placement of one of its finitely many
realizations; in the gauge `v17 = 0`, `v11 = -1` (a unit edge), a realization
lies in `L` up to isometry iff all its edge directions lie in `L` (translate by
an element of `L`, rotate by a unit of `L`).

## 1. The graph, the certified realization, and its structure

The House of Graphs adjacency list and the paper's edge list are isomorphic in
exactly two ways (`|Aut(EI17)| = 2`). The paper's Appendix A realization (24
digits) satisfies all 31 unit edges; Newton refinement in the full 30-equation
system converges from it to residual `1e-3051` with the 24 printed digits
unchanged, so it is a genuine faithful realization accurate to 3000 digits.

EI17 has nine unit 4-cycles. In any faithful realization a unit 4-cycle is a
parallelogram (four unit vectors summing to zero pair off as opposites; the
other pairing collapses two vertices). Hence the 31 edges carry only **13
distinct directions**, with the four edges `(2,10), (3,12), (6,16), (11,17)`
all equal to the gauge direction `1`, and the 15-dimensional cycle space yields
**6 independent integer relations** among the 13 unit vectors
(`u1 + u7 + u10 + u11 + u13 = 2u6` and five more). So a realization is six
free unit vectors `a = u7, b = u9, c = u10, d = u11, e = u12, f = u13`
with `u6 = 1`, subject to six equations `|L_k(a,...,f)| = 1` for the six
dependent directions, e.g. `|1 + b + e - c| = 1` and
`|2c + d - 2e - f - 2b| = 1`: six angles, six equations, **no sign branches**.
Every apex choice of the paper's two-angle construction is forced by a
parallelogram.

## 2. The published field certificate does not hold at the published realization

arXiv:2607.19995 states that the coordinate `x_v1` of this realization is a root
of an explicit irreducible integer polynomial `p` of degree 20 with leading
coefficient near `1e97` (Appendix B; "found by algdep at 2600 digits and
verified `|p(alpha)| < 1e-1900`"), with Galois group `S20`, and infers that
EI17's coordinates are not expressible by radicals.

Exact facts about `p` (reproduced here): it is primitive and irreducible over
`Q`; its Frobenius cycle types include `(20)` at `p = 3` and `(3,17)` at
`p = 71` (distinct-degree factorization, cross-checked by brute-force root
counts at 3, 5, 7, 29, 71, 109), so by Jordan's theorem and the odd 20-cycle,
`Gal(p) = S20`. That much is right.

But `p` **is not the minimal polynomial of `x_v1`**:

* `p` has four real roots, `-2.49880625682906..., -1.3496, -0.1783, 1.0626`;
  the certified coordinate is `x_v1 = -2.49880616862106307378...`. The nearest
  root differs from it by `8.8e-8`, and `|p(x_v1)|` is `3.8e-8` times the size
  of its terms, at a coordinate correct to 3000 digits.
* No coordinate of the realization, nor its negative, is a root of `p`.
* Integer LLL (python-flint) on `1, x_v1, ..., x_v1^d` at 2900 digits finds,
  for `d = 20, 24, 30, 40, 48, 60`, only vectors of the random height
  (`1e138, 1e116, 1e93, 1e71, 1e59, 1e48` against random expectations
  `1e131, 1e111, 1e90, 1e69, 1e58, 1e46`). In particular no degree-20 integer
  polynomial of height below about `1e128` vanishes at `x_v1`; the published
  `p` has height `1e98`.

So the number that was fed to algdep agreed with `x_v1` to about eight digits
and not beyond; `p` is an exact relation for some other number (its height is
far below the random expectation at 2600 digits, so it is not a random
artefact), and the certificate, the Laman number 20 read from it, and the
"not by radicals" conclusion for EI17 are unsupported. This is the paper's own
Remark 22 pitfall one level up. Its EI19 analysis was not examined.

## 3. Placement in Q(zeta55): excluded to height 1e63

With the 3000-digit realization, each of the 13 directions `u` was tested for
membership in `L` by integer LLL on the lattice spanned by the 40 power-basis
numbers `Re(zeta^k) + pi * Im(zeta^k)` and `Re(u) + pi * Im(u)` at 2900
digits (a relation with the transcendental mixing constant `pi` forces the real
and imaginary relations separately). The gauge direction `u6 = 1` returns its
exact relation; **each of the other twelve returns only vectors of height
`1e70`**, the random expectation being `1e69`. Hence no placement of the
certified realization in `L` has coordinates of height below about `1e63`
(LLL approximation factor `2^20` absorbed).

The same holds for the other real realizations (section 4) at 1000 digits in
the scratch runs; the committed script tests the certified one.

Consequently HYP-9247 survives its cheapest decisive test, and the published
exact obstruction cannot be cited for it. The honest status of "EI17 in `L`"
is **excluded numerically to height `1e63`, exactly OPEN**.

## 4. At least eleven real faithful realizations, not two

The two-angle construction with all `2^13` apex sign branches, seeded from a
grid over the two angles where both closure residuals are small, and
Newton-polished in the full system, finds non-congruent real faithful
realizations (all 17 points distinct, no extra unit distances; distinct
distance multisets, which are labelling-independent). With seeding threshold
`0.08`: six at grid `96`, nine at grid `160`; with threshold `0.12` (the
committed run): ten at each grid and **eleven in the union**, one of them
seeded only at grid 96 and one only at grid 160; the paper reports "more than
one". The count grows with the seeding and is only a lower bound; two of the
realizations have a vertex pair just `0.0014` apart. A second realization
shares the certified `x_v1 = -2.4988` (under the automorphic relabelling) and
differs from the certified one in exactly the four vertices `4, 6, 13, 16`
(the free angle at `6` and the rhombus `6-16` re-placed), the other 13
vertices coinciding to `1e-16`. A random-start search over the full
31-dimensional system finds, by contrast, thousands of **collapsed** solutions
(coincident non-adjacent vertices, many on positive-dimensional families:
Jacobian ranks 28-30), which is why the unit-length system is not
zero-dimensional and why an exact elimination must saturate the degenerate
components.

## 5. What would decide it exactly, and what the structure suggests

* **Exact route.** Eliminate the six-angle system of section 1 (six Laurent
  equations in `a,...,f`) after saturating the collapse components; the
  elimination polynomial of `a` (or of `x_v1`) mod two good primes gives the
  degree and Frobenius cycle types; a factor degree not dividing 40, or
  non-uniform cycle types on the orbit of the certified solution, proves
  `x_v1` outside every abelian field and settles the question. sympy's Gröbner
  engine did not finish the 12-variable version in 50 minutes; a dedicated
  engine (msolve, Singular, Magma) should.
* **Reading.** The LLL searches show that `x_v1` has either degree above 60
  or a minimal polynomial of height above `1e46` (above `1e128` if the degree
  were 20). An element of `L` has degree dividing 40 and would need height
  above `1e63` in the direction coordinates. Both are possible in principle
  (unit vectors of `L` have unbounded denominators), so the exclusion is
  numerical, not a proof.

## 6. Connection ledger

| Source object | Map and preserved predicate | Lost data / retained sidecar |
|---|---|---|
| HoG adjacency, paper edge list | two isomorphisms | automorphism swaps the gauge vertex `11` with `14` |
| Nine unit 4-cycles | parallelogram forcing | sign branches of the apex construction |
| 31 edges | 13 directions, 6 relations | the degenerate (collapsed) solution components |
| Published `p` | exact Galois facts (`S20`) | its root is not the published coordinate (gap `8.8e-8`) |
| `x_v1` at 3000 digits | LLL over degrees to 60 | heights only; no exact degree |
| 13 directions | LLL in the `zeta55` power basis at 2900 digits | no relation below height `1e69` |

## 7. Reproduction

    python3 -B 04-computation/experiments/ei17_placement_q_zeta55_20261012.py
    python3 -B 04-computation/experiments/ei17_placement_q_zeta55_20261012.py --census 96,160

[Script](../../04-computation/experiments/ei17_placement_q_zeta55_20261012.py);
[output](ei17_placement_q_zeta55_20261012.json) (fast part: isomorphisms,
3000-digit refinement, 4-cycles and relations, the exact facts about `p` and
its failure at `x_v1`, the LLL degree search, the placement LLL);
[census output](ei17_placement_q_zeta55_20261012_census.json) (adds the merged
two-grid census of real faithful realizations, about 35 minutes). The paper's coefficients were parsed from the PDF
text and agree with an independent transcription. Integer LLL uses sympy's
`DomainMatrix.lll` (python-flint backend); mpmath Newton in the 30-equation
system supplies the digits.
