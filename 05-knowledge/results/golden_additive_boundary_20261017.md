# Additive colourings free the golden boundary: a pre-screen, a dead v^2 candidate, and a forcing design loop

2026-10-17, opus-2026-10-17-S41. Owner prompt: golden-coordinate six-colour forcing, preserving sources, ports and
integer clocks; approach from the other end (design a system where the desired control is inherent, then merge).
Builds on clock_duality_oct15's golden contact:
- the complete unit-direction classification in `K + uK + vK` (`K = Q(zeta15, sqrt(-11))`, contact family
  `q(d, 1, 1)`, `d = 1 + zeta5`). This is lemma G1 of
  [marked_contacts_root_current_20261017e.md](marked_contacts_root_current_20261017e.md) section 5, with its hypotheses
  certified for the golden contact in [joint_contact_field_20261017f.md](joint_contact_field_20261017f.md) section 5;
- [shared_golden_contact_20261017g.md](shared_golden_contact_20261017g.md) (the 5400-vertex source `A + uC + vC`, all 570
  port patterns extend; proposed next step: `v^2` coordinates).

**Concurrent work (read after the computations below).**
[golden_square_contexts_20261017h.md](golden_square_contexts_20261017h.md) (clock_duality_oct15, committed at 20:25) uses
the same certificate type for its constructions.
- Its section 2 shows that the single translate `G ∪ (G + v^2)` is a Cartesian matching context: add 1 mod `k` on the
  copy. Section 3 below extends this to every `v^2 C` layer over the same coefficient sets.
- Its formula (3), `c(p + e l) = c(p) + c_L(l) mod k`, is an additive colouring.

Prior art for the method. Separate-coordinate additive (coordinate-sum) colouring CSPs, and "Cartesian products
retaining additive colourings" as a canonical hostile, already appear in clock_duality's
[graded_tail_carrier_20261017b.md](graded_tail_carrier_20261017b.md) (committed 14:10). The same note also gives the
divisibility argument that excludes additive homomorphisms. The lemma below is that method, stated for the golden
contact.

What is new here:
- the additive extension lemma for the non-Cartesian contact edges;
- the shared-coefficient result for all 570 patterns;
- the complete `v^2` census;
- the correspondence reformulation and the design loop.

Their degree-192 escape `E' = F(s, h)` with the direction `e` (their G1/G2) lies outside `K + uK + vK`, so HYP-9289 does
not cover it. The pre-screen is the cheap necessary test for their stated next target: a unit edge joining different
port and `e`-coefficient fibres.

**Typing.**
- PROVED (elementary): the additive extension lemma (clock_duality's coordinate-sum method), and the reduction of the
  additive problem to a `Z/6` colouring of the first coordinate with forbidden difference sets.
- PROVED: THM-4668 (clock_duality's degree-192 direction `e` is Cartesian on `K + uK + vK + eK`).
- FINITE-EXACT:
  - every pattern of clock_duality's receipt extends additively (verified on all 30,856 edges);
  - the `v^2` census (every unit difference with a nonzero `v^2` part is a pure `v^2`-axis unit);
  - the additive pre-screen on four larger linear golden sources, up to 2,297,178 mixed edges (a fifth, contact2, was
    built but its pre-screen did not finish in the session).
- OPEN: six-colour forcing; whether every finite golden source in `K + uK + vK` is additively free.

## 1. The other end: systems whose boundary is free for a structural reason

For a Cartesian sum `A □ C □ C` every colouring of the first factor extends: colour `(a, b, c)` by
`c_A(a) + c_1(b) + c_2(c) (mod 6)`. So the "inherently uncontrolled" (always-extending) system is the additive one. By
G1 (marked_contacts_root_current_20261017e.md section 5) a golden source differs from a Cartesian sum only in its
mixed edges
`(a, b, c) ~ (a + dq, b + q, c + q)`, `|q| = 1`.

**Lemma (additive extension; PROVED).** Let the source be `A' + uC' + vC'` with `A', C' ⊂ K` and ports `(p, 0, 0)`.
Suppose:
- `c_A` is a proper colouring of the unit graph of `A'` with the port pattern;
- `c_1` and `c_2` are proper colourings of the unit graph of `C'` with `c_1(0) = c_2(0) = 0`;
- on every mixed edge, `c_A(a + dq) - c_A(a) + [c_1(b + q) - c_1(b)] + [c_2(c + q) - c_2(c)] ≢ 0 (mod 6)`.

Then `c_A(a) + c_1(b) + c_2(c)` is a proper six-colouring of the whole source extending the pattern.
- *Proof.* An axis edge changes one coordinate along a unit edge of that factor. A mixed edge is covered by the last
  condition.
- The problem has `|A'| + 2|C'|` colour variables, independent of the size `|A'| |C'|^2` of the source.

**Corollary (pre-screen).** A golden source can force a port pattern only if its additive problem is unsatisfiable for
that pattern. This is a necessary test that runs in seconds before any full SAT run on a large graph.

**Correspondence form (PROVED).** Fix `c_1` and `c_2`. The constraints on `c_A` are then:
- `c_A(a) ≠ c_A(a')` on unit edges;
- `c_A(a + dq) - c_A(a) ∉ -S_q`, where `S_q = {Δ_q c_1(b) + Δ_q c_2(c)}` over the mixed edges at `(a, q)`.

So the additive extension problem is a `Z/6` colouring with forbidden difference sets (a correspondence/DP colouring
when every `S_q` is a single value) of the golden two-distance graph on
`A'`: distances 1 and `|d| = φ`, the golden ratio. The coefficient colourings choose the shifts. Two special cases:
- With `c_2 = -c_1` homomorphism-like, every shift is 0 and the problem is plain 6-colouring of the two-distance graph.
- That choice is too rigid. The golden pentagon's diagonals are `d`-steps (`P_(i+2) - P_i = d α^i`), so the ports
  span a `K_5`, and 465 of the 570 patterns fail it. The general additive problem escapes by using nonzero shifts.

## 2. clock_duality's source is additively free

For each of the 570 canonical patterns, `golden_additive_colourings_20261017.py` (input: clock_duality's receipt)
- finds an additive colouring with **`c_1 = c_2`**, i.e. one colouring of the shared golden coefficients used in both
  axes;
- checks the resulting 5400-vertex colouring against all 30,856 edges (146 x 100 A-axis, 15 x 540 on each coefficient
  axis, 56 mixed).

All 56 mixed edges change all three coordinates.
- This explains the receipt's 570/570 without its SAT witnesses.
- The doubled mixed-edge count (28 to 56) could not matter: the coefficient shifts absorb the contact constraints.

## 3. The proposed v^2 coordinate creates no cross edges

clock_duality proposed adjoining `v^2` coordinates. The quadratic relation `v^2 = -Zv - Z/conj(Z)` introduces a `uv`
coefficient outside the classified linear space, and they asked for a test of whether the new product directions create
cross edges.

`golden_vsquare_unit_census_20261017.py` enumerates every difference
`dA + u d1 + v d2 + v^2 d3`, with `dA ∈ D(A) ∪ {0}`, `d1, d2, d3 ∈ D(C) ∪ {0}` and `d3 ≠ 0`. That is
`1603 x 83 x 83 x 82 ≈ 9.1 x 10^8` candidates, with `|D(A)| = 1602` and `|D(C)| = 82`. The audit's independent float
census of 905,531,494 candidates confirms the result.
- **Sieve.** A modular sieve at the good primes 1021, 1171, 1741 tests the norm equation in `K(w)(s)`, using both
  embeddings. It leaves 1018, then 22, then 22 survivors.
- **Exact check.** All 22 survivors are exact units, and all 22 are the pure `v^2`-axis units
  (`dA = d1 = d2 = 0`, `d3` one of the 22 unit differences of `C`).

So `A + uC + vC + v^2 C` has the old edges plus Cartesian `v^2`-axis edges only. Adding `c_3 = c_C` extends the
additive colourings of section 2. **With these coefficient sets (`|A| = 54`, `|C| = 10`), the `v^2` candidate cannot
force.** Larger coefficient sets need a new census.

## 3b. clock_duality's degree-192 direction e is Cartesian on the whole linear space (THM-4668, PROVED)

clock_duality's escape `e = Zs/h` (`h^2 = 4 - N`, degree 192; golden_square_contexts_20261017h.md G2) makes `σ`-mates
unit-distance pairs. They named the next target: a unit edge joining different port and `e`-coefficient fibres.
- **Census.** `golden_e_direction_census_20261017.py` tests all `1603 x 83 x 83 x 130 ≈ 1.44 x 10^9` differences
  `dA + u d1 + v d2 + e d3` with
  `d3 ∈ D(C) ∪ D(L) \ {0}`, where `L` is their 13-point attachment set. Exactly the 26 pure `e`-axis units survive.
- **THM-4668** proves the general statement. In `K + uK + vK + eK` every unit with a nonzero `e`-coefficient is a pure
  `e`-axis unit `q e`, `|q| = 1`.
  - If the `v`-coefficient is nonzero, the norm equation forces `4 - N ∈ K`, which is impossible because `u` has degree
    3 over `K`.
  - If the `v`-coefficient is zero but the `u`-coefficient is not, the equation reduces modulo the exact cubic
    `u^3 = m2 u^2 + m1 u + 1` to an explicit nonzero element of `Q(zeta5, j)`
    (golden_e_direction_proof_check_20261017.py).

So every source `A + uC + vC + eL` over `K` is the old source multiplied Cartesianly along `e`; every boundary colouring
extends, and the named target edge does not exist in this space. Together with section 3 (`v^2`), the two natural
nonlinear moves of the golden program add no cross edges.

## 3c. Product coordinates: two new contact families, still additively free

`golden_product_direction_census_20261017.py` repeats the census for the product coordinates `Y = ue, ve, uv, e^2`. It
uses all four norm components in `F(s, h)`, with `d3 ∈ D(C) ∪ D(L)` (`1.44 x 10^9` differences each). Results
(FINITE-EXACT):

| `Y` | units with nonzero `Y`-part | pure `Y`-axis | cross units |
|---|---|---|---|
| `ve` | 26 | 26 | 0 |
| `ue` | 26 | 26 | 0 |
| `uv` | 36 | 26 | 10 |
| `e^2 = -Z/conj(Z)` (an element of `F`) | 36 | 26 | 10 |

`golden_product_relations_20261017.py` shows that each set of 10 cross units is one contact family, and both follow
from clock_duality's contact `|d + u + v| = 1`:
- **`uv`.** The family is `q (0, 1, 1, conj d)` in `(1, u, v, uv)`. Since `|u| = |v| = 1`,
  `|u + v + conj(d) uv| = |conj(d + u + v)| = 1`: the reflected contact.
- **`e^2`.** The family is `q (d, 1, 0, conj d)` in `(1, u, v, e^2)`. Indeed `d + u + conj(d) e^2 = Z conj(u)/conj(Z)`,
  a unit. This coordinate lies in the degree-48 base field `F`, so it needs neither `s` nor `h`.

The `e^2` family shares its `(A, u)`-step `(dq, q)` with the `v`-contact. That is the first genuinely coupled pair of
mixed families in the golden program.

`golden_forcing_search4_20261017.py` builds `A + uC + vC + e^2 C`: 54,000 addresses, 56 `v`-contact and 42
`e^2`-contact edge patterns per layer. For these coefficient sets the census plus G1 give the complete edge set. The
additive pre-screen still blocks **no** pattern, with or without shared coefficient colourings. Two coupled contact
families sharing an `(A, u)`-step are absorbed by independent coefficient shifts.

## 4. A design loop, and what four larger sources show

`golden_forcing_search_20261017.py` builds `A' + uC' + vC'` for `A', C' ⊂ K` directly from G1: unit edges by a
three-prime sieve plus exact norm in `K`; mixed edges combinatorially from the contact family. It then runs:
- the two-distance test;
- the additive pre-screen (shared and independent coefficient colourings);
- a full SAT run only for patterns the pre-screen blocks.

It reproduces clock_duality's source exactly: 146, 15 and 56 edges.

| recipe | `\|A'\|` | `\|C'\|` | A edges | C edges | mixed edges | two-distance blocked | additively blocked |
|---|---|---|---|---|---|---|---|
| base (clock_duality) | 54 | 10 | 146 | 15 | 56 | 465 | 0 |
| `dM`: `A' = P + M + dM` | 364 | 10 | 1139 | 15 | 2142 | 465 | 0 |
| contact1: `A' = A + d·U(C)` | 1099 | 10 | 3794 | 15 | 5302 | 465 | 0 |
| `dM_C2`: also `C' = C + U(C)` | 364 | 174 | 1139 | 532 | 2,297,178 | 465 | 0 |
| `MM`: `C' = M + M + pentagon` | 364 | 29 | 1139 | 73 | 72,358 | 465 | 0 |
| contact2: two contact steps | 10,597 | 10 | 44,745 | 15 | 87,846 | (skipped) | not finished (CSP of about 2 x 10^7 clauses) |

Here `U(C)` is the set of unit differences of `C`. With 4 x 10^4 times as many mixed edges, every pattern still extends
additively, even with `c_1 = c_2`.

**Reading.** Within the classified linear space `K + uK + vK`, enlarging the first coordinate or the coefficient sets
has not produced an additive obstruction. The coefficient colourings supply enough shifts to absorb every contact
constraint. This is evidence, not a theorem, that the linear golden space is boundary-free.

## 5. What a forcing construction must do next

1. **Defeat the additive problem first.** Forcing a pattern requires that no choice of coefficient colourings makes the
   forbidden-difference `Z/6` colouring of the golden two-distance graph on `A'` solvable. A candidate that passes the
   seconds-long pre-screen for every pattern is closed, as `v^2` now is for these coefficient sets.
2. **Make the shifts rigid.** Homomorphism-like coefficient colourings give one shift per direction. Coefficient sets
   whose unit graph forces many different shifts per direction (rigid, non-Cartesian coefficient graphs) are the
   natural targets.
   - The `MM` recipe (Eisenstein-rigid `M + M`) tests this; see section 4.
3. **Leave the linear space only where the census says edges appear.** The `v^2` census shows that the natural product
   coordinate adds no edges, and THM-4668 shows the same for the degree-192 direction `e` over `K`. Products such as
   `ue` and `ve`, or coefficients outside `K`, are the remaining candidates. Any future coordinate should be census-tested with `golden_vsquare_unit_census_20261017.py`
   adapted to it, before a graph is built.
4. **Gluing at literal shared vertices.** clock_duality's other route, certified incompatible extensions glued at shared
   ports, needs two sources whose *additive* boundary relations already disagree. The pre-screen also evaluates that
   cheaply.

**Connection to the Collatz lane (ANALOGY, no transfer claimed).** In both lanes the "inherently controlled" system is
a product/sum structure in which control is free: the random-merge martingale and the Cartesian additive colouring. In
both, the real object deviates from it only through a sparse coupling: sibling covariances and golden contact edges.
The remaining difficulty is whether the coupling can be made rigid (forcing) or kept tame (cancellation).

## Reproduction

```
python 04-computation/experiments/golden_additive_colourings_20261017.py
python 04-computation/experiments/golden_vsquare_unit_census_20261017.py
python 04-computation/experiments/golden_e_direction_census_20261017.py
python 04-computation/experiments/golden_e_direction_proof_check_20261017.py
python 04-computation/experiments/golden_forcing_search_20261017.py base     # also dM, contact1, dM_C2, MM (contact2 with --no-two)
```
