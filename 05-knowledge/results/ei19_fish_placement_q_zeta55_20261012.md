# EI19 and the Hochberg-O'Donnell fish in Q(zeta55): EI19 excluded by a forced angle (PROVED), the fish is a one-parameter family (OPEN)

2026-10-12, opus-2026-10-12-S24 (second addendum). **PROVED:** the 19-vertex
Exoo-Ismailescu graph EI19 has no faithful unit-distance placement in
`L = Q(zeta55)`, nor in any field lacking `sqrt(-2)` or lacking `sqrt(-35)`;
the direction-class structure of all three graphs (parallelogram lemma plus
exact integer linear algebra). **NUMERICAL:** the EI19 census (32 real faithful
realizations, all with the forced angles), the EI17 census revisited (12), the
fish family samples and the fish L-unit scan. **OPEN:** whether the fish's
one-parameter realization family has a point over `L`. **CITED:** House of
Graphs 51376 and 51404 (adjacency lists, API read 2026-10-12); MathWorld for
the fish's provenance (Hochberg-O'Donnell, Geombinatorics 5 (1996) 137-141);
arXiv:2607.19995 for EI19's base field (consistent with the forced angles found
here; its cap polynomial checked only as stated below). Self-audited; no
independent peer audit.

## The method: direction classes and forced triangles

In a faithful unit-distance realization every unit 4-cycle is a parallelogram
([EI17 note, section 1](ei17_placement_q_zeta55_20261012.md)). Identifying
opposite sides of all unit 4-cycles groups the edges into **direction classes**
`u_1, ..., u_m` (unit vectors), and the cycle space of the graph gives exact
integer linear relations among them. A relation of **support three**,
`n_a u_a + n_b u_b + n_c u_c = 0`, is a rigid triangle of side lengths
`|n_a|, |n_b|, |n_c|`: it forces the ratios `u_b/u_a`, `u_c/u_a` to specific
unit-modulus algebraic numbers. A placement in a field `F` (up to isometry)
needs every direction ratio in `F`, so a forced ratio outside `F` excludes the
graph from `F` **for every realization at once**, without any census.

| graph | vertices, edges | unit 4-cycles | classes | relations | free angles | smallest relation support |
|---|---|---:|---:|---:|---:|---|
| EI17 (HoG 51375) | 17, 31 | 9 | 13 | 6 | 6 | 5 (no forced triangle) |
| EI19 (HoG 51376) | 19, 35 | 11 | 13 | 6 | 6 | **3 (two forced triangles)** |
| fish (HoG 51404) | 23, 42 | 9 | 24 | 11 | 12 | 4 (one relation; a rectangle) |

The direction formulation also removes the sign branches of apex
constructions and most collapsed solutions: random starts over the free
angles converge almost always, which makes the realization censuses cheap
(EI17: 3,000 starts in 45 s give the same twelve realizations as 20,000).

## 1. EI19: two forced angles, hence no placement in Q(zeta55)

**Theorem.** Let `G` be EI19 with the House of Graphs labelling `0..18`. In
every faithful unit-distance realization of `G` in the plane,

    (v6 - v2) + (3/2) [(v13 - v6) + (v2 - v18)] = 0,

so the unit vectors `(v13 - v6)/(v6 - v2)` and `(v2 - v18)/(v6 - v2)` sum to
`-2/3` and therefore equal `(-1 +- 2 sqrt(-2))/3`; and, with the classes
`u4` (edge `(2,3)` and its parallels) and `u9` (edge `(3,16)` and parallels),
`u5 = 3(u4 + u9)`, so `u4/u5`, `u9/u5` sum to `1/3` and equal
`(1 +- sqrt(-35))/6`. Consequently a field containing all direction ratios of
any faithful realization of EI19 contains `sqrt(-2)` and `sqrt(-35)`. Since
`Q(zeta55)` contains neither (its quadratic subfields are `Q(sqrt5)`,
`Q(sqrt-11)`, `Q(sqrt-55)`; the split prime 661 has `-2` as a non-residue and
the script's witness prime does the same for `-35`), **EI19 has no faithful
placement in Q(zeta55)**, nor in any `Q(zeta_n)` with `8` not dividing `n`.

*Proof.* The eleven unit 4-cycles of EI19 are parallelograms, giving 13
direction classes; the cycle space yields six independent relations, among
which (after exact row reduction) the two of support three displayed; the
class of `(2,6)` is a single edge, the class of `(6,13)` contains `(6,13)`,
`(7,17)`, `(11,18)`, and the class of `(2,18)` contains `(2,18)`, `(3,17)`,
`(13,16)` with the stated orientation. Two unit vectors with a prescribed
sum of modulus `2/3` (resp. `1/3`) have real part `-1/3` (resp. `1/6`)
relative to the sum's direction. A square root of `-2` in `Q(zeta55)` would
reduce, at a prime above `661` (which splits completely), to a square root of
`-2` in `F_661`, which does not exist. ∎

*Census (NUMERICAL).* Random starts over the six free angles converge in
every case; 3,000 starts give **32 non-congruent real faithful realizations**
(hit 49-172 times each; all isolated; all 17 pairwise non-adjacent distances
above `0.044`), every one with the forced ratios `-1/3 +- 0.9428 i` and
`1/6 +- 0.9860 i`. Eight realizations re-solved in vertex space from
perturbed starts (35 unit-edge equations only, no direction machinery)
return to realizations satisfying the relation to `1e-15`. The seven classes
`{4,5,9,10,11,12,13}` are the paper's "exact base": their pairwise ratios are
identified by LLL in `Q(sqrt2, sqrt5, sqrt7)(i)` in all 32 realizations, in
agreement with arXiv:2607.19995's base field; the forced values
`(-1 +- 2 sqrt(-2))/3` and `(1 +- sqrt(-35))/6` are its `sqrt2` and
`sqrt35 = sqrt5 sqrt7`. The paper's degree-12 cap polynomial (irreducible, eight
real roots `0.0034, ..., 0.9196`) is corroborated here, unlike its EI17
counterpart: in all 32 realizations the cosine of the angle between the
directions of classes `8` and `10`, and of classes `8` and `11`, is one of
those eight roots (VERIFIED to `1e-9`); its Galois-group claims were not
audited.

## 2. The fish: a one-parameter family, no forced angle, no L-point found

The fish (23 vertices, 42 = 2*23 - 4 edges; hub `22` adjacent to `0..7`,
`20` and `21` of degree six, `18` and `19` of degree five) has nine unit
4-cycles, 24 direction classes and 11 relations: **12 free angles against 11
equations**, a one-dimensional realization family, as the edge count
predicts. The only relation of support at most four is
`u2 - u4 + u6 - u8 = 0`: the spokes to `0,1,2,3` satisfy `v0 + v2 = v1 + v3`
(hub at the origin), an inscribed parallelogram, hence a rectangle: spokes
`0,2` and `1,3` are antipodal in every faithful realization (verified on
every sample). There is no support-three relation, so no angle is forced and
the field is not fixed by the graph: along the family the coordinates are
transcendental in the parameter, and an `L`-placement exists iff the family's
algebraic curve has an `L`-point. That is a Diophantine question this note
does not settle.

*Evidence (NUMERICAL).* Fixing one free direction to each of 117 units of `L`
(the 110 roots of unity, five Mersenne directions `(2-zeta11^j)/(2-zeta11^-j)`,
the thirds unit and its conjugate) and solving the remaining eleven equations
from 60 random starts gives 0-10 faithful realizations per parameter; at
every one of them, LLL on the `zeta55` power basis at 200 digits finds no
further direction in `L` (shortest vectors at the random height `1e4`). So no
`L`-point of the family passes through those parameter values with small
height. A decisive answer needs the family's curve: its genus and
`L`-rational points, or a forced non-`L` invariant along the whole family.

## 3. EI17 revisited by the same method

EI17 has no relation of support at most four, so no angle is forced
combinatorially, which is why its field is hard (and why the published
"S20" certificate could be wrong without any contradiction from structure;
see the [EI17 note](ei17_placement_q_zeta55_20261012.md)). In direction space
the census is clean: 3,000 random starts converge 1,902 times to 13 distinct
solutions, 12 of them faithful, all isolated, each hit 64-312 times; 20,000
starts give the same twelve. The earlier branch-grid census (eleven) missed
one. The exact elimination of the six-angle Laurent system is the remaining
route to EI17's field.

## 4. Connection ledger

| Source object | Map and preserved predicate | Lost data / retained sidecar |
|---|---|---|
| Unit 4-cycles | parallelogram identification of directions | collapsed (bow-tie) solutions excluded by faithfulness |
| Cycle space | integer relations among classes | none (exact) |
| Support-three relation | rigid triangle of unit vectors | the sign of the forced imaginary part |
| Forced ratio | square-root field obstruction | realization-independent; needs no census |
| Fish relations | one-dimensional family | the field of generic points |
| L-unit scan | 117 parameters, LLL at 200 digits | heights above `1e4`, parameters outside the bank |

## 5. Reproduction

    python3 -B 04-computation/experiments/ei19_fish_placement_q_zeta55_20261012.py
    python3 -B 04-computation/experiments/ei19_fish_placement_q_zeta55_20261012.py --lscan

[Script](../../04-computation/experiments/ei19_fish_placement_q_zeta55_20261012.py),
[output](ei19_fish_placement_q_zeta55_20261012.json). The script recomputes
the structures exactly (sympy row reduction over `Q`), searches small-support
relations, verifies every reconstructed realization on all its edges,
confirms the forced ratios on every EI19 realization and in vertex space,
computes the non-residue witness primes, samples the fish family with its
local degree of freedom, and (with `--lscan`) repeats the L-unit parameter
scan; the scan's declared run is summarised above from the scratch logs.
