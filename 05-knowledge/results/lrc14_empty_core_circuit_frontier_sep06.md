# Complete carrier circuits: an arithmetic-progression hostile and two finite circuit types

**Status: FINITE-EXACT + independently verified controls; REFUTED additive
compression; universal three-direction classification now PROVED in the
linked colored-diamond note.** No theorem ID is
reserved. None of the height observations below proves an all-height bound.

## Source, inheritance, and consequence object

Use the primitive sorted distinct positive odd ternary-unit triples
`w=(a,b,c)` and the complete raw support of
[THM-4414 — six-separated contact capacity collapse](../../01-canon/theorems/THM-4414-lrc14-six-separated-contact-capacity-collapse.md):

```text
Lambda(w)={C in Z^3: C.w=0, C_i!=0 mod3,
                         14|C_i|<3(sum(w)-w_i) for every i},
E_i=sum_(C in Lambda) min(3/(7c),
                  [3(sum(w)-w_i)-14|C_i|]/(14w_jw_k)).
```

The target is `min_i E_i<=6/77`. The physical mass is the sum of the
pointwise minimum and can be strictly smaller. All scans here retain the
three separate network sums. A direction identifies the primitive integer
vector with its negative; raw multiples are retained as distinct carriers.

The closest proved mechanism is
[THM-4422 — projection deficit and Beatty-row reduction](../../01-canon/theorems/THM-4422-lrc14-projection-deficit-and-beatty-row-reduction.md).
Its sufficient count gate is `N<=2c/11`, and its first dense carrier circuit
is `(19,23,29)`. The incoming
[two-ray result](lrc14_two_ray_overnight_hexagon_sep05.md) makes every remaining
potential failure have at least three directions. The current session's
[three-through-six result](lrc14_empty_core_three_ray_sep06.md) moves that
boundary farther. These inherited results prompted the circuit probe, not
an assumption that every three-direction hull has the same additive shape.

The methodological source is the retained empty core in Heule–Scheucher,
[arXiv:2403.00737v1, Section 4.1](https://arxiv.org/pdf/2403.00737).
This paper concerns point sets in general position. Our centrally symmetric
integer support may have collinear multiples and omits one residue class.
No empty-polygon conclusion is transferred to this different object.

## Exact refutation and strongest survivor

The proposed shortcut was: if the complete support has at least three
primitive directions, it contains independent raw vectors `u,v` with `u+v`
also live. This is **REFUTED** by

```text
w=(41,47,49),
Lambda(w)=+/-{(11,5,-14),(14,-7,-5),(17,-19,4)}.
```

It is the first hostile in the complete height-499 universe, ordered by
`(c,a,b)`. A separate direct-box enumeration confirms completeness of this
six-point support and minimality through height 79. All choices of signed,
independent raw pairs fail the proposed additive test. Nevertheless,

```text
(11,5,-14)+(17,-19,4)=2(14,-7,-5),
(E_1,E_2,E_3)=(276/16121,376/14063,10/287).
```

The first failed implication is **three directions to unweighted circuit**.
The omitted coordinate is the coefficient two in the primitive dependence.
The strongest survivor is a weighted integer circuit, which is always
available for three directions in a plane. The repair is to retain its
primitive coefficient vector and every raw multiple. This witness refutes
the compression argument, not the network target.

## Complete height-499 observations

The finite universe contains **753,853** eligible triples, of which
**710,364** have at least two directions. In that multi-direction subset:

| Statistic | Exact maximum | Witness |
|---|---:|---|
| `max_i E_i` | `18/301` | `(5,37,43)` |
| `min_i E_i` | `12/343` | `(35,47,49)` |
| `(E_1+E_2+E_3)/3` | `1088/24087` | `(11,31,37)` |
| `N/c` | `6/29` | `(19,23,29)` |

Every projection is strictly below `6/77` in this subset. The count gate
fails only at `(19,23,29)`, whose six carriers are paid by the actual roof
deficits. Neither statement is extrapolated to arbitrary height.

There are **7,646** exactly-three-direction rows. Their primitive dependence
coefficient multisets are precisely:

| Absolute coefficients after common-gcd division | Rows | First witness |
|---|---:|---|
| `(1,1,1)` | 5,428 | `(19,23,29)` |
| `(1,1,2)` | 2,218 | `(41,47,49)` |

There are **41,760** rows with at least three directions but no independent
raw additive circuit. The two circuit types above classify only the
exactly-three-direction rows in this finite box. They are not a claim about
every triple selected from a larger support.

## A broader colored convex model

The new representation is an integer lattice with an invisible index-three
sublattice. For the cheap abstract test, take `Z^2`, make a point live when
its first coordinate is nonzero modulo three, and take the closed symmetric
convex hull of three primitive live vectors. Reject the hull if it contains
a fourth primitive live direction. Including the closed hull is legitimate
for testing a candidate subset of an open convex body: the hull of finitely
many interior points remains inside that body.

The complete abstract universe with representatives `1<=x<=5`,
`-5<=y<=5`, `gcd(x,y)=1`, `x!=0 mod3` has 29 surviving three-direction
hulls: ten of type `(1,1,1)` and nineteen of type `(1,1,2)`. Every possible
extra live point has a primitive representative in the same coordinate box,
so the finite rejection test is complete. This removes special speed-strip
geometry from the experiment, while retaining convexity, primitivity, and
the owner color.

The initial question is now **PROVED** by the independently audited
[colored-diamond theorem](empty_core_colored_diamond_sep06.md): a diamond on
two primitive live vectors of lattice determinant at least three contains
four live primitive directions. It implies exactly these two circuit types,
for any proper invisible subgroup. Concurrent THM-4431 owns the shared
classification namespace; this session adds an independent complete-fan proof.

Here is the exact carrier identification needed to apply that theorem.
Let `L={C in Z^3:C.w=0}` and

```text
Gamma={C in L:w_1 C_1=w_2 C_2=w_3 C_3 mod3}, H=3L.
```

Reduction of `L` modulo three is onto the two-dimensional kernel of `w`
over `F_3`: for a lift `z`, subtract `(w.z)b` where `w.b=1`; the correction
is divisible by three. The equal-owner condition is a one-dimensional line.
Thus `[L:Gamma]=[Gamma:H]=3`. The complete raw support is exactly
`(Gamma minus H) intersect K`, where `K` is the open symmetric convex roof
region. Content division of a live vector preserves this set: its content
is a unit modulo three and hence its primitive reduction stays on the same
nonzero owner line. Consequently primitive live vectors in `Gamma` are
exactly the ambient primitive carrier directions. The colored-diamond
theorem applies with no change in direction count. Its lattice determinant
is `det_xy/(3c)`, retaining the inherited owner factor three.

This connection preserves direction incidence and eligible residues. Forgetting the roof
lengths loses network weights; these must return as a sidecar before any
projection conclusion. The abstract body probe therefore cannot by itself
prove the `6/77` inequality.

## Reproduction and verification

```bash
c++ -std=c++17 -O2 04-computation/lrc14_empty_core_census_sep06.cpp -o /tmp/lrc14-empty-core-census-sep06
/tmp/lrc14-empty-core-census-sep06 499
python3 -B 04-computation/lrc14_empty_core_census_audit_sep06.py
python3 -B -O 04-computation/lrc14_empty_core_census_audit_sep06.py
```

The [primary C++ enumerator](../../04-computation/lrc14_empty_core_census_sep06.cpp)
solves one congruence per first-coordinate row and accumulates all roof
terms with integer denominator `14abc`. The height range is explicitly
restricted to `5..999`, within the integer arithmetic bound. Its
[height-499 output](lrc14_empty_core_census_sep06.out) is the finite evidence.

The [separate Python audit](../../04-computation/lrc14_empty_core_census_audit_sep06.py)
uses direct integer boxes and rational projection arithmetic on every one
of the 2,910 eligible triples through height 79. It matches all leaders,
direction and circuit counts, and all 266 low additive hostiles. Nine named
controls, including equality, a norm-four every-projection failure, the two
multi-ray hostiles, and a higher-height row, are checked against the explicitly
imported literal six-sheet interval engine from the incoming one-ray verifier.
The [audit output](lrc14_empty_core_census_audit_sep06.out) has 4,534 explicit
gates and is byte-identical under `python -O`. The full height-499 scan is
not claimed to have an independent full replay.

Frozen raw-byte SHA256:

```text
C++ source 167ee7e0219ba5f40ff32f62fbb5470e0ae8cf4b8f3dcff67722ad4261e8a0e7
C++ output ccbe3cc4aa95637d5534cce04e7bf7bafdcf09626addaab320e4914d76f84cec
audit source 9c6a8ced10ccc8c9342fe611d09a67051fc7ab5396c8f9eb09b697bc686a02c1
audit output 6f357a617762868af063705fdfb869e5f39b42ab782131554e4a9fa95d78f5d3
```

The session board compares these objects with the independent short-relation
slices and trinomial carries. The shared research mechanism is retaining the
coefficient that a tempting quotient discards; no mathematical equivalence
between those problems is asserted.
