# Sources and recovered connections

## Primary external inputs

Accessed2026-10-30. The author-maintained revision18 of
[Small Ramsey Numbers](https://www.cs.rit.edu/~spr/ElJC/ejcram18.pdf),
dated2026-04-24, gives lower bounds45 and61 in TableXIb and upper
bounds57 and91 in section6.1d for the two requested numbers.
The [revision page](https://www.cs.rit.edu/~spr/ElJC/eline.html) identifies
the version. Exact values are not supplied there.

[Luis Boza, New Upper Bounds for the Classical Ramsey Numbers
R(4,4,4), R(3,4,5) and R(3,3,6)](https://arxiv.org/html/2603.10851v1),
2026-03-11: introduction for the small exact inputs and recursive
upper bound; Theorem2.4 for91. The saturation and triangle-divisibility
method is credited to this source. Our single-port group restriction
and retained-variable calculations are derived in PROOFS.md, not
attributed to Boza's paper.

[Brendan McKay's complete Ramsey graph banks](https://users.cecs.anu.edu.au/~bdm/data/ramsey.html)
provide the completeness input and these downloaded factual data:

| File | Universe | Direct source |
|---|---|---|
| `inputs/r35_13.g6` | one `(3,5;13)` graph | [graph6](https://users.cecs.anu.edu.au/~bdm/data/r35_13.g6) |
| `inputs/r36_17.g6` | seven `(3,6;17)` graphs | [graph6](https://users.cecs.anu.edu.au/~bdm/data/r36_17.g6) |
| `inputs/r44_17.g6` | one `(4,4;17)` graph | [graph6](https://users.cecs.anu.edu.au/~bdm/data/r44_17.g6) |

The SHA256 digests are in verification.json and manifest.json. Only the
graph data, not copies of the articles, are stored in this packet.

[Wei Li's1987 thesis abstract](https://repository.rit.edu/theses/6642/)
records symmetry-assisted searches producing `R(3,3,5)>=45`, including
13 nonisomorphic cyclic44-vertex examples. This is relevant precedent
for orbit compression. It does not justify restricting an upper-bound
proof to symmetric colorings. We used the repository abstract, not an
unread full-text proof or the abstract's broad complexity wording.

A web search surfaced marketing text suggesting stronger lower bounds,
but the accessible text was internally inconsistent about graph order
versus Ramsey lower bound and supplied no checkable construction. Those
figures were not adopted. The bounds reported here are explicitly those
in the checked primary survey, not an assertion that every later claim
on the web has been adjudicated.

[Przybocki, Mackey, Heule and Subercaseaux, Doubly Saturated Ramsey Graphs](https://arxiv.org/abs/2604.21187),
2026-04-23, defines the edge-addition/edge-deletion saturation property
used in P3. Its abstract was checked; the general infinite-family proof
is not imported. The local H13 saturation counts here were computed
directly.

## Repository inheritance

The pinned starting revision is
`98f961a7aa142a739403de1ba343720371d6125c`. Recovery used bounded topic
searches, exact statements and selected proof files, not a chronological
log as authority. The current correction ledger was searched before
using historical claims. One overly broad grep reached generated JSON
indexes; that output was discarded and the search narrowed to Markdown.

| Source and status | Mechanism carried forward | Boundary |
|---|---|---|
| [Collatz Catalan Ramsey synthesis](../collatz_catalan_ramsey_20260927.md), historical | Paley17 and the distinction between algebraic constructions and count coincidences | Its unguarded carry/value-order claims are not inherited; see the [repair](../affine_tree_recurrences_20261028/AUDIT.md) |
| [Binary dependence and free bits](../unifying_binary_dependence_free_bits_7_15_56_20261028/README.md), proved finite mechanisms | `R(3,3)=6`, binary cycle/cut carriers, and marked quotient losses | `56=A000568(6)` is not a map to56 Ramsey vertices |
| [57/27/16 atlas](../number_bridges_57_27_16_20261028/ATLAS.md), source map |57 as a Singer plane size and16 as a binary field carrier | Norm/aliquot identities for57 do not construct a Ramsey coloring |
| [Singer Heisenberg E7 obstruction](../session35_dimensions_dualities_20261023/SINGER_HEISENBERG_E7_OBSTRUCTION.md), proved scoped obstruction | Compare the actual restricted group action and its fixed points | A56-dimensional augmentation representation is not a56-point coloring |
| [THM-4624 finite unit-circle affine rigidity](../../../01-canon/theorems/THM-4624-finite-beckman-quarles-for-unit-circle-cayley-graphs-affine-rigidity-except-clebsch-hamming-and-ternary-frame.md), scoped canon |GF(16), the binary five-frame, Clebsch symmetry1920 | General all-field affine rigidity remains outside the finite claim |
| [THM-4701 Kuratowski Petersen unit-circle graphs](../../../01-canon/theorems/THM-4701-kuratowski-graphs-the-petersen-graph-and-the-finite-unit-circle-graphs.md), scoped canon |Petersen as the Clebsch distance-two layer | No Euclidean unit-distance theorem follows from this finite graph bridge |
| [Regular parent symmetry](../regular_parent_symmetry_20261029/PROOFS.md), proved and exact |Preserve actual attachment actions;81/27 hostile control |Cyclic triple sets alone do not determine the group |
| [Odd-eighth owner families](../odd_eighth_owner_families_20261030/PROOFS.md), proved scoped families |A smaller carrier must be faithful to the operations that remain legal |No Ramsey-to-LRC implication is asserted |
| [Criterion C and Ramsey saturation reflection](../../../07-reflections/criterion-C-formalized-part-A-reduced-to-realization-plus-paley-QR-resonance-klein-S204.md), historical section3 only |Recover the double-saturation operation and test it directly on H13 |Its Paley13 has degree6; the current H13 has degree4. No historical claim about completing LRC is inherited |

## Explicit bridge tests

**GF(16) to Ramsey17.** Source: additive GF(16) and its multiplicative
cube cosets. Target: a three-coloring of K16. Map: color xy by the coset
of x+y. Preserved predicate: no zero sum of three elements in one coset,
equivalently no monochromatic triangle. Forgetting two colors loses
the global partition, enlarging symmetry160 to1920. Sidecar: the three
named cosets. Cheapest decisive test: all560 vertex triples, plus exact
color-preserving automorphisms.

**Saturated56 to one rooted port.** Source: a full avoiding coloring
and its root stabilizer. Target: Aut(H13). Map: restriction. Preserved
predicate: equality of group elements, by P4 injectivity. Destroyed data:
the rest of the coloring is not reconstructed from H13. Sidecar: the
full coloring is still needed to decide which local actions extend.
Decisive obstruction:13-cycle restrictions force three blue fixed points.

**Singer57 to the new symmetry method.** The prior source has a
C57-action on a56-dimensional augmentation module; the Ramsey target
has a permutation action on56 vertices. No object-identifying map is
claimed. What transfers is the proof operation: restrict the proposed
action, compute fixed points, and reject an extension when those fixed
points violate the target predicate. The numerical56 alone is insufficient.

**Tournament enumeration to Ramsey growth.** Source and target both
admit a new-vertex attachment operation, but their edge alphabets and
legal predicates differ. The shared map is deletion of the marked new
vertex in the corresponding groupoids. It preserves automorphism-weighted
mass and yields the rooted orbit identity in RECURRENCES.md. Losing the
mark loses vertex-orbit multiplicity. The81/27 pair remains a control
against replacing the action by a cheaper triple statistic.

**Collatz/LRC analogy.** The precise survivor is the warning that an
average does not control each designated target. The nonuniform local
matrix supplies a finite, exact instance of that failure. No map from
these Ramsey matrices to a Collatz descent or an LRC witness has been
constructed, and neither open problem is advanced by a claimed implication.

## Incoming checkpoints

During close-out, `origin/main` advanced through043b760966. The checkpoint
summaries for [Cayley-Dickson integral carriers](../cayley_dickson_integral_20261030/README.md)
at7f70dcd601 and [paid returns and core transplants](../return_kernel_core_transplants_20261030/README.md)
at043b760966 were read as new research signal. They are not proof
dependencies of P1-P7. The former distinguishes a minimal carrier under
named future operations from a current-response quotient; accordingly,
our67/1065 orbit inventories are not labelled minimal closed carriers.
The latter keeps entry-dependent returns in a changing-context argument;
the related lesson here is to preserve shared-neighbor incidence rather
than replace each count by its mean. The source manifest pins these two
read summaries separately from the original inheritance revision.
