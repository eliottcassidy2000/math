# Publication audit and reproducibility scope

2026-10-26. Proof status belongs to the reviewed mathematical statements;
finite receipts certify only their named universes. No Lean formalization,
six-colour obstruction, or representative atlas beyond the saved bound is asserted.

## Six-colour boundaries: PASS

The parent reviewed global colour equivariance, the difference-kernel slice,
fresh-colour construction, independent-port hypotheses, and the literal
union/private-interior gluing condition. Counts are separately verified by
[an independent colour-class algorithm](audits/check_six_color_independent.py):
it contracts port blocks, adds inequalities between distinct blocks, evaluates
the chromatic polynomial using independent-set colour classes, and divides by
the falling-factorial colour orbit. It imports no producer functions.
[Receipt](audits/six_color_independent.json).

A second proof review and exhaustive small-graph probes are recorded in
[information/SIX_COLOR_AUDIT.md](information/SIX_COLOR_AUDIT.md). The four
NAE partition orbits are not the four endpoint-frame states. Only the stated
symmetry-compatible tables are classified, not all graph-realizable tables.
The exact conditional unit-distance example remains three-colourable.

## Metric transport: PASS

The independent [metric audit](enumeration/METRIC_AUDIT.md) reviews the
cofactor signs, covariance, arbitrary-network sufficiency and hostile
holonomies. Its Gaussian/barycentric and metric implementation tests unequal
as well as equal lengths and two independent cycle constraints. It realizes
the incompatible maps in three-dimensional octahedral geometry, demonstrating
that planarity is essential. The parent's110592 triangle-network census is a
separate control.

## Tournament information: PASS

[six_color/INFORMATION_AUDIT.md](six_color/INFORMATION_AUDIT.md) independently
reviews I1--I5 and the signed coefficients, checks coefficient word expansions
through8, actual full-response arrays through5, and the alignment hostile.
A wording ambiguity was repaired: one orientation bit is added to the full
response/correlation code; curl alone plus one bit need not determine T.

The parent's [independent audit](audits/INFORMATION_AUDIT.md) uses path-subcube
accumulation and integer Walsh transforms. The producer uses subset dynamic
programming. Both cover all labelled tournaments through6 and agree on the
variance/Fisher counterexamples. The historical source is preserved below a
correction banner, with the mechanism entered in MISTAKES.md. The binary code
length is not an oracle-query complexity claim.

## Complete family enumeration: PASS

The parent reviewed both directions of the alternating criterion, the weight
and phase stabilizer actions, and the deletion owner. The refinement tree
orders degree cells, splits them by invariant neighbour-count vectors, and
individualizes the first nonsingleton cell. Relabelling transports the entire
tree. All minimum-code leaves are retained, so their first vertices form
exactly one automorphism orbit. Keeping that orbit removes root multiplicity
without deduplicating accepted children. This is an isomorphism normal form,
not the minimum over all unrestricted labellings.

An [independent exhaustive audit](audits/check_enumeration_independent.py)
partitions the entire labelled cube through6 using adjacent-vertex
transpositions, without the producer canonicalizer. It checks representatives,
stabilizer mass, and the core by enumerating all subsets and finding maximal
transitive modules. [Receipt](audits/enumeration_independent.json).

The producer materializes191536 classes through9, with147095 new cores at9,
and checks labelled orbit mass at every stage. Normal and optimized order-nine
compressed banks agree byte for byte. At7, all32768 fixed-path tilings and
their literal lines are decoded independently and matched to family path
fibres. The metagraph order is7, not9. The first degree-only implementation
timed out at its five-minute probe; no incomplete census was retained. The
final refinement implementation completes that order.

The coverage inequality comes from a surjection of tournament inflations,
not constant fibre sizes. Its budgets use the inherited reverse-valuation
engine and exact reference counts through128; these are count-only checks
beyond9.

## Reproduction and publication checks

`verify.py` checks a frozen manifest before replaying programs into temporary
outputs and comparing with pre-existing receipts. It never regenerates the
expected result first. `--replay` repeats bounded checks in both Python modes;
`--full` also regenerates the complete order-nine bank in both modes. Source
pins identify inherited inputs.

The separately published norm-seven arithmetic/F1 commit is integrated before
the final navigation check. Documentation checks, local link checks, a clean
worktree and fetch/push target verification complete publication.
