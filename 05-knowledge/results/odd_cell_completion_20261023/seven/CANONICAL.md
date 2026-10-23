# Source-only canonical charts and literal Hamiltonian ports

23 October 2026. **PROVED**, for a supplied named equitable partition
whose cells have size1,3,5 or7, with exact owned block data. **FINITE-EXACT**
controls below. No candidate target or target-class inventory is supplied
to the canonicalization algorithm.

This is the next step after the incoming
[partial-bijection comparison rule](../../odd_cell_gluing_20261023/README.md).
That rule compares TWO supplied systems. Here the new one-sided kernel
theorems let us build a complete canonical code directly from ONE source.
The internal templates are only the individual at-most-seven-vertex
cells, obtained by direct permutation, not an inventory of full-size
tournaments. This does not compress the input to a scalar observer:
actual block entries and their shared owner sets remain supplied.

## 1. Fix the small templates, then synchronize one tree

Keep the ordered cell names. For each cell, enumerate all internal vertex
orders and choose the lexicographically least tournament matrix as its
template. Retain every actual chart to that template. Their number is its
internal automorphism order, at most21.

Form the support graph of nonuniform cross blocks. By equitability and
coprime sizes, each connected component has only one non-singleton cell
size. Choose its least named cell as root and a deterministic breadth-first
spanning tree with neighbours in cell-name order. These choices are
invariant under any relabelling preserving the cell names.

For EACH root chart, proceed outward along the tree. Once the parent
chart is fixed, consider every internal chart of the child and serialize
the actual oriented parent-to-child block in those coordinates. Choose
the chart whose block word is least.

**Lemma.** This child minimizer is unique.

If two distinct child charts gave the same block word, their quotient
would be a nonidentity automorphism of the child fixing its entire
one-sided block. The [seven-cell kernel theorem](PROOF.md) and the
[five/prime-rotation theorem](../five/README.md) exclude that for every
nonuniform equitable block in this carrier. The C3 case is immediate.
Thus ordinary lexicographic comparison replaces a choice of unrelated
pairwise witnesses; each shared cell obtains exactly ONE chart.

After all tree charts are set, serialize the full component tournament,
including every NON-TREE block, in cell-name order. Minimize this full code
over the at-most21 root charts. The greedy tree code is not asserted to
be the globally lexicographically smallest matrix over all conceivable
cell charts. What matters is that the deterministic construction commutes
with every actual cell-preserving isomorphism, and the final code retains
the whole tournament.

Distinct support components interact only through uniform blocks. Hence
their chart choices are independent; canonicalize each separately and
then assemble them with those retained uniform quotient bits. There is
no product search over all component root charts.

## 2. Completeness and the exact residual symmetry fibre

**Theorem.** Two tournaments with the supplied named equitable cells
receive the same code if and only if there is a cell-preserving
isomorphism between them.

An isomorphism carries the small template-chart sets to each other and
preserves the support, its deterministic tree, and every oriented block
word. It therefore carries every greedy root trial to the corresponding
root trial in the other source. Their minima agree.

Conversely, equality of the final full matrix codes provides an explicit
vertex correspondence, matching the selected actual owners at each
canonical coordinate. Every tournament arc agrees, including the
non-tree and inter-component blocks. This is the required isomorphism.

The number of root trials attaining the minimum in one support component
is exactly the order of its cell-preserving automorphism group G.
Composing a minimizing chart with an automorphism gives another minimizing
trial. The action is free on root charts, because the root restriction is
faithful and the other charts propagate uniquely. Conversely, any two
minimizing full chart systems compose to an actual automorphism. Thus the
minimizing trials form one G-torsor, not merely a list of matching scalar
statistics.

The procedure may choose one representative chart by input labels when
G is nontrivial; this is a gauge choice. It does NOT intrinsically name a
vertex moved by G. The code itself is canonical, and the entire minimizing
chart fibre remains available. Across independent components, the whole
fibre is a factored product, which need not be expanded.

This gives a polynomial procedure in the input size for the SPECIFIED
named bounded-cell carrier: at most7! internal trials per cell, at most21
root charts per component, and at most21 child-chart comparisons per tree
step. Full matrix output itself has quadratic size. It is not a polynomial
canonicalization claim for unrestricted tournament input, where finding
such a partition with bounded cells remains an additional problem.

## 3. Gcd savings for compatible attachments

Restriction of G to EVERY cell in a connected informative component is
injective, because a fixed chart propagates uniquely across the support.
Therefore

    |G| divides gcd_i |Aut(T[D_i])|.                         (1)

In particular a seven-cell component containing both an internal
group-of-three type and a group-of-seven type is automatically rigid,
regardless of the particular nonuniform blocks. The actual matching
attachment in the checker gives a concrete rigid order14 example.
Components with some non-Paley seven-cell have |G| at most3 or7 and need
at most one suitably chosen moved owner for pointwise chart uniqueness.
Only a component consisting entirely of Paley-seven cells can retain
|G|=21 and genuinely require two owners. It may still have a smaller
global group because of its shared attachments.

Once the canonicalizer fibre has been found, its exact root-owner base
cost is0,1,1,2 for |G|=1,3,7,21 respectively, as proved in PROOF.md.
The calculation no longer needs a supplied candidate target. This remains
distinct from ordinary refinement: a rigid unmarked component can have
nondiscrete equitable cells, and the canonicalizer uses complete blocks
that refinement had compressed.

## 4. Retain the actual Hamiltonian port without expanding all components

Suppose an ACTUAL directed Hamiltonian path P is supplied. For each
support component, take its retained canonicalizer fibre and list, in
the order of occurrence along P, the canonical coordinates of that
component's vertices. Minimize this coordinate word within that fibre.
Because every vertex occurs exactly once, different charts give different
words; the minimizing chart is unique.

This minimization can be done independently in each support component.
The positions at which P visits each named component are already fixed.
Comparing two complete encoded paths first differs at a position belonging
to one component, so the lexicographic choice there is precisely the
component's own restricted-word choice. No enumeration of the factored
product of all component chart fibres is necessary.

The resulting pair

    (canonical source matrix, canonical actual path word)

is a complete cell-preserving source/port code. Matching two such codes
constructs an isomorphism carrying the SAME actual path to the other one.
Entrywise arc inspection then transports its actual complement partner:

    f(K_P(T)) = K_{f(P)}(f(T)).

This is not an endpoint-only decoder, and it never constructs a new
Hamiltonian path from a quotient. Native converse and integer mod8 guards
remain separate acceptance coordinates. Under converse the path must be
reversed and retained owner roles transported; no claim identifies a
generic matrix-word reversal with that operation.

A full actual Hamiltonian path already gives a trivial vertex order in
which its complete labelled tournament can be encoded. The new part here
is the SOURCE-CLASS chart canonicalizer without a full-size target
inventory or product search across components. The path tie-break is a
lossless shared-port transport check, not a new hardness or speed claim
for canonicalizing a fully supplied labelled path.

## 5. Order31 tests and precise universes

The source is the same literal D(D(P7)) as the inherited endpoint packet.
The canonicalizer is run on all930 ordered distinct two-mark partitions.
It receives only that source, its cell names and actual vertices. It
recovers the inherited group census through minimum-fibre sizes:

| Largest stable cell | Minimizing chart count | Ordered mark pairs |
|---:|---:|---:|
|1|1|504|
|3|1|126|
|3|3|168|
|7|1|42|
|7|3|84|
|7|21|6|

For each of the930 cases, a separately generated global relabelling of
all31 vertices changes both owner names and the actual cell subsets.
The canonical code and fibre counts stay equal, and matching canonical
coordinates constructs an entrywise verified tournament isomorphism.

Two different path universes are retained explicitly:

- In all132 seven-cell marked cases, an actual insertion path is transported
  with the two arbitrary marks retained. Those marks are NOT asserted to
  be that path's endpoints.
- The inherited933 literal path records, comprising932 distinct actual
  Hamiltonian paths, are separately replayed with THEIR ACTUAL ENDPOINTS
  as marks. Every canonical port is invariant under a fresh relabelling;
  its constructed map carries the full ordered path and its complement.
  Every computed complement matches the frozen literal target-outmask
  record. The inherited native mod8 receipt fingerprint is retained.

The930 mark pairs are not promoted to930 path-feasible endpoint pairs.
The933 path records are finite controls, not a full Hamiltonian-path census
of the order31 source. The canonicalization theorem itself does not depend
on either finite count.

## Reproduction

    python -B 05-knowledge/results/odd_cell_completion_20261023/seven/canonical.py
    python -O -B 05-knowledge/results/odd_cell_completion_20261023/seven/canonical.py

Use `--write` to regenerate `canonical.json`. The explicit failure checks
survive optimized Python. Source hashes normalize CRLF to LF. No git
mutation or candidate-target inventory is required by this checker.
