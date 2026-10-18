# Coupled residue branches and the support of the 714-point construction

**2026-10-18. PROVED:** the source-marked residue compiler, its complete
residual colouring criterion, and the alias-matching criterion below.
**FINITE-EXACT:** both supplied coordinate graphs retain all 570 original
six-colour boundary patterns. **OPEN:** a geometric equality fragment
with an equal colouring and no separating colouring, and the plane six/seven
question. A failed restricted repair chart is never reported as full UNSAT.

## 1. Inheritance, the map, and the retained coordinates

The closest mechanisms are G1--G4 in
[source orientation and six-coordinate contacts](source_orientation_boundary_contexts_20261018k.md),
the [exact shade fibres](snark_boundary_transfer_20261018a.md), and the
[source-preserving palette/alias compiler](shade_difference_support_20261018b.md).
The canonical hostile is a triangle whose three addresses inherit one
two-colour pair: all such charts fail although a full three-colouring exists.
The corrected near miss is to replace a joint edge condition by its two
modular marginals. The least-used sidecar is the original directed pair
behind every imposed difference, including aliases of physical vertices.

The live board is: literal source ports, modular cycle sums, joint residual
edge clauses, physical aliases, and complete boundary support. The map here
changes the representation of colour variables, not the input geometry or
the source. It preserves exactly the branch equations, pins, and proper
six-colourings. Dropping the residual clauses or a source alias destroys
that equivalence.

This packet uses arithmetic colour residues

\[
c\longleftrightarrow(c\bmod2,c\bmod3),\qquad
c=3b+4t\pmod6.                                      \tag{1}
\]

These are not the three nonzero labels of an \(\mathbb F_2^2\) snark flow,
and do not automatically identify with an arbitrary palette matching.
The coordinate \(c\bmod3\) has the specific matching
\(\{\{0,3\},\{1,4\},\{2,5\}\}\); other chosen charts require their retained
relabelling. For example, the consecutive-pair label \(\lfloor c/2\rfloor\)
is not additive: the labels of 1 and 1 sum to 0 modulo 3, but the label of
their sum 2 is 1. The chart in section 5 retains its pairing sidecar explicitly.

## 2. Exact branch reduction, with the remaining obstruction

Let \(G=(V,E)\) be a finite simple graph. A branch consists of actual equations
\(c_v-c_u=a\pmod p\), for \(p=2\) or 3, and literal pins \(c_v=c_0\).
Introduce a new ground vertex with value zero; each pin supplies one edge
from ground to its vertex in both modular systems.

For each modulus, propagate differences along a spanning forest. Either an
actual closed walk has nonzero sum, certifying inconsistency, or every
solution has the unique form

\[
c_v\bmod p=o_p(v)+z_{p,C_p(v)}.                       \tag{2}
\]

Here \(o_p(v)\) is the retained path sum, \(C_p(v)\) its component, each
ungrounded component has one free phase in \(\mathbb Z/p\), and the grounded
component's phase is zero. Isolated vertices are included. Thus (1)--(2)
are a bijection between phase assignments and all solutions of the branch
equations and pins. Different component partitions in the two moduli are
allowed.

For every *original* edge \(uv\), retain

\[
\neg\bigl[(o_2(u)+z_{2,C_2(u)}=o_2(v)+z_{2,C_2(v)})
\ \land\ (o_3(u)+z_{3,C_3(u)}=o_3(v)+z_{3,C_3(v)})\bigr]. \tag{3}
\]

**Residual criterion (PROVED).** A branch extends to a proper six-colouring
iff both modular systems are consistent and their free phases satisfy
every clause (3). The equivalence follows directly from the bijection and
the fact that equality modulo both 2 and 3 is equality modulo 6. If either
coordinate already forces inequality, that edge needs no clause; if both
force equality, the branch is impossible.

`compile_branch` exports one-hot phase variables and the complete forbidden
assignments for (3). `audit_compiled` reconstructs them from the exact source.
`solve` returns a checked colouring, an actual nonzero-voltage closed walk,
or a complete Boolean decision tree. Negative tree leaves must be literal
unit-propagation conflicts; every internal node includes both truth values.
`audit_unsat_tree` checks that tree independently. Exhausting its explicit
node budget raises `ProofLimit`, without exporting UNSAT. Solver output alone
is not a negative certificate.

A useful residual hostile is \(K_4\) with its binary coordinate fixed
constant. Both modular systems are consistent, but (3) requires a proper
three-colouring of \(K_4\), which is impossible. The saved three-node decision
tree certifies this. The underlying graph nevertheless has a proper
six-colouring. This is the distinction between a rejected branch and a
rejected graph.

## 3. What the 375 branches do and do not cover

In inherited G4, there are six independent direction factors in an additive
Cartesian contact. A nonempty difference set with no consecutive residues
is contained in one coset modulo 2 or modulo 3. Five factors containing
consecutive residues have full sumset \(\mathbb Z/6\); hence a valid contact
has at least two non-rich factors. Selecting those two factors and one of
their five cosets gives \(\binom62\,5^2=375\) overlapping local branches.

`contact_branch` imposes the selected residue equation on **every retained
directed pair** of each selected factor. A chosen branch can then be tested
jointly with all other graph edges, source pins, and already imposed branches
by section 2. Shared vertices identify the same variable, rather than fresh
copies. This is a reusable conditional branch consumer, not a claim of only
375 global cases across many contacts.

The independent positive control for this logic is an artificial pinned
18-vertex example. Its first five direction sets equal \(\{1,2\}\), and the
sixth equals \(\{1\}\). All 375 selected pairs contain a rich factor and
produce a checked voltage contradiction. The underlying edgeless graph is
six-colourable: its contradiction applies only if the additional additive
contact premise is imposed. No unit embedding is claimed.

This does **not** eliminate any branch in the old geometric construction.
Incoming [carry descent and transported branches, section 4](carry_descent_transport_branches_20261018l.md)
at `b252ba37bc` independently exhibits full-context witnesses for all 375
branches in each of its eleven groups, 4,125 surviving branch cases. Those
are different, actual constraints. Furthermore, the geometric inputs below
have three address factors, or four if their triangle factor is unpacked;
G4 does not automatically supply a 375-case cover of their full colourings.
An additive contact structure and its source-identifying equations would
have to be supplied first.

## 4. Complete boundary support on the new literal graphs

The separately authenticated input is
[coupled fractional geometry](coupled_fractional_geometry_20261018c.md)
and its JSON. The consumer retains all physical vertices, all actual edges,
the coordinate-content hash, and original ports 0 through 7. The induced
source is the one-point union of \(C_5\) and \(C_4\), not their graph join.
Its complete restricted-growth-string bank has 570 proper patterns using
at most six colours. Every labelled source colouring is a permutation of
exactly one such pattern.

For ordinary six-colour support, a non-port vertex with at most five live
neighbours can be deleted. Every colouring of the remainder extends back
by a greedy choice, so iterated relative peeling preserves support exactly.
It does not preserve extension counts. This reduction is not used to discard
vertices subject to additional branch equations. The remaining graph uses
direct six-state one-hot variables, the exact unconstrained specialization
of (1). Every exported positive colour word is checked against all original
edges and all eight pins after reverse lifting.

| Authenticated graph | Vertices / edges | Retained vertices / edges | Mixed edges retained in core | Boundary support |
|---|---:|---:|---:|---:|
| `candidate` | 316 / 1,004 | 64 / 190 | 0 of 45 | 570 of 570 |
| `lifted_candidate` | 714 / 2,914 | 472 / 1,903 | 94 of 180 | 570 of 570 |

Both rows have **123 equal and 447 separating** patterns at terminals 2,4.
Each of the 1,140 patterns has a saved full colour word, not merely a SAT
status. An earlier independent encoding using the unreduced binary/ternary
phase clauses also found all 570 witnesses in each graph. The final packet
uses the faster direct encoding and independently checks its full witnesses.
The 714-point result therefore survives genuine mixed edges in the relative
core, but removes no source pattern.

Terminals 2,4 are not a unit-separated pair in the literal pentagon source.
Even a future equality result at that pair would require an authenticated
terminal connector before it implied a seven-chromatic unit graph. There
is no such inference here. The separate abstract \(K_7-e\) control has a
checked equal colouring and a checked separating UNSAT tree; it has no
asserted unit embedding.

## 5. All palette matchings, aliases, and the simpler triangle obstruction

Fix one source pattern \(\phi\) and a perfect matching of the six colour
names. The restricted model assigns to an address \((a,m,t)\) the pair
containing \(\phi(a)\); only its shade may vary. Every address of one
physical vertex must agree on that pair. Source shades remain pinned.

**Alias criterion (PROVED).** Form a graph on the six colour names, joining
two names whenever aliases or a literal source pin require them to belong
to the same pair. Such a matching exists iff every connected component has
size at most two. A two-vertex component forces its pair. If there are \(k\)
such components, the remaining names have
\((6-2k-1)!!\) matchings, with \((-1)!!=1\). This follows because a matching
pair cannot contain three names, while isolated names may be paired freely.
The receipt retains the physical vertex and address indices behind each
colour-name edge.

After alias consistency, retain precisely the actual graph edges whose pair
labels agree. They require opposite binary shades. An odd cycle certifies
failure; a path between pinned vertices certifies failure if its parity
disagrees with the two pinned shades. Otherwise the chart has exactly
\(2^f\) lifts, where \(f\) counts unpinned components, including isolated
vertices. These are complete fibres for this fixed chart. Different charts
overlap and their counts must not be added.

For the 714-point graph, all 570 patterns have exactly three alias-compatible
matchings. The other twelve fail on actual aliases. Thus the sweep has
6,840 alias conflicts and 1,710 odd-cycle failures, and no positive chart.

This failure is already explained by a **seed triangle**, independently of
the 94 retained mixed edges. Addresses \((0,0,0),(0,1,0),(0,2,0)\), at indices
0,6,12, give distinct physical vertices 0,1,9, with all three actual edges.
They inherit the same pair in every chart. Every surviving chart's recorded
odd cycle is \((1,0,9,1)\). Thus the sweep is an exact classification of a
repair model that was already too narrow; it is not new six-colour forcing.
The saved full witnesses necessarily permit pair-label changes across this
triangle. A useful next constrained model must allow that change, for example
through factor-dependent pair labels, and still retain all aliases and
joint edge clauses.

## 6. Reproduction and validity boundaries

The exact finite universes include all 75 simple graphs on one through four
vertices with vertex 0 pinned, 225 coupled two-equation cases on a triangle,
all 63 nonempty residue supports, all 375 artificial contact branches,
the abstract equality control, both complete 570-pattern geometric banks,
and all 8,550 physical-address palette charts on the larger graph.
Direct colour enumeration checks the small CRT counts independently.
Malformed exact types, wrong compiled sources, and incomplete negative
receipts are rejected. Checks remain enabled under Python optimization.

```text
python -B 04-computation/experiments/coupled_support_branches_20261018c.py
python -B -O 04-computation/experiments/coupled_support_branches_20261018c.py
```

The JSON stores source hashes, literal negative receipts, complete positive
colour words (one digit 0--5 per physical vertex), and every palette-chart
obstruction. These words certify the unchanged finite graphs; neither a
finite support result nor the inherited G4 cover settles the whole plane.
Normal and optimized runs both pass **12,458 checks** and reproduce the
saved JSON exactly after LF normalization. Its SHA256 is
`3fc76dde4a41e7eb099a720dc1f58da9cbcd2a1fddf473900a924c9ceaa69d11`.

**Independent proof/API and solver-free receipt audit: PASS.** The reviewer
independently rebuilt all 570 source patterns, authenticated both coordinate
and edge hashes, replayed all 1,140 full words and all 8,550 palette charts,
and checked the source-preserving peeling data (2,253,721 gates). A separate
unit-propagation tree consumer checked freshly reconstructed K4 and
\(K_7-e\) CNFs and their 3-node and 47-node negative trees; all 375 marked
voltage cycles were also authenticated (1,580 further gates). No correction
remained at freeze.
