# Six colours through signed boundary parity

2026-10-11. **PROVED:** the elementary reformulations, fixed-label interface,
gluing theorem, and witness lifts below. **FINITE-EXACT:** the stated complete
small-graph comparisons. **OPEN:** a universal six-colouring of the Euclidean
unit-distance graph, or a finite seven-chromatic unit-distance witness.
This packet supplies neither. Its signed-graph principles are classical;
the concrete contribution here is a checked boundary compiler with retained
path witnesses.

## Inheritance and connection contract

The closest repository mechanism is
[boundary palette compression](boundary_palette_compression_20261011.md):
a fragment exports its complete relation on named boundary ports, with
colouring lifts. The hostile is replacing that relation by a colour count.
The present specialization retains binary parity after fixing a ternary
label. The least-used sidecar is the actual path that proves each boundary
equation. The live board is **ternary assignment / signed forest / named
ports / odd-cycle witness / geometric realization / complete quantifier**.

The source is a finite simple graph with a fixed ternary assignment. The
target is a signed relation on its boundary. The map removes interior
vertices while retaining boundary connectivity and relative binary parity.
It preserves exactly the extendible boundary bit assignments. It forgets
interior choices and coordinates; a spanning forest restores colouring
witnesses and lifts seam failures to actual cycles. Geometric coordinates
must be supplied separately. The cheapest decisive test is an even path
glued to an odd path with the same endpoints.

Harary balance and its all-negative/bipartite specialization are classical;
see Zaslavsky's [Reference Notes for Lectures on Signed Graphs and Geometry,
Theorem 2.3 and Corollary 2.4](https://people.math.binghamton.edu/zaslav/Tpapers/pala-notes.pdf).
The proofs needed here are given directly. No signed-colouring convention
or geometric representation theorem is imported.

## 1. The exact six-colour reformulation

Write a colour as
\[
 c(v)=2q(v)+b(v),\qquad q(v)\in\{0,1,2\},\quad b(v)\in\mathbf F_2.
\]
For a fixed \(q\), let \(G_q\) retain precisely those edges \(uv\) for which
\(q(u)=q(v)\), keeping all vertices. Then
\[
 c\text{ proper}\quad\Longleftrightarrow\quad
 b(u)+b(v)=1\quad\text{for every }uv\in E(G_q).                 \tag{1}
\]
Edges between different ternary labels impose no binary constraint.
Consequently the following statements are equivalent:

1. \(G\) has a proper six-colouring.
2. Some ternary assignment \(q\) makes \(G_q\) bipartite.
3. The vertices partition into three, possibly empty, classes each inducing
   a bipartite graph.
4. The hypergraph whose hyperedges are vertex sets of odd cycles of \(G\)
   has a three-colouring with no monochromatic hyperedge.

Indeed, a graph is bipartite exactly when it has no odd cycle. One may
restrict (4) to chordless odd cycles: a chord splits a monochromatic odd
cycle into a shorter odd cycle and an even cycle, and repeated shortening
ends at a chordless odd cycle. This is an exact reformulation, not a bound
on the number or length of cycles that must be considered.

Giving every edge of \(G_q\) sign \(-1\), balance means that every cycle has
positive sign. Its cycle sign is \((-1)^{\text{length}}\), so balance is
precisely bipartiteness. A spanning-forest parity potential proves balance,
or a violated edge together with its forest path proves an odd cycle.

The graph-product distinction matters:
\[
 K_3[K_2]\cong K_6,\qquad K_3\times K_2\cong C_6.              \tag{2}
\]
The first is the lexicographic product: different ternary coordinates are
already adjacent. The second is the categorical product, where both
coordinates must differ. A triangle maps to \(K_6\) and cannot map to the
bipartite graph \(C_6\). Requiring both coordinates to change is therefore
a strictly stronger, incorrect replacement for (1).

## 2. A complete boundary signature for a fixed label

Let \(F\) be a finite simple fragment and \(B=(b_1,\ldots,b_s)\) its ordered
list of distinct named boundary vertices. Fix \(q\) on every vertex.
If \(F_q\) is nonbipartite, export a monochromatic simple odd cycle and the
empty extension relation. This obstruction remains relevant even if
\(B\) is empty.

Otherwise choose a spanning forest of \(F_q\) and parity potential
\(p(v)\), zero at each forest root and changing by one on every forest
edge. For every residual component meeting \(B\), choose its first port
\(r\) in the stated boundary order. Export, for every port \(b_i\),
\[
 \bigl(r_i,\epsilon_i\bigr),\qquad
 \epsilon_i=p(b_i)+p(r_i)\pmod2.                             \tag{3}
\]
Retain \(q|_B\), the boundary names and order, and the forest paths from
\(r_i\) to \(b_i\). The complete extension predicate is
\[
 b(b_i)+b(r_i)=\epsilon_i\quad\text{for every }i.             \tag{4}
\]

Necessity follows by summing (1) along a path. For sufficiency, in a
component meeting \(B\), flip the forest potential by the unique constant
that gives the assigned bit at its first boundary vertex. Equation (4)
then matches every other port. In a component avoiding \(B\), choose either
flip. All residual edges are satisfied; combine with \(q\) using (1).
Thus the interface is exact and the lift is constructive.

Changing a component's global binary gauge does not change (3).
Connectivity without relative parity is insufficient, as the next section
shows. A component count loses still more information: it omits which
named ports are related. For \(s\) boundary ports there are at most
\(\sum_{j=0}^s S(s,j)2^{s-j}\) consistent signatures before taking ternary
labels into account, plus an inconsistent state. This is an interface
state bound, not an efficient algorithm for finding the correct ternary
assignment.

When the interior ternary assignment may vary, the correct interface is
the **union** of the corresponding relations, indexed by compatible
boundary ternary labels. It need not be a single system of binary
equations. A table of valid positive lifts is not automatically complete.

## 3. Exact gluing and a lifted obstruction

Suppose \(F\) and \(H\) intersect exactly in their same ordered boundary
\(B\), their ternary labels agree there, and their union has exactly the
declared edges. Then their fixed-label colourings glue if and only if
neither fragment has an internal odd-cycle obstruction and the union of
their boundary equations (4) is consistent. Equivalently, intersect their
complete extension relations, with an inconsistent packet representing
the empty relation even when its boundary is empty.

Necessity is restriction. For sufficiency, solve the binary equations,
lift the same seam bits separately in \(F\) and \(H\), and merge the two
colourings. Disjoint interiors and equal boundary values make this a
well-defined colouring of the whole union.

The implementation solves these boundary equations as a signed graph.
If it finds a contradictory seam cycle, replace each seam edge by its
retained path in its original fragment. The resulting closed walk has
odd total length: its constraint parities sum to one. Every edge is an
actual edge of the union, and every vertex on the walk has the same
ternary label. An odd closed walk contains a simple odd cycle. Construct
one by retaining any repeated-vertex subloop of odd length, or deleting
that subloop when its length is even. Each operation strictly shortens
the walk and preserves oddness until only the starting vertex repeats.

The returned cycle is independently checked against the original graph.
This proves failure of the **supplied ternary assignment**. It proves
six-uncolourability only if every ternary assignment is eliminated by a
complete argument. Reassigning its ternary labels can remove this local
obstruction.

The minimal parity hostile has boundary \((0,1)\). The path \(0,2,1\)
forces equal bits; the edge \(0,1\) forces unequal bits. Each fragment is
bipartite and has one connected boundary component, but their union is a
triangle. Both path witnesses and the assembled odd cycle survive in the
receipt. Exact coordinates
\[
 (0,0),\quad (1,0),\quad (1/2,\sqrt3/2)
\]
make every triangle edge unit length. These are subgraph fragments; they
are not claimed to be induced graphs of their respective coordinate sets.

For geometric assembly, shared vertices must be identified at the same
point, distinct vertices must remain distinct, and any desired additional
cross-unit edges must be included. Pure boundary parity does not certify
those conditions.

## 4. What a bounded obstruction bank loses

For every fixed length bound, an odd cycle longer than that bound has no
shorter odd cycle. Giving its vertices one ternary label therefore passes
all tests that look only for bounded-length monochromatic odd cycles, but
fails (1). This hostile is geometrically available: a regular odd polygon
scaled to side length one realizes the cycle as a unit-distance subgraph.
For odd order at least five, its nonconsecutive chords are longer than
one, so it is also the induced unit-distance graph on those points.

The compiler has no fixed cycle-length cutoff. Its forest paths can
return the entire long cycle. Conversely, no bounded boundary width,
finite set of fragments, or universal ability to find a good \(q\) is
proved here. Long odd cycles themselves have chromatic number three;
their role is to refute a bounded test for a fixed \(q\), not to establish
any plane lower bound beyond what is already known.

A useful next instrument is an adaptive ternary search. Whenever a
candidate \(q\) fails, the returned odd cycle supplies a valid constraint:
its vertices may not all receive the same ternary label. Retain that
cycle with its actual edges and continue. Every such learned constraint
preserves every proper six-colouring. A finite graph gives a finite
complete search if all remaining assignments are eventually tested; no
small bound or plane-wide termination follows. This is a precise way to
combine the parent counting formula with this constructive rejection
sidecar.

## 5. The number 1099 and the finite controls

The structural identity is
\[
 \sum_{n=1}^{5}2^{\binom n2}=1+2+8+64+1024=1099.             \tag{5}
\]
Each unordered vertex pair independently chooses whether to be an edge,
so these are all labelled simple graphs of orders one through five.
The same bit count encodes labelled tournaments after choosing a vertex
order and using each bit to orient its pair. That encoding preserves a
bitstring, not graph colouring or a tournament certificate. The empty
graph would add one to (5).

The literal decimal expression supplied in discussion satisfies
\[
 10^3+10^2-10=1090,\qquad 10^3+10^2-1=1099.
\]
The latter numerical equality supplies no general graph-count formula.
Both equalities, including the failed 1099 assertion with subtraction
by ten, are explicit arithmetic controls.

The checker uses all 1099 graphs, and all ternary assignments modulo
global permutation of their three names, encoded by restricted-growth
strings. This is a complete symmetry reduction because (1) depends only
on equality of ternary labels. Independently enumerating all binary
assignments supplies each reference extension relation.

- 42,925 graph/ternary cases: 37,821 extendible and 5,104 obstructed.
- 128,730 boundary relations, using empty, first-three, and full
  boundaries with duplicates removed; 814,385 checked colouring lifts.
- All 5,184 two-port gluings of fragments on \((0,1,2)\) and
  \((0,1,3)\), with every ternary label choice: 4,977 compatible and
  207 obstructed.
- Exact unit-triangle, graph-product, long-cycle, empty-boundary,
  path-simplification, and 18 malformed-type/provenance controls.

The graphs through order five are all six-colourable; their obstructed
cases refer only to their specified ternary assignment. These controls
validate the compiler and gluing theorem, not a difficult six-colour
instance. Exported packets are reauthenticated, including exact integer
types, before equality comparisons; bool/float aliases cannot impersonate
vertices, parities, or labels. The implementation rejects hidden shared
interiors, reordered ports, incompatible boundary ternary labels, and
forged forests or signatures.

Reproduction:

```text
python -B -X utf8 04-computation/experiments/plane_six_reframes_parity_20261011b.py
python -B -O -X utf8 04-computation/experiments/plane_six_reframes_parity_20261011b.py
```

The companion [script](../../04-computation/experiments/plane_six_reframes_parity_20261011b.py)
exports `compile_boundary`, `audit_packet`, `accepts`, `lift`, and `glue`.
The [saved output](plane_six_reframes_parity_20261011b.out) records
2,040,717 explicit checks. Normal and optimized output match the saved
output after LF normalization, with SHA256
`ba29bc5a4a1ddcfd555667404fc0c4b0c2d37cc169d6de0e205ea06997798c2c`.
