# A graph of path and domination fibres isolates the missing exchange lemma

2026-11-11. **PROVED:** an exact defect decomposition, a constructive
matching criterion, and the transitive cut-leaf lemma. **FINITE-EXACT:**
the candidate component bound on all old tournaments through order five.
**OPEN:** that bound in general, and the resulting quantitative domination
inequality. No finite survival is promoted to a theorem.

## The inherited target and the new carrier

Use Session260's
[positive cut law](../../session260_quasiperfect_chanpak_20261111/bridges/DEFECT_QUANTITATIVE.md)
and [domination injection](../../session260_quasiperfect_chanpak_20261111/bridges/BIPARTITE_DEFECT.md).
For U and a newborn z with z->A and B->z, let G be the bipartite graph
whose cross edges are the old arcs a->b, a in A and b in B. The remaining
quantitative question is

    H(U+z) >= H(U)+d(G)-1.                              (1)

The earlier injection chooses a newborn path for each dominating cut,
but does not reserve one distinct insertion for every old Hamilton path.
The faithful new object keeps both sorts of fibres simultaneously.

Make a bipartite multigraph with internal vertices of two types:

* every Hamilton path P of U;
* every dominating cut L of U, where the corresponding dominating set is
  `(A\L) union (B intersection L)`.

Each Hamilton path Q of U+z has a unique cut L (vertices before z).
If deleting z leaves an old Hamilton path P and L is dominating, draw
an edge P--L. If exactly one of these predicates holds, draw a stub at
that one internal vertex. If neither holds, retain Q as an unused object.
Distinct Q remain distinct edges or stubs. Both properties can be tested
from the literal path; no scalar count is used as a substitute for it.

Every old path has at least one insertion slot for z, so its vertex has
positive degree counting stubs. Every dominating cut has a positive
endpoint-path fibre by the inherited theorem. Thus no internal vertex
is silently isolated. A component is called closed if it has no stub.
Let tau be the number of closed components that are trees, and let u be
the number of unused newborn paths.

## Exact defect and a sufficient injection theorem

For each component C let n_C be its internal vertex count, m_C its edge
count, and s_C its stub count. Connectedness gives m_C>=n_C-1.
If C is open, s_C>=1; if C is closed and not a tree, m_C>=n_C.
Consequently the quantity

    E = sum_(C not a closed tree) (m_C+s_C-n_C)

is a nonnegative integer, and exact counting gives

    H(U+z)-H(U)-d(G)+1 = u+1-tau+E.                    (2)

This is an all-order identity. It localizes the possible negative terms
to closed tree components. In particular **tau<=1 implies (1)**.

There is also a literal injection under that condition. In an open
component choose one stub and a spanning tree. Assign the stub to its
incident root and each other vertex the tree edge toward that root. In a
closed component containing a cycle, assign cycle edges cyclically to
their vertices and tree edges toward the cycle to the remaining vertices.
Any extra edges are unused. In the possible single closed tree omit one
cut vertex and assign every other vertex the edge toward it. A closed
tree has a cut vertex because it contains no stubs and every internal
vertex has positive degree. These assignments use distinct newborn paths.
They inject all old paths and all but at most one dominating cut into
newborn paths, exactly the desired simultaneous reservation.

The all-order statement tau<=1 is still **OPEN**. Equation (2) also
shows it is sufficient, not logically necessary for (1): unused paths
and surplus in other components could compensate for multiple trees.

## A rigid leaf lemma

In a closed component, an old path P has degree exactly `1+#AB(P)`.
This is the inherited binary insertion identity. A path-side leaf
therefore has no AB transition: its colour word is all B followed by all
A, and its unique neighbouring cut is L=B.

A cut-side leaf has an endpoint fibre of size exactly one. Its two
induced subtournaments U[L] and U[V\L] must both be transitive.

To prove the assertion for L, use its domination condition: every A
vertex in L has an outgoing neighbour in B intersection L. The last
strong component of U[L] meets B. If that component is nontrivial and
contains just one B vertex b, all its A vertices must point to b: they
cannot point to any earlier strong component and have no other B choice
inside their component. This contradicts strong connectivity. It therefore
contains at least two B vertices, giving at least two B-ending paths by
the Hamilton-cycle consequence used in the inherited injection proof.
Any earlier nontrivial component likewise contributes at least three
path choices. Endpoint count one thus forces every component to be a
singleton. The argument for the A-starting right side is dual. Empty
sides cause no exception.

For a cut leaf adjacent to P, this says all backward arcs relative to
the order P must cross that cut. This gives a concrete place to seek an
exchange: a second legal insertion slot not crossed by such an arc must
produce additional fibre structure. Turning this observation into a
cycle-or-stub lemma remains open.

## Hostile probes and recovery from older work

The candidate that every closed tree contains L=B is false: take old
U=(a->b), A={a}, B={b}. One closed tree is
`empty cut -- (a,b) -- full cut`; the cut B has a deletion-invalid stub.
Even the weaker candidate restricting this phenomenon to transitive U
is false. With rows `[6,12,8,1]`, A mask10 and B mask5, H(U)=5; a closed
tree consists of path `(0,1,2,3)` and cuts1 and7, excluding cut B.
Both are local failures of proposed matching arguments, not refutations
of (1) or the tau<=1 conjecture.

The route was suggested by two older side branches: the July cut/cycle
reflection
`07-reflections/the-half-tiling-is-an-abelian-square-complex-the-metagraph-is-an-anti-LTC-and-cut-cycle-is-alexander-duality-kps.md`
and the connected-cluster/excess organization in **THM-4757**,
`01-canon/theorems/THM-4757-linear-forest-loop-expansion-of-hamiltonian-paths.md`.
Only the method of retaining incidence and counting component excess is
transferred. No historical topological or spectral assertion is needed:
(2) is proved directly for the present multigraph.

Source->target: literal newborn paths -> their two retained fibres.
Preserved: every path count and the exact defect. Lost by keeping only
degrees: connectedness, closed tree components, and matching feasibility.
The necessary extra data are the actual incidences and stubs. The two
hostiles above are the cheap tests against a premature canonical-root
argument.

## Reproduction

    python3 -B 05-knowledge/results/multiplicand_saturation_chanpak_20261111/exchange/probe.py --write

The receipt exhausts 33,866 pairs (old labelled tournament, newborn cut),
with old orders1 through5. It finds at most one closed tree component in
that finite universe. It records the nontransitive hostile separately.
The program enumerates literal paths, endpoint fibres, domination
conditions and connected components; it does not infer an all-order law.
