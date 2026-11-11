# Closed-tree exchange components: branching restrictions and a unique-star reduction

Date: 2026-11-11. **PROVED:** the local branching restrictions and uniqueness
of closed star components below. **OPEN:** every closed tree is a star.
Consequently this note does not prove the full quantitative inequality.

Use the incidence graph defined by the sibling `exchange/probe.py`.
Its left vertices are Hamilton paths P of U, its right vertices are
dominating cuts L for the apex z, and a Hamilton path of U+z is an edge
when deletion gives P and its z-cut is L. A path with only one such
incidence is a stub. Write A=N+(z), B=N-(z). A component is closed when
it has no stubs. All edges mentioned below are literal newborn paths.

The endpoint existence and SCC facts used below are inherited from
[BIPARTITE_DEFECT.md](../../session260_quasiperfect_chanpak_20261111/bridges/BIPARTITE_DEFECT.md).
The target injection and cut-fibre notation are in
[DEFECT_QUANTITATIVE.md](../../session260_quasiperfect_chanpak_20261111/bridges/DEFECT_QUANTITATIVE.md).

## 1. Exact degrees and the first hostile example

For an old path P, its full incidence degree including stubs is

    1 + number of A-to-B adjacencies in its colour word.

In a closed component all these insertions are edges. Thus a path-side
leaf has colour word B* A*, and its unique cut is L=B. Every closed
star centered at a cut therefore has center B, so there is at most one
such star.

A closed tree need not contain B, and its old tournament need not be
transitive. An exact order-four witness is

    rows=[6,12,8,1],  A=10, B=5,  H(U)=5.

Here a bit j in rows[i] means i->j. The component containing the old
path `(0,1,2,3)` has just the two cut leaves `L=1` and `L=7`, two
edges, and no stubs. It does not contain cut B=5. This refutes the
proposed strengthening that every closed tree in a nontransitive old
tournament contains B. The even smaller transitive example U=(a->b),
A={a}, B={b} has the tree `empty -- (a,b) -- V`.

## 2. A cut leaf has two transitive sides

At a dominating cut L, R=V\L, the full fibre size is
`e_B(L)*s_A(R)`. If it is1, then U[L] and U[R] are both transitive.

For the left side, domination says every A vertex in L has an
out-neighbour in B∩L. Its terminal SCC therefore meets B. If that SCC
were nontrivial and met B in only one vertex b, every other vertex in
the SCC would be in A and would have to point to b: it cannot point
to a B outside the terminal SCC. Then b has no outgoing arc inside
the SCC, a contradiction. A nontrivial terminal SCC therefore has at
least two B vertices, either of which can end a Hamilton path, giving
e_B>=2. Any nontrivial earlier SCC similarly multiplies the endpoint
count by at least3. Hence e_B=1 forces all SCCs singleton. The right
side is dual. Empty sides are transitive with endpoint count1.

## 3. The four-cycle restriction on every branching cut

**Lemma.** In a closed tree, at any cut L with complement R,

    e_B(L)>=2  implies R⊆A,
    s_A(R)>=2  implies L⊆B.                              (1)

For the first implication, suppose R also contains a B. Fix one
A-starting Hamilton path q of R and two different B-ending Hamilton
paths p1,p2 of L. Closedness ensures both concatenations pi q are old
Hamilton paths; otherwise their newborn paths would be stubs at L.
Since q starts in A and contains B, its colour word has a legal z-slot
strictly after the original cut: after a BA transition, or at the end
if it ends in B. That slot has the same cut L' in p1q and p2q because
their left vertex set L and their right path q agree. Closedness makes
L' dominating; otherwise the corresponding newborn paths would be
stubs on the old-path side. The four edges between p1q,p2q and L,L'
form a cycle, contrary to the tree hypothesis. The other implication
is obtained by reversing all paths and exchanging A,B.

Consequently every cut with both `L∩A` and `R∩B` nonempty is a leaf.
Every branching cut lies either below B (`L⊆B`) or above B (`L⊇B`).
If both endpoint factors exceed1, the cut must equal B.

This uses the full Cartesian product of endpoint paths and the fact
that closedness retains every incident fibre. Choosing one canonical
path per cut would lose the four-cycle argument.

## 4. At most one closed star, and why a star theorem would suffice

A path-centered closed star with at least two cut leaves has a central
old path P with at least one AB adjacency. Split its order as

    P = X Y Z,

where X is its initial all-B run, Z its terminal all-A run, and Y
starts in A and ends in B. The earliest legal z-slot is after X and
the latest is before Z. Their cut fibres have size1, so Section2 gives

    U[X union Y] and U[Y union Z] are transitive in P-order. (2)

Thus both A and B are transitive, their orders in P are fixed by U,
and every backward arc relative to P runs from Z⊆A to X⊆B. In
particular no B-to-A arc can run backward in P.

**Uniqueness.** Suppose P,Q are two different path-centered closed stars
with at least two leaves. Their orders on each shore agree, so some
A vertex a and B vertex b occur as a-before-b in P and b-before-a in
Q, after interchanging the names P,Q if necessary. The arc cannot be
b->a, since that would be a backward B-to-A arc in P. Thus a->b is
backward in Q, forcing b into X_Q and a into Z_Q. Let a0 be the first
A vertex of Y_Q. Then b->a0 by (2), and a0 precedes a in the fixed
A-order. In P we have `a0 before a before b`, making b->a0 a backward
B-to-A arc, a contradiction. Therefore such a star is unique.

It cannot coexist with a closed star centered at cut B. Indeed the
first A vertex a0 and last B vertex b0 both lie in Y_P, where a0->b0.
The shores are transitive, so their unique B-then-A concatenation
would need the opposite seam b0->a0. The newborn path formed from
those two shore paths is therefore a deletion-invalid stub at B.

A one-edge component has both centers and is already covered by the
cut-B case. We conclude:

    If every closed tree is a star, then there is at most one
    closed tree in the entire exchange graph.                 (3)

The finite evidence so far supports the antecedent, but (1) does not
prove it: a hypothetical nonstar tree could alternate branching cuts
below and above B. Excluding that configuration, or producing a
counterexample, is the exact remaining structural task.

## Finite scout, not a theorem

An instrumented read of the sibling probe, retaining each component's
path and cut vertices, tested500 pseudorandom cases at each old order3
through7 using Python `random.Random(125)`. In all2,500 cases, every
closed tree was a star. The tested shapes included both a single old
path with up to3 cut leaves and cut-B stars with many old-path leaves.
An earlier independent160-case scout (120 order6,40 order7; seed20261111)
also found at most one closed tree. These are finite exploratory scopes,
not exhaustive universes or evidence of an all-order proof by themselves.

The saved `exchange_shape.py` reproduces the2,500-case shape scout and
the named order-four hostile witness without modifying the root probe.
