# Marked sources, two endpoint marks, and a cyclic-three-cell decoder

2026-10-23. **PROVED:** the single-mark DRT refinement obstruction and the
scoped three-cell attachment decoder. **FINITE-EXACT:** all1,860 ordered
two-mark tests on the two incoming order31 examples, and2,626 cyclic-core
decoder checks through order5. No unrestricted small-signature theorem.

The incoming [minimal-signature packet](../metagraph_minimal_signatures_20261023/README.md)
already proves the unrooted DRT obstruction and separates its order31 pair
by a transitive-five motif. We retain that stronger motif result. The new
question here is different: what does retaining the ACTUAL port source,
or both endpoint identities, buy a neighbour-refinement decoder?

## 1. One distinguished source still cannot resolve any DRT pair

Let T be doubly regular of order4t+3. Start with the uniform vertex seed
retained by the incoming score/triangle/strong-four observer and mark a
single vertex r. Put U=N+(r), W=N-(r), both of size2t+1.

After one refinement the cells are {r},U,W. Their outgoing counts into
these three cells are

| Source cell | into {r} | into U | into W |
|---|---:|---:|---:|
| {r} | 0 | 2t+1 | 0 |
| U | 0 | t | t+1 |
| W | 1 | t | t |

For u in U, common outneighbors of u and r number t, giving the U entry;
the total outdegree gives the W entry. For w in W the same common-neighbor
law gives its U entry, and its arc to r gives the last entry. Every count
depends only on t, not on the vertex or the tournament.

Thus this partition is stable. Every further directed neighbour-refinement
round, including both in- and out-counts, agrees for every singly marked
DRT of the same order. This is stronger than just observing a remaining
automorphism: even distinct rooted DRTs have the same whole sequence of
these records. The preserved source is necessary information, but this
particular observer still discards its higher neighborhood structure.

This theorem explains the exact order31 result[1,15,15] for all31 roots
of both P31 and D(D(P7)). It prevents promoting the successful order13
source-mark repair from README.md into a universal one-mark rule.

## 2. Two actual marks break the two examples differently

Use two distinct, ordered vertex marks, with every other vertex initially
uncoloured. Because the inherited DRT tags are uniform, appending those
tags produces the same partition refinement. Exhaust all31*30=930 ordered
pairs in each incoming graph, constructed independently from Paley7 and
the inherited skew-Hadamard doubling formula.

| Graph | Stable cell-size multiset | Ordered mark pairs |
|---|---|---:|
| P31 | 31 singletons | 930 |
| D(D(P7)) | 31 singletons | 504 |
| D(D(P7)) | 7 singletons and8 triples | 294 |
| D(D(P7)) | 3 singletons and4 seven-vertex cells | 132 |

For P31, marks(0,1) give

    [1,1,29]
    [1,1,7,7,7,8]
    [1,1,2,2,2,2,2,2,3,3,3,4,4]
    [1,...,1].

For D(D(P7)), the same numerical marks stabilize at[1,1,1,7,7,7,7].
Every stable nontrivial cell is odd, as the odd-cell lemma predicts.
The tables refer to ALL ordered vertex pairs. They do not assert every
pair is the endpoint pair of a Hamiltonian path. For actual ports, retain
the literal source/path and apply the test only to its actual endpoints.
A discrete partition plus the directed outgoing-colour record recovers
every labelled arc. A nondiscrete partition is not itself proof of a
different target class: further constraints or the existing five-motif
observer can still decide the class.

## 3. A local decoder for a surviving cyclic triangle

Suppose the outside vertices have retained identities, a specified three
vertex cell D is known to induce a directed3-cycle, and all outside arcs
and the per-outside-vertex margins

    m_j = #{u in D : u->v_j}

are retained. For m_j=0 or3 the attachment is fixed. For m_j=1 or2 define
a normalized unit column z_j on D:

    z_j(u)=A[u,v_j]        if m_j=1,
    z_j(u)=1-A[u,v_j]      if m_j=2.

Each mixed column now selects one actual core vertex. Retain the pair
overlaps

    C_ij=sum_{u in D} z_i(u)z_j(u).

They are0/1 and say whether the SAME core vertex owns both selected
incidences. Their equality classes recover every attachment up to a
permutation of the three core vertices. If at least two selected core
vertices occur, retain one directed arc between representatives of the
first two different equality classes. This is one chirality bit. It fixes
the directed3-cycle, and hence every arc, up to an isomorphism fixing all
outside identities. If zero or one core vertex is selected, no chirality
bit is needed: exchanging two unused core vertices identifies the two
cycle orientations while fixing every attachment.

This is a complete constructive proof: name the at most three equality
classes, fill the unit/co-unit columns from m_j, and choose the cyclic
orientation from the bit. Each missing internal arc is then forced. For
untrusted data, reject unless the alleged overlaps define an equivalence
relation with at most three classes and all retained score/profile guards
hold. The checker tests the reconstruction up to all six core permutations
while fixing every outside vertex.

There is a real change from the two-vertex decoder. Equality tests on a
spanning tree no longer suffice: labels(0,1,0) and(0,1,2) have the same
two adjacent inequalities, but different endpoint equality. A forest
alone cannot select between them. Full equality data suffices; adaptively,
testing each new column against at most two representative classes also
suffices, since a third class is the last available one.

This is a LOCAL attachment theorem. If several outside cells are themselves
unresolved, their vertex identities cannot silently be chosen separately.
Their alignment constraints must be retained and solved jointly. That is
the next phase-synchronization target suggested by the eight triple cells
in the doubled order31 example.

## 4. Relation to the older triangle grammars

The [closed deficit-six classification](../linked_tail_and_deficit_six_closure_20261022.md)
already proves R6(z)=z^5+z^6+z^7. Its older quadratic-cutoff and incomplete
center notes are superseded as frontiers; this packet does not reopen
their settled coefficients. Those results count actual path-orbit deficits
using synchronized cyclic shores. Here the three-cell object is a boundary
attachment representation with named outside vertices, not a new center
count or proof that a local triangle move preserves all Hamiltonian ports.

Source -> target: a specified cyclic core with its literal boundary ->
unit/co-unit selector classes plus one orientation bit. Preserved: all
arcs and every outside port identity, up to internal core isomorphism.
Lost without the bit: cyclic orientation relative to two independently
selected core vertices. Lost by independent decoding of multiple cells:
their relative alignments. Cheapest hostile: the ternary tree example;
cheapest positive control: compare decoded arrays under all six core
permutations while fixing the boundary.

Run `python3 -B 05-knowledge/results/incidence_wiring_20261023/marked_cells.py`
and repeat with `-O`. Both compare marked_cells.json without rewriting it.
