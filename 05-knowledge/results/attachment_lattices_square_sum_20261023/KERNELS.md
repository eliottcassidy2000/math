# Larger cells: factor the invisible symmetries, retain their owners

23 October 2026. **PROVED scoped.** These are exact group and marking rules
for supplied named cells and their actual cross blocks. They do not discover
a compact cellular description of an arbitrary tournament.

## 1. Inheritance and the new operation

The previous [odd-cell completion](../odd_cell_completion_20261023/NETWORKS.md)
proved trivial one-sided kernels for regular cells of sizes3,5,7, and
exhibited a regular9-cell with a kernel of size27. The incoming
[whole-frame compiler](../session41_owned_attachments_fibonacci_20261023/residual-cycles/README.md)
already handles arbitrary finite cell groups without promising efficiency.
The new step is to factor simultaneous one-sided kernels before compilation,
and to prove exact budgets for a nested family with arbitrarily large cells.
The elementary group-extension mechanism is not claimed as a new general
group-theory theorem.

## 2. Exact kernel factorization

Let G_i be the full allowed internal automorphism group of named cell i,
including every existing owner mark. For each pair ij, let

    H_ij = {(g,h) in G_i x G_j : B_ij(g u,h v)=B_ij(u,v)}.

Uniform blocks may be omitted, with their relation understood to be the
whole product. Put P_i equal to the intersection of all projections of
incident H_ij onto G_i. Put

    K_i = intersection_j {g : (g,1) belongs to H_ij}.

For an isolated cell, P_i=K_i=G_i. K_i consists of precisely the internal
symmetries that can move independently while fixing every external owner.
It is not the set of arbitrary permutations of equal binary rows.

For a subgroup H of A x B, its one-sided kernel is normal in its projection:
conjugate (k,1) by any (a,b) in H. Applying this to each edge shows
K_i is normal in P_i, even when it is not normal in all G_i.

Let G be the full cell-preserving automorphism group of the supplied system.
Then K=product_i K_i is a normal subgroup of G. Moreover G/K is exactly the
compatible relation system on the quotient domains P_i/K_i: on each edge
retain the image of H_ij intersect (P_i x P_j). Every representative of an
accepted quotient tuple is a valid full tuple, because multiplying its
coordinates independently by K_i leaves every incident block unchanged.
Thus each quotient tuple has exactly product_i |K_i| lifts.

The statement extends to comparing two systems. Their nonempty local
relations are left cosets R_ij=r_ij H_ij of the source stabilizers. Intersect
their incident projections at cell i. Reject if this domain is empty;
otherwise choose an actual representative a_i, so the intersection is
a_i P_i. Translate by a_i to obtain domains P_i, then quotient by K_i as
above. Local relations may still be empty after restriction and global
compatibility must still be checked. Independent lifts preserve all the
constraints, so the same constant multiplicity holds. No guessed compatible
global frame is required to form this reduction.

### Marking and observation are additional obligations

Write b(K_i) for the minimum number of actual vertices whose pointwise
stabilizer in K_i is trivial. Any base for G needs at least

    sum_i b(K_i)

owners: a nonidentity element in an unbroken K_i extends by identity on
every other cell. Choose local bases B_i and let B be their union. Then
G_B intersects K trivially, so G_B injects into G/K. Further owners only
need to kill this residual group. In particular an upper bound is

    sum_i b(K_i) + floor(log_3 |G/K|).

For the last claim, tournament automorphism groups have odd order: an
involution would exchange some pair and reverse its arrow, and Cauchy's
theorem excludes even order. Every nontrivial residual group therefore
has a vertex orbit of size at least3. Marking a vertex in that orbit reduces
the group order by at least3. The actual residual after B can be smaller
than G/K; use it for the sharper bound.

The quotient alone preserves feasible comparison frames modulo K and their
counts. It preserves an observation P only if P(f k)=P(f) for every k in K.
Literal path owners usually fail this condition. Either retain the actual
lifts, fix the necessary owners before forming G_i, or retain the exact
fibre enumerator sum_(k in K) z^(P(f k)). Replacing that enumerator by |K|
is legitimate for unweighted counting, not for arbitrary port queries.

## 3. A nested tournament family with an exact all-size budget

Let R_d be the lexicographic power of the directed triangle, on words
x in {0,1,2}^d. At the first position where x,y differ, orient x->y exactly
when y's digit is x's digit plus1 modulo3. Thus R_1=C3 and R_2 is the prior
regular9-cell.

The prefix hierarchy is intrinsic. If two vertices have common-prefix
length r, their number of common out-neighbours is

    (3^d - 2*3^(d-r-1) - 1)/2.

To count, common outgoing branches above their split contribute
(3^d-3^(d-r))/2. If x->y at the split, their common outgoing vertices there
lie inside the child containing y, contributing
(3^(d-r-1)-1)/2. These values are distinct for r=0,...,d-1, so the hierarchy
is recovered from the tournament itself. Each internal tree node permits
exactly a cyclic permutation of its three children. Consequently the full
automorphism group consists of one independent C3 rotation at every
internal node, with the usual composition of these address-dependent maps.
Its order is 3^((3^d-1)/2); it is generally nonabelian.

Now take m named copies of R_d. An edge has a depth label ell in {1,...,d}
and cross block

    B_ell(x,y)=1 iff the first ell digits of x and y agree.

Give the support edges fixed orientations. Other blocks are fixed uniform
ones. All blocks are equitable, with row and column degree3^(d-ell).
For k=1,...,d, let H_k be the graph consisting of edges with ell>=k and
let c_k count all its connected components, including isolated cells.

**Theorem.** The full cell-preserving automorphism group has order

    3^(c_1 + 3*c_2 + ... + 3^(d-1)*c_d).

The minimum number of actual marked vertices needed to kill it is exactly

    3^(d-1)*c_d.

That many marks also suffice for ordinary directed equitable refinement
to become discrete, retaining the named original cells.

**Proof.** Equality of a B_ell block says exactly that the two cell actions
induce the same permutation of length-ell prefixes. In the source-prefix
rotation normal form, the rotation at a node of depth r must therefore
agree across every edge of H_(r+1). For each of the3^r source prefixes there
is one independent ternary choice per component of H_(r+1). Counting these
normal-form choices proves the group order. There is no assertion that
the group law on these coordinates is coordinatewise addition.

At depth d-1, every prefix and every component of H_d admits an independent
rotation of its bottom triple, simultaneously in precisely the cells of
that component. It fixes all other owners. Each such triple union needs
a marked vertex, proving the lower bound.

For the upper bound, in one representative cell of each H_d component mark
one actual vertex in every bottom triple. These vertices form an internal
base. They also distinguish all remaining vertices by adjacency profiles:
vertices in one bottom triple are distinguished by its mark; for vertices
in different bottom triples, a mark in the third child at their lowest
common ancestor distinguishes them. H_d edges are identity-matching blocks,
so a discrete representative cell propagates discrete refinement throughout
its component. Every cell belongs to such a component. QED.

At d=2 this gives the particularly useful rule **3c_2 marks**. A connected
coarse-only network needs3m marks, but fine matching edges synchronize whole
cells and reduce the cost to three times the number of fine-edge components.
For a three-cell path with depths(1,2), the exact cost is6, not9. This is a
native larger-cell improvement from interaction, not merely a scalar analogy.

The c_k form a nested connectivity profile. They give different information:
all levels determine the number of frames, while the deepest level controls
the actual-owner base. Distinct groups can therefore have the same marking
cost. Retaining only that cost loses the intermediate symmetries.

## 4. Boundaries and next target

These blocks compare prefixes in one retained coordinate system. Twisted
prefix selectors require the actual transports and their compatibility;
the c_k counts alone need not classify them. General larger cells need not
have this prefix hierarchy. The two-congruence [lattice compiler](lattice/README.md)
applies to cyclic chart subgroups, whereas the full groups above require
their nonabelian frame relations or the kernel-factorized finite compiler.

This is the same exact-elimination discipline as the Fibonacci zero slot:
an eliminated coordinate contributes a relation and, when needed, a
multiplicity or observation polynomial. Its disappearance from geometric
length is not permission to discard its shared-address information.
