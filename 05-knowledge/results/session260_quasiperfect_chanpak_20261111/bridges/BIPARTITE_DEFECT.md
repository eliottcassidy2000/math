# Domination injects into Hamilton paths: an exact nonnegative defect

**Status: PROVED for every finite tournament and every marked vertex.** This closes the inequality and equality question left OPEN in DOMINATION_TRANSFER.md. That accepted earlier note is preserved unchanged as the discovery record. The proof does not rely on its finite scout.

## 1. Main theorem and graph interface

Let T be any tournament with marked vertex z. Put A=N_T^+(z), B=N_T^-(z), and V=A union B. Define a simple undirected bipartite graph G_z on V by

    {a,b} in E(G_z) iff a in A, b in B, and a->b in T.

There is no restriction on the tournaments induced by A or B.

**Theorem.**

    H(T) >= d(G_z),

with equality **if and only if T-z is transitive**. In particular, if T-z contains a directed cycle, then

    H(T) >= d(G_z)+2.                                      (1)

Thus the chain-graph adapter of DOMINATION_TRANSFER.md is exactly the zero-defect locus of a general comparison, not merely a convenient example.

Classical inputs used below are Redei's existence of Hamilton paths, the linearly ordered strong-component decomposition of a tournament, and Camion's Hamilton-cycle theorem for nontrivial strong tournaments. Their usual consequences are: a Hamilton path can start at any chosen vertex of the first strong component, and can end at any chosen vertex of the last strong component. Singleton components require no cycle theorem. The oddness of d(G) is inherited from the Brouwer-Csorba-Schrijver theorem recorded in THM-4996; oddness of H is Redei/OCF, THM-002.

## 2. An explicit injection, indexed by the cut at z

For a dominating set S of G_z define

    L=(A minus S) union (B intersect S),
    R=(A intersect S) union (B minus S).                    (2)

These form a partition of V. Domination says exactly:

- every vertex a in L intersect A has an out-neighbour in L intersect B;
- every vertex b in R intersect B has an in-neighbour in R intersect A.

If L is nonempty, the last strong component of T[L] must meet B. Otherwise it consists of A vertices, each of which would have an outgoing edge to B outside that terminal component, contradicting terminality. Therefore T[L] has a Hamilton path ending in B. Similarly, if R is nonempty, the first strong component of T[R] meets A, and T[R] has a Hamilton path beginning in A.

Choose such paths canonically, for example lexicographically with respect to fixed vertex labels, and concatenate

    (path on L), z, (path on R).                            (3)

The edges into and out of z are legal by the definitions of A and B. If a side is empty, omit its path. Empty sides cause no difficulty: a nonempty L with no B would contradict its domination condition, and the dual statement holds for R.

The path determines L as the vertices before z, hence recovers

    S=(A minus L) union (B intersect L).

Thus (3) is an injection from dominating sets into actual Hamilton paths. No equality of unrelated scalar counts and no conjectural conflict-graph realization is used.

## 3. Exact positive defect formula

For L subset V define

- e_B(L): number of Hamilton paths of T[L] ending in B intersect L, with e_B(empty)=1;
- s_A(R): number of Hamilton paths of T[R] beginning in A intersect R, with s_A(empty)=1;
- epsilon(L):1 exactly when the two domination conditions following (2) hold, and0 otherwise.

Splitting every Hamilton path at z and using the bijection (2) between sets S and cuts L gives

    H(T) = sum_{L subset V} e_B(L) s_A(V minus L),
    d(G_z) = sum_{L subset V} epsilon(L),

and hence

    H(T)-d(G_z)
      =sum_{L subset V} [e_B(L)s_A(V minus L)-epsilon(L)].   (4)

**Every summand in (4) is a nonnegative integer.** If epsilon=0 this is immediate; if epsilon=1 it is exactly the endpoint-path existence proof above.

Formula (4) separates two genuine sources of excess:

1. a dominating cut supports several endpoint-path pairs, giving multiplicity excess;
2. a nondominating cut nevertheless supports an endpoint-path pair, giving existence without direct domination.

This is the positive sidecar that the earlier scalar comparison lacked. The endpoint counts, not d alone, determine the excess.

## 4. Equality forces both shores to be transitive

Assume equality in the theorem. Since all summands in (4) are nonnegative, each vanishes.

The dominating set S=V corresponds to L=B,R=A. Its fibre has size

    H(T[B]) H(T[A]).

It must equal1. Since H is a positive integer, both factors equal1, and OCF implies both induced tournaments are transitive: any cyclic triangle would already contribute2 beyond the empty packing.

This argument also handles an empty shore, using H(empty)=1.

## 5. A mixed triangle forces a positive defect

Suppose the shores are transitive but T-z is not. Then there is a mixed directed triangle. Consider first a triangle with two vertices in A and one vertex b in B.

The set S=V minus {b} dominates G_z: all other vertices are selected and b has an A-neighbour from the triangle. Its cut is

    L=B minus {b},     R=A union {b}.

The left endpoint-path count is1 because L is transitive and consists of B vertices. We show that T[R] has at least two Hamilton paths beginning in A.

The strong component C containing the triangle contains b and at least two A vertices. It is the only possible nontrivial strong component, since deleting b leaves the transitive tournament T[A].

- If C is the first strong component, take a Hamilton cycle of C and cut it at either of two distinct A vertices. These give two different A-starting Hamilton paths of C, and concatenate with paths of the remaining components.
- If C is not first, the first component is an A singleton, so every Hamilton path obtained by independently choosing a Hamilton path through C starts in A. The same cycle gives at least |C|>=3 choices.

Thus this dominating fibre has size at least2. At least one summand of (4) is positive, so H(T)>d(G_z).

For a mixed triangle with one A vertex a and two B vertices, use S=V minus {a}; the identical argument with terminal rather than initial strong components produces at least two B-ending paths in L=B union {a}.

Every nontransitive tournament contains a directed triangle. After excluding triangles wholly inside a shore, these two cases exhaust all possibilities. Therefore equality implies T-z is transitive.

Conversely, if T-z is transitive, take its linear order and colour its vertices A or B. The graph G_z and tournament T are exactly the word constructors of DOMINATION_TRANSFER.md, whose proved identity gives H(T)=d(G_z). This proves the iff. Finally both numbers are odd, so strict inequality implies the gap at least2 in (1).

## 6. Exact cancellation criterion for a fixed bipartite input

Now start from an arbitrary bipartite graph G with fixed shores A,B. Keep the cross orientations a->b for graph edges and b->a for nonedges, and attach z->A, B->z. Allow **arbitrary choices of tournaments inside the two shores**.

**Corollary.** Some such choice has H(T)=d(G) if and only if G is a bipartite chain graph.

- If equality is attained, the theorem forces T-z transitive. Its total order represents all edges of G as A-before-B pairs, hence gives nested A-neighbourhoods.
- If G is a chain graph, choose its interleaving order and use the exact adapter.

Therefore for a non-chain bipartite graph, reordering a shore or even replacing its internal orientation by a nontransitive tournament cannot cancel the defect while retaining the input graph and apex interface. Every such orientation has gap at least2. For a chain graph, compatible sorting really does cancel it.

The minimal structural obstruction is an induced2K2. With two independent graph edges and two vertices in each shore, d(G)=9; the receipt checks all four possible shore orientations and records their H-values. This is a graph-class obstruction rather than an inference from unequal finite spectra.

To escape it, a constructor must change the cross-incidence graph, change the marked-vertex interface, or replace the source graph by a different graph of the same dominating count. None of those steps is supplied merely by the solved even-graph domination theorem.

## 7. Consequences for Chan-Pak and what is not proved

This theorem genuinely transports graph domination into tournament counting as an injection and identifies the exact equality locus. It also gives a nonnegative, locally indexed defect formula. It does **not** say every odd d can be attained with zero defect; the graph class required for that is precisely the restricted chain class, which already misses29 as shown in the earlier note.

A promising further interface is to retain the two endpoint-count arrays in (4) and seek a controlled operation changing the excess by a specified amount. The strict fixed-input obstruction above means that an arbitrary even-graph hub move cannot simply be followed by a shore reorientation and expected to preserve d. A new graph replacement or multi-apex interface is needed. No such universal repair is claimed here.

## 8. Controls and status

`bipartite_defect_check.py` enumerates every labelled tournament on orders1 through5 and every marked vertex. For each root it independently enumerates Hamilton paths, computes every cut fibre, enumerates dominating sets of G_z, checks nonnegativity cut by cut, and checks the structural equality criterion. It also checks100 seeded order-six tournaments and every root, plus the complete2K2 shore-orientation case. The receipt states exact counts. These finite controls support the all-order proof; they do not replace it.

The earlier chain adapter and register formulas remain unchanged. This file closes the specifically labelled OPEN inequality/equality question in its section6; it does not modify the other open questions or assert Chan-Pak completeness.
