# Complement pairs need an observer; their full incidence can retain much more

23 October 2026. **PROVED elementary:** the all-size adjacent-swap obstruction,
its exact twin boundary, the constant-transfer obstruction, incidence balance
and triangle expectation. **FINITE-EXACT:** the full census through order seven.
**OPEN:** whether the distribution defined below reconstructs every graph up to
complement. No literature-priority claim or general tournament-wiring theorem.

## Inheritance and the repaired objects

The closest mechanism is the observer chart in
[the trijections and tree-gauges packet, sections 1--4](../../tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md).
For a simple graph G on V and an actual total order P, let R_P point from earlier
to later and define T(G,P) by reversing exactly the pairs in E(G). Thus an edge
of G points backwards and a missing edge points forwards. This is a lossless
same-size correspondence **with P retained**. Complementing G or reversing P
each gives the converse tournament; doing both leaves T literally unchanged.

There is a separate inherited correspondence between graphs on ordered n
vertices and tournaments on n+1 vertices with a marked Hamiltonian path, in
[the marked-path trijection](../../trijection_reset_20261020/graph_audit.md).
Its complement is the off-path operation K_P, not full converse. The distinction
is also explicit in [THM-4733, divisor-plane/reflection dictionary](../../../../01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md).
Neither construction is a graph-complement-pair/tournament-type bijection.

The canonical hostile is a graph with its sole edge uv and an isolated w.
The orders (u,v,w) and (u,w,v) produce, respectively, a transitive tournament and
a directed triangle. Their graph/complement pair is identical. Both tournaments
are self-converse, so forgetting which member of the pair was called G does not
repair the ambiguity. At n=4 there are already six graph-complement pairs but
only four tournament types (three after merging converses). An arbitrary
assignment of a tournament type to each pair is possible; reversible transfer
of the objects' theories does not follow from such an assignment.

There is also a direct naturality obstruction for n>=2. The named pair
{empty graph, complete graph} is fixed by every vertex permutation. An intrinsic
equivariant orientation on those same vertices would therefore have to be fixed
by every transposition. A transposition reverses the arc joining its exchanged
vertices, so no tournament can have that symmetry. Selecting an unlabelled type
by an arbitrary canonical convention is a different operation.

The least-used sidecar is the actual observer order, rather than merely an
orientation-parity bit. The live concept board is observer order / graph twins /
directed triangles / shared owners / incidence distributions. Anchor: lower the
sidecar honestly. Niche: replace an impossible single-valued transfer by exact
weighted incidence. Wildcard: test whether that incidence reconstructs the
complement pair.

## 1. Every nontrivial graph has an order-sensitive tournament

For distinct u,v define

    d_G(u,v) = #{w outside {u,v}: [uw in E(G)] != [vw in E(G)]}.

**Theorem 1.** There are two orders differing only by swapping adjacent u,v
for which the directed-triangle counts of T(G,P) differ by exactly d_G(u,v).
Consequently the triangle-count range over all orders has width at least
max d_G(u,v).

**Proof.** Put e=[uv in E(G)]. For each w with unequal incidences a=[uw in E(G)]
and b=[vw in E(G)], put w before the adjacent pair if a!=e, and after it if a=e.
Put all other vertices arbitrarily outside the pair. Compare P with u before v
to Q with v before u. Exactly the tournament arc uv changes. Each triple uvw
with unequal incidences changes from transitive to cyclic, by the stated choice
of side; equal incidences give no cycle in either order. Every other triple is
unchanged. The increase is exactly d_G(u,v). QED.

**Corollary 1.** The following are equivalent for any graph on n>=2 vertices:

1. T(G,P) has the same tournament isomorphism class for every P.
2. It has the same converse-merged isomorphism class for every P.
3. Its directed-triangle count is independent of P.
4. G is complete or empty.

For (3)->(4), Theorem 1 makes every d_G(u,v) zero. For every three vertices their
three graph bits must then agree. This forces all graph bits to agree; the
n=2 case is automatic. Complete and empty graphs produce transitive tournaments
under every order. Converse preserves the number of directed triangles, proving
the remaining implications. Thus **every** nontrivial graph/complement pair,
not merely a small hostile, fails the proposed order-free single-class reading.

### Exactly which adjacent swaps are universally safe?

Two graph vertices are twins if their incidences to all other vertices agree.
Swapping adjacent twins in P changes only uv; the vertex transposition (u v)
is a tournament isomorphism between the two results, because their other
tournament incidences agree. Theorem 1 proves the converse with the necessary
quantifier: a non-twin pair has at least one context where its adjacent swap
changes the tournament class. It may still be harmless in a special context.

This permits compression of a word P to its word of twin classes, forgetting
the order inside each occurrence class. For twin class sizes m_i there are
n!/prod(m_i!) such words. This is a safe upper bound on the retained order
states, not a claim that different words produce nonisomorphic tournaments.
Automorphisms may make further identifications.

If actual ports are marked, the implementing transposition must fix them.
Unmarked twin vertices can be exchanged; a port cannot be silently exchanged
with its twin. More generally the number of **ordered-graph** representations
modulo allowed graph automorphisms is exactly n!/|Aut(G; ports)|, because the
action on total orders is free. This is not the number of unmarked tournament
outputs. On named vertices, G and T together recover every pair comparison in
P, so P has n! distinct possibilities; a fixed-length lossless code for that
literal observer requires at least ceil(log2(n!)) bits.

## 2. A universal unqualified graph-property transfer is necessarily constant

**Theorem 2.** Suppose f is a graph statistic and F is a tournament statistic
such that f(G_P(T))=F(T) for every tournament T and every observer P, where
G_P(T) is its backward graph. Then f is constant on all graphs of the fixed
order n. The same conclusion applies if f is a statistic of complement pairs.

**Proof.** Take any G and pair uv. Choose P with u,v adjacent and set T=T(G,P).
Swapping them to obtain Q changes the backward graph of this fixed T from G to
G xor {uv}. Therefore f(G)=f(G xor {uv}). Single-edge changes connect the entire
Boolean graph cube. Its complement quotient is connected as well. QED.

This does not prohibit useful translated predicates: the inherited cyclic
triangle formula, path-cut coverage and attachment constraints retain their
observer coordinates. It says that dropping those coordinates cannot silently
turn *every* graph theorem into a theorem about a single tournament type.

## 3. The positive replacement: an exact incidence kernel

Let G and T now denote unlabelled graph and tournament types on n vertices.
For a fixed representative of G put

    b(G,T) = #{orders P: T(G,P) is isomorphic to T} / |Aut(G)|.

This is a nonnegative integer. It counts isomorphism classes of pairs consisting
of a graph and a total order whose encoded tournament has type T. Counting the
same pairs from the other side proves the exact balance law

    b(G,T) = #{orders Q of T: G_Q(T) is isomorphic to G} / |Aut(T)|.

Consequently

    sum_T b(G,T) = n! / |Aut(G)|,
    sum_G b(G,T) = n! / |Aut(T)|,
    b(complement G,T) = b(G,converse T).

All these statements follow from simultaneous isomorphism of the **ordered**
objects; the free action on total orders explains the divisors. This is not a
bijection between the two unmarked class sets.

For the unordered complement pair C={G,complement G}, define mu_C as the
distribution of the converse-merged type [T(G,P)] under a uniformly chosen P.
The last symmetry law makes mu_C independent of which G was selected.
Converse-invariant statistics can be averaged directly against mu_C. For an
arbitrary tournament-class statistic, first lift each non-self-converse orbit's
mass equally to its two types. This is canonical and recovers the full observer
distribution: reversing P pairs every occurrence of T with an occurrence of
its converse. With that lift understood, a tournament property holds for every
compatible order precisely when its indicator has expectation one; existence
corresponds to positive mass.
It does not select an actual P or shared ports, so it cannot be used as an
attachment certificate by itself.

**FINITE-EXACT.** The distributions distinguish every graph-complement pair
through n=7, even though there are more such pairs than target tournament types:

| n | Graph types | Graph/complement pairs | Tournament types | Converse-merged types | Distinct mu_C |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 1 | 1 |
| 2 | 2 | 1 | 1 | 1 | 1 |
| 3 | 4 | 2 | 2 | 2 | 2 |
| 4 | 11 | 6 | 4 | 3 | 6 |
| 5 | 34 | 18 | 12 | 10 | 18 |
| 6 | 156 | 78 | 56 | 34 | 78 |
| 7 | 1044 | 522 | 456 | 272 | 522 |

The cheaper statistic consisting only of the **histogram of directed-triangle
counts** over observer orders already separates exactly the same 522 pairs
through n=7. Thus full tournament-type information is unnecessary in this
tested range. This gives two precise next targets: does the triangle histogram,
or failing that mu_C, determine C in every order? Both are **OPEN**, not inferred
from the census. A single tournament type fails immediately, while the weighted
ensemble has retained considerably more information in every tested order.
The [incidence-flow extension](INCIDENCE_FLOW.md) combines this packet with
incoming all-size connectedness to give an exact reversible chain and an
electrical conditional-variance identity. It also finds a sharp limitation:
the star's observer distribution equals a one-quarter empty-source plus
three-quarter matching-source mixture. Distinct individual rows therefore
do not make arbitrary source mixtures reversible.

### A cheaply computable statistic of this incidence

Let d_v be the graph degrees, and c_3(T) the directed-triangle count. Then

    E_P c_3(T(G,P)) = (1/6) sum_v d_v (n-1-d_v).

Indeed, an induced triple with zero or three graph edges is always transitive.
A triple with one or two edges is cyclic in exactly two of its six orders.
The number of mixed triples is (1/2) sum_v d_v(n-1-d_v), since each mixed triple
has exactly two vertices incident with one edge and one missing edge. This
proves the formula. It is complement-invariant and supplies an immediate
necessary test for equality of mu_C, but is not a reconstruction theorem.

### Vertex insertion is an exact cut-energy response

Write Z_G(z)=sum_P z^{c_3(T(G,P))}, the observer triangle polynomial. If an
isolated new graph vertex is inserted at a cut A|B of P, then

    c_3(new encoded tournament) = c_3(T(G,P)) + e_G(A,B).

The only newly cyclic triples have the new vertex between the endpoints of
an edge of G. Every other new triple is transitive. Thus

    Z_(G+isolated vertex)(z)
      = sum_P z^{c_3(T(G,P))} sum_(cuts of P) z^{e_G(A,B)}.

For a universal new graph vertex, replace e_G(A,B) by |A||B|-e_G(A,B), the
number of crossing nonedges. This is a literal link to cut/electrical energy:
if s is the 0/1 shore indicator and L_G the ordinary graph Laplacian, the
increment for an isolated vertex is s^T L_G s. An arbitrary new neighborhood
requires the corresponding marked-cut response rather than this special
edge/nonedge formula.

The right-hand side retains the correlation between the old triangle count
and the actual cut; a histogram of cuts taken independently would lose it.
The complementary two-vertex graphs both have Z=2, but adjoining an isolate
gives Z=6 for the empty graph and Z=4+2z for the one-edge graph. Thus an
unoriented complement pair does not make the *isolated-vertex growth operation*
well-defined: complement sends isolated insertion to universal insertion.
This is an explicit operation that the proposed identification must retain,
and an exact place for the user's newly arriving square-sum vertices to keep
their actual neighborhood/cut data.

### Arbitrary neighborhoods have an addressed score-vector transition

The parent lane supplied the following generalization. Let T=T(G,P), insert
a new graph vertex w at the cut L|R of P, and prescribe its actual graph
neighborhood N. Its outneighbor set in the new tournament is exactly

    O = (L intersect N) union (R minus N).

Every new cyclic triangle has the form w->u->v->w with u in O and v outside O.
Therefore its increment is

    Delta c_3 = e_T(O,V minus O)
              = sum_(u in O) d_T^+(u) - binom(|O|,2).             (2)

The second equality subtracts the one internal outarc of each unordered pair
inside O from the total outdegrees of vertices in O. Simultaneously the exact
score update is

    d_new^+(u) = d_T^+(u) + [u outside O] for old owners u,
    d_new^+(w) = |O|.                                          (3)

Thus the **owner-indexed score vector**, together with the actual cut and
neighborhood, is sufficient for the next triangle count and score vector.
No enumeration of all old pair orientations is needed for this observable.
It is not sufficient to reconstruct the full tournament or its wiring.

The smallest address hostile has two old vertices with arc a->b and scores
(1,0). Insert w at the end and give it the singleton graph neighborhood {a}
or {b}; then O is respectively {a} or {b}. In both cases the old score multiset
and |O|=1 agree, but (2) gives increments 1 and 0. Therefore an unindexed score
multiset and a neighborhood cardinality cannot replace the owner membership.
This is exactly the same kind of address retention required by the Fibonacci
and shared-port constructions. The checker tests (2)--(3) for every graph,
observer, insertion cut and new neighborhood through four old vertices.

## 4. Connection contract and next use

| Source | Target and map | Preserved | Lost / required sidecar |
|---|---|---|---|
| G with actual order P | T(G,P), reverse the graph pairs in R_P | All pair incidences; exact complement/converse and relabelling laws | Nothing if P retained; unmarked graph alone loses the observer |
| Graph plus order modulo allowed automorphisms | Ordered tournament representation | Exact isomorphism and fixed-port constraints | Actual labels need the implementing permutation |
| Graph/complement pair | Distribution mu_C | Averaged, existential and universal compatible tournament-type predicates | The sampled order and joint alignment across attachments |
| Twin-class word | Tournament type | Safe swaps of adjacent unmarked twins | Actual marked owners cannot be forgotten |

This matches the Fibonacci address lesson: equal total lengths do not identify
internal addresses. Here equal graph/complement membership does not identify
the observer that directs the pair. It also matches the attachment lesson:
marginal distributions on individual blocks do not synchronize their shared
cell. Joint owner/order variables must be retained before projection. No
claim that the distribution itself resolves general larger-cell kernels follows.

## Reproduction and validity gate

From the repository root:

```text
python -B 05-knowledge/results/attachment_lattices_square_sum_20261023/complement/check.py
python -B -O 05-knowledge/results/attachment_lattices_square_sum_20261023/complement/check.py
python -B 05-knowledge/results/attachment_lattices_square_sum_20261023/complement/cut_response.py
```

The checker canonizes every labelled graph and tournament through n=7 by the
entire symmetric group, then examines all orders of every graph isotype. It
checks the constructive delta for every pair in every graph isotype, twins as
positive controls, complete/empty extremes, complement/converse, the exact
triangle expectation, and incidence by independent counting from both sides.
The saved receipt gives the bounded universe and small exact incidence matrices.
No assertion disappears under optimized Python. The all-size proofs above do
not depend on extrapolating these computations.
