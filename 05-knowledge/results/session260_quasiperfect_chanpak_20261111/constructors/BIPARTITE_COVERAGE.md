# Bipartite domination: the extra obstruction43 and lawful boundary-state arithmetic

**Session260 third wave, 2026-11-11. PROVED WITH FINITE EXACT BASE: the all-order spectrum below225. PROVED: leaf, false-twin and two-shore-hub identities. OPEN: the full bipartite spectrum and any general tournament adapter.**

## 1. Inherited results and what is actually different

We use THM-4996 (`dominating-set-counts-of-even-graphs-miss-exactly-3-7-13-below-512...`), THM-5040 (`sharp-even-graph-chan-pak`), THM-5055 (`even-graph-chan-pak-theorem...`) and THM-5056 (`odd-graph-chan-pak...`) as inputs. They settle the graph domination spectra in the stated even/odd automorphism classes. The unrestricted all-odd threshold-graph construction is the cited Brouwer–Csorba–Schrijver theorem, not a new result here.

**Bipartite is a different restriction.** The cones and joins used in those
proofs need not preserve a bipartition. Their scalar arithmetic cannot be
imported as a bipartite constructor. Also “even graph” is the automorphism
edge-reversal notion, not even order, even degree, or bipartiteness.

The initiating conjectural guess “all odd numbers except7,13,21” is false:
**43 is also impossible.** The exact new bounded-range statement is

> Over bipartite graphs of arbitrary order, the dominating-set counts
> below225 are exactly the positive odd numbers other than7,13,21,43.

This does not assert that those four are the only exceptions at larger
counts. The finite evidence is compatible with that question, but supplies
no uniform induction.

Research board: bipartite support / rooted domination states / false twins /
shore hubs / chain-graph equality / positive tournament defect. The useful
change is to retain domination boundary states, not recycle the threshold
binary recurrence under a new class name.

## 2. Why a finite bipartite census gives an all-order theorem below225

For a rooted tree T with root v let

- A count sets selecting v and dominating every vertex;
- B count sets not selecting v but dominating every vertex;
- C count sets not selecting or dominating v, while dominating all other vertices.

The singleton state is(1,0,1). Join two rooted trees by an edge, retain the
first root, and let their states be(A,B,C) and(a,b,c). The exact new state is

    ( A(a+b+c),   B(a+b)+Ca,   Cb ).                       (1)

Each formula follows by splitting on whether the child root is selected.
For example, an undominated first root requires it to have been in state C
and the child root to be unselected but internally dominated (state b).
Every rooted tree with at least two vertices decomposes this way by cutting
an edge incident to the root. Hence taking all splits of n, and deduplicating
only equal triples, is a complete finite state compiler for every rooted
tree of order n. It does not rely on an isomorphism census.

The exact results through order10 are:

| n | Distinct rooted triples | Minimum A+B |
|---|---:|---:|
|1|1|1|
|2|1|3|
|3|2|5|
|4|4|9|
|5|8|15|
|6|17|25|
|7|36|45|
|8|77|75|
|9|172|125|
|10|380|225|

These minima are finite exact results; no extrapolation of their apparent
powers-of-five pattern is needed.

Adding edges never destroys a dominating set. Every connected ten-vertex
graph contains a spanning tree, so it has at least225 dominating sets.
Every larger connected graph has a connected induced ten-vertex prefix.
Adding a vertex with a neighbour cannot decrease d: map each old dominating
set S to S plus the new vertex. Thus every connected graph of order at least10
has d≥225.

It follows that each nontrivial connected component of a graph with d<225
has at most9 vertices. Isolated vertices contribute factor1. Domination
counts multiply over disjoint components, so the entire arbitrary-order
problem below225 reduces to component products of the bipartite census
through order9.

`bipartite_spectrum.cpp` enumerates every0/1 adjacency matrix for shore
sizes p and n−p, 1≤p≤floor(n/2), n=2,...,9. It includes disconnected graphs
and isolates; duplicate colourings or isomorphisms are retained. Every
bipartite graph of the relevant order has at least one such presentation.
For each graph it counts dominating sets directly by closed-neighbourhood
unions. Every attained value has a literal witness encoded as(n,p,code),
where bit i(n−p)+j is the cross edge i—(p+j).

Taking component products below225 gives exactly the complement of
{7,13,21,43} among the odd numbers in that range. In fact every attained positive odd
value in that range is already present on at most9 vertices, so the receipt
provides direct witnesses without needing products.

### The smallest new obstruction can be isolated much more cheaply

43 is prime, so a graph with d=43 has exactly one nontrivial connected
component. The seven-vertex tree minimum is45, which excludes connected
components of every order at least7 by the same spanning-tree/prefix
argument. The complete bipartite census through order6 has no43. This
already proves its all-order exclusion. Thus43 is not merely an unfilled
entry in a nine-vertex search.

The analogous low-count reasoning excludes7,13,21. The inherited general
graph uniqueness results identify7 and13 with a triangle and a diamond,
respectively, neither bipartite. The exact component-product computation
also handles21; it cannot be rescued by3·7 since7 is excluded.

## 3. Leaf growth preserves bipartition but blocks the naive +1 digit

Attach one new leaf to a root v of any graph, retaining v as root. The same
three boundary counts give

    (A,B,C) -> (2A,B+C,0),
    d_new = 2d + C−B.                                    (2)

If v has a neighbour u, adding u to a C-set injects it into a B-set. The
original C-set contains no neighbour of v, so u was absent and can be
uniquely removed afterward. All other vertices were already dominated;
adding u also dominates v. Hence C≤B. Since dominating-set counts are odd,

    d_new ≤ 2d−1                                         (3)

for every nonisolated root. In particular, a leaf birth at such a root
**never** implements d→2d+1. At an isolated root A=C=d and B=0, so the same
operation gives3d instead.

After adding k≥1 leaves at the same root, the exact state is

    (2^k A, B+C, 0),     d_k=2^k A+B+C.

This gives positive infinite families but keeps a nonnegative offset. It
cannot be mistaken for a freely chosen binary digit. Moving the root to
the newly attached leaf instead gives

    (A,B,C) -> (A+B+C,A,B),

which retains the boundary memory rather than erasing it.

All these moves preserve bipartiteness by placing the new leaf on the
opposite shore from its parent.

## 4. False twins supply a genuine +1 family, but keep a rigid coefficient

Replace vertex v by m≥1 independent false twins with its original
neighbourhood. This preserves bipartiteness by placing all twins on v's
shore. On the exterior V\{v}, define

- α: sets S meeting N(v) for which S∪{v} dominates the original graph;
- A0: such sets S avoiding N(v);
- B1: sets S dominating the original graph without selecting v.

Then the exact count is

    d_m=(2^m−1)α+A0+B1 = 2^m α+β,
    β=A0+B1−α.                                         (4)

If S meets N(v), every nonempty selection of clones is allowed, giving
(2^m−1)α. If it does not, every clone must be selected, giving A0.
Selecting no clone gives B1. These cases are disjoint and exhaustive.
In particular

    d_(m+1)=2d_m−β.

β may be negative, unlike the fixed-root leaf offset. For complete
bipartite graphs, direct counting gives

    d(K_(p,m))=(2^p−1)(2^m−1)+2
              =(2^p−1)2^m+3−2^p,       p,m≥1.         (5)

The p=1 family is2^m+1 and follows2d−1. The p=2 family is

    3·2^m−1 = 5,11,23,47,95,...,

and follows **2d+1** indefinitely within bipartite graphs. These elementary
complete-bipartite families are not claimed as new in themselves. Their
role is to distinguish a lawful positive digit operation from an illicit
leaf/cone inference: the +1 chain fixes oddpart(d+1)=3 and does not yield
arbitrary binary words or permit an unproved return to a minus state.

At d=5 both a leaf move to9 and a twin move to11 exist. This local branching
must not be extrapolated to every subsequent target:11→21 cannot occur
within the bipartite class, by the exact obstruction above.

## 5. Shore hubs have a six-coordinate closed arithmetic interface

Let G have fixed **nonempty** shores L,R, possibly with isolates. Define

    F_L=sum_(S⊆L) 2^|N(S)|,
    F_R=sum_(S⊆R) 2^|N(S)|,
    P=2^|V(G)|,
    ε_L=1 iff no vertex of L is isolated,
    ε_R=1 iff no vertex of R is isolated.

A new vertex on L adjacent to every old R vertex has the exact update

    d'   = d+F_R−ε_R,
    F_L' = F_L+P,
    F_R' = 2F_R−1,
    P'   = 2P,
    ε_L' = ε_L,       ε_R'=1.                            (6)

For d', sets selecting the new hub contribute F_R: choose the old R subset,
then every old L vertex not hit is forced into the dominating set, and the
other old L vertices are free. Sets not selecting the new hub must have a
nonempty R selection. Among old dominating sets, exactly ε_R of them have
empty R selection, namely all of L. The F updates follow by whether the
summed subset is empty or includes the new hub. The R-hub rule is symmetric.

Adding an isolated vertex on L instead gives

    (d,F_L,F_R,P,ε_L,ε_R)
       -> (d,2F_L,F_R,2P,0,ε_R).                         (7)

Thus isolated padding plus shore hubs is a lawful finite-dimensional
arithmetic machine. It is not the one-register threshold cone recursion.

### Exact arbitrary two-shore family

Add a≥1 new L hubs and b≥1 new R hubs, each adjacent to every vertex on the
opposite shore, including the new vertices there. Then

    d_(a,b)=P(2^a−1)(2^b−1)
             +(F_R−1)(2^a−1)+(F_L−1)(2^b−1)
             +d+2−ε_L−ε_R.                             (8)

This is a genuine infinite bipartite constructor for every starting G.
It follows either by iterating (6), or by splitting according to whether
new vertices on each shore are selected. For example with no new R
vertices selected, a nonempty old R selection permits any nonempty choice
of L clones, but an empty old R selection forces **all** L clones to be
selected. That forced-all boundary is the source of the −1 corrections.
Dropping it produces an incorrect count and can even falsely manufacture21.

For a single L batch the corresponding formula is

    d_(a,0)=d+(2^a−1)(F_R−1)+1−ε_R,       a≥1.           (9)

The constant and correction differ from simply multiplying the one-hub
increment by2^a−1. The full boundary state matters.

Equations(6)–(9) give explicit recurrences, integer sidecars, and positive
families. No interval-coverage or all-odd surjectivity theorem follows from
their dimension or positivity alone.

## 6. What this contributes to the tournament route, and what it cannot do

The companion bridge lane proves H(T)≥d(G_z), with equality exactly when
the root-deleted tournament is transitive. Its fixed bipartite graph has a
zero-defect orientation exactly in the chain-graph case. We inherit that
scope, rather than claim every bipartite graph has a count-preserving
orientation.

Consequently even proving a four-exception bipartite spectrum would leave
an adapter problem: non-chain witnesses may have a strictly positive
Hamiltonian surplus at every such orientation. The extra exclusion43 also
rules out the initial three-exception guess immediately. These are distinct
obligations: support coverage, availability of a marked boundary state,
and vanishing of the tournament defect.

The concrete remaining question is whether leaf/twin/shore-hub operations,
possibly after graph replacement, admit a finite set of boundary types with
an exact terminating arithmetic descent for every odd target outside a
finite exceptional set. Equations(2),(4),(6) retain the needed state rather
than assume the desired return. A second, independent obligation is a
constructor that controls the positive tournament defect. Neither is solved
here; no conjectural universal spectrum or adapter is promoted to a theorem.

## Reproduction and evidence scope

    g++ -O3 -std=c++17 05-knowledge/results/session260_quasiperfect_chanpak_20261111/constructors/bipartite_spectrum.cpp -o /tmp/bipartite_spectrum260
    /tmp/bipartite_spectrum260 9
    python3 -B 05-knowledge/results/session260_quasiperfect_chanpak_20261111/constructors/bipartite_boundary.py

The committed spectrum receipt contains every attained value by order and
one explicit bipartite matrix witness per value. The boundary receipt
contains the complete rooted-tree state counts/minima, the component-product
closure below225, and literal controls on the new move interfaces. The
all-order cutoff is proved above, while the order≤9 matrix enumeration and
order≤10 tree-state enumeration are finite exact computations. No external
package, approximate arithmetic, or assumed isomorphism classification is
used. The program does not search higher orders or claim an unbounded census.
