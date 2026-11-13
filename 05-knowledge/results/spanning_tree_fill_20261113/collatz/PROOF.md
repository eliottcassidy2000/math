# Genuine reverse edges and the extra even-quotient shortcut

2026-11-13. **PROVED:** model A is exactly the classical basin question;
model B has an explicit logarithmic-length descent, a linear height buffer,
and a spanning arborescence toward1. The finite boundary identity is exact
without a Collatz assumption. **FINITE-EXACT:** the attached controls.
No conclusion about classical Collatz convergence is inferred from model B.

## 1. The two edge sets

On positive integers let the ordinary, unshortened Collatz map be

    C(n)=n/2 for even n, and3n+1 for odd n.

**A** permits every edge n→C(n) and its reverse. The genuine preimages
of n are exactly

    2n always;
    (n−1)/3 when n=4 mod6.                            (1)

The second preimage is odd. Thus the ordinary preimage graph is at most
binary, not a full binary tree at every vertex.

**B** retains all A edges and additionally permits

    n→(n−1)/3 when n>1 is odd and n=1 mod3.           (2)

Here n=1 mod6 and the quotient is **even**. Equation(2) is not a reverse
edge of C: C on its even quotient halves instead of multiplying by3.
The extra edge is permitted only in its displayed direction; the proof
does not assume its reverse is also available. At n=1 the quotient0 is
outside the positive vertex set, so no such edge exists.

## 2. A preserves the classical reachability question

Let Gamma be the classical basin of the cycle1→4→2→1. Then

    n has a finite A-path to1 iff n belongs to Gamma.          (3)

Indeed x∈Gamma iff C(x)∈Gamma. Membership therefore propagates in both
directions across every A edge. A finite A-path to1 puts its starting
vertex in Gamma. Conversely the classical forward orbit is itself an
A-path. This uses no size or stopping-time assumption.

The underlying undirected graph on Gamma has exactly the cycle{1,4,2}.
Removing its edge{1,4} leaves a tree spanning Gamma, directed toward1
by C with C(1) deleted. Any finite undirected cycle in a functional graph
must be a directed cycle, since each vertex can supply at most one of its
edges. There can be no second directed cycle in Gamma. Consequently this
is a spanning tree of **all** positive integers iff Collatz holds.
Adding genuine reverse edges alone does not weaken that condition.

These are the ordinary-edge versions of the guarded inverse-closure
equivalences in the old [reverse-tree pieces, §4](../../collatz_mod6_20260917_reverse_tree_pieces.md)
and [inverse-tree transport, Theorem1](../../collatz_procgen_20260924_inverse_tree_mod192.md).

## 3. The additional edge gives a complete descent policy

Define a macro-parent D for n>1 by

    n even:               D(n)=n/2;
    n odd, n=0 mod3:      D(n)=n/3;
    n odd, n=1 mod3:      D(n)=(n−1)/3;
    n odd, n=2 mod3:      D(n)=(2n−1)/3.              (4)

Every D(n) is a positive integer and D(n)≤2n/3<n. Each macro edge is
realized by the following literal B walk:

| Case | Walk | Edge count |
|---|---|---:|
| n even | n,n/2 |1|
| odd n=1 mod3 | n,(n−1)/3 |1|
| odd n=2 mod3 | n,2n,(2n−1)/3 |2|
| odd3 divides n | n,3n+1,6n+2,12n+4,4n+1,4n/3,2n/3,n/3 |7|

For odd n=2 mod3,2n=4 mod6 and the second edge is the genuine odd
preimage from(1). In the last row the first edge is forward C, the next
two are reverse even edges, and the fourth is a genuine reverse odd edge
because C(4n+1)=12n+4. The fifth is exactly(2); its input4n+1 is1 mod6
and its quotient4n/3 is positive and even. The final two edges halve.
Every parity and divisibility guard is thus explicit.

The first four edges in that row are the old sibling transport
n↔4n+1, from [inverse-tree transport, Proposition11](../../collatz_procgen_20260924_inverse_tree_mod192.md).
Combining it with the formerly discarded even inverse branch changes the
closure problem: it supplies the missing descent for odd multiples of3.

Iterating(4) reaches1 by strict descent. Hence **every positive integer
has a directed B-path to1**, using the additional edges only as stated.
The last row handles multiples of3 rather than omitting them from the
vertex set.

## 4. Exact growth, an explicit buffer, and a literal spanning tree

The graph with edges n→D(n), n>1, is a canonical rooted tree. Its restriction
to[1,N] is already a spanning tree, and adjoining N+1 adds exactly one
leaf edge to its earlier parent. This is a **macro tree**. Its edges are
walks in B, not necessarily single native edges. It is not binary:

    D(5)=D(6)=D(9)=3.

Each macro walk has at most7 edges and peak at most12n+4. Every later
macro starts below the preceding one. Thus the concatenated route from n
to1 has

    peak≤12n+4,
    length≤7 ceil(log_(3/2) n).                       (5)

The length bound follows because after j macro steps the current value is
at most(2/3)^j n. For n=1 the length is0. Consequently

    every n≤N reaches1 in B induced on[1,12N+4].       (6)

This is a finite spatial-fill certificate with a named boundary budget,
not a density assertion or a conjectural compactness argument. It does
not say the native B prefix[1,N] is connected: at N=3, vertex3 is isolated
while1 and2 are connected.

There is also an actual native-edge spanning arborescence of B toward1.
Let d(n) be the shortest directed B distance to1, finite by(5). For n>1
select the least-valued allowed successor p(n) with d(p(n))=d(n)−1.
Every selected edge is native, and every selected path reaches1. The
underlying selected graph is therefore a tree spanning all positive
integers. This is different from the macro tree and its one-leaf-per-label
growth rule. The union of arbitrarily chosen expanded macro walks need
not itself be a tree.

## 5. Finite induced boundary accounting for the ordinary graph

Let S be any finite set of positive integers. Retain the original C
orientation on edges with both endpoints in S; reverse duplicates are
not counted again. Put

    b(S)=#{n∈S:C(n)∉S},
    c(S)=number of directed C cycles wholly contained in S.

The induced graph is a pseudoforest: every weak component contains
exactly one sink or exactly one directed cycle. Therefore

    components(S)=b(S)+c(S),                         (7)
    edges(S)=|S|−b(S),    cycle rank(S)=c(S).

For a proof, follow the unique internal successor until exiting S or
repeating a vertex. An internal edge joins vertices with the same final
sink or cycle, so two different terminal objects cannot occupy one weak
component. This proves the first identity; the edge identity counts
outgoing arcs, and Euler's formula gives the last. Standard C has no
positive fixed point or two-cycle, so no undirected edge is counted twice.
Deleting one edge of each contained cycle leaves a forest with b+c roots.
This argument does not assume that{1,4,2} is the only positive cycle.

For S=[1,N], even successors always stay inside; exactly the odd n with
3n+1>N are boundary sinks. Hence

    b_N=ceil(N/2)−floor((N+2)/6),
    components([1,N])=b_N+c_N.                       (8)

The cycle count c_N is retained, not silently replaced by1. At N=3,
the two sinks are1 and3, whose genuine future successors are4 and10.
The induced components are{1,2} and{3}. At N=4, the known cycle is
present and3 remains a boundary sink, again giving two components.
Thus native integer-by-integer growth can create an isolated component
whose attachment arrives at a later, larger label. The macro-parent D
has the smaller-parent property precisely because it compresses a walk
that may leave the current prefix.

## Reproduction and remaining scope

Run `python3 -B check.py` and `python3 -O -B check.py` here. The checker
validates every primitive edge of all constructed routes through20000,
tests(7) on all subsets of[1,10], tests(8) through500, and independently
checks the bounded B reachability in(6) for selected cutoffs through1000.

Model B is solved for the specified reachability question. Model A is
still exactly classical Collatz. No theorem about all classical cycles,
classical peak bounds, or the old directed E-SCC conjecture follows.
