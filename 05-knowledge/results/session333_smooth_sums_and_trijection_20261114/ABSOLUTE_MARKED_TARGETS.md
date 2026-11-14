# From a smooth antichain to an exact marked target, not just a residue

The exact unmarked Wave 6 compiler is obstructed by the forced
[paired-frame debt](TWO_HUB_FRAME_DEBT.md). Here are two positive
variants. Every positive integer is the exact count of a marked
Hamilton chamber in a small strong tournament; the second construction
also gives a bipartite graph with one fewer edge than vertex.
Neither construction removes the marking for free.

## 1. Elementary smooth-antichain normal form

Every N>0 has a positive expansion

    N=Σ 2^a_i 3^b_i,    a_1>...>a_t≥0,    0≤b_1<...<b_t.

Here is a greedy proof. Strip any factor of 3, recording its exponent.
For the nonzero residue m modulo 3, choose the largest 2^a≤m whose
exponent parity makes 2^a≡m mod3. Then m<2^(a+2), so
m'=(m−2^a)/3<2^a. If m'=0 stop. Otherwise recur, increasing
the accumulated exponent of 3 by at least one. The next 2-exponent
is strictly smaller by the strict inequality; further stripping
of 3 cannot spoil it. This gives an exact constructive antichain,
with a_1≤log_2 N and b_t≤log_3 N.

This is part of the earlier 3-smooth representations programme,
not a priority claim; compare [Blecksmith–McCallum–Selfridge
(1998)](https://www.tandfonline.com/doi/abs/10.1080/00029890.1998.12004921),
also cited for the more restrictive representation questions in
repository THM-4640. The present elementary greedy procedure imposes
neither a minimum term size nor consecutive b_i. Its results here are

    34 = 16 + 18,
    45 = 18 + 27,
    69 = 24 + 18 + 27 = 24 + 45.

The two latter expansions share an *actual* terminal part in this
normal form; [the unmarked lower bound](TWO_HUB_FRAME_DEBT.md)
still forbids an unmarked realization using the weight-27 frame
(three ternary sites) at those absolute targets.

## 2. Literal three-run paths with a branch mark

Take a forward linear order with exactly two crossing reversed arcs,
0<i<a<b, namely a→0 and b→i. As in THM-4748, restrict to Hamilton
paths using both descents in this order:

    R0 ends at a, R1 starts at 0 and ends at b, R2 starts at i.

Every vertex strictly between 0 and i has two freely available runs,
R0 and R1. Every vertex strictly between i and a has all three
available runs. Take no vertices between a and b. Increasingly sort
each run; the fixed descents supply the joins and the unwanted
reversed pairs never occur inside an increasing run.

For the expansion with t terms, give the first pool size
a_1+ceil(log2 t) and the second pool size b_t. Reserve the first
ceil(log2 t) binary choices as a branch tag. At tag i permit
the first a_i remaining binary sites and first b_i ternary sites
their full choice sets; force other sites into one specified legal
run. Distinct tags select disjoint literal HP subsets; their counts
are exactly 2^a_i 3^b_i. The marked count is therefore N, with only
4+a_1+b_t+ceil(log2 t)=O(log N+1) vertices. The host is strong:
the forward chain and two backward intervals give a path from b
back to 0 and thence to all vertices. The arbitrary marks are
a genuine selection; the unmarked HP count also contains other
tags, assignments, and descent types.

## 3. Trijection edge budget via a permanent circuit

There is also a balanced bipartite graph G with |E|=|V|−1 and
exactly N perfect matchings, on O(log N+1) vertices.

Build a sparse directed acyclic path-count circuit with exactly N
paths from s to t. First a shared doubling ladder supplies all
2^a_i: two fresh intermediate branches from the previous state
and then a join give a new state with twice its path count. A
fresh accumulator, initially fed by the last binary tap a_t,
then runs backwards through the b_i gaps. To multiply by three,
use three such branches in a cell. For each term, after the
required b-gap multiply, add one arrow from the corresponding
binary tap into a fresh accumulator node (and one from the
old accumulator). Finally multiply by 3^b_1. This implements

    3^b_1 (2^a_1 + 3^(b_2−b_1)(2^a_2 + ...)).

All vertices of the DAG are reachable from s and s≠t. Its order
and size are O(a_1+b_t+t+1)=O(log N+1).

For each of its d vertices introduce left and right copies, and
interpret each directed edge as an edge from the left of its tail
to the right of its head. Add the diagonal edges for all labels
except s and add a sole return edge t→s. A perfect matching is a
permutation. It must move s, because s has no diagonal. Every
nontrivial permutation cycle must use the sole return edge, so
there is exactly one and it consists of one s→t DAG path followed
by the return. Thus #perfect matchings = N exactly. The other
labels use their unique diagonal edges.

The resulting graph has 2d vertices and d+|A| edges, where A
is the DAG arc set. Since all labels are reached from s,
e−v+1=|A|−d+1≥0. Add exactly that many disjoint K2 components.
They contribute one matching each, one edge, and two vertices,
so the count is unchanged and |E|=|V|−1 is exact. The graph
remains balanced, and the padding is O(log N+1).

This is a literal new use of the old trijection's "one less edge
than nodes need not mean connected" principle. In fact if the
graph has c connected components, its binary cycle rank is c−1.
Relative to any fixed perfect matching M0, the map M→M Δ M0
injects matchings into that cycle space. Hence N≤2^(c−1);
a connected tree carries at most one matching and disconnectedness
is essential for N>1. Powers of two attain equality: take a
disjoint C4 for each factor 2 and one K2. The ternary atom is
K3,3 minus two disjoint edges: 6−2−2+1=3 matchings, with cycle
rank 2. A disjoint union of a binary atoms and b ternary atoms
has exactly 2^a3^b perfect matchings; adding b+1 neutral K2
components restores |E|=|V|−1. The three choices at a ternary
atom are genuine matching branches, not three independent
binary cycles masquerading as a fourth choice.

## 4. Exact absolute matched paths in a strong tournament

Let the balanced crossing graph have d left vertices O, d right
vertices I. Take a transitive module U of d+1 vertices, set
U→O and I→U, orient an O/I pair O→I precisely when the crossing
graph has that edge (otherwise I→O), and complete O and I each
as transitive blocks. The host is strong: a perfect matching
supplies a return through every outside vertex. A marked
Hamilton word is required to have d+1 U visits (all singletons),
the prescribed order of U points, and the prescribed order of
O points. Such a word is necessarily

    u_0,o_0,i_sigma(0),u_1,...,o_(d−1),i_sigma(d−1),u_d,

one per perfect matching, so this strong tournament again has
exactly N marked Hamilton paths on O(log N+1) vertices.
Unmarking the visit orders already introduces factorial factors,
and paths with fewer U visits survive. The previous Matching
Sections theorem shows exactly why the normalized maximal
coefficient is N; neither theorem proves that this marked
whole count is the ordinary Hamilton count.
