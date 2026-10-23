# Palindromic square clocks: independent structural audit

**PASS. The first nontrivial palindromic square-sum Hamiltonian clock occurs at
N=48.** The empty clock at N=1 is the explicitly excluded trivial case.

This proof uses the previously audited adjacent-sum reconstruction theorem:
odd N has one Hamiltonian realization; even N has at most two, and a double
fibre is exactly a divisor-block twin. No unbounded Hamiltonian census is
claimed.

## 1. Quotient reduction with a retained central mark

For N=2m and d dividing m, let A_d be the alternating d-block lower halves,
and let tau_d swap each vertex a in A_d with a+d. Define the simple quotient
graph J_(N,d) on A_d by

    a~b iff a!=b and a+b+d is a square.

Also mark a precisely when 2a+d is a square. These marks are not ordinary
quotient edges: they encode the possible middle edge between the two members
of the same tau_d pair.

There is a Hamiltonian path of the square-sum graph with palindromic ordered
adjacent-sum word and twin displacement d **if and only if** J_(N,d) has a
Hamiltonian path ending at a marked vertex.

Proof. If the sum word is palindromic, reversal is a second realization, so
it must equal tau_d applied to the original path. The first m vertices
therefore contain exactly one member of each tau_d pair. Replace each by its
representative a in A_d. This gives a permutation of A_d. The first-half
positions alternate block sides, so an adjacent pair sums to a+b+d regardless
of the chosen initial side. Its final vertex a is followed in the full path
by its partner, giving middle sum 2a+d. These are exactly the quotient edges
and endpoint mark.

Conversely choose a quotient Hamiltonian path and either initial block side.
Lift its representatives to alternating sides to obtain a first half H.
Append tau_d(reverse(H)). All vertices occur once, all adjacent sums are
squares, and the sum word is palindromic.

Each oriented quotient path ending at a mark yields two reversed full paths,
with one common palindromic word. The at-most-two reconstruction theorem makes
this correspondence injective on words. Thus the exact number of square
palindromic words of order 2m is the sum, over d|m, of the numbers of oriented
quotient Hamiltonian paths ending at a marked vertex. This is an exact
reduction, not a closed-form square-subgraph enumeration.

## 2. Minimality below 48 needs only small structural certificates

Odd orders N>1 cannot have palindromic words: the distinct reversed path would
contradict unique reconstruction. For every even N<48 and every d|N/2, first
retain only cross-block edges in the original square-sum graph. All but two
of the 76 candidate (N,d) pairs have an isolated vertex, more than two leaves,
or a proper connected component. The remaining candidates are (36,6) and
(40,4).

For (40,4), the quotient has **no marked vertex**, so a middle square edge
cannot exist.

For (36,6), the quotient vertex set is

    {1,2,3,4,5,6,13,14,15,16,17,18,25,26,27,28,29,30}.

Deleting vertices {1,14} leaves exactly four components:

    {2,4,6,13,15,17,26,28,30},
    {3,16,27},
    {5,18,25},
    {29}.

A Hamiltonian path, after removal of any two vertices, has at most three
nonempty path segments. Every remaining graph component must contain a
segment, so four components are impossible. This excludes ALL quotient
Hamiltonian paths, not merely those ending at its marked vertices {5,15,29}.

The independent memoized subset recurrence also confirms zero marked-endpoint
paths at (36,6) in 442 states, but this dynamic program is redundant: the
two-vertex separator is a complete structural certificate.

## 3. Positive order-48 witness

For N=48,d=8 the independently checked full path is

    3,46,18,31,5,44,37,27,22,42,7,9,40,41,8,28,21,43,6,10,39,25,24,12,
    4,32,17,47,2,14,35,29,20,16,33,48,1,15,34,30,19,45,36,13,23,26,38,11.

It is a permutation of {1,...,48}. Its ordered square clock is

    7,8,7,6,7,9,8,7,8,7,4,7,9,7,6,7,8,7,4,7,8,7,6,4,
    6,7,8,7,4,7,8,7,6,7,9,7,4,7,8,7,8,9,7,6,7,8,7.

This is palindromic. Reversal is exactly tau_8 applied to the path. The first
half descends to a quotient Hamiltonian path whose endpoint representative
is 4, and 2*4+8=16 supplies the middle square edge.

## 4. Audit evidence and boundaries

The checker reconstructs all 76 lower-order candidate graphs, verifies their
obstructions, checks the stronger (36,6) separator, verifies the literal
order-48 witness, and independently checks the quotient lifting equivalence
for every pair ordering and both initial sides through N=8. Its normal and
optimized outputs byte-match, using explicit failures instead of assertions.

The first general square-clock twin occurs at 36, while the first nontrivial
palindromic clock occurs at 48. These are different questions. Neither is the
threshold for existence of a square-sum Hamiltonian path, and neither implies
that every larger order has the corresponding feature.
