# Square-sum path overlaps: common blocks and the first cubic response

**Status.** The general block normal form and source-preserving blow-up
lemma below are elementary proofs. Square-guard spectra through order 30
are **FINITE-EXACT**, checked by two separate counting algorithms. The
order-25 formula, general compatible-subset expansion, and its two-block
swap example are **RECOVERED, not new**: they appear in the pinned
`joint_incidence_clocked_frames_20261102/PATH_KERNEL.md` source. No new
all-order arithmetic Hamiltonicity or ordinary tournament-spectrum theorem
is claimed.

## 1. Exactly what is counted

On the labelled vertex set `[n]={1,...,n}`, let `G_n` contain the unordered
pair `{u,v}` precisely when `u != v` and `u+v` is an integer square. A
**directed square-sum path** is a permutation `p=(p_1,...,p_n)` whose
successive unordered pairs belong to `G_n`. Its arc set is
`A(p)={(p_i,p_(i+1)):1<=i<n}`. Reversal is a different path when `n>1`.

For a full labelled tournament `T`, retain the graph guard and put

    H_G(T) = #{p : A(p) is a subset of the arcs of T},
    Z_G(z) = sum_T z^(H_G(T)).

This is NOT the unrestricted tournament Hamilton-path count unless `G` is
complete. In particular even values of `H_G` are allowed. Write
`N=binom(n,2)` and `m=|E(G)|`. Orientations of `G` have exactly `2^(N-m)`
full-tournament extensions; all have the same guarded response.

The inherited expansion is

    Z_G(z) = sum_(S compatible) 2^(N-|U(S)|) (z-1)^|S|,          (1)
    U(S) = union_(p in S) A(p).

Here `S` runs through subsets of the full directed path bank; the empty
term is `2^N`. Compatible means that `U(S)` never contains both `(u,v)`
and `(v,u)`. Necessity and sufficiency are literal: inconsistent arcs
cannot belong to a tournament, while every consistent union extends by
freely orienting all unprescribed pairs. Expanding
`product_p [1+(z-1) 1_(A(p) subset T)]` proves (1). The coefficient of
`(z-1)^j` is `sum_T binom(H_G(T),j)`, not the ordinary power moment.
Pairwise compatibility suffices, but the union-cardinality weight is joint.
The largest compatible subset equals the largest possible `H_G`: extend
its union to obtain the lower bound; the legal paths in any tournament
are themselves compatible, giving the upper bound.

## 2. Canonical common-block normal form

**Theorem 1 (all finite banks).** Let `S` be a nonempty family of distinct
directed Hamilton paths on the same vertex set, not assumed compatible.
Intersect their directed arc sets. Its components, including isolated
vertices, are directed paths. Call these maximal common blocks
`B_1,...,B_b`, each with its inherited internal direction. Then:

1. Every member of `S` is uniquely a permutation of these whole blocks.
2. Let `J` be the set of ordered block pairs occurring as consecutive
   blocks in any of those permutations. Then

       |U(S)| = n-b+|J|.                                      (2)

3. The family is compatible **if and only if** there are no two singleton
   blocks `B_i,B_j` with both `(i,j)` and `(j,i)` in `J`.
4. Consequently the term of `S` in (1), when compatible, is exactly

       2^(N-n+b-|J|) (z-1)^|S|.                               (3)

**Proof.** The intersection is a subgraph of any one directed simple
spanning path, hence a directed linear forest. Every shared arc is a
required consecutive step in every member, so each whole component must
be traversed contiguously and in the same direction. This proves the
unique block permutations. There are `n-b` common internal arcs. Each
ordered block pair `(i,j)` contributes the distinct endpoint arc
`last(B_i)->first(B_j)`, giving (2).

An internal arc cannot conflict with another member, because it belongs
to every member. Opposite connector arcs would require
`last(B_i)=first(B_i)` and `last(B_j)=first(B_j)`: both blocks must be
singletons, and the block orders must use both directed adjacencies.
Conversely those two singleton joins really are opposite orientations of
the same original pair. This proves both directions of assertion 3 and
then (3). QED.

This explains precisely why the endpoint quotient is not generally a
tournament. With nontrivial blocks, `i->j` and `j->i` use DIFFERENT pairs
of original vertices and may coexist. Indeed any family of permutations
of a fixed partition into directed blocks of length at least two is
orientation-compatible, provided every selected join is allowed by the
guard. The theorem does not say those joins exist in a square-sum graph.

**Pair refinement.** For distinct compatible `p,q`, put
`d=|A(q) minus A(p)|=|A(p) minus A(q)|`. Their common-block count is `d+1`,
and their union has `n-1+d` arcs. If `d=1`, the two block orders are `AB`
and `BA`, and their union is a directed Hamilton cycle. Thus in any guard
with no undirected Hamilton cycle, compatible distinct paths require
`d>=2`, equivalently at least `n+1` prescribed union arcs. This bound is
sharp at order 25. A singleton pair would instead conflict, not produce a
valid two-vertex directed cycle. The bound is not asserted for arbitrary
Hamiltonian-cycle guards.

For `n=18,...,30`, vertex 18 is a leaf of `G_n` (its only neighbour is 7),
so those guards have no Hamilton cycle. At 25 the compatible pairs attain
the lower bound `n+1=26`.

## 3. Recovering all four order-25 overlaps, not just one

All-start DFS gives `m=32` and 20 directed paths, ten up to reversal. The
compatibility graph is `4 K2` disjoint union `12 K1`; in particular it
contains no triple. Each of the four edges has union size 26 and shares
22 of the 24 arcs. The two pairs beginning at 18 are `UAB, UBA`, with

    U = (18,7,9,16,20,5,11,25,24).

For the first pair,

    A = (1,3,22,14,2,23,13),
    B = (12,4,21,15,10,6,19,17,8),

and the four joins `U->A,U->B,A->B,B->A` have sums `25,36,25,9`.
For the second pair,

    A = (1,8,17,19,6,10,15,21,4),
    B = (12,13,23,2,14,22,3),

and the same four endpoint sums are `25,36,16,4`. The other two pairs
are simultaneous reversals of these pairs, not four additional pairs.
Every listed block is nontrivial, so Theorem 1 explains compatibility
without collapsing opposite block directions to contradictory arcs.
The exhaustive bank in `normal.json` proves there are no other pairs.

The guard polynomial is therefore

    2^32 + 20*2^8 (z-1) + 4*2^6 (z-1)^2
      = 4294962432 + 4608 z + 256 z^2.

Multiplying by `2^268`, since `N=300`, recovers exactly the quoted full
labelled-tournament polynomial. Its coefficients sum to `2^300`; its
first moment is `20*2^(300-24)`; its second binomial moment is
`4*2^(300-26)=2^276`. Equivalently the whole polynomial is

    2^300 + 20*2^276 (z-1) + 2^276 (z-1)^2.

The coefficient of `z^2` equals the pair moment ONLY because no triple
exists. This equality must not be extrapolated to order 26.

## 4. Order 26: the first compatible triples and a cubic polynomial

At 26 there are `m=34` guard edges and 24 directed square-sum paths.
The compatible-subset union census is

| Subset size | Union arcs : number of subsets |
|---:|---|
|1|25 : 24|
|2|27 : 8; 29 : 2; 30 : 2|
|3|29 : 2|
|4 or more|none|

Thus the guard binomial moments are `2^34,12288,1120,64`, since

    M1 = 24*2^9,
    M2 = 8*2^7 + 2*2^5 + 2*2^4,
    M3 = 2*2^5.

Inverting gives the new finite-exact response

    Z_(G26)(z) = 2^291 (17179857952 + 10240 z + 928 z^2 + 64 z^3)
               = 2^296 (536870561 + 320 z + 29 z^2 + 2 z^3).    (4)

The normalized coefficients in the second line sum to `2^29`; total mass
is `2^325`. The first moment is `24*2^300`. The two triples are reversals
of one another, each fixing 29 arcs. Their intersections have 32 guard
orientations each. Consequently exactly `2^297` full tournaments have
three guarded paths. At orders 2 through 25 no compatible triple exists,
so 26 is the first triple in this independently checked finite range.

A canonical triple is carried by five common blocks

    B0 = (18,7,9,16,20,5,11,25,24),
    B1 = (1,8,17,19,6),
    B2 = (3,13),
    B3 = (12,4,21,15),
    B4 = (10,26,23,2,14,22).

Its three block orders are `01234`, `01423`, `03142`. Their eight
connector arcs are

    0->1, 0->3, 1->2, 1->4, 2->3, 3->1, 3->4, 4->2.

In actual labels these are
`24->1,24->12,6->3,6->10,13->12,15->1,15->10,22->3`, with square sums
`25,36,9,16,25,16,25,25`. There are 21 common internal arcs and eight
joins, giving union 29 exactly. All five blocks are nontrivial, so the
singleton criterion makes their simultaneous compatibility transparent.
This is a three-way endpoint relation, NOT two independent Boolean swaps
(which would have supplied four choices).

**Theorem 2 (source-preserving block blow-up).** Let `D` be a finite
loopless directed graph with a specified vertex `s` of indegree zero.
Replace every vertex `i` by a nonempty internally directed path `B_i`,
and replace every arc `i->j` by `last(B_i)->first(B_j)`. Then Hamilton
paths of the expanded directed graph are in bijection with Hamilton paths
of `D`, by concatenating whole blocks.

**Proof.** The first vertex of `B_s` still has indegree zero, so any
Hamilton path must start there. Every noninitial vertex of every block
has the unique incoming arc from its predecessor within that block, and
cannot be the global start. Those internal arcs therefore belong to any
Hamilton path. Contracting them gives exactly a Hamilton path of `D`;
expanding is the inverse. If there is no Hamilton path both sides are
empty. QED.

For the eight-arc port graph above, vertex 0 is a source, and its ONLY
Hamilton orders are `01234,01423,03142` (enumerating the possible outgoing
choices proves this, also checked literally). Hence the oriented union
of the displayed 26-vertex triple is itself a closed three-path directed
carrier. Arbitrary nonempty subdivisions of its five vertices preserve
that count. Only the displayed arithmetic labelling is proved square-sum;
the abstract blow-up theorem does not give all-size square-sum embeddings.
The source hypothesis matters: expanding a vertex of a directed triangle
makes a directed four-cycle, increasing the number of Hamilton orders
from three to four through cuts inside a block.

## 5. Small-order boundaries and finite continuation

For orders 2 through 30, the only nonempty path banks occur at
`15,16,17,23,25,26,27,28,29,30`. Their respective directed counts are
`2,2,2,6,20,24,70,104,38,40`. Order 1 has its single empty-arc path.
In particular the guarded response at 24 is identically zero:
`Z_(G24)=2^276`. This is not a claim that ordinary tournaments at 24 lack
Hamilton paths.

The table gives the complete GUARD-orientation histogram. To count full
labelled tournaments, multiply the indicated row by `2^(N-m)`.

| n | m | N-m | coefficients for H_G=0,1,2,... |
|---:|---:|---:|---|
|15|15|90|32764, 4|
|16|16|104|65532, 4|
|17|17|119|131068, 4|
|23|28|225|268435072, 384|
|24|30|246|1073741824|
|25|32|268|4294962432, 4608, 256|
|26|34|291|17179857952, 10240, 928, 64|
|27|36|315|68719413784, 54976, 7360, 488, 120, 8|
|28|38|340|274877723488, 158096, 21840, 2928, 528, 64|
|29|40|366|1099511486464, 129024, 10240, 2048|
|30|42|393|4398046212096, 274432, 20480, 4096|

Maximum guarded response is therefore nonmonotone here: `2,3,5,5,3,3`
at orders 25 through 30. This is a finite observation, not an eventual law.
The n25 matching-based complement carrier in the pinned prior source
cannot simply be applied at n26, whose compatibility graph has triangles.

Finally, 20 directed square-sum paths is a cardinality of a labelled path
bank. An unrelated occurrence of 20 proper divisors, or the difference
`704-694=10`, supplies no map between the counted sets and is not evidence
for this polynomial. Ten here has the precise meaning 'paths up to
reversal at order 25'; all counts in (1) use BOTH directions.
