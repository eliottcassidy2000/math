# Joint path incidence: an exact square-sum response spectrum at 25

**PROVED:** the elementary intersection compiler and marked counting formulas.
**FINITE-EXACT:** the square-sum responses below, by two distinct algorithms.
**INHERITED:** the small square-sum path census and tree-edge trijection.
**OPEN:** a faster unmarked A000568 recurrence and general strong-atom closure.

## 1. Keep the allowed-edge guard

For a specified labelled simple graph G on n vertices, let P(G) be its
directed Hamiltonian orders (each undirected path has two orientations).
For a tournament T on those same vertices, define

    H_G(T)=#{P in P(G): every successive arc of P belongs to T}.

This is an actual tournament response with G as a retained guard. It is
ordinary H(T) only when G=K_n. A square-sum graph having a Hamiltonian path
does not turn it into a tournament, or force every tournament to contain
one of its allowed paths. The independent relation G is indispensable.

Write N=binom(n,2). Each specified P prescribes n-1 distinct arcs, so its
completion cube contains2^(N-n+1) tournaments. For any finite set Q of
paths, let U(Q) be the union of their prescribed arcs. Their cubes intersect
iff no unordered pair is prescribed in both directions; if compatible,
the intersection has size

    2^(N-|U(Q)|).                                      (14)

Compatibility is pairwise, but the union size is genuinely joint incidence.
Pairwise intersection cardinalities alone do not in general recover it.

## 2. The all-finite-configuration compiler

Define the integer spectrum polynomial Z_G(z)=sum_T z^H_G(T). Then

    Z_G(z)=sum_(Q subset P(G), Q compatible)
                          (z-1)^|Q| 2^(N-|U(Q)|).      (15)

**Proof.** For each T, expand
z^H_G(T)=product_(P in P(G))[1+(z-1)1_(P subset T)].
Sum the expansion over T and apply (14). QED.

Equivalently the jth factorial moment sum_T binom(H_G(T),j) is the sum
in (14) over compatible j-path subsets. This is an elementary indicator
expansion, not a novelty claim for inclusion-exclusion. Its use here
keeps source paths, shared arcs, and directions in the same object.

Algorithm1 enumerates cliques of the path-compatibility graph, retaining
the signed union for each clique. Algorithm2 performs Shannon elimination
on edge orientations and memoizes the residual path constraints, INCLUDING
their multiplicities. They are distinct computations of the same polynomial.
Both use exact integers; reducing the resulting coefficients modulo any
positive modulus is safe. No division by a nonunit is introduced.

The parameter is the number and intersection structure of allowed paths,
not n alone. This can be dramatically smaller than2^N, but is exponential
in general. It is not a uniform efficient Hamilton-path algorithm.

## 3. The first square-sum overlap

Take G_n on {1,...,n}, joining distinct i,j when i+j is a square. The
[prior exact census](../stability_priority_moduli_20261102/AUDIT.md)
has1,1,1,3,10 paths up to reversal at15,16,17,23,25; all other n in2..25
have none. We enumerate the actual paths, retaining their orientations.
For n<=23 all the corresponding directed-path cubes are mutually disjoint.
At25 the20 directed paths have compatibility graph

    4 K2 disjoint union 12 K1.

Each compatible pair uses26 distinct arcs. The square-sum graph itself
has32 edges. Therefore (15), now counting orientations of G_25 alone, is

    Z_guard(z)=2^32+20*2^8(z-1)+4*2^6(z-1)^2
              =4294962432+4608 z+256 z^2.              (16)

Every orientation of the guard extends in2^268 ways to a tournament on25
vertices. Thus the full labelled tournament polynomial is2^268 times(16).
The two independent algorithms agree exactly; the clique compiler visits
25 nodes and the Shannon compiler598 states AFTER the path bank is known.
These are measured implementation counts, not minimal-state claims.

| n | Guard edges | Directed paths | Guard orientations with H_G=1 | With H_G=2 |
|---:|---:|---:|---:|---:|
|15|15|2|4|0|
|16|16|2|4|0|
|17|17|2|4|0|
|23|28|6|384|0|
|25|32|20|4608|256|

All other orientations have H_G=0. This makes25 the FIRST order in this
finite range with a positive even guarded response. Ordinary tournament
H remains odd: the guard has changed the counted path set. The occurrence
of overlap at the same endpoint as the cited eventual square-sum onset
is an exact local observation, not a proof explaining that all-size onset.

One compatible pair factors as UAB and UBA, where

    U=(18,7,9,16,20,5,11,25,24),
    A=(1,3,22,14,2,23,13),
    B=(12,4,21,15,10,6,19,17,8).

Both use only square sums. The four connecting edges have totals
24+1=25,24+12=36,13+12=25,8+1=9. Their directions impose no contradiction.

**General two-block swap criterion (PROVED).** For three disjoint directed
paths U,A,B, both UAB and UBA exist exactly when the end-to-start arcs
U->A, U->B, A->B and B->A all exist at their specified endpoint pairs.
This is immediate by checking the four joins. These are port relations:
collapsing A and B to single vertices of an ordinary tournament would
incorrectly demand opposite directions on the same pair. The blocks here
are not asserted to be tournament modules. This is a concrete benefit of
retaining joint ports instead of forcing every quotient into a tournament.

## 4. What survives the tree trijection and unlabelled counting

The [Session90 outer-endpoint theorem](../session90_even_graphs_odd_descent_20261029/trijection/TREE_COMPRESSION.md)
already gives an equivariant bijection

    tournaments extending an oriented spanning tree X
                    <-> tournaments on E(X).

It transports every non-tree arc through the first and last tree edges
on its unique path. LCA/path data compute this transport; dropping X loses
the decoding frame. Hamiltonian marking is the path-tree special case.
The graph cycle/component-syndrome boundary is independently recorded in
[the older trijection packet](../tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md).

A000055 counts unmarked tree shapes; A000568 counts unmarked tournament
types. The coinciding first four converse-merged counts stop at tournament
order5, where the merged tournament count is10 and the tree count at6 is6
([older exact census](../trees5_tournaments4_erdos_sos_20261005.md)). The actual
bridge is the marked completion map, not equality of those two sequences.

For a group Gamma preserving G, Burnside can be applied to (15) while
retaining the guard. If g has even order, it fixes no tournament. For odd
order g, close U(Q) under g. If two constraints in this closure prescribe
opposite orientations, the fixed intersection is empty. Otherwise it has

    2^(o_2(g)-c_g(Q))

members, where o_2(g) is the number of unordered-pair orbits and c_g(Q)
the number constrained by Q. Summing these terms with(z-1)^|Q| and then
averaging over Gamma gives the spectrum on Gamma-orbits of tournaments.
This follows by making one orientation choice per pair orbit; odd order
prevents internal pair reversal. It is an exact formula, not a claim of
an improved asymptotic algorithm.

For G=K_n,Gamma=S_n,z=1 it reduces to the existing A000568 Burnside formula.
For square-sum G, Gamma must preserve the arithmetic guard; its orbits are
different objects. Integer Burnside sums must be divided BEFORE reduction
modulo a non-coprime modulus, unless an integral replacement is supplied.

An analogous sum over invariant spanning trees recovers the Session90
marked count; forgetting marks has a nonconstant fibre. A nontrivial
odd-order vertex permutation with no fixed vertex fixes no finite tree,
because it must fix the tree's center, and cannot interchange the ends
of a central edge. This elementary pruning applies to tree-marked terms,
not to all tournament terms. It explains one genuine way tree structure
can help a marked enumeration without equating A000055 and A000568.

## 5. Algorithmic target left open

The promising next object is a shared-arc decision diagram for actual
attachment paths, with coefficients in the port-response ring and explicit
complement transport. The square-sum example compresses a bank because its
compatible subsets are exceptionally sparse. General strong atoms can
revisit children and change their complemented relation; those features
must remain in the state. The81/27 control must still separate before
claiming an unmarked recurrence or a smaller universally closed carrier.
