# Collatz inverse growth cannot fit a uniform planar lattice

**PROVED elementary geometric obstruction; FINITE-EXACT small certificate.**
The inverse binary tree is recovered from SF5 in
[Collatz source floors that survive refinement](../collatz_source_refinement_floor_20261005.md).
Concurrent packets already prove the same polynomial-growth obstruction:
[clocks/carriers, section3](../clocks_carriers_genus_20261102/COLLATZ_GEOMETRY.md)
(66cf347030, root85, macro-clock<=7) and
[joint incidence](../joint_incidence_clocked_frames_20261102/COLLATZ_LATTICE.md)
(252a039419, the same SF5 root5). This session credits those proofs and adds
the radius21/26 certificates, exact lattice capacities, finite-fibre bounds,
and [native joint-port reconstruction](JOINT_PORTS.md).
None proves or disproves Collatz convergence.

## 1. Specify the clock and the meaning of a hexagonal grid

Use the ordinary positive-integer map

    C(n)=n/2 if n is even; 3n+1 if n is odd.

Every vertex has exactly one outgoing edge. The incoming vertices at m
are2m and, exactly when m=4 mod6, (m-1)/3. Thus the incoming degree is
one or two, as suggested. The shortcut map T replaces the odd step by
(3n+1)/2; the odd-only map U(n)=oddpart(3n+1) has infinitely many
incoming vertices. These are different clocks and different graphs.

The vertices of a hexagonal tiling form the degree-three honeycomb graph
H; its cell centres form the degree-six triangular graph L. They have
different adjacency and colouring. In particular H is bipartite whereas
L has triangles. The native cycle1->4->2->1 already rules out an injective
nearest-neighbour map from the entire C graph to H.

The following obstruction remains after removing all cycles: the inverse
tree rooted at5 alone cannot map injectively to either lattice with every
native edge joining lattice neighbours. Nor can a uniform bound on edge
length, or any fixed finite number of sheets, repair the infinite model.

## 2. Exact lattice ball capacities

Use H vertices(a,b,s), a,b integers, s in{0,1}; a black(a,b,0) is adjacent
to white(a,b,1),(a-1,b,1),(a,b-1,1). The triangular lattice on Z^2 has
steps±(1,0),±(0,1),±(1,-1). With graph metric,

    |B_H(r)|=1+3r(r+1)/2,
    |B_L(r)|=1+3r(r+1).                               (1)

For L the radius-k ball is |a|,|b|,|a+b|<=k; its shells have6k points.
In H, black vertices at distance<=2k form that triangular ball. White
vertices at distance<=2k+1 form its union with shifts(-1,0),(0,-1).
This union has

    -k-1<=a,b<=k,    -k-1<=a+b<=k,

and3(k+1)^2 points. One verifies the union description by shifting any
point on a new lower boundary back into the triangular ball. To count,
put K=k+1 and shift both coordinates by K: a 2K by2K box remains, with
sum between K and3K-1. The two removed triangular corners have altogether
K^2 points. The even-radius white count is3k^2. Adding the two colours
gives(1) at both parities. This is a proof, not extrapolation from BFS.

## 3. A finite obstruction inside the known convergent basin

Let I_r be the set of positive integers whose first C-hit of5 occurs in
at most r steps, including5. Generate successive levels using the two
exact predecessor formulas above. Distinct levels are disjoint:5 is not
periodic, as its forward orbit is5,16,8,4,2,1,4,... .

The complete exact computation gives

    |I_21|=704 >694=|B_H(21)|,
    |I_26|=2263>2107=|B_L(26)|.                       (2)

Every source and its first-hit clock are retained in
inverse_root5_depth26.tsv. Every edge and every full path back to5 are
replayed by the checker; lattice capacities are independently generated
by BFS as well as(1). These are the first radius-count violations for
this root in the respective tested prefixes. No minimal obstruction over
all roots or embeddings is claimed.

An edge-preserving injective map would send I_r inside the radius-r ball
around the image of5, contradicting(2). This finite subgraph is a tree;
the argument is independent of odd-cycle parity or unknown trajectories.
It applies to arbitrary placements, not just a particular spiral drawing.

## 4. All-size distortion and multiplicity bounds

Recover the SF5 tree rooted at5. At an odd parent n>1 not divisible by3,
use candidate exponents{2,4,6} when n=1 mod3, and{3,5,7} when n=2 mod3.
The associated predecessors are

    j_a(n)=(2^a n-1)/3.

Exactly one of the three candidates is divisible by3: modulo9, the
eligible powers of2 cycle through the three lifts of the required
nonzero residue modulo3. Keep the other two. They are odd, positive,
not divisible by3, and strictly greater than n. Their U-image is n,
with exact valuation a. Unique forward parents and strict increase give
a full binary tree without collisions or repeated depths.

Each inverse macro-edge costs a+1<=8 ordinary C steps. Consequently
the2^(h+1)-1 distinct odd nodes through depth h all lie within native
distance8h of5. This recovers an explicit exponential-growth subfamily
without invoking the Collatz conjecture. At h13 it already gives

    16383>16381=|B_H(104)|.

Let f map the native graph to H, with every native edge's images at
graph distance<=L, and at most K source vertices over any image vertex.
For every h,

    2^(h+1)-1 <= K[1+3(8hL)(8hL+1)/2].                (3)

For integer graph metrics take L an integer upper bound. Fixed finite
K,L contradict(3) as h grows. More generally a target with ball size
at most C(1+r)^d requires

    2^(h+1)-1 <= KC(1+8hL)^d.                        (4)

Thus an injective fixed-dimensional lattice drawing must allow edge
stretch growing exponentially in h up to the polynomial h denominator.
Keeping bounded stretch requires unbounded, quantitatively growing fibres.
Calling those fibres sheets does not remove the necessary growth.

For a planar Euclidean drawing whose distinct image positions have mutual
separation at least delta>0, disjoint radius-delta/2 disks give instead

    2^(h+1)-1 <= K(1+16hL/delta)^2.                   (5)

Here L bounds each native edge's Euclidean length. This does NOT exclude
arbitrary unit-distance drawings with no uniform separation, arbitrary
curved planar edges, or subdivisions with unbounded route length.

## 5. A faithful geometric repair on the certified subtree

Retain the binary address w of each SF5 node and its native edge clocks.
At depth h, let j be the binary value of w and draw the node at

    x=(j+1/2)/2^h,       y=1-2^(-h).                  (6)

This is injective. Each child is drawn in its parent's horizontal interval;
distinct branches cannot cross. A macro-edge from depth h to h+1 has
length sqrt(5)*2^(-h-2). Thus the picture accommodates exponential growth
by shrinking its scale, exactly as(5) requires. The limiting boundary is
not a set of additional positive-integer nodes. Independent clock-marked
subdivisions duplicate shared EVEN tails: the native graph requires their
identification, as the concurrent root85 packet already warns. Injectivity
or planarity of this macro drawing does not automatically survive that
gluing. JOINT_PORTS.md gives its exact identifications and count correction.

The representation retains source labels, branch addresses, exponents,
and the root stop. Erasing those labels leaves an abstract binary tree
shared by many unrelated maps. LCA in this tree is the longest common
address prefix; native first-hit times are sums of a+1, not address length.
The first native meeting of two branches can precede their odd LCA; its
exact correction uses the two final exponent labels (JOINT_PORTS,J4).

There is no single forward Collatz trajectory visiting every positive
integer: after it reaches1 it is periodic, and branching predecessors
cannot both immediately precede the same vertex in a nonrepeating path.
The construction is therefore a directed network, not a single
space-filling successor curve. This distinguishes a faithful network
picture from a proposed convergence argument.

## 6. The five-times-plus-one hostile and a prime-clock generalization

Let p be an odd prime for which2 generates(Z/p^2 Z)^*. Put
A=ceil(log_2(p+2)). For any odd n coprime to p, consider

    a in[A,A+p(p-1)-1],    j_a(n)=(2^a n-1)/p.

Exactly p exponents make j_a integral, and exactly one of those makes it
divisible by p. The other p-1 give strictly larger positive odd unit
predecessors of n for U_p=oddpart(pn+1). Hence they generate a collision-free
(p-1)-ary inverse tree with native macro-clock at most A+p(p-1).
The proof uses the orders p-1 and p(p-1) of2 modulo p and p^2, followed
by unique forward parents. It is conditional on that order hypothesis;
no assertion that it holds for every prime is made.

For p5 the hypothesis is elementary and yields four branches with a in3..22.
But U_5 has the positive cycle13->33->83->13. Thus abundant inverse trees,
planar drawings at shrinking scale, and the degree restrictions do not
distinguish convergence to the principal cycle. The missing assertion for
Collatz is universal forward entry into its known basin.

Applegate--Lagarias give substantially stronger quantitative inverse-tree
growth estimates in [Density Bounds for the3x+1 Problem I, Theorem1.1](https://websites.umich.edu/~lagarias/doc/applegateI.pdf).
Their growth theory is classical. We use the weaker repository SF5 tree
because its short proof supplies a self-contained geometric certificate.
