# Sharp order-seven weighted stability from a finite graph certificate

Status: PROVED, independently audited. The reduction is analytic; its last
graph inequality has a complete finite exact certificate and an independent
full-universe edge-toggle replay.

## Statement

For every weighted tournament of order seven,

    Delta >= (37/48)V + (7/24)tau + C4/4 + C5 + C6 + 3M33.       (1)

In particular the optimal unreserved variance coefficient is 37/48.
Equality with V>0 in Delta>=37V/48 occurs exactly for ordinal sums of three
fair 2-blocks and one singleton.

Definitions H,u,v,Delta,V,C_j,M_ab and q_e are as in
SHARP_LOW_ORDER_STABILITY.md. M_223 is the sum over a directed triangle
and two disjoint complementary 2-cycles, counted once per cycle collection.

## Exact cycle expansions

At order seven the determinant and odd-cycle formulas are

    H=1+2tau+2C5+2C7+4M33,
    u=1-V+tau-C4+C5-C6+C7+M22-M23+M24-M25
         +M33-M34-M222+M223,
    v=-C7+M25+M34-M223.

Therefore, with P=V-M22+M222,

    G_-=P+tau+C4+C5+C6+2C7+M23-M24+3M33,
    G_+=P+tau+C4+C5+C6+M23-M24+2M25+3M33+2M34-2M223.

For each directed 4-cycle, its complementary triple has total edge
variance at most3/4, so M24<=3C4/4. For each directed triangle, its
complementary four-set has three 2-edge matchings, each of weight at most
1/16. Hence M223<=3tau/16. All other displayed cycle weights are
nonnegative. Both signed margins consequently satisfy

    G_+,G_- >= P+(5/8)tau+M23+C4/4+C5+C6+3M33.       (2)

No even-cycle cancellation or sign of det(A) has been assumed.

## A small local variance bound and its disjoint lift

TRIANGLE_VARIANCE.md proves, over the entire weighted triangle cube,

    tau_R >= sum_{unordered {e,f} of distinct edges in R} q_e q_f,

where unordered pairs of its edges are counted once. Summing gives

    tau >= W(q):=sum_{incident unordered edge pairs e,f} q_e q_f.

Multiplying the local bound by the variance of a disjoint edge and summing
also gives

    M23 >= P21(q),

where P21(q) sums q_e q_f q_g over every two-edge wedge plus a vertex-disjoint
edge. The isolated edge is unique in such a three-edge pattern; no factor
of two or three is inserted.

## Finite graph inequality, continuously extended

The polynomial

    F(q)=P+(1/3)W(q)+P21(q)-(37/48)V                 (3)

is affine in every individual q_e, since every matching or wedge pattern
uses distinct edges. Thus nonnegativity on the 21-dimensional box
0<=q_e<=1/4 is equivalent to nonnegativity at all its endpoint vertices.
An endpoint vertex is an arbitrary simple undirected graph G on seven
labeled vertices, with q_e=1/4 on its edges and zero elsewhere.

Let m be its edge count, m2 its number of unordered 2-edge matchings,
m3 its number of unordered 3-edge matchings, w its number of two-edge
wedges, and p21 its number of wedge-plus-disjoint-edge patterns. Direct
scaling gives

    384 F = 22m-24m2+6m3+8w+6p21.                    (4)

`sharp7_graph_certificate.cpp` checks every one of the 2^21=2,097,152
labeled simple graphs with no filters or isomorphism assumptions. It counts
matchings by increasing edge indices and disjointness masks. Wedges are
w=binom(m,2)-m2. To count p21, for each candidate isolated edge it computes
sum_v binom(deg(v),2) in the graph on its complementary five vertices.
All arithmetic is integral and far below int32 limits. The minimum in (4)
is zero; the receipt records the exact zero count and a witness.

Because (3) is separately affine, the finite graph certificate proves (3)
for every real vector q in the closed box, not just graph endpoints.
Combining (2),(3) with tau>=W and M23>=P21 gives

    P+(5/8)tau+M23
      >= (37/48)V + (5/8-1/3)tau
      = (37/48)V+(7/24)tau,

which is (1).

## Sharpness and equality

Three fair pairs and a singleton give H=1,u=27/64,v=0,V=3/4, so
Delta=37/64=(37/48)V. If equality holds in the unreserved inequality with
V>0, (1)'s positive triangle coefficient forces tau=0. Every supported
binary completion then has no directed triangle and is transitive, hence
H=1. The all-order H=1 classification and scalar minimization in
H1_STABILITY.md force exactly three fair pairs and one singleton.

## Reproduction and controls

    g++ -O3 -std=c++17 sharp7_graph_certificate.cpp -o sharp7_graph_certificate
    ./sharp7_graph_certificate > sharp7_graph_certificate.json

The weaker graph reduction without p21 is false: C4 disjoint-union K3
has residual -41 in the earlier 15w formula. This hostile is documented
in TRIANGLE_VARIANCE.md. The p21 coordinate is necessary to this proof
route; dropping it is not a harmless simplification.

At fixed graph target37/48, the minimum coefficient alpha for
P+alpha W+P21>=37V/48 is1/3. The graph P3 disjoint-union two edges witnesses
necessity, and the same full enumeration verifies sufficiency. This does
not assert optimality of the triangle reserve coefficient7/24 for actual
weighted tournaments, since tau>=W loses information.
