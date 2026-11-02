# Sharp order-nine stability through cycle-complement grouping

Status: PROVED, independently audited. This separate fourth-checkpoint
packet is outside earlier frozen publication manifests. No claim for
orders10..1499 follows from this note.

## Theorem

Every weighted tournament of order nine satisfies

    Delta >= (175/256)V+(105/256)tau+C4/16+(5/8)C5+C6/4+C7+C8.   (1)

Thus the optimal variance coefficient is175/256. Equality with V>0 in the
unreserved Delta>=175V/256 is exactly four fair ordinal pairs and one
singleton.

The notation is as in the reviewed lower-order notes: C_j counts directed
j-cycles once up to cyclic rotation, M_ab... counts unordered disjoint
cycle collections, and a2-cycle has weight q_e=Aij Aji. Put

    P=V-M22+M222-M2222.

## 1. Full order-nine cycle-cover ledger

The exact expansions are

    H=1+2tau+2C5+2C7+2C9+4M33+4M35+8M333,

    u=1-V+tau-C4+C5-C6+C7-C8+C9
       +M22-M23+M24-M25+M26-M27
       +M33-M34+M35-M36+M44-M45
       -M222+M223-M224+M225-M233+M234+M333
       +M2222-M2223,

    v=-C9+M27+M36+M45-M225-M234-M333+M2223.

Let

    X=P+tau+C4+C5+C6+C7+C8
      +M23-M24+M25-M26+3M33+M34+3M35-M44
      -M223+M224+M233.

Then

    G_-=X+2C9+8M333,
    G_+=X+2M27+2M36+2M45-2M225-2M234+6M333+2M2223.

Every partition of nine into parts>=2 appears in this ledger. In
particular the full covers3+3+3 and2+2+2+3 have their required multiplicities
and signs; no cancellation by a guessed sign of det(A) is used.

## 2. Group directed4-cycles with their five-vertex complements

For every directed4-cycle C, let R be its complementary five-set. At order5,

    u_R+v_R=1-V_R+tau_R-C4_R+M22_R,
    u_R-v_R=1-V_R+tau_R-C4_R+2C5_R+M22_R-2M23_R.

Therefore the G_- combination is

    C4-M24+M34-M44+M224
       =sum_C w(C)[u_R+v_R+C4_R/2],                 (2-)

and the G_+ combination is

    C4-M24+M34-M44+M224+2M45-2M234
       =sum_C w(C)[u_R-v_R+C4_R/2].                 (2+)

The factor1/2 is essential: summing complementary4-cycles counts M44 twice.
All other displayed collections contain a unique selected4-cycle and have
no such factor. The inherited uniform margin bounds on a five-vertex
remainder give u_R+v_R>=1/16 and u_R-v_R>=5/16. Thus both (2-) and (2+)
are at least C4/16.

## 3. Other negative even-cycle couplings

A directed6-cycle has a complementary triple with total variance<=3/4,
so C6-M26>=C6/4. A directed5-cycle has a complementary four-set with three
2-edge matchings, each at most1/16. Thus the G_+ term

    C5-2M225>=5C5/8.

The G_- term C5 is of course also at least5C5/8. The remaining C7,C8 and
all unused cycle products in the ledgers are nonnegative.

## 4. Triangles and a six-vertex variance estimate

For any weighted six-vertex tournament R, write V_R and W_R for variance
and incident-pair variance. Since

    M22_R=(V_R^2-sum_e q_e^2)/2-W_R <=V_R^2/2-W_R,
    tau_R>=W_R>=V_R^2/3-V_R/4,

we have

    M22_R-(3/2)tau_R
       <=V_R^2/2-(5/2)W_R
       <=-V_R^2/3+5V_R/8
       =75/256-(V_R-15/16)^2/3
       <=75/256.                                  (3)

For each directed triangle C in the nine-vertex matrix, its complement
has six vertices. Summing (3) against w(C) gives

    tau+3M33-M223
       =sum_C w(C)[1+(3/2)tau_R-M22_R]
       >=(181/256)tau.                             (4)

Here every unordered pair of directed triangles is counted twice, giving
the coefficient3 from the local coefficient3/2. The M223 collection has a
unique triangle, so its coefficient is one.

Combining (2)--(4) with the complete signed ledgers proves

    Delta>=P+(181/256)tau+M23+C4/16+(5/8)C5+C6/4+C7+C8.        (5)

This reduction is analytic and uses only already-proved lower-order
identities and the local triangle-variance lemma.

## 5. Multiaffine graph reduction

The triangle-variance lemma yields tau>=W(q) and M23>=P21(q), where P21
is the weighted count of a two-edge wedge and a vertex-disjoint edge.
Consider

    F(q)=P+(19/64)W(q)+P21(q)-(175/256)V.

Every monomial uses distinct edges, so F is affine in each q_e on the box
0<=q_e<=1/4. At an endpoint undirected graph G on nine vertices, let
m,m2,m3,m4 count edges and unordered matchings of the indicated size,
w count two-edge wedges, and p21 count a wedge plus disjoint edge. Then

    1024F=81m-64m2+16m3-4m4+19w+16p21.              (6)

`graph9_certificate.cpp` verifies nonnegativity on every nine-vertex graph
up to relabeling through a self-contained augmentation certificate:

1. Partition all2^21 labeled seven-vertex graphs into permutation orbits,
   explicitly verifying total coverage2,097,152 and obtaining1044orbits.
2. Attach an eighth vertex with all128 neighbor subsets to every orbit
   representative, producing133,632augments that cover every eight-vertex
   graph up to relabeling.
3. Canonicalize each eight-vertex augment by sorting vertex degree classes
   and enumerating EVERY permutation within each equal-degree class.
   This is an exact isomorphism-invariant canonicalization, not a hash or
   incomplete color-refinement surrogate. It yields12346representatives.
4. Attach a ninth vertex with all256neighbor subsets to each representative.
   These3,160,576augments cover every nine-vertex graph up to relabeling.
5. Count1260three-edge matching patterns,945four-edge matching patterns,
   wedges, and p21 by exact pattern masks and complementary degrees.

The minimum of (6) is zero. The14zero augmentations may include duplicate
isomorphism types and are not presented as a14-orbit classification.
No graph-atlas file, external data, or random sample is used. All arrays
are generated in memory; the receipt is a small JSON object.

Because F is separately affine and (6) is relabeling-invariant, this finite
certificate proves F>=0 over the full real q-box. Substituting into (5)
gives

    P+(181/256)tau+M23
       >=(175/256)V+(181/256-19/64)tau
       =(175/256)V+(105/256)tau,

which is (1).

The coefficient19/64 is optimal for this graph reduction: P3 disjoint-union
three edges has m=5,m2=9,m3=7,m4=2,w=1,p21=3, and the residual without the
w term is-19. This does NOT claim that105/256 is the optimal actual
triangle reserve coefficient for weighted tournaments.

## 6. Sharpness and equality

Four fair ordinal pairs and one singleton give

    H=1, u=81/256, v=0, V=1, Delta=175/256.

Conversely, equality with V>0 in the unreserved variance inequality forces
tau=0 by (1). All supported binary completions are therefore transitive,
so H=1. The reviewed H=1 classification and scalar minimization yield exactly
four fair pairs and a singleton.

## Reproduction

    g++ -O3 -std=c++17 graph9_certificate.cpp -o /tmp/graph9
    /tmp/graph9 > /tmp/graph9.replay.json
    cmp graph9_certificate.json /tmp/graph9.replay.json

The script also records the weaker coefficient alpha=1/2 and the exact
maximum required alpha=19/64. Those corroborative outputs are not additional
assumptions. The independent audit reconstructed both the grouping and the graph
certificate using a distinct BFS/equitable-canonicalization implementation,
vertex-addition pattern counts, and full-DP cross-check snapshots.
