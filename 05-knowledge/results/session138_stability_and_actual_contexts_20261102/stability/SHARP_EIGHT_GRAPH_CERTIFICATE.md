# Sharp order-eight stability from a finite graph certificate

Status: PROVED, independently audited. This is a separate second-checkpoint
packet, not part of the frozen first Session138 publication manifest.

For every weighted tournament of order eight, the theorem is

    Delta >= (87/128)V + tau/16 + C4/8 + C5 + C6/2 + C7 + 3M33.  (1)

Notation is as in SHARP_LOW_ORDER_STABILITY.md; in particular M_ab...
counts unordered disjoint directed-cycle collections, and 2-cycles have
weight q_e=Aij Aji. The proof reduction and zero-wedge graph certificate have both been
independently audited. A distinct BFS-orbit/subset-DP implementation verifies
the entire eight-vertex scope, with matching normal and UBSan receipts.

## Complete cycle-cover ledger

Let P=V-M22+M222-2M2222. At order eight,

    H=1+2tau+2C5+2C7+4M33+4M35,

    u=1-V+tau-C4+C5-C6+C7-C8
       +M22-M23+M24-M25+M26+M33-M34+M35+M44
       -M222+M223-M224-M233+M2222,

    v=-C8+M26+M35+M44-M224-M233+M2222.

For G_+=H-u+v all full-cover determinant terms cancel, leaving their H
coefficients. For G_-=H-u-v they double with negative sign. Thus

    G_+=P+2M2222+tau+C4+C5+C6+C7
         +M23-M24+M25+3M33+M34+4M35-M223,

    G_-=P+tau+C4+C5+C6+C7+2C8
         +M23-M24+M25-2M26+3M33+M34+2M35-2M44
         -M223+2M224+2M233.

This ledger explicitly includes all partitions of eight with parts>=2,
including the full covers4+4,2+6,3+5,2+2+4,2+3+3,2+2+2+2.

## Grouping the potentially negative cycle terms

For a directed4-cycle C, let R be its complementary4-set. Then

    C4-M24+M34-2M44
       =sum_C w(C)[1-V_R+tau_R-C4_R]
       =sum_C w(C)[u_R-M22_R].                      (2)

The factor2 on M44 is essential: each unordered pair of disjoint4-cycles
is counted twice on the right. The inherited universal determinant bound
gives u_R>=5/16 at order4, whereas M22_R<=3/16 because there are three
2-edge matchings, each at most1/16. Hence (2)>=C4/8.
For G_+ the -2M44 term is absent, so inserting it only weakens a lower bound.

Every6-cycle has a complementarypair with variance<=1/4, so

    C6-2M26>=C6/2.

Every directedtriangle has a complementary5-set with15two-edge matchings,
each of weight<=1/16. Consequently

    M223<=15tau/16, and tau-M223>=tau/16.

All other omitted terms in the two signed ledgers are nonnegative.
Therefore both margins, and hence their minimum, satisfy

    Delta >= P+tau/16+M23+C4/8+C5+C6/2+C7+3M33.       (3)

## The graph certificate

The audited triangle-variance lemma implies M23>=P21(q), where P21 is the
weighted count of a two-edge wedge and a vertex-disjoint edge. Consider

    F(q)=V-M22+M222-2M2222+P21(q)-(87/128)V.

It is affine in each q_e in [0,1/4]. At a box endpoint, let G be the
undirected graph of q_e=1/4 edges. Write m,m2,m3,m4 for its edge/matching
counts, and p21 for wedge-plus-disjoint-edge count. Then

    512 F=41m-32m2+8m3-4m4+8p21.                    (4)

`graph8_zero_alpha.cpp` certifies nonnegativity of (4) on all graphs up to
relabeling without relying on a graph-atlas data file:

1. Enumerate all2^21 labeled seven-vertex graphs. Whenever an unmarked
   graph is found, mark all5040vertex-permutation images. This produces
   1044representatives and explicitly verifies coverage2,097,152.
2. For each representative, attach vertex8 with each of its128possible
   neighbor subsets. Every eight-vertex graph is isomorphic to one of
   these133,632augmentations, because deletion of any fixed vertex gives
   one of the covered seven-vertex graphs.
3. Count all420three-edge matchings,105perfectmatchings, and p21 directly
   from precomputed pattern masks and complementary degrees. The minimum
   of (4) is zero, with two zero augmentations (empty and perfectmatching).

The graph polynomial is relabeling-invariant. Multiaffine interpolation
therefore extends the finite graph certificate to every real q-box point.
Combining F>=0 with (3) proves (1).

## Sharpness and equality

Four fair ordinal pairs have H=1,u=(3/4)^4=81/256,|v|=1/256,V=1.
Thus Delta=174/256=87/128, proving sharpness. Equality with V>0 in the
unreserved bound Delta>=87V/128 forces tau=0 by (1), hence H=1. The H=1
classification then forces four fair pairs.

## Reproduction

    g++ -O3 -std=c++17 graph8_zero_alpha.cpp -o /tmp/graph8
    /tmp/graph8 > /tmp/graph8.replay.json
    cmp graph8_zero_alpha.json /tmp/graph8.replay.json

An earlier weaker graph certificate (`graph8_certificate.cpp`) retained
an additional positive2w term; it is not needed for the proved result.
No claim for ordernine or arbitraryorder follows from this finite scope.
