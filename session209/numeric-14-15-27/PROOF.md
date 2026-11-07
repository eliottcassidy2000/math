# Numeric 14, 15 and 27: counts, quotient degrees and an efficient atom

Status: PROVED consequences of the cited strong-tournament bound and inherited complete small-class atlas; independently recomputed finite inventories and VERIFIED order14 witness. This packet explores a numeric interpretation of “14, 15 and 27”; it does not assert that this was the intended reference rather than prior session numbers. No priority claim is made.

## 1. Three different numerical roles

Let H(T) count directed Hamilton vertex words. Let N_n(H) count ordinary tournament isomorphism classes and S_n(H) the self-converse ones. Converse-merged classes number (N_n(H)+S_n(H))/2. Let m(H) be the least realizing order.

The following statements concern different objects:

- At order n=14 there are exactly110 ordinary classes with H=15 and210 with H=27. None is self-converse, so the merged counts are55 and105.
- m(15)=5 and m(27)=6. A raw Hamiltonian-path count14 is impossible by Rédei parity.
- An explicit strong tournament of order14 has H=15^6=11390625. Its repeated dominance joins imply the construction rate upper bound7/3 for powers of15. This is not a claim that either m(15^6)=14 or that7/3 is optimal.

The first statements are exact specializations of earlier SCC multiplicity formulas. The further contribution here is an exact degree-resolved comparison of their order14 strata, with a standalone replay and explicit separation from native metagraph degrees.

## 2. Complete small-atom data and why it is global

For a nontrivial strong tournament C, write a(C)=|Aut(C)|. The relevant complete inventory is the following multiset of (order,H,automorphism order,self-converse):

- H3: (3,3,3,yes).
- H5: (4,5,1,yes).
- H9: two (5,9,1,yes).
- H15: (5,15,3,yes), (5,15,5,yes), two (6,15,1,no).
- H27: two (6,27,3,no), two (7,27,1,no), one (7,27,1,yes).

The two non-self-converse classes in each equal-sized pair are converses. RESULT.json supplies a checked adjacency code for every atom. Its encoding is bit k on the kth lexicographic pair i<j, equal to1 exactly when i beats j.

Completeness of the456 ordinary order7 representatives is inherited from the audited Session24 census, not claimed to be a new exhaustive enumeration of all labelled order7 tournaments. Taking every induced subset recovers all ordinary classes of each smaller order: any smaller tournament extends to order7 by adding sources. The script canonicalizes every such subset by all degree-preserving permutations, recovering exactly2,4,12,56,456 classes at orders3 through7. It independently recalculates H, strong connectivity, automorphism orders and self-converseness.

The cited bound of Arthur H. Busch, *A Note on the Number of Hamiltonian Paths in Strong Tournaments*, Electronic Journal of Combinatorics13 (2006), N3, https://doi.org/10.37236/1141, is

    H(C) >= 5^((|C|-1)/3).

Since5^7>27^3, no strong atom with H<=27 has order>=8. Thus the small inventory is global, not a fixed-order extrapolation. Its use and primary-source provenance are inherited from Session24. The current packet does not newly certify Busch's proof.

## 3. Exact ordinary and self-converse counts

Every tournament has a unique ordered strong-component word. H and automorphism orders multiply across components. For a core word with r nontrivial components and total order s, insert n-s singleton components in its r+1 gaps. This gives

    binom(n-s+r,r)

ordinary classes when n>=s, and zero otherwise. Different ordered SCC words cannot become isomorphic.

For H15 the possible core types (r,s,multiplicity) are

    (1,5,2), (1,6,2), (2,7,2).

For H27 they are

    (1,6,2), (1,7,3), (2,8,4), (3,9,1).

These lists exhaust multiplicative factorizations15=15=3*5 and27=27=3*9=3*3*3. The product9 has its two strong order5 atoms. Hence

    N_n(15) = 2(n-4)+2(n-5)+2 binom(n-5,2)
            = (n-4)(n-3),                                  n>=5,

    N_n(27) = 2(n-5)+3(n-6)+4 binom(n-6,2)+binom(n-6,3)
            = n(n-5)(n-4)/6,                              n>=6.

The displayed polynomial expressions include vanishing terms at their low boundary; the literal core contribution is always zero below its size. They give N_14(15)=110 and N_14(27)=210.

Converse reverses the SCC word and each component. A self-converse tournament therefore has form A -> Z -> A^op with empty or strong self-converse center Z, and H=H(A)^2 H(Z). For15 only H(A)=1 is possible; the center is one of the two odd order5 atoms. Thus S_n(15)=2 for odd n>=5, and0 for even n. For27 either H(A)=1 with the unique self-converse order7 H27 center, or H(A)=3 with center C3. The former contributes1 at each odd n>=7; the latter contributes the number of H3 ordinary classes of order (n-3)/2, namely (n-7)/2 when n>=9 is odd. Therefore

    S_n(27)=(n-5)/2 for odd n>=7, and0 for even n.

In particular both order14 strata consist entirely of non-self-converse pairs. The two minimum orders follow from the explicit strong order5 H15 and order6 H27 witnesses and the complete lower-order spectra (maximum H at order4 is5, and at order5 is15). Source/sink extension makes these spectra nested.

## 4. Exact order14 paired-black degrees

Use the fundamental reflection-paired black converse-merged multigraph Q of THM-4731, Statement3. It is obtained by identifying sigma-twin black lines; loops contribute two to degree. At a non-self-converse pair its degree is

    deg_Q = H(T)/|Aut(T)|.

At a self-converse class the separate formula is (H/|Aut|-B)/2, with B its blue path-orbit count. No order14 class considered here is self-converse, so the first formula always applies. Degree is taken in the full Q_14 graph, not in the subgraph induced by a fixed-H stratum.

For H15:

- The order5 atom with automorphism order5 contributes10 ordinary padded classes, hence5 merged nodes of degree3.
- The other order5 atom contributes10 ordinary classes with automorphism order3. The two ordered C3/S4 cores contribute2 binom(9,2)=72 ordinary classes, also with automorphism order3. Together these give41 merged nodes of degree5.
- The two asymmetric order6 atoms contribute18 ordinary classes, hence9 merged nodes of degree15.

Thus the exact H15 histogram at order14 is

    degree3:5 nodes; degree5:41 nodes; degree15:9 nodes.

For H27:

- The three-C3 core gives binom(8,3)=56 ordinary classes, each with automorphism order27. Hence28 merged nodes have degree1.
- The two order6 atoms give18 ordinary classes with automorphism order3. The four ordered C3/H9 core types give4 binom(8,2)=112 more with that automorphism order. Hence65 merged nodes have degree9.
- The three asymmetric order7 atoms give24 ordinary classes. Hence12 merged nodes have degree27.

Thus the exact H27 histogram at order14 is

    degree1:28 nodes; degree9:65 nodes; degree27:12 nodes.

The28 degree-one nodes are genuine leaves: a loop contributes two, so cannot account for degree1. Their fraction is28/105=4/15, NOT a majority. There are no H15 leaves at order14.

More generally, for even n>=10 the respective histograms are

    H15: degree3 -> (n-4)/2,
         degree5 -> (n-4)/2 + binom(n-5,2),
         degree15 -> n-5;

    H27: degree1 -> binom(n-6,3)/2,
         degree9 -> n-5+2 binom(n-6,2),
         degree27 -> 3(n-6)/2.

This gives a precise finite precursor to the inherited typical-degree result: at fixed H15 and growing n, almost all merged nodes have degree5; at fixed H27, almost all have degree1. The latter conclusion holds also along odd n after removing the lower-order palindromic core words. It is an asymptotic statement, not an order14 majority claim.

Session14's rigid-word leaf classification identifies the degree-one mechanism: nonpalindromic words of singleton and C3 components. Session15 closes the other leaf mechanism as A -> S4 -> A^op with rigid A. Their raw counts are respectively3^a and5*9^a. Hence H15 has no Q-leaves at ANY order, while all H27 leaves arise from the nonpalindromic three-triangle words. These stronger classification conclusions are inherited, not re-proved by this finite packet.

## 5. The latest order14 atom is a different use of14

Session208 supplied the strong adjacency code

    292187853813044375413825338

on14 vertices, with H=11390625=15^6. The standalone checker independently reconstructs its arcs, verifies all-source reachability and recounts paths by endpoint subset DP. It does not invoke the original search/oracle.

Dominance joins multiply H and add orders. For k=6q+s, 0<=s<6, use q copies of this atom and s copies of an order5 H15 atom:

    m(15^k) <= 14q+5s.

Consequently the existing subadditivity limit satisfies lim m(15^k)/k<=7/3. The witness is not a minimum-order proof. Larger strong atoms might improve it; a positive linear lower bound remains open. This atom does not belong to either H15 or H27 stratum above: its H is15^6.

The mixed numeric construction

    m(15^(6q) * 27^r) <= 14q+6r

also follows immediately from the checked H27 order6 atom. It is merely a valid dominance-join bound, not an optimal mixed-prime solution.

## 6. Boundaries and reproduction

The labelled native metagraph retains a Hamilton path and reverses all other arcs. Its labelled degree is raw H, so degrees15 and27 there must NOT be replaced by the histograms above. The single-arc search graph is a third distinct graph. Degree histograms prove neither neighbors, connectivity, nonleaf bridges nor a global spectral threshold.

Run `python3 verify.py > RESULT.json` and `python3 -O verify.py > OPTIMIZED.json`; the outputs must agree. The checker uses only the standard library and bundled atlas input. Its326 active checks survive optimized Python. The complete source atlas is an inherited input, and the Busch inequality is cited; no claim of a fresh all-labelled order14 census is made. Source paths, hashes and the exact current publication commit are in DEPENDENCIES.json.
