# The Fano quotient and its three-port lift

This is about the exact 21-vertex Cayley tournament U used in Session319's arc-shadow collision. Its underlying group consists of pairs (a,b) with a modulo 7, b modulo 3, multiplication (a,b)(c,d)=(a+2^b c,b+d). The generating set has three slices:
```
S_0={1,2,3},   S_1={0,1,2,4},   S_2={1,2,4}.
```
The new ingredient is a canonical triangular partition, its Fano quotient, and a precise incidence observer for U in the merged path-tiling metagraph. The use of {1,2,4} as the Paley/Fano difference set is inherited (e.g. collatz_mod6_20260922_paley_fano_octonion_design.md and THM-2523); what follows concerns the specific Cayley lift.

## Canonical triangles and Fano lines

Use the number of common out-neighbors of each unordered pair, an orientation-dependent pair weight. At (0,0), grouping other vertices by their second coordinate, these values are:
```
b=0, a=1..6: 5,4,3,3,4,5
b=1, a=0..6: 7,5,5,3,4,5,4
b=2, a=0..6: 7,5,3,5,4,4,5
```
These are direct intersections using the three S slices, replayed from literal adjacency. By left transitivity the weight-7 pairs are exactly the pairs inside each triangle F_a={(a,0),(a,1),(a,2)}. Each triangle is directed in increasing b modulo 3. This partition is canonical under automorphisms *and anti-isomorphisms*: U is regular of outdegree 10, hence every pair has equally many common in- and out-neighbors (subtract the internal arc and cross counts from the two degrees, or use complementation of the external in/out sets). Thus the weight remains the same under global reversal. This pair weight is not the Hamilton arc-use weight from Session319.

Put D={1,2,4}. For d=c−a in D, between F_a and F_c the arcs with target b'=b+1 or b+2 always point from a to c, since D is invariant under multiplying by 2. Of the three same-b arcs two point from a to c. Their exceptional b is determined by d=2^(b+2) mod 7: the same-b arc allows offsets {2^b,2^(b+1),3·2^b}, exactly two of D. So there are eight forward arcs and exactly one dissent, joining the same port b at both ends. For d∈−D the directions are reversed by interchanging a,c. The resulting majority quotient is the Paley tournament on 7 points. In each triangle its three outgoing quotient arcs select its three dissent ports bijectively.

For the field interpretation, choose alpha with alpha^3+alpha+1=0 in F8. The three nonzero field elements alpha, alpha², alpha⁴ sum to zero, and alpha has order 7. Thus the three-point translates a+D on exponents are precisely the Fano lines on nonzero field vectors (seven distinct sets, each a line). The affine exponent action a→c+2^d a is multiplication by alpha^c followed by a Frobenius power. The dissent formula literally encodes the three Frobenius powers as ports; this is a finite field construction, not a structural implication of octonion/Bott periodicity.

## Automorphisms and why no converse lifts

An automorphism induces an automorphism of the quotient, and one fixing the quotient fixes all dissent ports by the outgoing bijection. Thus the induced map is faithful. The quotient automorphism group has size at most 21: after fixing 0, its out-neighbors {1,2,4} form a directed triangle and permit at most three rotations; each of the three in-neighbors beats a distinct single out-neighbor, so these rotations determine the action on them too. Left translations of the Cayley group already induce all 21 affine maps a→c+2^d a, hence Aut(U) has order exactly 21.

Every quotient anti-map is an automorphism composed with inversion a→−a. If U had an anti-map, compose it with the available lift of that automorphism's inverse and assume it induces inversion. A quotient arc from a to c has offset d∈D, and its reversed image from −c to −a has the same offset d, hence the same dissent port b. Therefore the anti-map must take (a,b) to (−a,b) for all a,b. This preserves the cyclic orientation *within* each triangle rather than reversing it, a contradiction. The quotient is self-converse; the three-port lift is not. An exact new obstruction to lifting a quotient converse is thus a conflict between forced port motion and the internal orientation.

## Marked merged-metagraph incidence

In the usual line convention kappa reverses all arcs except those of the distinguished Hamilton path. The path-tiling presentations of U identify with Aut(U)-orbits of Hamilton paths, and left regularity gives one representative per orbit beginning at (0,0). Since U is not self-converse, none is blue; their sigma twins give exactly the same edges in the merged graph and are counted once in the fundamental black graph. Thus the degree of its merged node in that graph is H(U)/21. Every neighbor is near-regular: on a path p0→...→p20, the kappa partner's degrees are 11 at p0, 9 at p20, and 10 elsewhere. So no loop is possible at this node, and every endpoint in the partner is intrinsically marked by its degree (converse interchanges the 9 and 11 labels). This degree-boundary argument holds for every regular source, not just U.

A further source marking is the displacement of the first Hamilton arc after normalizing p0=(0,0). These are the ten directed-arc orbit types because the automorphism group acts regularly. Their exact incident dart counts from the subset DP in replay_fano_count.cpp are:
```
(1,0)   585723718963
(2,0)   592255137742
(3,0)   594981636439
(0,1)   571650274131
(1,1)   584204812473
(2,1)   585280890671
(4,1)   590182955148
(1,2)   585182864820
(2,2)   595482271370
(4,2)   590121718978
sum    5875066280735
```
The underlying H is 123376391895435, so v_3(H)=v_7(H)=1, not a deeper 3- or 7-adic factor. This incidence counts weighted lines with multiplicity, not distinct adjacent classes. In the full merged black multigraph each incidence is doubled. The DP fixes the start, computes on the remaining set a suffix count with prescribed first vertex, and sums over its 10 out-neighbors. No 64-bit overflow is possible: even for the whole graph a crude upper bound is 21 times the product over remaining choices r=20,...,1 of min(10,r), below 2^64. The computation and the port quotient do not imply that H(U)/21 itself is attained as a tournament path count.

## Six merged nodes from a complete Fano twist family

Reverse the cycle inside F_a independently at each of the seven points a, recording the chosen set as a mask m. All 128 resulting tournaments remain regular (one internal out-arc at each vertex) and strong (contract the seven strong triangles and use the strong Paley majority quotient). The canonical partition survives: for a pair inside one triangle its common external neighbors are unchanged and its common internal neighbors number zero. For cross pairs only the orientations of the two triangles containing its endpoints can affect the intersection. Use the original affine symmetry to normalize the first endpoint to (0,0); for second endpoint (a,b) the maxima over the two independent flips are:
```
b=0, a=1..6: 5,4,3,3,4,5
b=1, a=1..6: 5,5,4,4,5,5
b=2, a=1..6: 5,4,5,5,4,5
```
These 72 local cases are simply the two membership tests per endpoint from the S slices, with the two internal cyclic bits independently toggled (literal replay included). No cross pair can acquire the signature 7.

Dissent ports remain unchanged. Thus an isomorphism between two twists is exactly an affine Paley permutation carrying one mask to the other; every such permutation has its original unique lift. A converse isomorphism is an affine permutation composed with inversion on a, with the mask also complemented, because the normalized inversion preserves port indices. This converse action changes the size of a mask from k to 7−k, so no twist is self-converse. Ordinary Burnside yields
```
(2^7 + 6·2 + 14·2^3)/21 = 12,
```
using one identity, six order-seven rotations (one cycle), and fourteen order-three rotations (cycle lengths 1,3,3). There are exactly six merged classes. They have unique types for mask sizes 0,1,2 and three types at size 3: transitive triple, Fano-line directed triangle, nonline directed triangle. For instance the respective masks 7,11,13 represent those three types. Their orbit sizes are respectively 21,7,7 and cover all 35 triples; this can be checked by the affine action or by directly partitioning triples by their oriented and Fano types.

All six merged vertices form an independent set with no loops in the actual metagraph, since every one is regular and every neighbor has one exceptional degree at either end. They have the following fundamental black degrees, from H/|Aut| with |Aut|=21/(mask orbit size):
```
mask    |Aut|     H                    black degree
  0       21      123376391895435        5875066280735
  1        3      123381234928431       41127078309477
  3        1      123384468287435      123384468287435
  7        1      123386086601285      123386086601285
 11        3      123385527546237       41128509182079
 13        3      123386644838787       41128881612929
```
The six H values are distinct, so this tiny orbit classification yields an exact nontrivial family of six distinct spectra even though cross-triangle incidences and degree sequence do not change. Run replay_fano_family.cpp on these six masks for an independent full-start DP. No claim is made about arbitrary tournaments with a seven-triangle quotient.
