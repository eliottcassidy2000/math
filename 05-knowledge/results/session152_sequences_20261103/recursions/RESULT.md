# Ear cut reciprocity and path-cover representation

Status: elementary identities proved for every finite tournament. The five-vertex data are exact exhaustive checks, not a resolution of the spectrum conjecture. The response polynomial and curvature substantially overlap THM-4114, THM-4115, THM-4118; their convention takes S to be the OUT-neighbors of the new vertex, whereas S below is its IN-neighbor set. The example's ten-value image is already documented in THM-4118. Session 151's incidence packet already observes that an anti-automorphism fixes at most one vertex. No claim of literature priority is made for the reciprocity observation.

Let h=H(T), S subset V(T), and F(S)=H(T+x_S), with i->x iff i in S. For i != j let M_ij count pairs of nonempty directed paths whose disjoint vertex sets partition V(T), the first ending at i and the second starting at j. Let start_i and end_i count old Hamiltonian paths with the indicated endpoints. Splitting a new Hamiltonian path at x gives

  F(S)=h+sum_{i in S}(end_i-start_i)+sum_{i in S,j outside S} M_ij.  (1)

Thus F is quadratic, with mixed difference at i,j equal to -(M_ij+M_ji). In particular it is submodular. Redei's parity theorem says each F(S) is odd, so beta_ij=(M_ij+M_ji)/2 is an integer. Define

  d_i=end_i-start_i + (sum_j M_ij - sum_j M_ji)/2.

Then d_i are integers, sum_i d_i=0, and exactly

  F(S)=h+d(S)+beta(delta S),
  F(S)+F(V-S)=2h+2 beta(delta S).                              (2)

Here each edge crossing the cut is counted once. The integer claim follows because row-minus-column and row-plus-column have the same parity. Adding one vertex cannot decrease H: for every old path, either x can precede its first vertex, follow its last vertex, or fit between a consecutive in-neighbor and out-neighbor (scan the path). A fixed rule chooses one insertion for each old path, giving an injection. Consequently (2) implies the nontrivial universal inequalities |d(S)| <= beta(delta S), for every S. These are exact constraints, not a claim of interval coverage.

## Converse reciprocity

Suppose T is self-converse and f is any anti-automorphism. Reversing every arc and relabeling by f proves

  F(S)=F(f(V-S)).                                               (3)

Hence the bivariate enumerator E_T(u,z)=sum_{S subset V} u^|S| z^F(S) obeys the exact polynomial reciprocity

  E_T(u,z)=u^n E_T(u^{-1},z).                                  (4)

For odd n, f has exactly one fixed vertex. An odd cycle of f with length >1 would return two of its vertices after an odd number of arc reversals, a contradiction. Two fixed vertices would also reverse their mutual arc. Since the other cycles are even, odd n forces exactly one fixed vertex. The permutation g(S)=f(V-S) toggles membership of this vertex. Every orbit of g therefore has even length, and (3) implies that every ear-count multiplicity is even, including after excluding empty and full cuts. No claim that g is an involution is needed: g^2(S)=f^2(S), which may differ from S.

## Exact example

Take lexicographic upper-pair code 8 on vertices 0,...,4, where bit 1 means lower->higher. This has h=9. The anti-automorphism f=(0 4)(1 3)(2) is unique among the 120 permutations. Ordered by lexicographic pairs, the ten beta values are

  (3,2,2,4,5,5,2,5,2,3),  d=(3,9,0,-9,-3).

Thus (2) provides a 10-edge weighted-cut certificate for each of its 32 responses. Explicitly,

 E(u,z)=(1+u^5)z^9
       +(u+u^4)(z^15+z^17+2z^23+z^33)
       +(u^2+u^3)(z^17+z^19+z^23+z^25+z^27+z^29+z^33+2z^37+z^41).

The proper nonempty-cut image is {15,17,19,23,25,27,29,33,37,41}; 21,31,35,39 are absent between its extremes. Even multiplicities do not establish these holes or determine the image.

## Reproduction

Run `python3 check.py`. It constructs all old permutations and ordered path covers, forms all 32 children, and independently counts their paths by direct permutations. It checks (1), (2), all anti-automorphisms, and reciprocity under each such map, and prints the exact multiplicities. This finite computation supports the example; the general statements above are proved separately. The code is intentionally small and does not claim independent novelty over the indicated records.
