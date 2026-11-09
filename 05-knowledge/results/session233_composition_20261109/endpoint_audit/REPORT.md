# Independent audit: binary endpoint universality

Verdict: PASS for the stated endpoint theorem and exact response identity. No claim of complete Hamilton-path-count spectrum, global quotient universality, or Chan–Pak completeness follows.

## Proof audit

1. Conventions: H counts directed ordered spanning vertex words, with singleton count 1; the response definition allows the complementary path to be empty with H(empty)=1. The endpoint is the LAST vertex. Response B is rooted at the FIRST vertex of its distinguished path. These different orientations are used consistently.
2. The base B_n contains the directed spanning cycle, hence is strong. Every path either has no descent (the unique increasing word), or has the unique descent n-1 -> 0. The partition argument gives precisely e=(1,1,2,...,2^(n-3),1); its sum is 1+2^(n-2).
3. For n>=4, the degree multiplicities are at most two. Degree preservation keeps each class invariant, and a permutation cannot exchange the two vertices of an oriented pair. All vertices are therefore fixed. No claim of base rigidity at n=3 is needed (and it would be false there).
4. Deleting the last x gives exactly the disjoint union of base paths whose last vertex beats x. The chosen interior incoming set selects binary weights without using either endpoint's extra weight 1. Every integer from 1 through 2^(n-2)-1 occurs, even values included.
5. Any positive selected set supplies an incoming edge to x; x beats both 0 and n-1. Both directions to the strong base imply strongness of the extension. An automorphism fixing x induces a base automorphism, hence is trivial. The entire automorphism group of U need not be assumed trivial.
6. For positive integer N, ceil(log2(N+1)) equals N.bit_length(). Thus n=max(4,N.bit_length()+2) gives exactly the stated bound n+1=max(5,ceil(log2(N+1))+3), with no power-of-two ambiguity. At N=1 the construction has order 5, e_x=1 and H=9. The n=4 selectors N=1,2,3 are all covered.
7. The endpoint fibre is the orbit set under Aut(U)_x, not division by Aut(U). Trivial root stabilizer gives fibre N. The free action is immediate because an automorphism fixing an ordered spanning word fixes every vertex. This also directly justifies the invoked THM-4988 step.
8. Prepending ordered C3 strong components multiplies both endpoint count and root stabilizer by 3 per component. SCCs are individually fixed because the condensation is a finite linear order. These towers are generally NOT strong; the text does not say otherwise. The strong witnesses are the original U_N.

## Exact sidecar

The formula B_x(A+x)=H(A)+sum_(x->v) B_v(A) is valid for every tournament A, regardless of strongness. If the rooted path is just x, the other path is a full path of A. Otherwise remove x from the rooted path and retain the second path unchanged; the next vertex v must be an outneighbor. The map and its inverse preserve every ordered pair of paths, including an empty second path.

For the binary selector this yields exactly the displayed complementary-digit response formula. These response weights are not the endpoint powers of two. The n=6 H and B profiles in PROOF.md are reproduced exactly by literal permutation-and-cut enumeration. Nonmonotonicity is a property of this finite profile, not an impossibility theorem for other compositions. In particular the endpoint theorem gives no independent choice of (H,e,B), and the one-exit and two-base identities still require the other coordinates.

## Independent checks

`audit.py` does not import or call producer code. It builds arcs as ordered-pair sets, enumerates permutation words, and counts ordered two-path covers by cutting each permutation at every permitted location. This directly counts B without the producer's subset start-path convolution. A separate memoized deletion recurrence checks endpoint counts. Strongness uses Floyd–Warshall closure. Automorphisms are literally enumerated and arc-tested.

- Base orders 4 through 7: exact endpoint vector and actual rigidity, plus strongness.
- All 25 selectors from base orders 4,5,6: full literal endpoint vector, independently recursive vector, actual root stabilizer, strongness, and literal response identity.
- All 246 selectors from base orders 4 through 9: recursive endpoint equality and strongness.
- All 1,098 one-vertex extensions of all labelled tournaments of base orders 1 through 4: general response identity by literal cover counts. This includes empty and full outgoing sets and non-strong controls.
- Named witnesses (N,order,H): (1,5,9), (7,6,23), (21,8,187), (37,9,461), (101,10,1573), (257,12,9899).
- `python3 audit.py` and `python3 -O audit.py` pass with byte-identical receipts. Checks use explicit RuntimeError guards, never assert. `--fail` deliberately fails under both modes; receipts record the exception.

The auditor initially supplied an incorrect guessed expected H=7 for the added N=1 control, and the guard rejected it. This was corrected to the independently computed H=9; it was not an error in the producer's claims. The subsequent explicit failure sentinels are separate controls.

## Dependencies and novelty scope

Read and checked THM-4988, THM-250, and THM-4865 in the designated worktree. THM-4988 supplies the exact endpoint orbit fibre and isolated 7/21 witnesses. THM-250 supplies the near-transitive family and total count; its erroneous inclusion of 2 in its corollary is correctly excluded by the producer. This packet's endpoint vector is independently established rather than relying on that corollary. THM-4865(f) explicitly credits HYP-9457(c) for the two-base law, as stated in this packet.

The meaningful new result relative to those checked sources is the universal binary endpoint selector with strongness, trivial rooted stabilizer, and logarithmic order. This is a scoped lineage comparison, not an exhaustive repository or literature novelty certification. The general response identity is a natural multi-outneighbor generalization of the inherited sole-exit identity, with its own valid elementary proof. No global spectrum completion is certified.

No producer files or repository files were edited; no commits were made.


Session233 close-out note: THM-250's erroneous corollary entry2 and its adjacent-flip proof sentence have now been repaired in canon; the formula used above is unchanged. See MISTAKE-744.
