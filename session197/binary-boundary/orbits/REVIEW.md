# Independent review of the quadratic boundary's converse orbits

**PASS, no correction.** This review binds the frozen `session197-orbit-wildcard/boundary-orbits/PROOF.md` SHA256 `1fa1e66a16168cc28c74b6306b820aca84677d96e4f2a7b8d7f8fc174fd0b2ce`, standalone verifier `6db533e31ccd662787e87fbad9905e659333d666e1d1e2a804acf8c0c067bdc5`, and optimized/debug receipt `c81b622613033e680f4c6e70adc028204b58ffe4920d15b81d3132e2945cdedb`. The source files were left unchanged.

## Proof and normalization

The fixed-word decomposition is complete. The exceptional O_a->I_a edge is fixed by reversal followed by the designated anti-involution. If used in a fixed Hamilton word, it must occur centrally, since a noncentral occurrence would force the identical edge to occur again. The first half then consists of exactly the q complementary-path arrangements. Otherwise the first half is entirely in I and its last vertex is not I_a, giving h-e choices. The cross joins cannot use the missing reversed arc elsewhere. This proves h+q-e, including the singleton boundary. The proof concerns the designated anti-involution for a general side graph and does not silently assume rigidity there.

For the quadratic family the degree list is exact: all degrees from 2 through n-3 are unique, with the two endpoint degrees each repeated twice. A degree-preserving automorphism cannot swap a tied pair, because doing so reverses their mutual directed edge. Hence every vertex is fixed. Composing two anti-isomorphisms would give an automorphism, so the known anti-involution is unique. The assertion of one ordinary class and n! labelled versions is correctly restricted to this one specified family member at each even order, not all self-converse C=1 classes.

The fixed count 2h-1 and Burnside formulas follow with the correct factors of two. The normalized fibre weight is distinct from both the ordinary class count and the unrestricted Hamilton-word mass. I read the named `CONVERSE_MASS.md` dependency and verified its SHA256 `07f4b20dde9f9d23b7711f65bd484d7468bdbef3cc6da0cc2888aa91d402a620`. In this rigid case, the normalization can also be justified directly: distinct fixed Hamilton words relabel to distinct normalized adjacency matrices, since equality would give an automorphism. The independent-bit normalized universe has exponent t(t-1) at n=2t. Thus the fraction (2^(t-1)+1)/2^(t(t-1)) has the right numerator, denominator and interpretation. There is no conflation with uniform sampling of ordinary tournament classes.

## Independent controls

The independent `check.cpp` imports no source checker and does not enumerate pair-transversal half words. Instead it enumerates complete Hamilton words by literal DFS and applies the anti-involution to each full word.

- It checks the general fixed-word formula on all 285 rooted labelled side tournaments through order four and 40 deterministic order-five controls, for 325 side cases.
- For all six family graphs, of orders six through sixteen, it explicitly forms the set of converse-orbit representatives and the set of normalized adjacency matrices obtained from fixed words. The fixed counts are 5,9,17,33,65,129; the orbit counts are 11,29,89,305,1121,4289; and the normalized-set sizes equal the fixed counts exactly.
- It independently enumerates every degree-compatible permutation, checking one automorphism and one anti-isomorphism for each family graph. It also checks the full displayed degree multiset.
- It counts the free reflected edge orbits after forcing the normalized Hamilton path, obtaining exponents 6,12,20,30,42,56.
- The displayed 14-vertex word is checked for distinct vertices, every directed arc and the full reflection equation.

Independent optimized and unoptimized assertion-enabled runs agree. The supplied standalone verifier was also compiled and replayed under both modes; both outputs exactly match its frozen receipt. These finite controls support the implementation, while the fixed-word decomposition and degree argument establish the unbounded claims.

This review is separate from the earlier binary deletion-boundary review and from any later full-bin mirror operator. It makes no additional arithmetic-coverage or literature-priority claim.
