# Completeness and nonduplication at degree thirteen

## 1. Inherited small transitive images

The predecessor [complete census through ten](../odd_mixed_prime_census_20261020/README.md) proves that the odd transitive images of degrees 3,5,7 are C3, C5, and either C7 or F21. Its invariant-tournament long-cycle argument shows that an odd transitive image of degree nine is a 3-group. The [degree-eleven/twelve extension](../mixed_prime_census_extension_20261020/README.md) proves that the odd transitive degree-eleven images are C11 and F55, including an elementary neighborhood certificate.

Here is the degree-nine exclusion in the form needed. An odd permutation group preserves a tournament, because no element exchanges the endpoints of a pair. Transitivity makes that tournament regular. If its order were divisible by 5 or 7, Cauchy's theorem would give a single p-cycle with 9/2<p<9. A vertex outside the p-cycle has either all p points among its outneighbors or none, so its outdegree is at least p or at most 8−p, never four. Thus only the prime three divides its order. This also handles a nine-orbit accompanied by a three-orbit at degree thirteen: the whole subgroup is a 3-group.

## 2. An explicit degree-thirteen certificate

**Lemma.** Every odd transitive subgroup of S13 is conjugate to C13 or F39=C13⋊C3.

An odd transitive H contains a 13-cycle. Label it by translations of F13. An H-invariant tournament must be a cyclic tournament T(S), where S selects one element from each of the six pairs {i,−i}, 1≤i≤6. This is exactly 64 sets, not a sample of tournaments.

For each S the independently executable certificate `prime13_certificate.py` does the following.

1. Form the entire 13×13 adjacency matrix A.
2. Let R=F13\({0}∪S), the six in-neighbors of zero. Verify that the six signatures (A[v,s])_(s∈S), v∈R, are pairwise distinct.
3. Inspect all 6! bijections of the out-neighbor set S. An automorphism fixing zero must preserve this subtournament. For each surviving bijection, the signatures force at most one image for each point of R.
4. Check the proposed extension on **all ordered vertex pairs**. The accepted extensions are precisely all automorphisms fixing zero.
5. Compare their actual vertex permutations with all multiplicative maps x↦a x preserving S. They agree exactly, not merely in group order.

The output is 60 tournaments with full automorphism group C13 and four with full group F39. Since translations are transitive, multiplying the exact root-stabilizer order by thirteen gives the full group order. A transitive subgroup of C13 or F39 contains their unique Sylow-13 subgroup and is C13 or F39. This proves the lemma with an explicitly bounded finite certificate. It uses neither a subgroup database nor a presumed affine classification.

## 3. Reduction to vertex orbits

Every vertex orbit of an odd-order subgroup has odd size. Through degree thirteen the possible nontrivial sizes are 3,5,7,9,11,13. The subgroup is faithfully a subdirect subgroup of the product of its actual transitive images. A subdirect subgroup of two groups of coprime orders is their full direct product: its order is divisible by both image orders and therefore their product.

The degree-thirteen mixed-prime orbit patterns are exactly:

| Vertex-orbit pattern | Number of classes |
|---|---:|
| 13 | 1 |
| 11+1² | 1 |
| 7+5+1 | 2 |
| 7+3+3 | 7 |
| 7+3+1³ | 3 |
| 7+1⁶ | 1 |
| 5+5+3 | 2 |
| 5+3+3+1² | 2 |
| 5+3+1⁵ | 1 |
| **Total** | **20** |

All remaining partitions produce prime-power groups. In particular an orbit of size nine leaves room for at most one triple and yields only a 3-group; two five-orbits without a triple yield only a 5-group; only triples and singleton orbits yield 3-groups.

The patterns with thirteen or eleven points give F39 or F55. For 7+5 the coprime-product argument gives C35 or F21×C5. For 7+3, a cyclic seven-image gives C21. A Frobenius seven-image has only the trivial common quotient or its quotient C3 in common with the triple-image. It gives F21×C3 or the graph of F21→C3. Reversing the triple makes the two possible surjections conjugate. On the lone seven-orbit the only mixed-prime image is F21.

For 5+3+3 the actual three-image is either a diagonal C3 or C3². All nonzero slopes are conjugate by independently reversing a triple. Taking its direct product with C5 gives two classes. For 5+5+3 the identical argument over F5 gives diagonal C5 or C5², independently times C3. There is only one class on 5+3+1⁵.

## 4. Exactly seven classes on 7+3+3

Write A for the seven-orbit, B,C for the two triples. If the A-image is C7, coprimality separates C7 from the three-image on B∪C. The latter is either a diagonal C3 or C3². These are two classes.

Suppose instead that the A-image is F21. The subgroup contains the full translation C7: its intersection with C7 has order divisible by seven because the other orbit-image orders are powers of three. Quotienting by C7 leaves a vector subspace

V ≤ F3³, with coordinates (u,b,c),

surjective onto each coordinate. Here u is the quotient F21/C7 and b,c are the actions on the triples. Ambient conjugation can independently negate b and c and can exchange them. It does **not** negate u: the normalizer of F21 on its seven-orbit is AGL(1,7), whose conjugation induces the identity on its abelian multiplier quotient.

The coordinate-surjective subspaces have precisely five orbits under these allowed transformations:

1. V=F3³: F21×C3², order 189.
2. V={b=c}: F21×diag(C3), order 63.
3. V={u=b}, c free: F21diag×C3, order 63.
4. V={u=b+c}: the plane-kernel-diagonal fibre product, called “F21diag×diagC3” in the incoming list, order 63.
5. V=span(1,1,1): the doubly diagonal F21, order 21.

For completeness, a one-dimensional subdirect space has all three coordinates nonzero; independent triple negations normalize it to the displayed line. A two-dimensional space is the kernel of one nonzero linear form αu+βb+γc. Surjectivity excludes forms supported on just one coordinate. If its support has size two, it either omits u, giving b=c, or involves u and one triple coordinate, giving u=b after a possible triple exchange. If its support has size three, sign changes and multiplication of the defining equation by a nonzero scalar give u=b+c. Full dimension gives the entire space. This is exhaustive.

The executable `subdirect_certificate.py` independently enumerates all subspaces and all allowed coordinate transformations. There are 15 raw subdirect subspaces of F3³ and five conjugacy classes; the C7 case has three raw subdirect subspaces of F3² and two classes. It extracts the quotient coordinates from the actual seven proposed permutation groups, checks that their conjugacy classes are all distinct, and checks exact equality with the exhaustive subspace list.

Nonduplication can also be read directly. The order and seven-image distinguish the one- and three-dimensional cases and the two C7 cases. Among the three order-63 groups with seven-image F21:

- V={b=c} has image of order three on B∪C and coincident pointwise triple kernels;
- V={u=b} has image of order nine on B∪C, with one triple kernel C7×C3 and the other F21;
- V={u=b+c} also has image of order nine on B∪C, but both triple kernels are F21.

These are invariants under arbitrary ambient conjugacy, even when the two triples are exchanged. Groups in different vertex-orbit patterns cannot be conjugate. Consequently all twenty degree-thirteen classes are distinct as well as exhaustive.

## 5. Promotion boundary

Combining the 20 classes above with the inherited 28 degree/class pairs through twelve gives 48. The classification does not rely on numerical equality surviving in a table. Conversely the classification alone does not prove any full-isomorphism count; those are separately certified in the census and its independent checks.

No argument above classifies degree fourteen. The all-degree inequality tau_T≥tau_E remains an open assertion, not a consequence of this finite census.
