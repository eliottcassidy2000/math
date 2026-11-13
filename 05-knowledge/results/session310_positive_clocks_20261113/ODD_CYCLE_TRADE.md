# The collision trades six triangle pairs for twelve pentagons

The labeled strong graphs A and B in [NEUTRAL_TAILS](NEUTRAL_TAILS.md) satisfy Q_B(t)=tQ_A(t). The identity has a direct combinatorial mechanism, not just an interpolation certificate. Enumerating directed odd cycles up to rotation and disjoint pairs gives:

| block | triangles | pentagons | heptagons | disjoint triangle pairs |
|---|---:|---:|---:|---:|
| A (seven vertices) | 10 | 18 | 7 | 6 |
| B (eight vertices) | 10 | 30 | 7 | 0 |

No other nonempty odd-cycle packings can fit these orders. Padding A by a final singleton leaves the cycle list unchanged. By [the positive packing expansion](POSITIVE_CLOCK.md), a single pentagon on an eight-vertex graph contributes 2t⁴ and a disjoint triangle pair contributes 4t⁴. Therefore the twelve additional pentagons of B replace precisely the weight of six triangle pairs of A→TT1, with the other degrees unchanged. Both graphs have

    Q(t)=t⁸+20t⁶+60t⁴+14t².

The mechanism does not assert a particular pairing between individual cycles or a local arc rewrite. More generally the observer remembers only the sum of (cycle length − 1) and the *weight*, not the decomposition of that even contraction among cycle lengths. A packing of k cycles with contraction 2j has the same monomial weight as 2^(k−1) individual (2j+1)-cycles. The algebraic identity is general; whether a directed graph can carry out any given trade remains a separate construction problem. Here both graphs are literal, the neutral substitution is exact in every uniform context, and their triangle-packing numbers differ, so the finer hypergraph can detect a physical change that the full slot polynomial does not.

## The hidden packing-rank gap grows in globally strong graphs

The maximum number of disjoint directed triangles in A is two, and in B is one (both statements follow from the complete packing table). Consequently the ordinal joins A^r→TT_r and B^r have packing ranks 2r and r. Wrap both modules into C3[block,TT1,TT1]. This context is strong even when the internal block is not: routes through the other two cyclic slots join every pair of internal vertices in both directions. Every cross-triangle must use both external singleton vertices, so at most one can accompany an internal packing. Conversely each maximum internal packing leaves at least one vertex unused, and that vertex together with the two singletons is a cross-triangle. The packing ranks of the resulting strong graphs are exactly 2r+1 and r+1. Their Hamiltonian path counts agree for all r by uniform substitution. In fact their *full* slot polynomials agree: on each color separately, the same itinerary bijection uses an ordered cover of the block across all colors, so equality of its F-profile controls every color coefficient. Thus even a growing additive triangle-packing-rank gap can be invisible to the full slot polynomial in strong tournaments.
