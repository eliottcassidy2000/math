# Five-count windows with arbitrary binary slopes

The [same proof](PROOF.md) allows larger constant-width seed
banks. Here is a separate bank of size five, all rooted at
vertex 8 in zero-based outgoing masks.

| Rows | h,g,R,S | Intercept I |
|:--|:--|--:|
| (2,260,273,7,43,15,319,383,57) | 83,15,88,94 | −75 |
| (0,257,267,35,271,87,415,319,41) | 79,15,91,91 | −73 |
| (0,257,323,23,39,271,315,383,25) | 81,15,91,91 | −71 |
| (272,5,9,259,14,31,319,383,54) | 89,15,88,94 | −69 |
| (0,257,11,275,71,287,431,319,21) | 85,15,91,91 | −67 |

They share g=15 and Z=167. A common near-sink append makes
all their continuation profiles identical (g,R,S)=(15,167,182).
The arbitrary binary address b of length n therefore gives
Z=334·2^n+15b, with all five offsets preserved under any common
further appends. The proof for distinct absolute windows changes
only the constants in the collision argument: 82 becomes 334,
g=9 becomes 15, and a disjoint selection loses at most a factor
of five rather than three. There are thus Ω(X^(2/3)) disjoint
**five-count** windows at all large scales. The new seeds can
have neutral singleton components; this claim is about actual
Hamiltonian counts, not the spectrum restricted to strong
tournaments. All matched reader swaps preserve every append arc
and can be made on the same vertex set.

For discovery, append a new root to each small-order atlas
parent. Its outgoing set Q gives A by summing the weighted
old root-path totals of the vertices in Q. If W_ij counts
weighted even-order old paths from i to j, its intercept is
exactly g−2 Σ_(i,j∈Q) W_ij: the new root count gains the
Q-to-complement cut, and the subtracted rooted odd term
uses all destinations. Hence a table keyed by (g,A)
finds consecutive intercepts, rather than merely consecutive
counts. seed_search.cpp replays the discovery with the complete
order-eight CSV and g≤121. The five-row statement follows
from direct profile replay and does not require atlas
completeness or a claim of maximal bank width.
