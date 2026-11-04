# Classify every minimum-increment marked interval chain

**PROVED consequence of INTERVAL_CHAINS.md; independently audited.**
This classifies equality inside Moon's interval chains for every order n>=9.
It does not assert that arbitrary strong minimizers are interval chains.

Let F denote the interval-chain family, and retain the marked vertex in all
isomorphisms. Chains on at least four vertices are rigid and their spine
composition is recovered from the score order and the orientations within
the two tied end pairs, as proved in inherited THM-4744. Thus a composition
and marked interior position specify one rooted class without ambiguity.

For k>=0, the complete minimizing rooted classes are:

1. **n=9+3k.** One part2 and k+2 parts3; the marked vertex is the middle
   of the part2, and that part is not an end part. There are k+1 classes.
2. **n=10+3k.** One part2, one part4 and k+1 parts3; the marked part2 is
   internal, with immediate neighbors3 and4 in either order. There are
   2(k+1) classes.
3. **n=11+3k.** Two parts2 and k+2 parts3; the marked part2 is internal
   with both immediate neighbors3. There are (k+1)(k+2) classes.

The increments are respectively26*5^k,46*5^k,78*5^k.

## Why these are all the equality cases

The eight terms in the preceding minimum formula have coefficient/cost
pairs (2,3),(4,4),(6,5),(10,6),(16,7),(26,8),(46,9),(82,10).
After normalizing by5^k, their values in the three residue classes are:

| coefficient/cost | n=9+3k | n=10+3k | n=11+3k |
|---|---:|---:|---:|
| 2/3 | 30 | 50 | 90 |
| 4/4 | 36 | 60 | 100 |
| 6/5 | 30 | 54 | 90 |
| 10/6 | 30 | 50 | 90 |
| 16/7 | 144/5 | 48 | 80 |
| 26/8 | 26 | 234/5 | 78 |
| 46/9 | 138/5 | 46 | 414/5 |
| 82/10 | 738/25 | 246/5 | 82 |

An entry with an invalid residual size (negative or1) is infinity instead:
in column1, row5 needs k>=1, row7 needs k>=1 and row8 needs k>=2;
in column2, rows6 and8 need k>=1; in column3 row7 needs k>=1.
The other entries apply already at k=0. Hence each column has its stated
unique winning coefficient/cost:26/8,46/9,26/8 respectively.

Every reduction that eliminates an unmarked part>=5, a marked part>=5 or
an interior neighbor>=5 is strict. The only non-strict neighbor reduction
was the end-neighbor5 case, which reduces to an end-triangle term; all those
terms are strictly above the winning values. A part4 among outside blocks
is equivalent in product to two2s, but the optimal outside residual is
0 modulo3 in cases1/2 and2 modulo3 in case3. Its optimal compositions are
respectively all3s, or one2 and otherwise3s. Therefore no outside4 or
other hidden equality case survives.

The winning26/8 is exactly a marked triangle between blocks3,3, and the
winning46/9 is exactly a marked triangle between blocks3,4 or4,3. Placing
the optimal outside blocks gives precisely the three families above.

## Count the rooted classes

In case1 there are k+3 blocks total, and the sole part2 can occupy any of
the k+1 internal positions.

In case2, choose the order of neighbors3/4 (two choices) and distribute the
k remaining part3s before and after this three-block segment (k+1 choices).

In case3 there are L=k+4 blocks. Choose the marked part2 at one of the
L-2 internal positions. The other part2 can lie anywhere except the marked
position or either adjacent position, giving L-3 choices. This counts each
rooted word once: (L-2)(L-3)=(k+2)(k+1). Even when both triangles are
eligible marks, they are distinct rooted objects because the chain is rigid.

The first values by n=9,...,20 are

    1,2,2,2,4,6,3,6,12,4,8,20.

Converse reverses the composition and transports the mark. Case1 has one
self-converse rooted minimum exactly when k is even; its marked triangle
is the central block. Cases2/3 have none: case2's neighboring lengths differ,
while in case3 a root-fixed reversal would fix one part2 centrally but would
leave the sole other part2 unpaired. This is rooted converse symmetry, not
an assertion about ordinary automorphism groups (which are already trivial).

`verify_chain_equalities.py` checks the classification, the class count and
the rooted converse statement directly across every composition through
order24, using the exact proven deletion formulas. These checks supplement
the combinatorial proof; they are not a substitute for its all-order scope.
