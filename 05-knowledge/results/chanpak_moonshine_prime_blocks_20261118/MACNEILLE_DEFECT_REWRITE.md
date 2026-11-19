# MacNeille cuts and a physical rewrite of the entire defect-six family

**Proved.** Every canonical odd target of defect six, $N_j=14\,2^j+1$
for $j\ge0$, has a tournament of exactly its canonical order and
Hamiltonian count $N_j$, obtained by reversing old arcs of the
canonical source itself. For all $j\ge4$, one fixed six-arc patch
works. It changes no edge incident to the growing transitive tail.
Combining with [the previous descent](LOW_DEFECT_DESCENT.md) closes
defects 0, 2, 4 and 6, except for the two proved exclusions 7 and 21.
Any further counterexample has canonical defect at least 8.

## The lower cut

Let $F_X(m)$ be the Hamiltonian count of marked tournament $X$
after replacing its mark by a transitive tournament of order $m\ge1$.
Let $c_r(X)$ count legal compressed Hamilton words with $r$ runs
of the marked block. Order the finite signed profile vectors
$$
p(X)=(c_1,-c_1,c_2,-c_2,\ldots,c_R,-c_R)
$$
coordinatewise, retaining only vectors with physical preimages for the
final selection. Here $R$ bounds the runs in the finite candidate bank,
with zero padding. A requirement $L\le c_1\le U$ is the intersection
of two principal downsets in the ambient product of extended integers,
one for each sign. Each equality $c_r=d_r$ likewise gives two.
Their intersection is a [MacNeille lower cut](https://pmc.ncbi.nlm.nih.gov/articles/PMC1076704/) in that product: a set of
common lower bounds, not a coordinatewise combination of witnesses.
A repair exists in the bank exactly when this cut contains **one
physical profile**. In particular intersection with the physical bank
must not be replaced by intersecting independently projected ranges.

There is an exact reason to keep these coordinates. For any single
transitive replacement the response is
$$
F_X(m)=\sum_{r\ge1} c_r(X)\;r!\,{m\brace r}.
$$
A compressed word fixes the order of its runs; partition the $m$
distinct vertices into $r$ labeled nonempty sets and traverse each
in increasing order. By inclusion-exclusion $r!{m\brace r}$ is
$r^m$ plus a linear combination of smaller exponential bases.
Consequently a difference $F_Y-F_X$ is bounded over all widths iff
all coefficients at $r\ge2$ match: the largest differing run index
would otherwise dominate. When they match the difference is exactly
$c_1(Y)-c_1(X)$. So for a ray with fixed initial defect $D>0$,
**every** uniformly improving, nonovershooting replacement in this
model comes from the finite cut
$$
c_r(Y)=c_r(X)\quad (r\ge2),\qquad
c_1(X)+1\le c_1(Y)\le c_1(X)+D.
$$
For an exact repair, require equality at the upper end. This is an
all-width criterion for single marked transitive inflation, not a
universal guarantee of a physical point. The defect here is the
original target minus the current ordinary count; the terminal reader
is the new ordinary count. There is no assertion that the old rooted
binary decoder continues to have its old internal state after a patch.

## The fixed patch

Use canonical labels in birth order, with 0 the retained root. Start
from the edge 0→1, append two sinks, then append $j+1$ near-sinks
(each new vertex beats only 0 at birth). The canonical word is
$001^{j+1}$. At $j=4$ its outgoing bit masks are
```
14 508 504 496 481 449 385 257 1
```
and mark vertex 8. Reverse the arcs on the six unordered pairs
```
(0,1) (0,2) (0,5) (1,7) (3,4) (3,7)
```
giving
```
40 381 505 352 489 448 385 267 1
```
In both quotients the mark beats only 0. Any run except the last
must exit to 0, so there can be at most two runs. Their complete
profiles are respectively $(219,112)$ and $(225,112)$.
The included replay computes these directly from permutations of
the eight unmarked vertices and placements of one or two copies
of the mark. Two copies use distinct gaps, with 0 immediately
after the first copy.

At general $j\ge4$, replace vertex 8 by $TT_{j-3}$.
This is literally the original larger canonical source; the patch
commutes with the substitution because it does not touch the mark.
The two counts are
$$
219+112(2^{j-3}-2)=14\,2^j-5,\qquad
225+112(2^{j-3}-2)=14\,2^j+1.
$$
Thus the patch raises the count by six, preserves order and the
full growing-block response, and reaches the original absolute target.
It is a simultaneous patch, not a sequence of neutral intermediate
switches: among the 64 subsets of its support, only the empty and
full sets preserve $c_2=112$. The intersection with $c_1=225$
is the full patch alone **in that cube**, with no global minimum
edit-distance claim.

At the four smaller heights, the following exact finite patches apply
directly to the word $001^{j+1}$.

| $j$ | Target | Unordered pairs reversed |
|---:|---:|---|
| 0 | 15 | (1,3) |
| 1 | 29 | (2,4) |
| 2 | 57 | (0,2), (4,6) |
| 3 | 113 | (0,5), (5,7) |

These are checked by subset DP counting paths with a prescribed last
vertex. The same independent DP checks the six-arc cube and the first
six widths; the compressed word proof handles every width.

## What is left

The cut rules out any attempt to trade a nonconstant response for a
bounded gain over an entire ray. It also tells us which coordinate can
move: the one-run count, with the higher-run responses locked. For the
next defects, the constructive question is whether changing the edit
support or the marked quotient supplies a physical point in the cut.
The present proof closes one entire infinite layer, not all such cuts
and not Chan–Pak itself.

Run `python3 macneille_defect_replay.py`.
