# Which cells do 625, 600, and 576 count?

**PROVED:** the explicit set maps and all-length ordinal-product formulas
below. **FINITE-EXACT:** the direct controls in `square_counts.py/json`.
The arithmetic `625+125+25+1=776` is correct. The separate equality
`576=(25-1)^2` deletes a distinguished row and column, not a diagonal.
This distinction can be realized on actual tournament Hamiltonian paths.

## Inheritance

The literal square/blue-diagonal mechanism and the balanced-family deletion
map are already proved in
[microcosm.md, sections 2 and 4](../zero_heavy_ports_20261020g/microcosm.md),
building on [double_return_counts.md](../strong_blue_games_20261020f/double_return_counts.md).
Their family has q=2^(k-1)+1, H=q^2 and blue count q. In particular q=25
is not a member of that binary family. Its arithmetic must not silently
change q into H. The ordinal-product mechanism below supplies a different,
explicit source for the user's numbers.

## 1. A four-vertex source for the powers of five

Let Q have vertices 0,1,2,3, every arc forward in this order except
2->0 and 3->1. Its five Hamiltonian paths are

    0123, 0312, 1203, 2031, 3120.

The last four are cuts of the directed four-cycle 03120;
the remaining path is 0123. Exhaustion follows by considering the starting
vertex: from 0 there are 0123 and 0312; from 1 only 1203; from 2 only 2031;
from 3 only 3120. Thus H(Q)=5.

The reflection sigma(v)=3-v is an anti-automorphism. Q is rigid: its two
score classes each have two vertices, and the internal arc prevents an
automorphism from exchanging either pair. Its unique anti-automorphism
therefore gives the path involution R(P)=sigma(reverse(P)). The three
fixed paths are 0123,2031,3120; the other two are exchanged.

Write Q^r for the ordinal sum of r copies: every vertex in an earlier
block beats every vertex in a later block. Every Hamiltonian path visits
the blocks in order, each as one consecutive run. Consequently

    Ham(Q^r) <-> Ham(Q)^r,        H(Q^r)=5^r.                 (1)

The empty tournament has its one empty path. Reflection of all 4r labels
is the unique anti-automorphism for r>0: strong components are the rigid
Q blocks, whose total order forces their reversal. The induced path action
is `(P_1,...,P_r) -> (R(P_r),...,R(P_1))`. Hence

    blue(Q^(2s))=5^s,     blue(Q^(2s+1))=3*5^s.              (2)

These are products of actual paths, not class-count identities.

## 2. The exact square and the two different deletions

Let A=Q^2, X=Ham(A), and p0=01234567. Then |X|=25, R(p0)=p0,
and the bijection

    X x X -> Ham(A ordinal-sum A),
    (P,K) -> P followed by the shifted path R(K)              (3)

turns positional reversal into `(P,K)->(K,P)`. The 625 paths have exactly
25 blue paths on the diagonal. Therefore

    |X^2 minus diagonal| = 625-25 = 600,
    |(X minus {p0})^2|  = 625-2*25+1 = 576.                  (4)

The second set still contains 24 diagonal/blue paths and 552 black paths.
Their quotients under swapping have respectively 300 black orbits in the
full square, and 276 black orbits plus 24 fixed paths in the pointed
subsquare. Thus 576 is not the black count, either raw or after reversal.

The source, map, and preserved predicate in (3) are literal: the full two
path words reconstruct every arc of the Hamiltonian path, and equality of
the two coordinates is exactly blue status. Keeping only |X| destroys the
words and edge incidence. Keeping only endpoints also loses blue status.
Already in this Q^4 source, the standard path 0,...,15 is blue, whereas
`(0,1,2,3,4,7,5,6,8,9,10,11,12,13,14,15)` is black; both have endpoints
(0,15). Their four blocks are legal Q paths, but the latter fails the
reverse/reflection test. This is independently checked in the receipt.

## 3. What 776 can count, and what this does not force

The disjoint union

    Ham(Q^4) disjoint-union Ham(Q^3)
      disjoint-union Ham(Q^2) disjoint-union Ham(Q^0)

has 625+125+25+1=776 elements, equivalently words in a five-letter
alphabet of lengths 4,3,2,0. This is an explicit realization. The choice
to omit length one is additional data; the complete length-0-through-4
union has 781 elements. Neither a square nor a diagonal operation forces
that omission. No recurrence, single-tournament path count, divisor
classification, or Collatz descent is inferred from this disjoint union.

## 4. Recovering the balanced-family deletion fibres

In the inherited balanced family, X_k is the subsets of a (k-1)-element
set together with one exceptional state infinity. The one-step map forgets
the last subset bit and fixes infinity; it has two lifts per finite state
and one for infinity. On X_(k+1)^2 -> X_k^2 the fibres are therefore 4,2,1,
according to whether zero, one, or two coordinates are infinity. This is
realized by an actual two-vertex deletion, not an arbitrary set surjection.

With q=|X_k| and H=q^2,

    q'=2q-1,       H'=4H-4q+1,
    (H'-q')=4(H-q)-2q+2.                                  (5)

A black parent has only black lifts. A finite diagonal parent has two
blue and two black lifts; the exceptional diagonal parent has one blue
lift. This accounts for every correction term. The source state infinity
is the necessary sidecar; “all fibres are four” already fails at the
standard path. Equations (4) and (5) are different operations on different
sources, even though both involve squares and a distinguished state.

## Verification

From the repository root:

    python3 05-knowledge/results/trijection_reset_20261020/square_counts.py
    python3 -O 05-knowledge/results/trijection_reset_20261020/square_counts.py

The checker compares direct legal-path DFS in all five full tournaments
Q^r, 0<=r<=4, with independently formed product words, verifies the
anti-automorphism and every blue count, and checks the complete square map,
both deletion sets, and the selected disjoint union. It imports no producer
code, uses explicit exceptions, and has no symmetry/pruning filters.

An independent agent checked the written proof using all 24 permutations
and a separately hand-listed arc set. A prepublication draft had accidentally
listed path words from a different labelling of the strong four-tournament;
the direct legal-arc check caught and repaired that list before publication.
The corrected words and symmetry are now explicit in section 1.

## Incoming comparison: another literal source of 576

Concurrent mainline work at commit70cb917742 supplies
[a different24-by24 product](../charge_palettes_boundary_20261020g/palettes.md).
For a fixed allowed output colour of its gate, the remaining graph is two
triangles joined by one edge. Its three-colour assignments are decoded by
3 choices for one bridge endpoint,2 for the other, and2 orders within each
triangle:24 assignments. Two copies with only the output in common give
576 independent extension pairs. This source retains a named shared colour;
our path square retains full words and positional reversal. Their products
have the same size, but no colouring-to-Hamiltonian-path map is supplied.
In particular a24-element set cannot carry an involution with just one fixed
point, so the odd pointed-set involution from the balanced-family model
cannot be transferred to that fibre without changing its carrier.
