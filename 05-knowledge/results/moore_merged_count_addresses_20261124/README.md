# Moore graphs: exact heights, marked merged lines, and a degree-57 fork

2026-11-24. Extends the [native spine ledger](../moore57_graph_spine_tiles_20261124/),
[earlier exact counts](../moore_chanpak_counts_20261124/),
[compact chart](../moore57_compact_blue_routes_20261124/), and
[Hoffman–Singleton symmetry obstruction](../five_thresholds_hoffman_singleton_20261124/).
Moore-graph bases: [Hafner](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v11i1r77).
A merged node means a tournament isomorphism class further identified
with its converse. A marked line here complements the nonspine tiles;
a physical single-arc flip is a different operation.

## 1. Exact native classifications

For an ordering along either an actual Hamilton spine of the graph
or a Hamilton spine of its complement, start with forward orientation
and reverse the nonspine edges of the Moore graph. In the former case,
recover the graph by putting the spine back; in the latter all actual
edges are already marked. The mate complements all optional tiles.

| Graph and spine | Distinct marked layouts | Merged source nodes | Merged mate nodes | Merged lines | Source H | Mate H |
|---|---:|---:|---:|---:|---|---|
| C5, graph | 1, blue | 1 | 1 | 1 | 9 | 11 |
| C5, complement | 1, blue | 1 | 1 | same line | 11 | 9 |
| Petersen, graph | 2, black | 1 | 1 | 1 | 5213 | 6961 |
| Petersen, complement | 796, black | 366 | 398 | 398 | 333 values, 4927–13837 | 362 values, 5571–15571 |

Source and mate sets are disjoint in each row. The last row is an
exact total of labels/values within this particular marked family,
not an interval coverage assertion. There are no loops. Its unique
marked black line of height difference 2 is **5625 ↔ 5627**,
from the spine 0,2,4,1,9,8,7,3,5,6 in the standard pentagon
0..4, pentagram 5..9 presentation (forward spine, reverse actual
edges). None of the 45 physical single arc flips of the marked
source yields the mate up to converse or relabeling.
The equal-distance analogy to Chan–Pak must therefore retain which
operation was used. The source and mate *families*, all marked black,
each contain one underlying self-converse node (height 11793 and 15549,
respectively); their full blue supports after rerooting are 165 and
181 symmetric tilings. Thus the native color is not the color of
every other tiling in its node.

The exact replay enumerates all permutations, deduplicates native
tile masks, and uses equitable refinement followed by complete
individualization. Each key is the smaller of the two canonical
tournament keys. Hamilton counts are independently obtained by
endpoint subset DP. A vertex ordering satisfying the spine condition
and its backward-edge mask specify the marked line, not just its
height. This is line incidence, not a claim of shortest distances
for the 398 marked lines.

Cyclic C5 (15 paths) and compact Petersen (8265 paths) remain
different decoders. The full compact/native 50-site heights and
56-site twin-blowup height are not evaluated here. The previous
block-contiguous contributions to the 50-site chart remain lower
bounds, not full totals. The 56-site chart is on a blue line and
its even self-converse height is 1 mod 4 by the odd-cycle collection
formula; the native even charts cannot themselves be marked blue.

## 2. Exact count-addressed encoding

A different encoding makes the entire labelled Moore graph recoverable
from a single Hamilton-path count. It is larger than the geometry-only
chart; the tradeoff is exact arithmetic inversion and a reproducible
placement. Choose root c, order its k neighbors, and label their
k branches (each of size k−1). Between every pair of branches there
is a perfect matching: every branch vertex has precisely one
neighbor in each other branch by degree and absence of squares.
Normalize all matchings from branch 0 to the same labels. The
remaining C(k−1,2) permutations are the complete payload; the
root, spokes, and branch stars are implicit. Encode by Lehmer ranks
in base (k−1)!, as a W-bit integer C. If an unlabelled invariant is
wanted, take the minimum C over the finite choices of root,
neighbor order and first-branch labels. Reconstructing from C
and checking the Moore condition decides whether a payload is valid.
The replay gives specified labelled examples, not a minimization.

Take a bit word w with 1010101010, the W-bit payload, then a zero.
Let L=|w|, K=L+7, B the integer value of w. On a transitive chain X
of length L, add apex a with a→i iff bit i is 1. Call the count h(w).
Counting by the last vertex before a and first vertex after it gives

    h(w) = w_1 + 1−w_L + #adjacent 01 + Σ_(w_i=0) value(w_1…w_(i−1)).

Append K ones and a final zero. In the new tournament the full
displayed chain then a is a Hamilton spine. Its exact height is
\(\mathbf{2^K(B+1)+h(w)}\). Since 1≤h(w)<2^L (partition choices
before a give an upper bound 2^L, and the count is odd), one
can recover B by quotient and hence recover the entire graph.
The chart is rigid: a is the unique vertex whose deletion leaves
a transitive tournament, by two triangles intersecting only at a;
the remaining order has no automorphisms. Reversal/converse would
reverse-complement the contact word, but the sentinel 101 contradicts
the first three bits 100 of the reverse complement (K≥2).
Thus the ordinary class has H tilings, the merged node 2H, all black.

### Exact marked mate in linear-many transfer steps

After flipping the nonspine tiles, the chain arrows descend except
for consecutive ascending steps. The apex still has final bit 0
and its other contacts are inverted. Divide a Hamilton path into
two lanes on each side of the apex. Every ascending run in either
lane is an interval of consecutive integers; the runs occur in
decreasing interval order. Scan the integer sites increasingly,
assigning each to one lane. Either extend that lane's highest run
(if the preceding integer is in the same lane), or create a new
highest singleton. Keep its run count clipped at 2, whether its
highest run is a singleton, and a pending-invalid flag if a new
singleton lies immediately above a singleton, which would make
a forbidden descent by one. Extending the highest run repairs
that flag; starting yet another run cannot repair it. Also retain
the endpoint contact bit: final vertex of the lowest run in the
left lane, first vertex of the highest in the right lane. This
is a bounded-state exact DP in O(L+K) arithmetic operations,
not a unit-cost assertion for arbitrarily long integers.
The replay checks both height formulas on 255 small tournaments
against subset DP.

| Input | Payload width W | Exact counted black line order | H and mate |
|---|---:|---:|---|
| C5 | 0 | 31 | 357827701 ↔ 1353711331 |
| Petersen | 1 | 33 | 1431832353 ↔ 5482528947 |
| Hoffman–Singleton | 143 | 317 | full integers in word_code.out |
| conditional degree 57 | 382925 | 765881 | formulas above, no witness graph |

The computed mate heights differ in all known examples. Neither
this count-address representation nor native height is an invariant
of the graph *across* encodings.

## 3. A graph-preserving physical +2 fork

The same code admits **one literal arc flip that adds 2 for every
payload**. Insert a 6-site tournament Y after the prefix X
and before the guard run, as an ordinal block relative to them.
Keep the same apex and give it outcontacts {0,1,4} in Y.
In the lexicographic upper-triangle arc convention, Y has
backward mask 177 or 181, a flip of (0,3). The common Y spine
is 0,2,1,3,4,5. The following are checked by endpoint subset DP:

* H(Y)=29; Q(Y)=320 counts ordered partitions into two directed
  paths, permitting empty paths.
* R(Y)=158 counts the same covers with the right path empty or
  beginning at an apex outcontact.
* J(Y)=162 counts the covers with the left path nonempty and
  ending at an apex incontact.
* H(Y+a) is 111 or 113, respectively.

For the entire prefix M=X→Y, R(M)=320B+158: if the right
path touches X it supplies B eligible splits and both Y lanes
are free, otherwise all of X lies left of a. In M+a, if X
lies entirely to the left then w_L=0 allows all H(Y+a) paths.
If the right path touches X and the left touches no Y, the
contribution is 29(h(w)−1); otherwise the left ends at a Y
incontact and the contribution is 162B. Now call the final
zero site z. It either ends the full path, in which case all
guards follow a and the count is H(M+a), or it immediately
precedes a. In the latter case guards are freely divided into
two sets around z,a and contribute 2^K R(M). Thus

    H_− = 2^K(320B+158) + 29(h(w)−1) + 162B + 111,
    H_+ = H_− + 2.

For fixed L,K, divide H−158·2^K−111 by 320·2^K+162
and take the quotient. The remainder is at most 29·2^L+2,
so the quotient is B at either end; the residual selects the
endpoint. The arc changes the count, **not the decoded graph**.
Both endpoints are explicit distinct merged nodes. They are rigid
and non-self-converse: the apex uniquely maximizes triangle
participation, since its triangles through z number its outdegree
and the sentinel supplies at least 15 more; every other site
participates in at most the apex outdegree plus 10 (guard width
makes outdegree ≥ indegree, and the only nonapex triangles lie
in Y). On deletion of a the unique nontrivial SCC is the rigid
Y; the singleton strings on its two sides have different lengths.
So the surrounding ordinary classes each have H tilings and the
merged nodes each have 2H, all black. This suffix changes the
simple chart's marked mate; no mate count for the six-site fork
is inferred from §2.

There are \(\gg\sqrt X\) distinct +2 seams with upper value ≤X
when allowing arbitrary binary payloads (L chosen with
2^(2L) a small fixed fraction of X). This goes beyond a finite
list of fixed-word rays because the contacts grow with X.
It still does not prove density one or an exact selector for all
odd targets. The actual degree-7 quotient produces a 323-site
pair with huge exact totals in word_code.out; degree 57 would
give order 765887, conditional on a valid payload.

## 4. A degree-57 encoding constraint before doing the tournament count

Let g be a nontrivial involution of a hypothetical (3250,57,0,1)
Moore graph. Write f for the number of fixed vertices, s for the
sum of their induced degrees, and A for the number of moved
vertices adjacent to their mate. As in the previous symmetry
session, those A vertices have no fixed neighbors, and between
the A pairs there are no edges (a square would result). Every
other moved vertex has exactly one fixed neighbor. Thus

    A=3250−58f+s, and A≤57f−s, in particular f≥30.

The fixed graph has diameter at most 2 (unique common neighbors
of nonadjacent fixed pairs are fixed), no triangles/squares,
and odd degrees. Such a graph is a star or a regular Moore
graph: its adjacency matrix D satisfies
D²+D=J+diag(deg−1); commutation with D forces nonadjacent
vertices to have equal degrees. If irregular, all different
degree groups see one another completely, giving two independent
groups by triangle-freeness and a singleton group by square-freeness.
By the Moore options, the only proper regular fixed graph big
enough is Hoffman–Singleton itself (50 vertices).

If it is a star, root at the center. Its f−1 fixed leaves are
some of the first-layer branch heads. Every second-layer
vertex in a moving branch belongs to an adjacent moved pair:
it has no fixed neighbor, since different branches cannot meet
the center through a second fixed leaf. There are perfect
matchings between every pair of moving branches, but A has
no edges between different adjacent pairs. Thus at most two
branches move: f=56 or 58. Finally use adjacency eigenvalues
57,7,−8. If x is the trace of g on the 7-eigenspace, then
trace(g)=f and trace(D_full g)=A give

    15x = 3185−50f+s.

At a star s=2(f−1), so f≡1 mod 5 and f=56. The two cases are:

| fixed graph | f | moved adjacent to mate | moved nonadjacent to mate |
|---|---:|---:|---:|
| 56-site star | 56 | 112 | 3082 |
| Hoffman–Singleton | 50 | 700 | 2500 |

### In the star case the raw substrate is a noncommuting K56 matching deck

Label the two moving width-56 branches U,V via their common
matching, and for each fixed branch i label its sites by
their partner in U. The involution in each fixed branch gives
a fixed-point-free permutation σ_i on the 56 labels. V_x
meets label σ_i(x) in that branch. For x≠y the pair U_x,V_y
has a unique common neighbor, so exactly one of the **55**
permutations sends x to y. The σ_i constitute a 1-factorization
of K56. Their largest pairwise commuting subfamily has size
at most **7**: a commuting set generates an elementary abelian
2-group whose orbit sizes are powers of 2, one of them ≤8
since 16 does not divide 56. Each generator takes a fixed
point to a distinct other site, so the smallest orbit has
size at least one more than the number of generators.
Thus no commuting translation model supplies the full deck.

For fixed branches i,j the unknown matching π_ij must
intertwine their involutions. For each fixed i and label y,
the exact two-step constraint from U into i says that
{y,σ_i(y)} together with the 54 values π_ij(y) (j≠i)
list all 56 labels once. This narrows each general permutation
to 2^28·28! equivariant choices before imposing the row and
remaining two-step constraints. The ambient rooted matching
conditions are exact: for every branch pair i,j and label x,
the direct image π_ij(x) and the through-k images
π_kj(π_ik(x)) over k different from i,j exhaust branch j.
These conditions ensure the Moore property on top of the
root and branch stars. The replay recovers the five matching
rows on six sites in the known Hoffman–Singleton star model.
The next constructive test is to fit the simultaneous row
conditions to a noncommuting K56 factorization, or treat the
alternative fixed Hoffman–Singleton directly. No degree-57
graph or native full Hamilton count is being asserted.

## Reproduce

Run python3 word_code.py (imports gadget.py); compile native.cpp
using C++17 with optimization and run it. The receipts are
word_code.out and native.out.
