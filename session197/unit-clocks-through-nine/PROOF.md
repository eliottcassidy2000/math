# Complete rooted unit ternary clocks through order nine

## Scope and result

A rooted tournament (T,r) is a unit ternary clone source when replacing r by a transitive module of m vertices gives

`H_m=A*3^m+beta`, with A>0 and beta in{-1,+1},

for every m>=1, or even just for all sufficiently large m. This packet classifies these sources only through order9. It is not an all-order amplitude classification and does not turn the observed residue pattern into a theorem for larger sources.

The complete amplitude/sign list through order9 is

- beta=-1: A in{60,100,156,180};
- beta=+1: A in{72,120,200}.

At order exactly9 there are16 rooted isomorphism classes, including the marked root as part of the isomorphism data. Their counts by profile are

- (60,-1):4;
- (100,-1):2;
- (156,-1):2;
- (180,-1):4;
- (72,+1):1;
- (120,+1):2;
- (200,+1):1.

Twelve are strong. The four remaining classes are neutral one-vertex paddings of the eight-vertex negative60 sources. The file catalog.json provides every actual outgoing-row array, with its root always labeled8, and its complete run vector. These are rooted source counts, not counts of all tournaments of a given order.

In particular a positive180 clock cannot have a source of order at most9. Larger or differently constructed sources remain open here.

## Exact reduction to matching number two

For a rooted tournament, let a_j count words using every nonroot vertex once and j indistinguishable root symbols, with no two root symbols consecutive. Cutting a Hamilton path into its transitive-clone runs gives the exact formula

`H_m=sum_j a_j j!S(m,j)`.

Let O and I be the root's outneighbor and inneighbor sets. In the bipartite graph of arcs O→I, let nu be the maximum matching size. Then the largest nonzero run index is exactly nu+1.

For the upper bound, each internal fragment between roots starts in O and ends in I and therefore contains a crossing edge. Different fragments give disjoint crossing edges.

For the lower bound, take a maximum matching of positive size and use its edges as internal fragments. Extend one matched edge to a directed path covering every vertex not used in the other matched edges, by successive vertex insertion; insertion preserves the old path order. Cut that extended path before its first O vertex and after its last I vertex. These positions occur in the correct order because the original matched edge remains an ordered subsequence. The middle fragment runs from O to I. Its discarded prefix contains only I vertices and fits before the first root; its discarded suffix contains only O vertices and fits after the last root. The other matching edges remain internal fragments. This constructs nu+1 root occurrences. When nu=0, an ordinary Hamilton path gives a_1>0.

Since every coefficient is nonnegative, an eventual positive pure ternary ray must have largest run index3. Thus nu=2, and comparison of the distinct bases1,2,3 gives the exact finite conditions

`a3=A`, `a2=3A`, `a1=3A+beta`.

Conversely these conditions with nu=2 prove the identity for every m>=1. This also explains why an eventual source is not missed by the finite conditions. No unsupported scalar interpolation is used.

## Complete physical generation at order nine

Fix the root as vertex8. The matching condition requires at least two outneighbors and two inneighbors. Taking graph converse preserves every clone count and exchanges their sizes, so it suffices to use the degree splits

`(|O|,|I|)=(2,6),(3,5),(4,4)`.

For a split (r,c), represent the crossing graph by r bit rows of width c. Every crossing graph can first have its rows sorted. Generate every nondecreasing r-tuple of bit rows, retain exactly matching number2, and minimize over all c! permutations of its columns, sorting its rows again. This gives a complete set of bipartite isomorphism representatives, with no probabilistic or heuristic pruning.

The resulting numbers of crossing patterns are respectively42,38,37. For each representative, enumerate every orientation inside O and inside I. All other arcs are fixed by the root partition and the crossing pattern. Therefore every rooted tournament of order9 with nu=2 is represented, up to rooted relabeling and, for the unbalanced splits, converse.

The exact census sizes are

- split2+6:42 patterns,2,752,512 orientations tested,0 unit hits;
- split3+5:38 patterns,311,296 orientations tested,40 unit hits;
- split4+4:37 patterns,151,552 orientations tested,72 unit hits.

The3,215,360 tested graphs are presentations, not isomorphism classes. The unit-hit counts include duplicates under remaining symmetries. Every graph is counted by an integer subset recurrence for its full root-skeleton coefficients up to index4; the proved matching bound ensures that all larger indices vanish.

Add the converses of the unbalanced hits, then canonicalize each rooted graph exactly. The canonicalization partitions nonroot vertices by their incidence to the root and their whole-graph outdegree, and minimizes over every permutation within those classes, keeping the root fixed. Every rooted isomorphism preserves those two invariants, so this is exact rather than a hash heuristic. It yields the16 graphs in catalog.json.

## Why the amplitude list covers all smaller sources too

Any rooted source of order below9 can be padded to order9 by universal source or sink vertices, without changing its marked root. The extra singleton strong components contribute a factor1 to every clone count. Thus its amplitude/sign pair must occur in the complete order9 catalog. Conversely every listed pair is realized by a graph in that catalog, so the displayed list is exact through order9.

The statement does not bound the order of arbitrary unit ternary sources. In particular it does not infer absence at order10 or above from this finite census. The complete classification of smaller-amplitude negative clocks from a separate strong-core argument is not used to enlarge this packet's scope.

## Independent checks and reproduction

Run `python3 verify.py` and `python3 -O verify.py`. A C++17 compiler is required. The Python verifier regenerates all crossing representatives, compiles the C++ census in an isolated temporary directory inside this packet, and runs the complete fixed-size enumeration. It checks exact case inventories and the complete canonical catalog, then recomputes every retained graph's run vector with a separate Python recurrence.

A further direct recursive word enumerator visits actual root skeletons rather than using subset dynamic programming. It verifies every coefficient of each of the16 retained sources through index4. Exact matching, strongness and rooted canonical forms are also checked. The receipt records hashes of the complete crossing-pattern input and census output, so normal and optimized runs can be compared byte for byte. An optional CLOCK_CPP_MODE=debug uses an unoptimized C++ build with standard-library assertions.

The result is a reproducible finite-order classification with actual witnesses. Its mathematical ingredients are the exact run-degree theorem, exhaustive crossing-pattern generation, complete side orientations, and exact rooted isomorphism reduction.
