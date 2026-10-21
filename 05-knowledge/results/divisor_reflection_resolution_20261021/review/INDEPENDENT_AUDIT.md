# Independent audit of the reflection resolution

Verdict: **PASS** for the exact scoped elementary claims in REFLECTION_RESOLUTION.md, SHA-256 eecbe71f9cd36aba84ed3836d3f24c81b62ee50b7303eb430b4db6ac7f5da80c. The producer files were not edited. This review does not claim literature novelty, unmarked tournament-class correspondence, a divisibility embedding or a Collatz result.

## Proof review

1. The exponent chart is injective after deleting the two endpoints: for g=0 the b=c=0 and b=c=1 ranges overlap only at f=0, precisely those endpoints; g=1 and g=-1 separately recover a. Q fixes b=c, giving2k surviving states. CQ fixes exactly the two a=k/2,b+c=1 states for even k and none for odd k. C fixes none.
2. The staircase graph/tournament dictionary is pairwise bijective. Reversal of vertices alone reverses the fixed directed path, so the converse in R is necessary. Its nonconsecutive pair image is exactly the reversal of the corresponding graph edge.
3. Every vertex involution on at least two vertices has a fixed unordered pair: a transposition orbit, or any pair if the involution is the identity. On this pair complement-reflection has no fixed edge bit. An equivariant map from a source CQ-fixed point would therefore have no possible image. This excludes many-to-one maps too. Interchanging source axis reflections still sends a source involution with fixed points to the forbidden product.
4. Splitting the two central fibres, preserving epsilon under C and reversing it under Q, yields commuting involutions and an equivariant projection to each original divisor. The printed twelve-row map was checked directly, not merely a copied expected map. Its two operations match graph complement and reversal; actual directed five-vertex tournaments have four blue states and six complementary pairs, exactly two blue pairs.
5. Minimality is source-fibre-local: CQ acts freely in each fixed source fibre, so each such nonempty fibre is even and at least two. The other eight source fibres are nonempty, giving12. A fixed Hamiltonian-path chart on at most four tournament vertices has at most8 states; five vertices are attained.
6. The resolved all-k action has k orbits with stabilizer<Q> and floor(k/2)+1 free orbits. On m-vertex graphs, reversal fixes floor(m/2) edge coordinates and has t=(binomial(m,2)+floor(m/2))/2 edge orbits. Neither complement nor complement-reversal has fixed graphs. Orbit-stabilizer therefore gives precisely the two capacities stated in Theorem2. Matching source and target orbits with equal stabilizers proves the sufficient direction; injections cannot merge source orbits, proving necessity. The construction is correctly described as noncanonical.

## Independent exact implementation

The accompanying independent_reflections.py uses sets of unordered edges and sets of directed arcs. It imports no producer implementation. It reads the literal Markdown table, translates each printed divisor label to its exponent/mark tuple, and tests the table's actual operations and colours. Further controls:

- original and resolved divisor actions and injective coordinate charts for every k=2,...,50;
- every graph on m=2,...,6 vertices and every vertex involution in those orders;
- explicit orbit-matched injective equivariant maps for each fitting pair k=2,...,12 and m=2,...,5;
- retained original-source projection and chart-size lower bounds.

Normal and optimized Python produce identical62,357-check receipts, SHA-25650299a869a64c001e8d4a20b1406002f6cf3f5e8c56bac482eca63d280a0c03d. Assertions are not used as load-bearing checks. Run from the published packet root:

    python3 -B review/independent_reflections.py
    python3 -B -O review/independent_reflections.py

Compare stdout with review/independent.normal.json. Uniform statements are justified by the proof above, not extrapolated from these finite controls.
