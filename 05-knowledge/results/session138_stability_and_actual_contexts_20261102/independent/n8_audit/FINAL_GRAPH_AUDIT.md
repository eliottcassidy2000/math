# Independent order-eight graph certificate audit: final no-wedge form

**Verdict: ACCEPT (complete finite certificate and full-cube consequence).**
This audits the graph residual and coverage, not an unreviewed translation
from tournament cycle sums to that residual.

For an undirected graph G on eight labelled vertices, let m be its number
of edges, m_j its number of j-edge matchings, w its number of unordered pairs
of incident distinct edges, and p21 the number of a wedge plus a disjoint
edge (not required to be induced). Then the complete certificate is

    41m - 32m_2 + 8m_3 - 4m_4 + 8p21 >= 0.

Equality holds exactly for the empty graph and the perfect matching.

## Independent orbit construction and coverage

`graph8_no_wedge.cpp` imports no producer code, data or representatives.
It constructs the complete seven-vertex graph orbits by breadth-first search
under the six adjacent vertex transpositions. These generate S7, so BFS
components are exactly isomorphism orbits. An exhaustive scan of the 2^21
labelled graph masks starts a new orbit only at an unseen mask. The result is
1044 disjoint orbits covering all2097152 masks. Every orbit size divides5040;
the receipt retains the entire orbit-size histogram.

For each representative, all128 possible neighborhoods of vertex7 are
adjoined. This covers all eight-vertex graphs up to isomorphism: delete any
fixed vertex of an arbitrary graph, relabel its seven-vertex card to its
representative, and extend that relabeling by fixing the deleted vertex.
Its seven neighbors become one of the128 tested subsets. The residual is
isomorphism invariant. Thus no full 2^28 scan is required.

## Independent counting method

The producer enumerates fixed matching patterns and sums wedges after each
edge deletion. This auditor instead computes the entire matching polynomial
of every vertex subset by the recurrence

    M_S(t) = M_(S-v)(t)
             + t sum_(u adjacent v, u in S-v) M_(S-v-u)(t),

where v is the least vertex of S. For p21 it first chooses the wedge center
and its unordered neighbor pair, then counts edges in the complementary
five-vertex induced graph using the same subset table. These are independent
counting directions. The identity m_2+w=binom(m,2) is also checked for every
augmentation. All checks throw exceptions and remain active in optimized
builds. The counts and all intermediate residual terms are small integers;
undefined-behavior sanitization provides a separate execution control.

## Exact output

- Seven-vertex representatives:1044.
- Seven-vertex labelled coverage:2097152.
- Eight-vertex augmentations:133632.
- Minimum residual:0; maximum residual:10808.
- Zero augmentations:2.
- Ordered residual FNV-style digest:9912251846543590531.

The zero witnesses in lexicographic edge-mask convention are:

1. mask0, the empty graph;
2. mask134226448, edges(0,5),(1,4),(2,3),(6,7), with
   `(m,m_2,m_3,m_4,w,p21)=(4,6,4,1,0,0)`.

The second graph is therefore explicitly a perfect matching, not a guessed
interpretation of the zero count. Orbit completeness implies that all
labelled zeros are the empty graph and the105 perfect matchings.

## Multiaffine zero-set consequence

If a separately affine polynomial on q_e in[0,1/4] has these scaled residuals
as its endpoint values, tensor linear interpolation proves nonnegativity on
the whole cube. Its zeros are only the endpoint configurations above. Indeed,
a zero with one interior coordinate would force both endpoint configurations
of an adjacent cube edge to have zero residual. But no two zero graphs differ
by a single edge: the empty graph has0 edges and every perfect matching has4.
Thus every q_e is an endpoint, and the zero support is empty or a fair perfect
matching. This conclusion concerns the graph polynomial; any additional
nonnegative tournament terms still have their own equality conditions.

## Reproduction

    g++ -std=c++17 -O2 -Wall -Wextra -Werror graph8_no_wedge.cpp -o graph8_no_wedge
    ./graph8_no_wedge > no_wedge.normal.json
    g++ -std=c++17 -O1 -g -fsanitize=undefined -fno-sanitize-recover=undefined -Wall -Wextra -Werror graph8_no_wedge.cpp -o graph8_no_wedge_ubsan
    ./graph8_no_wedge_ubsan > no_wedge.ubsan.json
    cmp no_wedge.normal.json no_wedge.ubsan.json

No Git operation or external-state mutation was performed.
