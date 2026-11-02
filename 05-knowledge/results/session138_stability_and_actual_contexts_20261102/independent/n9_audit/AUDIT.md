# Independent order-nine graph certificate audit

**Verdict: ACCEPT.** For every undirected graph on nine vertices,

    81m - 64m_2 + 16m_3 - 4m_4 + 19w + 16p21 >= 0,

where m_j counts j-edge matchings, w counts pairs of incident distinct
edges, and p21 counts a wedge together with a vertex-disjoint edge, without
an induced-subgraph requirement. This note audits the graph polynomial and
its complete universe; the tournament cycle-ledger translation is separate.

## Complete orbit coverage and canonicalization

The independent C checker does not import producer representatives or code.
It constructs seven-vertex orbits by BFS under the six adjacent vertex
transpositions, which generate S7. The scan covers every one of2097152
labelled graphs and produces1044 disjoint orbits. Every one of their128
vertex augmentations is retained before canonical deduplication, for133632
complete eight-vertex candidates.

The auditor's eight-vertex canonical map differs from the producer's plain
degree-partition map. It starts with degrees and repeatedly refines vertex
colors by the tuple (old color, numbers of neighbors in each old color).
Color names are assigned in lexicographic signature order. This is a
label-invariant refining partition and terminates after at most eight color
classes. The auditor then minimizes the adjacency encoding over every
permutation within each final color cell, placed in color order.

This canonicalization is exact, not a heuristic: any graph isomorphism
preserves every refinement color and maps the permitted orderings bijectively;
conversely equal canonical adjacency strings exhibit an isomorphism. Sorting
and deduplicating the133632 strings therefore gives precisely12346
isomorphism classes of eight-vertex graphs.

The producer's simpler degree-cell enumeration is also mathematically valid:
degree partitions are invariant, its within-cell permutations are exhaustive,
and equality of resulting adjacency strings exhibits an isomorphism. A full
refinement partition is an independent implementation, not a requirement for
the producer's method.

Finally all256 neighborhoods of the ninth vertex are appended to every
retained eight-vertex representative, giving3160576 cases. This covers every
nine-vertex graph up to isomorphism: delete any fixed vertex, canonically
relabel its eight-vertex card, and extend the relabeling to the deleted
vertex. Its neighborhood becomes one of the256 tested masks. No2^28 or
2^36 visited array, and no order-ten run, is used.

## Independent efficient counts

For each base graph G on eight vertices the checker computes all induced
matching polynomials by the least-vertex recurrence

    M_S(t)=M_(S-v)(t)+t sum_(u adjacent v in S) M_(S-v-u)(t).

For an added vertex with neighborhood S, matching counts are then obtained
from the independent vertex-addition identity

    M_(G+x)(t)=M_G(t)+t sum_(v in S) M_(G-v)(t).

This avoids the producer's enumeration of fixed matching masks.

The p21 evaluator also counts in a different direction. Let W=w(G), let
c_v be the number of old wedges containing v, let e_vu be the number of
edges of G-{v,u}, and let l_v=sum_(u adjacent v) e_vu. Then

    p21(G+x)=p21(G)
       +sum_(v in S)(W-c_v+l_v)
       +sum_({v,u} subset S) e_vu.

The terms respectively count an old wedge and the new disjoint edge x-v;
a new wedge having x as a leaf; and a new wedge centered at x. These cases
are disjoint and exhaustive. The base p21 count itself is wedge-first,
using matching-table edge counts on the complementary vertex subset.

Every augmentation checks m_2+w=binom(m,2). In344 fixed or zero-witness
snapshots a separate full nine-vertex matching DP and direct wedge-first
complement count reproduce all six counts. Every check is active and throws
an error on failure. All count intermediates are safely within signed32-bit
range; the unsigned64-bit rolling digest intentionally wraps.

## Exact results

- Seven-vertex representatives:1044; labelled coverage:2097152.
- Eight-vertex augmentations:133632; exact orbit count:12346.
- Nine-vertex augmentations:3160576.
- Minimum residual:0; maximum residual:60372.
- Zero augmentations:14.
- Exact maximum required wedge coefficient:19/64 in the convention
  `base+64 alpha*w`; the tested coefficient is sharp for this graph bound.
- Full-DP cross-checks:344.
- Ordered residual digest:4534999060215151363.

The14 zero augmentations consist of exactly three graph isomorphism types:

1. the empty graph (one augmentation);
2. four disjoint edges plus one isolated vertex (three augmentations);
3. a three-vertex path plus three disjoint edges (ten augmentations).

The receipt retains all14 masks and degree sequences. The latter uniquely
identify these types here: all zero degrees; one zero and eight ones; or
eight ones and one two. The last type has counts
`(m,m_2,m_3,m_4,w,p21)=(5,9,7,2,1,3)` and gives base=-19, proving the
required coefficient19/64 is not merely a numerical upper bound.

**Equality caution:** unlike the order-eight certificate, the order-nine
zero graph types can differ by a single edge. For example, adding an edge
from the isolated vertex of4K2+K1 to one matching endpoint gives P3+3K2.
Consequently the multiaffine continuous graph polynomial has genuine zero
segments. It is invalid to infer only endpoint zeros. Any tournament equality
classification must use its retained positive triangle reserve and/or other
nonnegative terms; this graph certificate alone does not force that result.

## Guarded execution and reproducibility

`graph9_independent.c` built and ran successfully in both unoptimized and
-O1 builds. Their complete JSON outputs are byte-identical:
`normal.json` and `optimized.json`.

Successful compilation used a24MiB data cap and40MiB RSS setting. Measured
peak RSS was18432KiB for the unoptimized build and24416KiB for -O1.
Both executables ran under a hard40MiB address-space cap; measured peak RSS
was6132KiB and6144KiB respectively. Exact resource receipts are retained.
The successful runs exited0 and fit the requested40MB RAM bound.

The first UBSan compilation attempts under the same24MiB data cap failed
allocation at -O1 and one -O0 retry (measured RSS34460KiB and32936KiB).
Those failures remain recorded in `compile_ubsan.json` and
`compile_ubsan_o0.json`; they are not retroactively claimed as passes.

The compiler data cap was then increased to 64 MiB while retaining the
40 MiB runtime target. The -O1 UBSan build succeeded with compiler peak
RSS 50,468 KiB. Its complete bounded n9
run then exited0 with runtime peak RSS6164KiB, below40MB. Runtime limits
were128MiB address space (to load sanitizer libraries),24MiB data and40MiB
RSS. The UBSan output `ubsan.json` is byte-identical to both successful
ordinary receipts. No undefined-behavior diagnostic occurred; stderr contains
only the expected representative-count progress line. Recovery compilation
and terminal runtime resource receipts are retained. No order-ten run was
started.

Reproduce the successful executable with:

    gcc -std=c11 -O1 -Wall -Wextra -Werror -Wno-misleading-indentation graph9_independent.c -o graph9_independent_optimized
    ./graph9_independent_optimized > optimized.json

Use the retained `guarded_run.py` to enforce the same caps and generate a
resource receipt. No Git mutation or expanded-order run was performed.
