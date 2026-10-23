# Three incidence colors force the order-eight reflected wiring

October 23, 2026. **FINITE-EXACT:** the complete 6,880-class quotient-alphabet
search, combined with the proved reflection-marginal forest theorem and
the exhaustive 2^21-presentation native ledger. **PROVED:** monotonicity of
forest failure under further alphabet quotienting. The minimum below is
for FOREST CERTIFICATION within quotients of one fixed nine-symbol
alphabet, not impossibility for all binary signatures or for models using
additional capacities.

## A simple successful ternary rule

For each actual arc u->v, color it by

    chi(u->v) = (c_v+c_uv) mod3,

where c_v is the cyclic-triangle count at its head and c_uv the count through
its unordered pair. At each root u retain only the three numbers of its
outgoing arcs of colors 0,1,2. Take the multiset of those root triples.
Call this class invariant K(T). It retains no vertex labels.

This is a quotient of E_33. If a=c_v mod3 and b=c_uv mod3, its literal color
array is:

| a \ b | 0 | 1 | 2 |
|---|---:|---:|---:|
| 0 | 0 | 1 | 2 |
| 1 | 1 | 2 | 0 |
| 2 | 2 | 0 | 1 |

The restricted-growth partition word is `012120201`. Exact row lists and
color-count triples carry the same information. The degree of each root is
the sum of its three counts, so no separate score tag is needed.

Across all 6,880 ordinary order-eight classes, K has 6,864 distinct labels.
Its joint ordered pair (K(T),K(T^op)) separates all classes, and its unordered
pair separates the 3,528 converse-merged nodes. More importantly, the actual
separate K-margins, together with AUTHENTICATED converse-reflection equality,
force every ordinary port cell by the forest criterion in
`../CONVERSE_MARGINAL_REPAIR.md`. This does not infer joint margins from
paired-label injectivity alone.

### A genuinely stronger use of the forest criterion

The signature-incidence graph has 13,728 vertices, 6,880 class edges, and
6,848 components:

* 6,820 isolated edges;
* 24 two-edge stars;
* four three-edge paths.

There are four classes whose own K-label and converse K-label are both
non-singletons. They are the middle edges of the four three-edge paths.
Thus the easier singleton-or-converse-singleton test fails, but forest
leaf elimination still recovers them after the endpoint edges. This is a
concrete improvement over the nine-symbol star-only repair.

The native order-eight table provides feasibility and the exact source,
endpoint, color and reflected-category identifications. Its 956,339
nonzero class-port entries and all colored merged multiplicities are forced
in this reflection-constrained model. No claim is made that K(T) alone
computes a partner without an actual source path.

## Exhaustive bounded search and its stopping boundary

Number the nine E_33 symbols by z=3a+b. An alphabet quotient is a set
partition of these nine symbols; color names themselves are irrelevant.
There are Bell(9)=21,147 such partitions. We start with the nine-symbol
identity partition, merge one pair of blocks at a time, and retain every
partition whose global signature-incidence multigraph is a forest.

| Colors | All partitions | Explicitly tested | Already excluded | Forest successes |
|---|---:|---:|---:|---:|
| 8 | 36 | 36 | 0 | 35 |
| 7 | 462 | 462 | 0 | 414 |
| 6 | 2,646 | 2,636 | 10 | 1,981 |
| 5 | 6,951 | 6,800 | 151 | 3,599 |
| 4 | 7,770 | 7,163 | 607 | 1,886 |
| 3 | 3,025 | 2,376 | 649 | 117 |
| 2 | 255 | 129 | 126 | 0 |

All excluded partitions have an ancestor already certified nonforest.
The search explicitly tests 19,602 partitions and stops at two colors,
where no success remains. One color is then excluded by the same monotonic
argument. In particular three is the minimum FOREST-CERTIFIED alphabet
among all quotients of E_33. The sum rule above is selected for its simple
formula, not because it has the largest number of ordinary class labels.
The complete successful and failed tested partitions remain in
`quotient_search.tsv`.

### Why pruning does not miss a successful quotient

A further symbol quotient induces a function on outgoing-root signatures,
and therefore identifies vertices within each shore of the signature-
incidence graph. It retains every distinct ordinary-class edge. If B is
the oriented incidence matrix before quotienting, the new matrix is RB,
where R merges vertex rows. Any nonzero old dependency Bx=0 remains a
nonzero dependency (RB)x=0. A multigraph is a forest exactly when its
oriented incidence columns are independent. Hence once a partition fails
the forest test, no further coarsening can restore it.

Conversely, every successful k-block partition has only successful
refinements. Some chain of single block merges from the identity therefore
reaches it through the retained parents. This proves the downward search
covers every possible success, including the partitions pruned without
an explicit full-class scan. Parallel edges are kept; deleting them would
invalidate the argument and the minimum claim.

## What remains unproved by this search

A cycle in the GLOBAL signature-incidence graph is not necessarily a
feasible port-table counterfeit. Endpoint exclusions, nonnegativity,
column capacities, currents or other additional constraints may still
force wiring. Therefore the zero binary-forest count does not establish
that every binary quotient fails the full constrained wiring problem.
Actual nonnegative reflected switches are a separate required certificate.
It also does not rule out binary signatures outside quotients of E_33.

The nine-symbol modular minimum in `CONVERSE_SIGNATURES.md` concerns only
product residues E_pq. The present three-symbol palette is a non-product
quotient and is a different, explicitly enlarged menu. Neither claim
contradicts the twelve-symbol one-direction ordinary-class minimum.

The all-size aligned decoder in `THREE_COLOR_DECODING.md` supplies a mechanism
for this plus palette, but requires row and edge associations that K's
outer multiset forgets. The order-fifteen circulant hostile still defeats
K and its converse: the rule is not an all-order class decoder.

## Reproduction and independent controls

    python -B incidence/quotient_full_check.py --classes wiring/classes_input_n8.tsv
    python -O -B incidence/quotient_full_check.py --classes wiring/classes_input_n8.tsv

The consumer regenerates the prepared class data, compiles the bundled
standard C++17 helper locally, and compares its exact transcript against
the frozen search. It verifies every candidate-child set, all level counts,
and the first zero-success stopping condition. It independently checks the
selected ternary forest and all 36 immediate eight-symbol quotients in
Python. The helper uses exact byte-string equality, not bare hash values,
and no solver. Failed partitions stop when a cycle appears; their recorded
`signatures_seen` is only the prefix count, not a full-class count. Successful
partitions process all 6,880 class edges.

`--read-existing` checks the saved transcript and Python controls without
rerunning compilation; it is not a substitute for the full replay commands.
The separate independent auditor checks the sum rule and all 255 binary
quotients without this failed-parent pruning.
