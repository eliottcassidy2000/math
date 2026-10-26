# Closed endpoint-four replacement for every order at least eight

2026-10-25. **PROVED by audited inheritance and complete finite
certificates:** for every tournament T of order n>=8 and every directed,
closed, blue Hamiltonian word P with endpoint outdegrees (n-5,4), there
are three directed black words Q1,Q2,Q3 at their actual successive
sources such that

```
K_Q3(K_Q2(K_Q1(T))) = K_P(T).
```

Equality holds on the original labelled vertices and every arc. There
is no target-strongness, mixedness, clone, module, disjoint-port,
isomorphism, or failed-guard filter. The threshold eight is sharp for
this unfiltered statement: a closed endpoint-four source at seven has
a pure-blue target and no black detour of any length.

This assembles three already specified ranges; it does not infer a
general theorem by extrapolating a finite census. The packet's original
[n>=12 frontier](PORT_FRONTIER.md) remains its own proved milestone.
Its old sentence leaving order eleven open is superseded by the
incoming theorem cited below. Full black-core connectivity remains
**OPEN**.

## 1. Exhaustive assembly and exact source universes

| Range | Status and mechanism | Full source scope |
|---|---|---|
| 8,9,10 | FINITE-EXACT: THM-4732 section 3, freshly rerun with a corrected decoder and independent literal verification | 640, 10240, 143360 unique marked sources, respectively; every one has a three-black certificate |
| 11 | FINITE-EXACT: incoming complete affine partition and seven exceptional-source repairs | All 224 nonempty dimension-14 source cells; 3,670,016 marked sources |
| 12,13,14 | PROVED by complete affine certificates and structural repairs | The complete base bank and all endpair residuals; the six failed old cells at 12/14 are repaired by source-level selectors |
| n>=15 | PROVED by the inherited structural/finite-affine bank and exhaustion of the remaining endpair family | Arbitrary three-port row I; all middle orientations and actual source arcs |

For the last two rows, the exact partition, thresholds and references
are the [n>=12 theorem](PORT_FRONTIER.md) and
[its independent review](PORT_FRONTIER_REVIEW.md). Its residual family
is exactly `I={h,n-3,n-2}` with closed middle; the h=2, h=3, h=4 and
h>=5 repairs exhaust it. No disjointness or prefix-guard condition is
introduced at the assembly step.

The order-eleven result is the incoming
[complete endpoint-four theorem](../session64_lca_reverse_valuation_20261026/session64-repairs/next11/README.md)
and its [independent audit](../session64_lca_reverse_valuation_20261026/session64-audit-eleven/AUDIT.md),
originally inspected at incoming commit
`9e79465c2cfae960c3b26762aa9f3e15340789ea` and then in the merged tree.
Its 212 base cells and 12 endpair cells are complete. The endpair cells
split into 25 nonempty, disjoint, exhaustive guarded leaves. All failure
loci are decoded exactly; their union is seven literal tournaments.
Seven independently replayed triples repair those sources. No theorem
whose structural threshold is twelve is simply assumed at eleven.
One repair really merges the ordered clone block (5,4) and later splits
it; full vertex rows are retained throughout.

For all rows, relabel the given P in path order to `(0,...,n-1)`.
Its source conditions are precisely the closed arc `(n-1)->0`, complete
positional anti-reversal, the directed base path, and
`N^-(0)={n-1} union I` for an arbitrary three-element
`I subset {2,...,n-2}`. Transport the obtained words back to the original
labels. This proves the all-word quantifier without quotienting the
source, intermediate tournaments, or target.

## 2. Finite-head validity gate, including the repaired decoder

The inherited source is THM-4732 section 3 in
`01-canon/theorems/THM-4732-black-detours-of-closed-blue-lines-order-9-census-order-10-detours-pure-blue-layers-and-the-klein-dictionary.md`.
Its k=4 raw receipts are
`04-computation/experiments/blue_line_black_detours_20261020_n8.out`,
`blue_line_black_detours_20261020_n9.out`, and
`blue_line_black_detours_20261020_n10.out`.
Inspection of the actual C++ iterator confirms that it includes every
closed standard-blue marked source in the requested degree layer.
The counts are recorded before the pure-blue-target test; pure-blue
cases would therefore be visible failures, not omitted sources.

Targeted correction review found that MISTAKE-646 preserves these
finite counts while correcting other endpoint thresholds and claims;
MISTAKE-648 repairs a pendant-node claim in another part of THM-4732.
Neither correction supplies an excuse to omit the pure-blue hostile at
seven, which is retained explicitly below.

The new source audit found a concrete issue in original line 61 of
`04-computation/experiments/blue_line_black_detours_20261020.cpp`:
the directed-path decoder could read an uninitialized second neighbor
before rejecting a premature endpoint. The [finite-head packet](finite_head/README.md)
records the original raw and LF hashes, the exact guard-order patch,
the preserved corrected copy, compiler command and minimal five-vertex
path-plus-triangle hostile. The original code's potential undefined
behavior was repaired before this threshold promotion; no false count
or false literal certificate was inferred from the bug itself.

The corrected full k=4 runs passed at 8,9,10 under both ordinary Python
and Python `-O`. An independent standard-library checker consumes every
certificate and recomputes source validity, actual-source directedness,
blackness, path-edge XOR, and the complete literal target. It rejects
duplicate sources and independently proves coverage by the exact
port/middle-orbit count

```
binomial(n-3,3) * 2^floor((n-3)^2/4).
```

The ordinary and optimized runs produce identical certificate-stream
hashes and mathematical results: see [normal](finite_head/normal.json),
[optimized](finite_head/optimized.json), and the
[minimal decoder control](finite_head/decoder_control.json).
No claim here depends on the conditional all-order classification of
pure-blue tournaments: every source in the three finite universes was
actually repaired. Compiled executables and raw certificate streams
are not committed; the runner builds in a temporary directory and
validates the entire stream before recording its digest.

## 3. Sharp unfiltered boundary at seven

For `P=(0,1,2,3,4,5,6)`, the source row tuple

```
T=(6,4,8,19,103,79,15)
```

is closed, P-blue, and has endpoint scores (2,4). Its exact target is

```
K_P(T)=(122,124,121,112,32,64,16).
```

This is the pure-blue join `C3 -> 1 -> C3`. The independent finite-head
checker enumerates all nine directed Hamiltonian words in that literal
target and verifies that every one is blue. A K move retains its word
and preserves whether its reversal is an anti-automorphism. Therefore
a black move entering this target would also provide a black word in
the target, a contradiction. There is no black detour of any length.
This is stronger than a failed three-word ledger and is separate from
the non-pure-blue F7 line whose shortest detour has length five.

## 4. The genuine next boundary

Closed endpoint-four replacement does not settle arbitrary open
endpoint-four words. If the source has endpoint scores `(n-1-k,k)`,
then K_P(T) has scores `(k+1,n-2-k)` and the closing pair changes
orientation. The layer therefore changes by

```
k -> n-2-k.
```

An open endpoint-four word becomes closed layer n-6, with scores
`(5,n-6)`, not generally closed layer four; it stays in layer four
only at n=10. If every closed layer were repaired, reversing the black
detour would also handle every open blue line. The present fixed-layer
result cannot replace that missing all-layer argument.

The next structural degree layer has four incoming interior ports
(endpoint scores `(n-6,5)`). The earlier three complementary row-pairs
plus one additional forced port only imply at least four ports, so
they no longer contradict its source budget. A repair must exploit
more positional information or a new operation. Even after all
admissible blue lines become black-redundant, full black-core
connectivity still needs the pure-blue-deleted line metagraph to be
connected, exactly as in THM-4732 section 2, Lemma 2. Neither condition
is supplied for all layers by the current endpoint-four theorem.

The inheritance board is: complete marked-source universes; literal
K targets; reversal/degree transport; pure-blue incidence; surviving
port budget. The anchor closes this one degree layer. The niche repairs
the exact finite-head decoder boundary. The wildcard is an all-layer
change-of-marked-path operation that preserves a route to the original
literal target. Its existence remains an open research question.
