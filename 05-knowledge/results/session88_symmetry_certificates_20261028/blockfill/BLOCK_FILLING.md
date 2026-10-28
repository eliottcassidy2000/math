# One complete named row fills one equivalence block

Session88, 2026-10-28. **PROVED** elementary finite-partition and transport
statements below. **FINITE-EXACT** tournament checks through order eight.
No method-novelty, asymptotic speedup, or new counting-sequence formula is claimed.

## Result first

For ordinary tournaments at order eight, block filling returns exactly the same
complete root-switch operator and tournament bank as the full-row reference,
while reducing classified rows from **54,256 to 6,880** and root-classification
requests from **434,048 to 55,040**. The equality check includes every pointed
state's complete, vertex-indexed transporter sets and full automorphism group,
not merely counts or anonymous row shapes.

Three alternating local runs with identical full outputs and cold per-stage
literal-matrix caches gave order-eight medians **12.255 s versus 3.323 s**
(3.69 times), including parent preparation. This is not a hardware-independent
or asymptotic factor. The reference here materializes more transport data than
the old Session87 benchmark, so its seconds must not be compared directly with
that old benchmark's class-only output.

A further, exact adaptive-root experiment reduced 55,040 first-row root queries
to 47,391. A single illustrative pass did **not** improve runtime: 3.056 s versus
3.075 s. Count reduction and actual acceleration are different conclusions.

## Inheritance, scope, and loss ledger

Pinned source: `7a0bab9bfa7bb87c5a10c7ba239592d90cfafd51`.
The hash-verified external Session87 baseline is
`05-knowledge/results/session87_four_vertices_enumeration_20261028/check_generation.py`,
SHA-256 `f5d3da961f34dc1a32eaa0c1b69fad9bda983154578ea0187c5fd6b2e2d10fae`.
`baseline_sources.py` loads those exact bytes from the repository; no inherited
source snapshot is vendored. A staging run can supply `--baseline-root PATH`.

Closest proved mechanism: Session79 `ROOT_SWITCH.md`, exact weighted root-orbit
blocks, inherited by Session87 `EXACT_GENERATION.md`. Corrected near miss:
parent canonical form plus an independently minimized attachment mask loses
alignment. Canonical hostile: identical anonymous rows need not mean identical
named pointed states. Least-used restored sidecar: all transports attaining the
**joint** parent-and-mask key, which can reconstruct the whole automorphism
action, rather than just root-orbit sizes.

The carrier is every ordinary pointed tournament class; converse is not merged.
All parent attachment masks are allowed. No Hamiltonian, prime, regularity, or
arithmetic filtering is used. Counts are consulted only after construction.
Every production canonical search has at most N-1 vertices; full-child searches
are isolated audit oracles and never supply the next bank.

Concept board: named equivalence fibers / weighted reset / complete root ports /
point stabilizers and transports / output-sensitive computation / cheap selectors.
The anchor is full-Q generation; the niche is adaptive first-row construction;
the wildcard is the finite measured-partition interpretation.

Source → target: pointed classes → ordinary children. Map: forget the root.
Preserved: ordinary child identity, full Q when indexed by named states, raw
representatives. Lost by anonymous weight profiles: the actual child, port
identities, stabilizers, and joint-port orbits. Sidecar: pointed-state dictionary
and aligned transporters. Cheapest decisive hostile: the two order-seven cyclic
tournaments below.

## 1. Finite weighted-fiber theorem

Let X be a finite set, pi:X→C a surjection, and assign each x a positive integer
weight w_x. Define

    Q[x,y] = w_y if pi(x)=pi(y), and 0 otherwise.

Let n_c=sum_{pi(y)=c} w_y. A complete row Q[x,-], **with its column names in X**,
has support exactly pi^{-1}(pi(x)). Every row indexed by that support is identical.
Therefore the following procedure computes the complete operator exactly:

1. Choose any as-yet-unfilled x.
2. Compute its complete row using the original exact row oracle.
3. Fill that row at every state in its positive support.
4. Continue until every state is filled.

**Proof.** Positivity makes the support exactly the entire fiber. The formula for
Q depends on x only through that fiber. Two fibers are disjoint; consequently
no earlier block can intersect a later unfilled block, and every state is filled
once. Exactly |C| oracle rows are used instead of |X|. QED.

No oddness, tournament structure, group action, or equal fiber mass is needed.
For incidence R[x,c]=1_{pi(x)=c} and weighted incidence
W[c,y]=1_{pi(y)=c}w_y, Q=RW and WR=diag(n_c). Thus

    Q^2 = DQ = QD,       D[x,x]=n_{pi(x)}.

When n_c=N is constant, Q²=NQ. Over the rationals or reals,
P[x,y]=Q[x,y]/n_{pi(x)} is idempotent. It replaces a function by its weighted
mean on each fiber: precisely conditional expectation onto the partition.
For arbitrary positive component constants lambda_c,
mu_x=lambda_{pi(x)}w_x satisfies mu_x P[x,y]=mu_y P[y,x], so P is reversible.
The map preserves functions constant on fibers and kills weighted-mean-zero
contrasts inside each fiber. This is a genuine operator identification, not
an inference that arbitrary dynamics have become a rank-one reset.

For a finite G-set of legal roots/ports, X is the orbit set over each object and
w_x is the orbit size |G|/|Stab(root)|. Tournament automorphism groups have odd
order, so their root weights are odd; that is an additional fact, irrelevant to
block filling but useful to the binary-rank theorem. For general groupoids one
must specify an actual finite marked carrier or relative measure: raw arrow
counts or arbitrary rational groupoid weights are not automatically orbit sizes.

**Boundary.** Positivity and completeness are essential. Reducing coefficients
modulo a prime can erase positive entries and whole blocks. Fill using exact
integer support first; do not use a reduced-characteristic support as though
it were the original equivalence relation. Knowing Q²=NQ alone does not certify
that a purported row has the right support or weights.

## 2. Transporting the full sidecar, not only Q

Represent a pointed class j by its canonical parent and canonical attachment
mask, with the distinguished new root last. For an anchor child U, classifying
root v returns its complete pointed key j and **all** parent orderings attaining
that joint key. Appending v gives all isomorphisms from the stored representative
of j into U sending its distinguished root to v. The code records their full
permutations, not just their number.

Write U_j for the **underlying unpointed tournament** of stored pointed class j.
For each j in the anchor row, choose one such map F_j:U_j→U. If B_k is the union
of all recorded maps U_k→U over all vertices of pointed type k, then

    Iso(U_k,U_j) = F_j^{-1} B_k.

This follows simply by composing an isomorphism with F_j or its inverse. In
particular Aut(U_j)=F_j^{-1}B_j is the full **unpointed** automorphism group,
not the stabilizer of j's distinguished root. Equivalently, with anchor group G,
all unpointed pairwise transports are F_j^{-1} G F_k. Restricting a record to
one target root v gives the corresponding pointed transporter fiber; restricting
to the distinguished target root gives its point stabilizer when the IDs match. Individual vertex indices are
transported by F_j^{-1}; therefore the record retains which map lands at which
root. Every target state's full root record can be filled without a new
classification request.

The first row thus supplies the **aligned** sidecar needed for every pointed
presentation in the block. It does not say that the numeric row itself contains
the group. In fact a known full point stabilizer plus one certified transporter
to each root in its orbit is a smaller sufficient carrier; retaining all maps
here deliberately makes the equality audit especially strong.

Implementation guardrails: all row entries are positive integers; source lies
in its support; every named target exists in the dictionary; total root mass is
N; no target is filled twice; no target remains uncovered; transport fixes the
correct distinguished port; all rows in a support agree; root records reproduce
weights. Independent direct replay verifies every recorded permutation preserves
every adjacency bit and checks sparse Q²=NQ at every state.

## 3. Matched full-output measurements

Both methods use exactly the same raw-mask state construction, parent
canonicalizer, literal-matrix cache, parent normalization, and full sidecars.
The reference classifies all N roots at every pointed state. Block filling
classifies all N roots at the first state in each block. Sorted state iteration
means both choose the same first representative, so the raw bank is equal
entry for entry, not only up to isomorphism.

| N | Classes | Pointed states / old rows | Filled rows | Old root queries | New root queries | Old / new searches after cache, including preparation |
|--:|--:|--:|--:|--:|--:|--:|
|1|1|1|1|1|1|1 / 1|
|2|1|2|1|4|2|1 / 1|
|3|2|4|2|12|6|2 / 2|
|4|4|12|4|48|16|7 / 4|
|5|12|48|12|240|60|28 / 15|
|6|56|296|56|1,776|336|203 / 73|
|7|456|3,040|456|21,280|3,192|2,386 / 643|
|8|6,880|54,256|6,880|434,048|55,040|48,326 / 11,074|

At order eight, raw-bank and parent preparation add 59,280 requests to both
methods: 58,368 raw masks plus 912 parent-map requests. Total requests are
493,328 versus 114,320. A request can be a cache hit; searches and tested
permutations are separately recorded in `receipt.json`.

| N | Full-Q median seconds | Block-fill median seconds | Ratio |
|--:|--:|--:|--:|
|7|0.408021|0.230785|1.76797|
|8|12.254597|3.323452|3.68731|

Timing includes parent preparation, actual full output construction, transport
reconstruction, and intrinsic validation. It excludes earlier levels, output
hash/equality comparison, and independent full-child/labelled audits. Each run
starts with a fresh per-stage cache; execution order alternates. The benchmark emits raw individual
times, including variability at order seven. No comparison to an
uncached historical executable or a C++ generator is asserted.

### Immutable rows, explicit-output size, and a factored API

The implementation stores each row as an immutable tuple of `(named_state,
integer_weight)` pairs, and maps every state in a block to that same tuple.
This has the same numeric operator semantics as separately stored equal rows.
It does **not** promise the object-identity semantics of distinct row objects;
sharing mutable row lists would be an API change and is not done here.

If a consumer insists on separately materializing every nonzero Q entry, the
output contains sum_c |X_c|² entries (or |X|² entries for a dense matrix).
Classification-call savings do not remove that output-size cost. This checker
also explicitly constructs all requested per-vertex transporter sets; their
output and reconstruction costs remain in the timing.

Alternatively retain one named row per block plus the state/transport dictionary.
The total number of factor entries is sum_c |X_c|=|X|. For a vector f, compute
one weighted sum s_c=sum_{j in X_c} w_j f_j and broadcast s_c to every state in
X_c. This applies Q in O(|X|) arithmetic operations and stores O(|X|) numeric
factor entries. Integer bit complexity depends on operand sizes. Q^k applies
n_c^(k-1)s_c on each block for k>=1; ranks and block counts follow separately
from the specified coefficient field. No named dictionary or joint-port sidecar
is erased. `check_factored_operator.py` checks exact Qf and Q²f against the
explicit sparse operator through eight.

**Different task from Session87's selector.** The neighbor-degree selector can
emit one raw class representative after only 693 competitor-key requests at
order eight, but does not produce full Q, all pointed roots, or full transport
sidecars. It remains a better choice when that smaller output suffices. The
present result is a fair same-output saving; it does not supersede the selector.

## 4. Adaptive first-row construction: further exact saving, no measured win

Start with the known automorphisms of the canonical parent that stabilize the
attachment mask, lifted to fix the original root. They form its **full** point
stabilizer H. The original root's key and all its fixed-root transports are
therefore already known without a new deletion classification.

Maintain a certified subgroup K≤Aut(U) and complete root records for some
vertices. If g∈K sends an already classified u to v, composition by g transports
its complete record to v, so v need not be classified. Otherwise classify an
uncovered v. Its returned stabilizer can be lifted into K; when its key equals
that of an earlier u, composing their two canonical transports supplies an
automorphism u→v and can enlarge K. Close the subgroup and propagate again.

Every inference is by a genuine automorphism and thus exact even while K is
proper. At termination K contains H and reaches every vertex in the true orbit
of the original root (classified matches add a transporter; inferred matches
already have one). Orbit-stabilizer then forces K=Aut(U). No assumption that
initial H acts transitively is used.

| N | Fixed first-row queries | Adaptive classified roots | Roots inferred by subgroup | Original roots already known |
|--:|--:|--:|--:|--:|
|1|1|0|0|1|
|2|2|1|0|1|
|3|6|3|1|2|
|4|16|9|3|4|
|5|60|38|10|12|
|6|336|244|36|56|
|7|3,192|2,593|143|456|
|8|55,040|47,391|769|6,880|

The last three columns sum to the second. Most of the extra saving is the
already-known source root: that alone would leave 48,160 queries at order eight.
Certified subgroup propagation saves a further 769. At eight there are 88
subgroup enlargements. `adaptive.json` compares all full outputs, and its optimized
Python replay checks the same equality. Timing is only one illustrative pass;
the subgroup closure overhead outweighed the small query saving at seven/eight.
A generator-based orbit data structure may reduce that overhead, but no such
implementation or acceleration is claimed here.

## 5. Hostiles: exactly which implication fails?

1. **Partial row.** Omitting a lawful root may omit a whole pointed type; then
   support is not the fiber. An unrenormalized omission fails total mass N.
   Renormalizing can conceal it: 3I on three states satisfies row sum three,
   odd positive nonzero entries and Q²=3Q, yet falsely splits a transitive
   three-root block. Repair: exhaust all roots or infer each through certified
   transport, never infer completeness from polynomial identities.
2. **Absent port.** A positive column name not in the pointed dictionary cannot
   be filled. Reject immediately. A correct scalar count does not repair a
   nonexistent named attachment state.
3. **Wrong orbit weights.** Both constant-row blocks `(1,3)` and `(3,1)` have
   positive odd weights, sum four, and Q²=4Q, but differ as named operators.
   A real source→C3 tournament has these two orbit sizes. Algebraic checks
   alone cannot tell which named state gets the weight three. Repair: count
   actual root occurrences and verify their common-frame transports.
4. **Anonymous row.** On Z/7, orient i→i+s for s in {1,2,3}, or for Paley steps
   {1,2,4}. Both are vertex-transitive and have the anonymous block `[7]`, but:

   | Steps | Named pointed key | H (directed Hamilton paths) | Aut order | Orbits on ordered distinct pairs |
   |---|---|--:|--:|---|
   |{1,2,3}|(280,56)|175|7|six of size 7|
   |{1,2,4}|(8850,56)|189|21|two of size 21|

   Full S7 permutation checks independently verify their automorphism groups.
   They have distinct canonical child codes. The first failed implication is
   “same row shape ⇒ same named block/object.” The strongest survivor is the
   weighted-fiber theorem **in a shared named dictionary**. Full joint ports
   require the transporter sidecar, not merely the one-port orbit profile.

## 6. Reproduction and verification boundary

```sh
python3 -B check_blockfill.py --max-order 8 --benchmark 3 > receipt.json
python3 -B -O check_blockfill.py --max-order 8 > receipt_optimized.json
python3 -B check_named_dictionary.py > named_dictionary.json
python3 -B check_adaptive.py --max-order 8 > adaptive.json
python3 -B -O check_adaptive.py --max-order 8 > adaptive_optimized.json
python3 -B check_transport_replay.py > transport_replay.json
python3 -B check_factored_operator.py > factored_operator.json
```

All meaningful checks use explicit exceptions, not removable Python assertions.
Normal/optimized receipt hashes agree for complete Q, states, raw bank, all
vertex-indexed transporter sets, all automorphism permutations and parent data.
The independent full-child oracle audits each ordinary class and its full
automorphism group through eight; the independent full labelled S_n-orbit census
covers all labelled tournaments through six (including all 32,768 at six).
Direct transport replay and Q-square tests cover every pointed state through
eight. Beyond order eight is untested, although the finite-fiber and transport
proofs are not bounded by eight.

The public checkers need no committed receipt files. Each reproduction command
emits its own local JSON. `verify.py` replays all five checker families and
checks normalized semantic results against compact pinned hashes; timing and
benchmark measurements are excluded from those hashes. Run both
`python3 -B verify.py` and `python3 -B -O verify.py`.
