# Independent audit of direct ordinary skip contraction

2026-10-26. **REVIEW PASS / FINITE-EXACT.** The separate
`contracted_ordinary.py`, SHA256
`63b45349cacb5eb9312c6f98a02e0e459a4e478b6319e6f87005b4c6211fd37a`,
correctly constructs the smaller ordinary graph directly. It does not first
materialize the full graph. Its scope is exactly deterministic unaffordable
stage skips and early completed-state skips; active defect distinctions remain.

## Proof of the precise contraction

A removed nonterminal vertex has S=0 or ell>S. It has precisely one legal
ordinary choice k=0. When S>0 and ell>S, the retained odd coordinates are
already capped by S; decreasing ell by two cannot pass an affordable odd
length until the active endpoint is reached. Hence S,D,u,b are unchanged
along every removed edge. For S=0, D is empty throughout. The exact formulas
therefore give T=1 and J=0 on these edges.

Following this unique chain ends at the terminal for S=0, or at
`(min(ell,largest_odd_at_most(S)),S,D)` for S>0. At a surviving source, retain
every original choice k with its T,J and redirect its target to that chain's
endpoint. Deleting the uniquely determined unit chains gives a weighted
complete-path bijection; reinserting them is its inverse. This proves both
exact ordinary counts and continuation compatibility, not just final numerical
agreement. Active stage still strictly decreases, so the resulting graph is
acyclic. No stage-blind merger or clone/excess substitution is involved.

V is constant on each unique zero-cost chain and consequently agrees at every
full vertex with its contracted image. Every retained reduced cost J' is
unchanged. P and accumulated Y/Z are different: a canonical target may receive
additional histories. They agree with the full recurrence **at surviving
vertices**, not at every vertex mapped to that survivor. The implementation
and its documentation respect this distinction.

## Independent verification

`audit_contraction.py` does not use the producer's `canonical_vertex` or
`audit_contraction` to decide the map. It identifies removed vertices from the
full graph's stage/remainder, then explicitly follows each actual outgoing
chain, verifies uniqueness and T=1,J=0, and finds the survivor. It compares
the complete survivor set and every labelled local transition against the
directly built graph, including edge stage, multiplicity, committed degree,
and exact exponent. V is compared on every original vertex; P and exact Y/Z
are compared on surviving aggregates only.

For N=0..32 it independently enumerates all odd cycle partitions. Each full
path projects to the same contracted path with precisely the forced edges
removed. The audit verifies distinct partitions yield distinct contracted
paths, then independently counts **all** contracted DAG paths and checks this
equals the number of partition paths. This establishes both directions of
the finite path bijection, not only an injection. Independently computed
partition valuations supply exact low-valuation edge-support unions for
K in {1,2,3,8}; these equal the contracted potential filter exactly.

Both normal and `python -O` independent runs pass with identical receipts:

| Check | Coverage |
|---|---:|
| Orders | N=0..50,60,80,100 |
| Surviving exact coordinates including terminals | 24,102 |
| Surviving labelled edges | 97,091 |
| Forced vertices removed | 8,773 |
| Modular full-coordinate comparisons, K=1,2,3,8,17 | 270 |
| Independent partition paths, N=0..32 | 2,765 |
| Exact low-valuation edge-support equalities | 132 |

The original direct-contraction producer program was additionally replayed
normally and under `-O`; both outputs match its supplied receipt.

## Hostiles and consequences

* N=6 has the stage-erasure self-loop hostile
  `(3,1,(0,)) -> (1,1,(0,))`. The specified chain contraction deletes that
  forced edge and retains an active-stage vertex, rather than keeping a loop.
* N=8 has the active source `(3,3,(0,2))` with two choices; its k=0 edge has
  J=3. It is retained. “Remove every k=0 edge” would be false.
* At N=6, early completed `(1,0,())` has P=1 but its terminal image has P=0.
  This is a minimal recorded counterexample to pointwise P descent. It does
  not invalidate V descent or equality at surviving aggregates.
* N=0 is retained as its one-vertex terminal graph, with count one.

At N=100 the directly built graph has 10,547 input states and 47,147 edges,
against 13,444 input states and 50,044 edges in the full graph. At K=2 the
combined precision-filtered update count falls from 568 to 64; at K=8 it
falls from 2,576 to 1,161. All figures are independently reproduced. This is
an actual administrative graph reduction, unlike precision filtering alone.
It is not a new equivalence of active defect states, a new Burnside identity,
an asymptotically smaller recurrence theorem, or a measured runtime claim.
The modular evaluator still allocates an array for all contracted vertices;
its edge-update counts exclude graph construction and potential passes.

## Reproduction

Place this audit beside `audit_reverse.py`, with the producer scripts in
the sibling `session64-counting` directory, then run:

    python audit_contraction.py > contraction-audit-normal.json
    python -O audit_contraction.py > contraction-audit-optimized.json
    cmp contraction-audit-normal.json contraction-audit-optimized.json

For other layouts, supply `--source /path/to/directory/containing/producer/scripts`.

The primary reverse audit remains separately sealed by `REVERSE_AUDIT.md`.
No repository source file was changed during either audit.
