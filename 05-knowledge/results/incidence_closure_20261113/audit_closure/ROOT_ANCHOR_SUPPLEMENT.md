# Source audit of the independent residual-antichain anchor solver

**Verdict: PASS on source logic, differential controls, and the completed normal
and optimized order-seven executions.** This supplement reviews `audit/anchors.py`
independently of the original optimizer. It does not change producer statuses.
The full receipts are byte-identical: both independently reject at most eight
anchors in each of the four palettes and verify nine-anchor witnesses. This
supplement is frozen after direct comparison of the completed receipts.

Reviewed source SHA256:
`1cadc0b12e072a0cecc18cdfa19e52aed438708bd3b5ac50295cde59395c2609`.
The shared, previously audited `closure/data.json` hash is
`c0484a99712edff88e8158b5f6000097f32f88874150d059da5aefa20f3de8f7`;
the nine-anchor witness input `closure/analysis.json` hash is
`41093c9b0b11e055ceaafcd3b6e247da1077cbe978a7addd9661fe1300e65023`.

## Independent data reconstruction

The new program imports no producer modules. It reconstructs ordinary ordered-
cover equality on the merged class space, obtains 184 collision pairs, builds
raw A and blue B from integer triples, and recomputes every source degree by
summing A. Source degrees are positive in the exact frozen data, as checked
separately by this lane's certificate audit and the root's rational-walk audit.
Normalized comparisons multiply by opposite **source** degrees, so use P=D^-1 A.

In the colored case it compares `(A,B)` rather than the producer's `(B,A−B)`.
These are equivalent output coordinates: equality of total and blue is exactly
equality of blue and black, by the invertible linear transformation
`(a,b)↦(b,a−b)`. The equivalence remains true after dividing both coordinates by
the same source degree. Thus this independent formula reconstructs the correct
coverage masks, including normalized and colored regimes.

The program reads nine-anchor lists from the producer's analysis, but directly
verifies their distinctness and coverage using its newly computed masks. The
lower-bound search does not use the claimed minima, search counts, packing,
pruning outcome, or failed budgets from the producer. Shared exact incidence
inputs are the scope of this independent optimization check; it is not a second
order-seven tournament census.

## Residual antichain and completeness

Let `left` be the remaining collision obligations and `remaining` the budget.
Every column becomes its intersection with `left`. Zero and equal residual
masks may be discarded. Sorting by decreasing cardinality ensures every strict
superset precedes its subset; the pass retaining masks not contained in an
already retained mask therefore returns exactly the residual maximal masks.
Any discarded column may be replaced by a retained containing residual column
without increasing the number selected or losing an obligation. This argument
holds afresh at every recursive node.

The retained residual masks are actual restrictions of original columns. A
previously selected column has zero intersection with the remaining obligations,
so cannot reappear as a useful retained mask. No omitted selected-column state
or availability condition is needed. If replacement merges two future choices,
that only lowers their number; the decision problem permits **at most** the budget.

The pruning rules are valid:

- If their union does not equal `left`, an obligation is uncovered by every column.
- Retained masks stay in nonincreasing cardinality order; their first mask gives
  a true maximum one-column coverage b. If `b*remaining < |left|`, no cover exists.
- For budget one, feasibility is exactly that a retained mask equals `left`.
- Each uncovered demand's candidate neighborhood is its set of retained masks.
  A collection of pairwise disjoint neighborhoods requires distinct anchors.
  The greedy packing is a lower bound even though it need not be maximum.

Finally select one uncovered demand and branch on **every** retained residual
column covering it. Every feasible cover can first be replaced by retained
maximal columns, and must then contain a column in this candidate neighborhood.
Each branch removes precisely the obligations covered by its selected column
and decreases the budget once. Thus the disjunction is exhaustive. The cache
key `(left,remaining)` is sufficient because all residual choices are determined
by `left` and the fixed global column list. `left=0` returns true even with budget
remaining, so an eight-budget failure excludes every cover of size at most eight.
Together with a verified nine-anchor witness, this proves minimum nine when the
complete execution returns false for budget eight.

This is a distinct implementation from the original static-antichain search:
residual dominance and candidate neighborhoods are recomputed at every node,
the branching demand is chosen from current rather than initial neighborhoods,
and it directly tests only budget eight. The shared set-cover problem itself
is necessarily the same.

## Independent differential controls

`check_root_anchors.py` loads only this callable module, never its complete n7
main. Its reference implementation enumerates literal subsets of original
columns. It tests every family of at most five distinct nonzero four-bit masks
(4,944 families) at every budget zero through four, plus eleven controls for
empty targets, uncovered targets, zero/duplicate columns, strict target subsets,
column permutations, dominance and noncontiguous bit indices.

All **24,775** comparisons pass. Normal and optimized receipts are byte-identical.
This is a differential test of the actual reviewed routine against an independent
reference, not a claim that loading the routine is itself an independent solver.
The all-finite completeness reasoning is given above.

Reproduction from the isolated worktree:

```sh
python3 05-knowledge/results/incidence_closure_20261113/audit_closure/check_root_anchors.py > 05-knowledge/results/incidence_closure_20261113/audit_closure/root_anchors.normal.json
python3 -O 05-knowledge/results/incidence_closure_20261113/audit_closure/check_root_anchors.py > 05-knowledge/results/incidence_closure_20261113/audit_closure/root_anchors.optimized.json
```

Checker SHA256:
`d6eed6c7f9cea63849f948ba8995a3f13f7cdbff3bc5e537f2d8d2c2419bfee5`.
Both control receipts SHA256:
`6f6d1fc939532bf25971c9a7a4507b85825ac9a89bbdfcb2734a235406aa1048`.

## Completed independent normal execution

The completed `audit/anchors.normal.json` reports PASS, 493 explicit checks and
480 internal small controls. Its reviewed source, incidence-data and witness-input
hashes match the three frozen hashes above exactly. Direct inspection confirms
`budget_eight_feasible=false` and `minimum=9` for every palette:

| Palette | Distinct memoized search states | Minimum |
|---|---:|---:|
| A | 38,757 | 9 |
| P | 36,606 | 9 |
| colored A | 38,620 | 9 |
| colored P | 35,943 | 9 |

These are search-state counts, not counts of eight-element subsets. Exhaustive
coverage follows from the audited recursion. The completed normal computation
therefore supplies the independent order-seven lower-bound execution that was
explicitly outside the earlier `REPORT.md` controls. Together with that report's
independent small-order enumeration, the fixed-anchor minimum sequence
`0,0,1,2,9` is independently verified across all four palettes in the declared
orders 3 through 7.

Normal execution receipt SHA256:
`2680658edcefaa86d353fbe6a45f38acfa956dbd919a9ee18db437142c3d7f99`.

**Optimized full replay: completed and byte-identical to the normal receipt.**
Both full receipt SHA256 values are
`2680658edcefaa86d353fbe6a45f38acfa956dbd919a9ee18db437142c3d7f99`.
The reviewed solver source hash remains unchanged. This complete replay is
separate from the small normal/optimized differential controls above. No source
repair is requested. This report is frozen; any later changes require a new
source/receipt binding.
