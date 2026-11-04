# Independent review of the three fixed q3 carriers

2026-11-03. **VERIFIED proof identities and FINITE-EXACT literal replay.**
All checks passed in normal and optimized Python. This audit edits no
producer files and imports no producer implementation. It establishes the
fixed-carrier claims below; a separate SCC classification is required to
turn them into an arbitrary-deletion-carrier theorem.

## Actual inherited inputs and independent counting

The inputs are the two actual named C2 children of counts1409 and1427 in
[`entry_selector/receipt.json`](../../c4_port_returns_20261103/entry_selector/receipt.json).
All13 names are retained. The audited apex is physical0. The three deleted
carriers are, in their exact physical order,

```
TT(99,10,11,2,1,3,4,12,200) -> C3(201->401->400->201),
TT(99,10,11,2,1,3,4,12) -> C3(201->401->400->201) ->200,
C3(201->401->400->201) -> TT(99,10,11,2,1,3,4,12,200).
```

The scripts reconstruct the second carrier by its three actual reversals
incident with200 and the third by its27 actual cap/prefix reversals. They
then reconstruct every one of4096 cuts around0 per carrier using an
independent adjacency-set implementation. They check the precise bit order,
completeness, uniqueness, all saved rows, and directed strongness.

[`exact_counter.cpp`](exact_counter.cpp) counts literal directed Hamilton
paths by subsets and endpoints on each graph and its converse. Its separate
subset-cover convolution supplies all E1/S1/E2/S2/J, P2, and every deletion
count. It does not use the producer's deleted-J quadratic to obtain H.
The underlying endpoint recurrence is standard; independence here is a
separate implementation applied to actual matrices, not a claim of a new
algorithm. No graph is accepted only through its target scalar count.

## Proof checks

The deleted-J identity in [`star/q3/PROOF.md`](../star/q3/PROOF.md) is exact:
putting the apex last, first, or between two nonempty directed paths gives
its three disjoint sums. The quadratic coefficient J is ordered and keeps
both boundary ports; identifying unordered covers would lose the formula.

For the two carriers with99 at the beginning, strongness forces0->99.
Vertex99 then has sole in-neighbor0. Removing99 leaves a deleted carrier
whose three Hamilton paths all start at10, giving exactly

```
E99(2)-H = 3[0->10].
```

This root guard is also verified directly on all1408 accepted matrices.
It is not assumed on an equal-count replacement.

For `TT_l -> C3 -> z`, the positive contribution J(z,r), where r is the
first prefix vertex, is `2^(l-1)*12`. The cap can occupy either or both
ordered path slots: its slot count12 includes the cases with an empty cap
piece. This agrees with the counter's P2 convention. Strongness supplies
the two apex boundary arcs, so these covers inject into actual Hamilton
paths and cannot cancel. At l=8 the resulting lower bound1536 excludes
both target counts even before any entry or root guard.

In the same relocated carrier, z has sole out-neighbor0. Its start defect
is the number of cap vertices beating0: an endpoint path in the deletion
of z must traverse the transitive prefix, then the triangle, then0. Each
cap endpoint contributes exactly one. Every guarded matrix independently
verifies this terminal law. Thus the relocation produces entries, but the
positive path contribution already makes its counts too large.

The [initial-cap proof](../star/q3/initial/PROOF.md) keeps its claim finite
and permits **any** possible end-zero root. Every exact-count target has
an end-zero vertex, but none has start defect1 or2 anywhere. Choosing a
different physical root cannot repair that missing positive entry.

## Complete replay results

| Carrier | Literal H checks | Strong graphs | Root99 end-zero graphs | Positive-entry graphs in that root chart | H1409 / H1427 hits |
|---|---:|---:|---:|---:|---:|
| TT9 -> C3 |4096|1792|896|0|6 /3|
| TT8 -> C3 ->200 |4096|1024|512|384|0 /0|
| C3 -> TT9 |4096|1792|not imposed|not censused outside target fibres|6 /3|

For the first two rows, all1408 guarded graphs receive full independent
profiles, including every named start/end defect and every ordered J
coordinate. Their nine exact-count witnesses have all rows, names, source
reversal pairs, full profiles and fresh F49/F343 observers checked. Four
additional profiles cover the deleted carriers and retained source children.

For the third row, all4096 H values are checked literally, followed by full
profiles of all nine strong target-count witnesses, the deleted carrier,
and both retained source children. All nine targets have some end-zero
vertex and **zero** start-defect1/2 vertices. All incoming/outgoing cuts and
fresh named phase charts are checked. A corrupted self-loop input is
rejected by the independent counter before it emits any accepted profile.

The total independent universe is12288 star matrices. Normal and `-O`
receipts match byte-for-byte. The required predicates use explicit checks,
not removable Python assertions.

## Scope and reproduction

These results do not assert native legality of the listed reversals,
strongness throughout a physical repair, an arbitrary-order q3 obstruction,
or reusable positive-digit closure. The all-order q1 exclusion is the
separate theorem in [`ONE_APEX.md`](../recovery/ONE_APEX.md). The all-order
terminal-cap no-entry theorem linked by the producer is not needed for
this finite replay and is not certified by the present note.

Run from the repository root:

```
python 05-knowledge/results/compensation_returns_20261103/audit/check_q3_review.py
python -O 05-knowledge/results/compensation_returns_20261103/audit/check_q3_review.py
python 05-knowledge/results/compensation_returns_20261103/audit/check_initial_review.py
python -O 05-knowledge/results/compensation_returns_20261103/audit/check_initial_review.py
```

The receipts `q3_review.json` / `q3_review_optimized.json` and
`initial_review.json` / `initial_review_optimized.json` pin the reviewed
producer proofs, actual data, independent checker sources and counter.
