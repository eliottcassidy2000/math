# V3 correction: a displayed Hamiltonian word, not a numerical theorem

**Status: corrected and source-aware rechecks PASS; approved by the final reviewer after source-aware replays of this replacement package.** The earlier V1/V2 manifest approval was withdrawn. The invalid displayed word was caught before the earlier package was published. V1/V2 stages and their original manifests remain preserved for provenance, not active publication approval.

## Exact defect and repair

In `regular_return_theorems.md`, Section 5.3, the second row of the displayed table incorrectly used

    P=(0,1,3,4,2).

The carousel R has arcs `i->i+1,i+2` modulo five. Thus `4->2` is absent, and the displayed word is not a Hamiltonian path of R. The corrected row is

    P=(0,2,4,1,3), with Q=P only.

This is the only change to the mathematical theorem file. The original SHA256 is

    ea9c2b020b1ed632682b023873806069153fd0eae8655038a23e0955d7925729

and the corrected V3 SHA256 is

    db2bd384ca2ad4251870723f2e8cea7a872c65ec487ac10fb0a35e84bbef001b.

Direct enumeration from the corrected printed table gives:

- `(0,1,2,3,4)` has intermediate code 733 and only its own matching-endpoint continuation.
- `(0,2,4,1,3)` has intermediate code 302 and only its own matching-endpoint continuation.
- `(0,2,3,4,1)` has intermediate code 654 and the three continuations printed in the table.

These are exactly all three Hamiltonian paths starting at zero. Therefore the conclusions H=15, M=25, the return weights, the two icosahedra, the Gram factorization, actual face incidence, and stationary-mass corollary are unchanged.

## How it was caught, and the audit correction

The final reviewer checked the displayed word against the carousel's literal arcs and found the forbidden adjacency `4->2`. The producer's executable receipt had already enumerated the correct paths. The initial independent audit verified the mathematical mechanism and finite output but failed to compare the prose table's literal words against that output. Its statement that the displayed representatives were correct was therefore wrong.

The amended `AUDIT.md` acknowledges this specific omission. The new `check_literal_examples.py` parses the current proof's table directly, verifies every printed word and continuation, and requires the table to equal the complete enumerated starting-at-zero path set. The original bad word is an explicit hostile control; the checker requires its only missing arc to be `4->2`. This closes the text-to-receipt gap rather than merely repeating the earlier generic computation.

## Other literal examples and figure

The same source-aware checker validates:

- The explicit path pair in equation (11).
- Every listed directed triangle and every listed one-vertex triangle pair.
- All twelve tournament codes and all rows of the icosahedral certificate, including their stated rotation/converse provenance.
- The Gram addendum's intermediate code 654, masks `(28,1,10,18,6)`, exceptional endpoints, interior directed cycle, and three regular neighbours.
- The audit's seven-vertex path pair and its claimed six-cycle support.
- The SVG's 24 actual regular-tournament labels and all 60 distinct six-edge return edges, with no omissions or extras.

The final reviewer performed the visual review of the figure. The independent check concerns its labels and graph edges. Its README caption states that it is a projected drawing of the derived labelled return graph, not the original merged metagraph or a unit-distance embedding.

A supplementary finite control also verifies the intermediate types: 120 private columns have H=11, automorphism order 1 and blue fibre 3; 120 private columns have H=13, automorphism order 1 and blue fibre 3; the forty face columns have H=15, automorphism order 3 and blue fibre 1. Here a blue Hamiltonian word means that its positional reversal is an anti-automorphism, and blue fibre counts are blue words divided by automorphism order. This refinement is recorded in the new receipt; no original theorem statement depends on it.

## Package preservation

`V2_TO_V3.patch` records every change to existing Markdown files exactly. The main theorem changes only the one table row; `AUDIT.md` and `README.md` are amended for the correction and the newly included figure. All original producer programs, receipts, the Gram addendum, the independent mathematical checks, and the stationary-mass proof/checker remain byte-identical.

Historical V1/V2 package checksum lists are not copied as current verification lists into V3, because the corrected files intentionally have different hashes. They remain with the preserved V1/V2 package. The original `FROZEN_SHA256SUMS` is retained as a clearly historical staging record. `PUBLICATION_SHA256SUMS_V3_1` is the current package-relative verification list; `READY_MANIFEST_V3_1.json` identifies the reviewed scientific replacement package.

## Receipt-only staging repair before approval

The first V3 proposal's new checker validated equation (11) correctly, but a later supplementary loop reused its local P variable before JSON serialization. Consequently that proposal's receipt printed the wrong P in its `literal_pair_11` field. This was caught immediately after staging and reported before approval. V3.1 snapshots the checked equation-(11) witness before that loop; its printed P is now exactly `(0,2,3,4,1)`. No theorem, validation outcome, or numerical result changed. The V3 proposal remains preserved separately, and `V3_TO_V3_1.patch` records the checker/provenance repair.
