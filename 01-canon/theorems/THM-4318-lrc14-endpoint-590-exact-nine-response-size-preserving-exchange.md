---
id: THM-4318
title: "LRC(14) endpoint-590 exact-nine response and size-preserving exchange"
status: >
  PROVED RELATIVE TO THM-4314 + FINITE-EXACT + INDEPENDENT O2/O3 COMPLETE
  REPLAYS. The unchanged 3,925-mask carrier has exactly 100 endpoint-590
  failures, all on row (210,590). The complete rank-8/rank-9 response
  hypergraph has cover number exactly nine. A certified nine-addition,
  nine-deletion exchange closes all 493 accumulated rows while retaining
  size 3,925 and every joint mask. The typed union has 2,113 rows and leaves
  20,534, maximum endpoint 589 on 28 rows. No physical entry or LRC(14)
  follows.
source: root + endpoint592_capacity + endpoint592_cover_opt + endpoint592_structure / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4314-lrc14-endpoint-591-complete-layer-closure-and-deletion-boundary
related:
  - THM-4313-lrc14-endpoint-592-fortythree-response-size-preserving-exchange
  - THM-4302-lrc14-endpoint-596-response-minimum-four-and-size-preserving-exchange
  - THM-4282-inactive-signature-deck-surgery-endpoint-663
artifact_root: 05-knowledge/results/lrc14_endpoint590_exact_nine_response_exchange_thm4318
artifact_manifest: 05-knowledge/results/lrc14_endpoint590_exact_nine_response_exchange_thm4318/SHA256SUMS
artifact_manifest_sha256: 37470a5228c8c1505756772c7a796e8219aeed3656a1385e7bb220a90f9289b5
primary_scripts:
  - 04-computation/lrc14_endpoint590_exact_nine_response_exchange_thm4318/endpoint590_final_carrier_audit.cpp
  - 04-computation/lrc14_endpoint590_exact_nine_response_exchange_thm4318/endpoint590_response_structure.cpp
  - 04-computation/lrc14_endpoint590_exact_nine_response_exchange_thm4318/endpoint590_exact_no8_search.cpp
  - 04-computation/lrc14_endpoint590_exact_nine_response_exchange_thm4318/endpoint590_cover_direct_audit.cpp
  - 04-computation/lrc14_endpoint590_exact_nine_response_exchange_thm4318/endpoint590_low_witness_census.cpp
  - 04-computation/lrc14_endpoint590_exact_nine_response_exchange_thm4318/endpoint590_protected_deletion_quotient.cpp
  - 04-computation/lrc14_endpoint590_exact_nine_response_exchange_thm4318/endpoint590_exchange493_audit.cpp
  - 04-computation/lrc14_endpoint590_exact_nine_response_exchange_thm4318/typed_endpoint590_consumer.py
  - 04-computation/lrc14_endpoint590_exact_nine_response_exchange_thm4318/verify_endpoint590_packet.py
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. Baseline, response-atlas,
  low-witness, quotient, direct-cover, exact-no-eight, full-exchange, and typed
  O2/O3 artifacts agree byte-for-byte. Two independent implementations rebuilt
  the response atlas; an independent adversarial audit checked every exact-DFS
  prune. The packet verifier has no assert statements and checks manifest
  closure under normal and optimized Python.
---

# THM-4318 -- LRC(14) endpoint-590 exact-nine response and size-preserving exchange

**PROVED RELATIVE TO THM-4314 + FINITE-EXACT + COMPLETE O2/O3 REPLAY.
LRC(14) REMAINS OPEN.**

## 1. The unchanged carrier has one hostile row

Let `C` be THM-4314's 3,925-mask carrier, with ranks 3,818/107, FNV
`a0d08a38c10bdab7`, and all 421 joint masks retained. Let `R_590` be its
13-row endpoint-590 typed frontier, ordered FNV `44aa8a793d162cf9`.

Independent O2 and O3 builds reconstruct `C` from canonical ledgers and test
all

~~~text
13*binom(30,9)=185,992,950
~~~

row-body obligations. Raw transcripts, pair ledgers, and failure ledgers agree
byte-for-byte. Exactly 100 bodies fail, all at `(210,590)`; their ordered FNV
is `8d19cba1e86e53b5`. Every other endpoint-590 row is already closed.

## 2. The complete response-cover optimum is nine

For each rank-eight or rank-nine mask absent from `C`, record the subset of
the 100 failed bodies that it covers while active at `(210,590)`. Complete
enumeration of every such mask gives

~~~text
rank 8 responders       36,285
rank 9 responders      568,812
distinct nonempty signatures 14,368.
~~~

Masks with empty response are irrelevant. Replacing any response signature by
an inclusion-maximal containing signature cannot increase cover size, so the
14,368-signature hypergraph reduces losslessly to 1,165 maximal signatures.

An exact integer dual of denominator three has total weight 22 and maximum
signature load three, proving only the lower bound eight. A deterministic
depth-eight search then branches on every maximal signature containing a
chosen uncovered pivot. Its sum and residual-dual prunes are necessary
conditions, and local dominated-gain removal is replacement by a containing
gain. Independent O2/O3 traversals agree on

~~~text
nodes          7,163,197
dead states    6,775,810
sum prunes       943,159
dual prunes    5,575,606
result        UNSAT at depth 8.
~~~

Thus no eight responses suffice. Direct arithmetic verifies that the nine
rank-nine masks

~~~text
20490236 22045017 29224016 0a439108 0b220096
0120403f 12844116 10686016 084a6016
~~~

cover every failure. Their ordered FNV is `d1cf49e4b811b958`; the 100 bodies
have hit distribution 69 once, 26 twice, and 5 three times. Therefore the
complete rank-8/rank-9 response-cover number for this fixed failure ledger is
exactly nine. The numerical MILP report is corroborative, not a proof
dependency.

## 3. Exact protected singleton quotient

Adjoin those nine responses to `C`. With all 421 joint masks protected, an
exact low-witness census classifies deletion of one old nonjoint mask across
the 493-row target consisting of the 467 THM-4313 rows, the 13 endpoint-591
rows, and `R_590`.

For inherited rows whose old witness minimum is at least two, adjoining masks
cannot create a singleton. The remaining 42 inherited rows are rescanned over
all bodies; the complete endpoint-591 and endpoint-590 low-witness ledgers
supply the other rows. After the additions, 1,600 obligations have a unique
old nonjoint witness, supported on exactly 425 old masks. Hence, among the
3,504 old nonjoint masks, exactly 425 are unsafe and exactly 3,079 are safe for
a single deletion from the augmented carrier. This quotient preserves the
single-witness predicate and forgets simultaneous interactions.

The nine lowest-activity safe masks selected for deletion are

~~~text
06021829 23222801 12444083 20827018 29c04082
02916180 13c00881 070c4840 2380408a
~~~

with ordered FNV `3546eb56552b4cde`. The singleton quotient alone does not
certify deleting all nine together.

## 4. Direct simultaneous exchange theorem

Delete those nine masks and add the nine responses. The resulting carrier
`C'` has size 3,925, ranks 3,809/116, FNV `eeae5518d84ccac5`, and retains all
421 joint masks.

Independent O2/O3 programs directly test every body on all 493 target rows:

~~~text
493*binom(30,9)=7,053,424,950.
~~~

Their complete pair and failure ledgers agree byte-for-byte and give zero
failures, pair-ledger FNV `1092fd57a8581a34`. This raw replay, rather than an
inference from singleton safety, proves the simultaneous nine-for-nine
exchange.

## 5. Typed consequence and scope

Consuming `R_590` gives an exact 2,113/20,534 partition of the frozen
22,647-row universe. The union FNV is `c806cce6b836fdff`; the residual FNV is
`11285b5a49f4150d`. The next residual endpoint is 589 on 28 rows, FNV
`5d9429c9f9971322`.

This theorem is finite-exact for one labelled carrier, one rank-8/rank-9
response universe, and one accumulated row family. It proves neither global
carrier optimality, safety outside the 493 rows, a physical entry/owner map,
terminating descent, nor LRC(14). The exact cover number nine concerns the
fixed 100-failure response hypergraph; it is not a universal response lower
bound.
