---
id: THM-4314
title: "LRC(14) endpoint-591 complete layer closure and deletion boundary"
status: >
  PROVED RELATIVE TO THM-4313 + FINITE-EXACT + INDEPENDENT O2/O3 COMPLETE
  REPLAYS. The unchanged 3,925-mask THM-4313 carrier closes all 13
  endpoint-591 rows. A complete all-witness census classifies every deletion
  set of size at most two: exactly 2 singleton deletions and 7,869 pair
  deletions are unsafe. The typed union has 2,100 rows and leaves 20,547,
  maximum endpoint 590 on 13 rows. No physical entry or LRC(14) follows.
source: root + endpoint592_capacity + endpoint592_structure / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4313-lrc14-endpoint-592-fortythree-response-size-preserving-exchange
related:
  - THM-4310-lrc14-endpoint-594-residual-layer-closure-and-single-deletion-boundary
  - THM-4302-lrc14-endpoint-596-response-minimum-four-and-size-preserving-exchange
  - THM-4282-inactive-signature-deck-surgery-endpoint-663
artifact_root: 05-knowledge/results/lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314
artifact_manifest: 05-knowledge/results/lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314/SHA256SUMS
artifact_manifest_sha256: cf5fef4ade136ec0669f26bda4d82e725049759842181daee90c52cb08667b96
primary_scripts:
  - 04-computation/lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314/endpoint591_final_carrier_audit.cpp
  - 04-computation/lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314/endpoint591_all_witness_census.cpp
  - 04-computation/lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314/endpoint591_twohit_census.cpp
  - 04-computation/lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314/derive_endpoint591_deletion_boundary.py
  - 04-computation/lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314/typed_endpoint591_consumer.py
  - 04-computation/lrc14_endpoint591_complete_layer_closure_deletion_boundary_thm4314/verify_endpoint591_packet.py
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. Capacity, all-witness, and
  protected-joint O2/O3 artifacts agree byte-for-byte. The packet verifier
  checks the complete low-multiplicity classification, typed partition,
  manifest closure, and canonical THM-4313 inputs under normal and optimized
  Python. The verifier contains no assert statements.
---

# THM-4314 -- LRC(14) endpoint-591 complete layer closure and deletion boundary

**PROVED RELATIVE TO THM-4313 + FINITE-EXACT + COMPLETE O2/O3 REPLAY AND
ALL-WITNESS CENSUS. LRC(14) REMAINS OPEN.**

## 1. Fixed carrier and layer closure

Let `C` be THM-4313's final carrier, of size 3,925, ranks 3,818/107, FNV
`a0d08a38c10bdab7`, retaining all 421 joint masks. Let `R_591` be the 13-row
endpoint-591 residual layer, ordered FNV `fc332c0697c671c7`.

Independent O2 and O3 programs reconstruct `C` from canonical ledgers and
test all

~~~text
13*binom(30,9)=185,992,950
~~~

row-body obligations. Their raw transcripts, pair ledgers, and failure
ledgers agree byte-for-byte. They find 28,791 joint-exposed bodies, 1,031,512
nonjoint hit incidences, and zero failures. The pair-ledger FNV is
`8191899a3e142a2c`. Hence `C` closes all of `R_591` without an exchange.

## 2. Exact deletion theorem through cardinality two

For an obligation `(q,591,B)`, define its witnesses to be all masks `m` in
`C` that are active at `(q,591)` and satisfy `m & B = 0`. A deletion set `D`
preserves `R_591` iff every obligation retains a witness outside `D`.

A separate complete O2/O3 all-witness census finds exactly zero obligations
with no witnesses, two with one witness, and 26 with two witnesses. The two
singleton-critical witnesses are the protected-joint rank-eight masks

~~~text
024085d0 for body 051c6208 at (96,591),
20a09640 for body 1d584001 at (96,591).
~~~

The 26 two-witness obligations induce 23 distinct witness pairs; exactly 22
avoid both critical masks. Therefore, for every `D subseteq C` with
`|D| <= 2`:

- `|D|=1` is unsafe exactly for the two critical masks, so 2 singletons are
  unsafe and 3,923 are safe;
- `|D|=2` is unsafe exactly when `D` contains a critical mask or is one of
  the 22 additional two-witness pairs, so 7,869 pairs are unsafe and
  7,692,981 are safe among all 7,700,850 pairs.

Indeed, an obligation with at least three distinct original witnesses cannot
be exhausted by at most two deletions, while each listed obstruction deletes
every witness of its displayed obligation. This proves both directions.

Retaining every joint mask therefore makes every single nonjoint deletion
safe. Exactly 12 nonjoint pairs are unsafe. Nothing here classifies deletion
sets of size at least three.

## 3. Typed consequence

Consuming `R_591` gives an exact 2,100/20,547 partition of the frozen
22,647-row universe. The union FNV is `3b2d991da091a7df`; the residual FNV is
`59ca49a11d140ec5`. The next residual endpoint is 590 on 13 rows, FNV
`44aa8a793d162cf9`.

## 4. Scope

This theorem is finite-exact for the fixed labelled carrier and fixed row
layer. It proves neither global carrier optimality, endpoint-590 closure,
larger-deletion robustness, a physical entry/owner map, terminating descent,
nor LRC(14).

Full ledgers, hashes, reproduction commands, and the typed successor are in
the artifact packet.
