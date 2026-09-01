---
id: THM-4323
title: "LRC(14) endpoint-587 direct literal closure"
status: >
  PROVED RELATIVE TO THM-4318/4322 + FINITE-EXACT + INDEPENDENTLY AUDITED.
  The unchanged 3,925-mask carrier leaves exactly 41 failures, all on q=50,
  across the ten endpoint-587 rows. Every failure has a strict direct literal
  Haar certificate. The layer closes without carrier surgery and the typed
  union advances 2,207->2,217 with next endpoint 586. No physical entry,
  terminating descent, or LRC(14) follows.
source: root + endpoint588_scout + endpoint587_independent / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4318-lrc14-endpoint-590-exact-nine-response-size-preserving-exchange
  - THM-4322-lrc14-endpoint-588-direct-literal-closure
artifact_root: 05-knowledge/results/lrc14_endpoint587_direct_literal_closure_thm4323
artifact_manifest: 05-knowledge/results/lrc14_endpoint587_direct_literal_closure_thm4323/SHA256SUMS
artifact_manifest_sha256: 1459e35062ac7b604f4eb069da6ac935914a2d5b9b8097192a8a63c3e871e237
independent_artifact_root: 05-knowledge/results/lrc14_endpoint587_independent_audit_thm4323
independent_artifact_manifest_sha256: d3462d90c317ef65fb11c3170c8ab0a732f4206bb7a945d8a80af26031be3df8
primary_scripts:
  - 04-computation/lrc14_endpoint587_direct_literal_closure_thm4323/endpoint587_exchanged_carrier_audit.cpp
  - 04-computation/lrc14_endpoint587_direct_literal_closure_thm4323/endpoint587_direct_literal_primary.cpp
  - 04-computation/lrc14_endpoint587_direct_literal_closure_thm4323/endpoint587_direct_literal_rawcell_independent.cpp
  - 04-computation/lrc14_endpoint587_direct_literal_closure_thm4323/verify_endpoint587_literal_closure_packet.py
  - 04-computation/lrc14_endpoint587_direct_literal_closure_thm4323/independent/endpoint587_cleanroom_carrier_audit.cpp
  - 04-computation/lrc14_endpoint587_direct_literal_closure_thm4323/independent/endpoint587_cleanroom_wall_audit.py
audit: >
  PASS / ACCEPT for both frozen manifests. Primary O2/O3 carrier and wall
  ledgers agree; aggregated-class and raw-cell masses agree on all 41 bodies.
  A no-import clean-room carrier scan, separate integer wall evaluator, and
  independent typed consumer reproduce all hashes and the successor. Both
  hardened verifiers pass normally and optimized.
---

# THM-4323 -- LRC(14) endpoint-587 direct literal closure

**PROVED RELATIVE TO THM-4318/4322 + FINITE-EXACT. LRC(14) REMAINS OPEN.**

## 1. Statement

Let `E_587` be the ten rows in THM-4322's frozen typed successor, ordered FNV
`f48ca5f1904d6f52`. For every `(q,587) in E_587` and every rank-nine body
`B` in the thirty-label pool,

~~~text
mu(G_(B union {q,587})) >= 4/63.                        (1)
~~~

The THM-4318 carrier remains unchanged.

## 2. Carrier exit and literal repair

Exact O2/O3 replay checks

~~~text
10 binom(30,9)=143,071,500                             (2)
~~~

row/body cases. Nine rows close outright. Only `(50,587)` has failures:
exactly 41 bodies, ordered body FNV `c76719ced1d5c52b`. Across all rows there
are 53,771 exposed bodies and 1,587,855 carrier-hit incidences; the global
failure FNV is `bae65dc3d3bd34d0` and pair-ledger FNV is
`0062cc50be726e54`.

On the exact q=50/r=587 wall grid, define the truncated and full masses as in
THM-4322:

~~~text
L(B)=sum_(F intersect B=empty, |F|<=9) w(F) <= M(B).    (3)
~~~

The grid is `53,537,802,887,368,800`; it has 8,385 open wall cells, 5,792
pair-safe cells, 2,420 full failure classes, and 2,365 low-rank classes. Every
one of the 41 bodies satisfies

~~~text
63 L(B)-4D > 0.                                        (4)
~~~

The minimum is `283,424,219,270,897,292` ticks at body `11186405`. Full mass
is positive for all 41; truncated and full masses agree on 31 and differ
strictly on ten. Aggregated-class and raw-cell outputs are byte-identical,
detail SHA `5a89b23b60fce1ec4be6d03bdd92c4a3cf3a861afe3d30226f8bfe831154c895`.
Thus `(4)` treats every carrier failure and proves `(1)`.

## 3. Independent replay and successor

The clean-room audit serializes the inherited carrier but imports no endpoint
587 audit source. Its separate carrier scan and integer wall implementation
reproduce `(2)`, the sole 41-body fibre, all ledger hashes, every detailed
mass, and the typed successor:

~~~text
typed union: 2,217, FNV e6592cbef9b616d8;
residual:   20,430, FNV 4710f750dfcf91ea;
next top: endpoint 586, 12 rows, FNV a1b617faa2e7f63f. (5)
~~~

The mechanism is the same literal exit as THM-4320/4322, but the failure
hypergraph has collapsed from 144,867 bodies on ten fibres to 41 on one. This
is a third consecutive positive layer, not a monotonicity theorem.

The result is finite-exact for the fixed pool, carrier, typed universe, and
endpoint layer. Equation `(5)` is a proof-graph ledger advance, not physical
entry. No arbitrary-pair theorem, terminating descent, owner map, or LRC(14)
follows. **QED.**
