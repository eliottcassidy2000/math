---
id: THM-4324
title: "LRC(14) endpoint-586 direct literal closure"
status: >
  PROVED RELATIVE TO THM-4318/4323 + FINITE-EXACT + INDEPENDENTLY AUDITED.
  The unchanged 3,925-mask carrier leaves exactly 4,090 failures, all on
  q=50, across the twelve endpoint-586 rows. Every failure has a strict
  direct literal Haar certificate. The layer closes without carrier surgery
  and the typed union advances 2,217->2,229 with next endpoint 585. No
  physical entry, terminating descent, or LRC(14) follows.
source: root + endpoint588_scout + endpoint586_independent / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4318-lrc14-endpoint-590-exact-nine-response-size-preserving-exchange
  - THM-4323-lrc14-endpoint-587-direct-literal-closure
artifact_root: 05-knowledge/results/lrc14_endpoint586_direct_literal_closure_thm4324
artifact_manifest: 05-knowledge/results/lrc14_endpoint586_direct_literal_closure_thm4324/SHA256SUMS
artifact_manifest_sha256: 1ada30d13656edd7ac1b59e11e23525219fd1a06aa205406eb6e87a24b0fc492
independent_artifact_root: 05-knowledge/results/lrc14_endpoint586_independent_audit_thm4324
independent_artifact_manifest_sha256: 0a617d34217c78adc96e26fdf9e18bd152407aaecc184671a5af2eec8627b340
primary_scripts:
  - 04-computation/lrc14_endpoint586_direct_literal_closure_thm4324/generic_endpoint_exchanged_carrier_audit.cpp
  - 04-computation/lrc14_endpoint586_direct_literal_closure_thm4324/generic_endpoint_direct_literal_primary.cpp
  - 04-computation/lrc14_endpoint586_direct_literal_closure_thm4324/generic_endpoint_direct_literal_rawcell_independent.cpp
  - 04-computation/lrc14_endpoint586_direct_literal_closure_thm4324/verify_endpoint586_literal_closure_packet.py
  - 04-computation/lrc14_endpoint586_direct_literal_closure_thm4324/independent/cleanroom_carrier_audit.cpp
  - 04-computation/lrc14_endpoint586_direct_literal_closure_thm4324/independent/cleanroom_literal_mass.py
audit: >
  PASS / ACCEPT for both frozen manifests. Primary O2/O3 carrier and wall
  ledgers agree; aggregated-class and raw-cell masses agree on all 4,090
  bodies. A no-import clean-room carrier scan, separate integer wall
  evaluator, and independent typed consumer reproduce every hash and the
  successor. Both hardened verifiers pass normally and optimized.
---

# THM-4324 -- LRC(14) endpoint-586 direct literal closure

**PROVED RELATIVE TO THM-4318/4323 + FINITE-EXACT. LRC(14) REMAINS OPEN.**

## 1. Statement

Let `P` be the thirty-label pool of THM-4323 and let `E_586` be the twelve
rows in its frozen typed successor, ordered FNV `a1b617faa2e7f63f`. For every
`(q,586) in E_586` and every rank-nine body `B in binom(P,9)`,

```text
mu(G_(B union {q,586})) >= 4/63.                        (1)
```

The THM-4318 carrier remains unchanged. The inheritance mechanism is the
literal exit of THM-4320/4322/4323; the hostile control is the surviving
`q=50` fibre, and the decisive sidecar is exact wall mass rather than the
projected carrier response.

## 2. Carrier exit and literal repair

Exact O2/O3 replay checks

```text
12 binom(30,9)=171,685,800                              (2)
```

row/body cases. Eleven rows close outright. Only `(50,586)` has failures:
exactly 4,090 bodies, ordered body FNV `14ce094f4ab4ba94`. Across all rows
there are 525,048 exposed bodies and 12,125,735 carrier-hit incidences; the
global failure FNV is `ffb884b2b17e6ef4` and pair-ledger FNV is
`46a2c17caecc55df`.

On the exact `q=50,r=586` wall grid, put

```text
L(B)=sum_(F intersect B=empty, |F|<=9) w(F) <= M(B),    (3)
```

where `F` is a pool-failure class and `M(B)` is the complete literal-safe
mass. The grid is `26,723,298,545,143,200`; it has 8,381 open wall cells,
5,802 pair-safe cells, 2,424 full failure classes, and 2,371 truncated
classes. Every failed body satisfies

```text
63 L(B)-4D > 0.                                        (4)
```

The minimum is `136,976,666,519,138,544` ticks at body `011c6405`. Full
mass is positive on all 4,090 bodies; truncated and full masses agree on
3,602 and differ strictly on 488. Aggregated-class and raw-cell evaluations
are byte-identical, detail SHA
`6dbfc6178c7da5c844abc95e980f2696cf02c3845a83bdfff95f881c12c9c4fd`.
Together with the carrier witnesses on every nonfailure, this proves `(1)`.

## 3. Independent replay and successor

The clean-room audit reconstructs the twelve-row frontier from the prior
partition and reads only a serialized inherited carrier. Its import-free
O2/O3 carrier scan, separate unaggregated integer wall evaluator, and typed
consumer reproduce `(2)`, the sole failure fibre, every ledger hash and mass,
and the successor:

```text
typed union: 2,229, FNV 035ebf12f02ecc62;
residual:   20,418, FNV 89b73e31224821c4;
next top: endpoint 585, 23 rows, FNV 8f1b7c8db8fd5e87. (5)
```

This is a fourth consecutive literal-exit layer, not an endpoint
monotonicity theorem. The result is finite-exact for the frozen pool,
carrier, typed universe, and endpoint layer. Equation `(5)` is a proof-graph
ledger advance, not physical entry. No arbitrary-pair theorem, terminating
descent, owner map, or LRC(14) follows. **QED.**
