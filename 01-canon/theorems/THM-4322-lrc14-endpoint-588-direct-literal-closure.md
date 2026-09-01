---
id: THM-4322
title: "LRC(14) endpoint-588 direct literal closure"
status: >
  PROVED RELATIVE TO THM-4318/4320 + FINITE-EXACT + INDEPENDENTLY AUDITED.
  The unchanged 3,925-mask carrier leaves 144,867 failures on ten of the
  sixty-six endpoint-588 rows, and every failure has a strict direct literal
  Haar certificate. Thus the whole layer closes without carrier surgery and
  the typed union advances 2,141->2,207 with next endpoint 587. No physical
  entry, terminating descent, or LRC(14) follows.
source: root + endpoint588_scout + endpoint588_independent / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4318-lrc14-endpoint-590-exact-nine-response-size-preserving-exchange
  - THM-4320-lrc14-endpoint-589-direct-literal-closure
related:
  - THM-4321-lrc14-endpoint-589-q50-response-nineteen-channel-obstruction
artifact_root: 05-knowledge/results/lrc14_endpoint588_direct_literal_closure_thm4322
artifact_manifest: 05-knowledge/results/lrc14_endpoint588_direct_literal_closure_thm4322/SHA256SUMS
artifact_manifest_sha256: a52f03acac9b178588be54b470ed9abf3e79f9cc19aaba34aee4c70554f88211
independent_artifact_root: 05-knowledge/results/lrc14_endpoint588_independent_audit_thm4322
independent_artifact_manifest_sha256: 8a23f23dcc4cceb09abd7a6a2813fe747141ce9512b627f028e7e91dc7c47b4a
primary_scripts:
  - 04-computation/lrc14_endpoint588_direct_literal_closure_thm4322/endpoint588_exchanged_carrier_audit.cpp
  - 04-computation/lrc14_endpoint588_direct_literal_closure_thm4322/endpoint588_direct_literal_primary.cpp
  - 04-computation/lrc14_endpoint588_direct_literal_closure_thm4322/endpoint588_direct_literal_rawcell_independent.cpp
  - 04-computation/lrc14_endpoint588_direct_literal_closure_thm4322/verify_endpoint588_literal_closure_packet.py
  - 04-computation/lrc14_endpoint588_direct_literal_closure_thm4322/independent/endpoint588_cleanroom_carrier_audit.cpp
  - 04-computation/lrc14_endpoint588_direct_literal_closure_thm4322/independent/endpoint588_cleanroom_literal_audit.cpp
audit: >
  PASS / ACCEPT for both frozen manifests. The primary O2/O3 scan checks
  944,271,900 cases and two separately written wall evaluators agree on all
  144,867 failure bodies. A no-import clean-room carrier/literal replay and
  independent typed consumer reproduce every count, row fibre, ledger hash,
  mass, and successor. Normal/optimized hardened verifiers pass.
---

# THM-4322 -- LRC(14) endpoint-588 direct literal closure

**PROVED RELATIVE TO THM-4318/4320 + FINITE-EXACT. LRC(14) REMAINS OPEN.**

## 1. Statement

Let `P` be the thirty-label pool of THM-4320 and let `E_588` be the sixty-six
rows in its frozen typed successor whose ordered FNV is `18cf9a572cf9a5be`.
For every `(q,588) in E_588` and every `B in binom(P,9)`,

~~~text
mu(G_(B union {q,588})) >= 4/63.                        (1)
~~~

The 3,925-mask THM-4318 carrier is unchanged. No response addition, deletion,
or exchange is required.

The closest mechanism is THM-4320's literal exit at endpoint 589. Its
canonical hostile is the large `(50,588)`/`(210,588)` carrier-failure pair:
the carrier misses almost 144 thousand bodies there despite every literal
target being positive. The corrected near miss is again that carrier failure
does not imply Haar failure. The decisive underused sidecar is the exact wall
mass before response projection.

## 2. Exact carrier exit

The inherited carrier has 3,809 rank-eight and 116 rank-nine masks, including
all 421 protected joint masks. Exact O2/O3 replay checks

~~~text
66 binom(30,9) = 944,271,900                           (2)
~~~

row/body cases. Fifty-six rows have no failure. The other ten fibres are:

| `q` | failures | ordered body FNV |
|---:|---:|:---|
| 25 | 9 | `8fe6c7431c54c6ff` |
| 50 | 99,836 | `4a8d6e89f7aa5d75` |
| 96 | 1,130 | `5766a653739cf191` |
| 100 | 14 | `f0fdda5ad80a46b7` |
| 105 | 60 | `a6f76a994fd08192` |
| 206 | 7 | `def96e7fae5db299` |
| 210 | 43,799 | `f756dc445e2790cd` |
| 256 | 7 | `392c4d67394b2e17` |
| 420 | 1 | `664e6def29193bd0` |
| 462 | 4 | `5f993deacddbf9d3` |

There are 144,867 failures in all, aggregate FNV `6f51fa88f3b09cdc`.
Aggregate exposed bodies are 5,001,257 and active carrier-hit incidences are
76,551,030. On every nonfailure the active response and monotonicity prove
`(1)`; only these ten fibres remain.

## 3. Direct literal inequality

For a fixed hostile pair, partition the circle at all walls of `P` and the
pair. On each pair-safe cell let `F subseteq P` be its pool-failure set and
`w(F)` its exact integer-grid width. Put

~~~text
L_q(B)=sum_(F intersect B=empty, |F|<=9) w(F),
M_q(B)=sum_(F intersect B=empty) w(F).                  (3)
~~~

Since omitted widths are nonnegative, `L_q(B)<=M_q(B)`. Exact evaluation gives

~~~text
63 L_q(B)-4D_q > 0                                     (4)
~~~

for all 144,867 failures, where `D_q` is the pair's wall grid. The row minima
are:

| `q` | minimum surplus in `(4)` | minimizing body |
|---:|---:|:---|
| 25 | 3,845,468,183,436,180 | `0d1cd000` |
| 50 | 3,359,314,066,788,540 | `003c6405` |
| 96 | 1,528,699,623,231,060 | `2d14c001` |
| 100 | 3,874,870,450,697,778 | `36704001` |
| 105 | 731,279,928,922,452 | `0458d082` |
| 206 | 81,246,765,687,627,708 | `14087c01` |
| 210 | 706,240,375,488,648 | `06b82090` |
| 256 | 13,754,370,514,849,146 | `1f106001` |
| 420 | 795,987,916,404,408 | `07584088` |
| 462 | 890,738,628,297,084 | `1e582001` |

The global minimum is the q=210 value. Full mass `M_q(B)` is also positive.
The aggregated-class and raw-cell implementations agree byte-for-byte on all
truncated and full values, detail SHA
`8c2652b9cba7cfdacd4d6e908334220eb08445306fbc09db816d0b0b9c48664d`.
Together with Section 2, this proves `(1)`.

## 4. Independent replay and typed successor

The clean-room audit imports maintained code only once to serialize the
inherited carrier. Its carrier scan, wall construction, literal evaluator, and
typed consumer are otherwise separate implementations. They reproduce all 66
rows, `(2)`, every failure fibre and ledger, every direct-mass detail row, and
the successor:

~~~text
typed union: 2,207, FNV 18d067b5614cf47f;
residual:   20,440, FNV 794bd808e92e27cd;
next top: endpoint 587, 10 rows, FNV f48ca5f1904d6f52. (5)
~~~

This is a second consecutive literal-exit layer, which is evidence for the
literal-first search policy but not a uniform theorem in the endpoint.

The result is finite-exact for the displayed pool, carrier, typed universe,
and endpoint layer. Equation `(5)` is a proof-graph ledger advance, not a
physical entry. No arbitrary-pair theorem, carrier termination, physical
owner map, or LRC(14) follows. **QED.**
