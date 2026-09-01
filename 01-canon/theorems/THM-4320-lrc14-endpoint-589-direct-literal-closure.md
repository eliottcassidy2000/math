---
id: THM-4320
title: "LRC(14) endpoint-589 direct literal closure"
status: >
  PROVED RELATIVE TO THM-4318/4242 + FINITE-EXACT + INDEPENDENTLY AUDITED.
  The unchanged 3,925-mask carrier leaves exactly 20,036 failures on the
  twenty-eight endpoint-589 rows, but every failure has a strict direct
  literal Haar certificate. Thus the whole endpoint-589 layer closes without
  carrier surgery, the fixed-fifty full-pool ray extends from r>=590 to
  r>=589, and the typed union advances 2,113->2,141 with next endpoint 588.
  No physical entry, terminating descent, or LRC(14) follows.
source: root + endpoint592_structure + endpoint590_structure_theory / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4318-lrc14-endpoint-590-exact-nine-response-size-preserving-exchange
  - THM-4242-fixed-fifty-direct-r590-tail-and-twenty-three-label-chart
related:
  - THM-4319-lrc14-endpoint-590-fortyone-body-response-critical-core
  - THM-4234-fixed-fifty-twenty-label-pair-haar-charts
artifact_root: 05-knowledge/results/lrc14_endpoint589_direct_literal_closure_thm4320
artifact_manifest: 05-knowledge/results/lrc14_endpoint589_direct_literal_closure_thm4320/SHA256SUMS
artifact_manifest_sha256: fb3f9158bbac86fedf0887d670b25d0c34ba6e893a877315f6afba24604be188
primary_scripts:
  - 04-computation/lrc14_endpoint589_direct_literal_closure_thm4320/endpoint589_exchanged_carrier_audit.cpp
  - 04-computation/lrc14_endpoint589_direct_literal_closure_thm4320/endpoint589_direct_literal_primary.cpp
  - 04-computation/lrc14_endpoint589_direct_literal_closure_thm4320/endpoint589_literal_lower_bound_independent.cpp
  - 04-computation/lrc14_endpoint589_direct_literal_closure_thm4320/verify_endpoint589_direct_literal_closure_packet.py
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. The 400,600,200-case
  carrier audit agrees byte-for-byte under O2/O3. Self-contained C++,
  clean-room Python, and a separate raw-cell/class-aggregation C++ replay
  agree on all 20,036 direct certificates; normal/optimized typed derivations
  and the hardened verifier are byte-identical. The independent structure
  audit also checks the labelled hub fibres and trivial automorphism group.
---

# THM-4320 -- LRC(14) endpoint-589 direct literal closure

**PROVED RELATIVE TO THM-4318/4242 + FINITE-EXACT. LRC(14) REMAINS OPEN.**

## 1. Statement

Let

~~~text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}
~~~

and let `E_589` be the twenty-eight rows in the frozen endpoint-589 ledger
whose ordered FNV64 is `5d9429c9f9971322`. For every `(q,589) in E_589` and
every `B in binom(P,9)`,

~~~text
mu(G_(B union {q,589})) >= 4/63.                       (1)
~~~

No change to THM-4318's 3,925-mask carrier is needed. In particular, combining
the row `q=50` with THM-4242 gives the complete fixed-fifty pool ray `(1)` for
every integer `r>=589`.

The closest proved mechanism is THM-4318's exact response exchange. Its
canonical hostile is the present carrier scan: a response carrier can fail
even when the literal target mass is positive. The corrected near miss is to
distinguish a missing rank-eight/rank-nine response from a failed Haar
inequality. The least-used decisive sidecar is the unprojected literal wall
mass, including cells whose pool-failure set has rank above nine.

The live concept board was

~~~text
carrier response | literal wall mass | failure hypergraph | fixed-50 petals
typed frontier | symmetry quotient.                              (2)
~~~

## 2. Exact carrier exit

The unchanged endpoint-590 exchange has 3,925 masks: 3,809 of rank eight and
116 of rank nine, with all 421 protected joint masks retained. The exact
O2/O3 replay checks

~~~text
28 binom(30,9) = 400,600,200
~~~

pair/body cases. It leaves 20,036 failures on exactly two rows:

~~~text
(50,589): 20,025, failure FNV ff421454f02d9099;
(96,589):     11, failure FNV 6f70a2d28ff6a957.       (3)
~~~

The other twenty-six rows have none. On the complement of `(3)`, the replay
provides a protected-joint or active-nonjoint mask disjoint from the body;
activity and monotonicity give `(1)`. It remains only to treat `(3)`.

## 3. Direct literal inequality

Fix one of the hostile pairs `(q,589)`. Partition the circle by all walls of
`P union {q,589}`. On a pair-safe open cell let `F subseteq P` be the labels
that fail there and let `w(F)` be its total width on the exact integer grid.
For a body `B`, its literal safe-set mass in grid units is

~~~text
M_q(B)=sum_(F intersect B=empty) w(F).                 (4)
~~~

Define the truncated mass

~~~text
L_q(B)=sum_(F intersect B=empty, |F|<=9) w(F).         (5)
~~~

Every discarded width is nonnegative, so `L_q(B)<=M_q(B)`. Endpoint-only
cells have Haar measure zero. Therefore

~~~text
63 L_q(B)-4D_q > 0                                    (6)
~~~

implies the strict version of `(1)`, where `D_q` is the wall grid.

Exact evaluation of every body in `(3)` gives

| row | grid `D_q` | low classes | minimum in `(6)` | minimizing body |
|:---|---:|---:|---:|:---|
| `(50,589)` | 2,827,379,709,554,400 | 2,383 | 14,566,818,763,788,984 | `013c6401` |
| `(96,589)` | 1,130,951,883,821,760 | 2,352 | 7,172,391,058,639,758 | `0d0c6401` |

Thus every one of the 20,036 carrier failures is strictly positive. Together
with Section 2 this proves `(1)` for the whole layer.

The independent raw-cell audit also computes the exact full mass `(4)`. At
`q=50`, full and truncated masses agree for 16,788 bodies and the full mass
strictly dominates for 3,237; at `q=96` the split is `4+7`. The minima of the
full masses occur at the same displayed bodies and equal the truncated
minima. These facts are checks of the lower-bound mechanism, not assumptions.

## 4. Fixed-fifty bridge and labelled structure

Use THM-4234's split `P=C union O`, with `|C|=18` and `|O|=12`. By
`k=|B intersect O|`, the 20,025 `q=50` carrier failures split as

~~~text
k:      0     1     2     3     4    5
count: 870  4519  7803  5502  1293   38.              (7)
~~~

The first 18,694 bodies lie in the inherited zero-through-three-petal chart
layers. The direct audit handles all 1,331 bodies with `k>=4`; their minimum
truncated surplus is `15,777,555,364,138,176` at `20744601`. This proves the
new `r=589` row needed to extend THM-4242's `r>=590` ray. It does not assert
casewise novelty against every selected chart in the earlier atlas.

The q=50 failure hypergraph has full support, empty common intersection, and
thirty distinct vertex degrees. Hence its label-permuting automorphism group
is trivial. Labels 95 and 193 form dominant hubs:

~~~text
neither 7; only 193: 1,347; only 95: 1,217; both: 17,454. (8)
~~~

The seven neither-hub bodies all contain labels 80 and 168. Exactly 5,651 of
the 20,025 bodies are inclusion-minimal transversals of the 1,398-mask active
q=50 carrier clutter. Finally, `gcd(95,589)=gcd(190,589)=19`, while their
failure degrees are 18,671 and 3,465: gcd data alone does not recover the
labelled structure. These are exact diagnostics, not dependencies of `(1)`.

## 5. Typed successor and scope

Consuming the twenty-eight rows gives

~~~text
typed union: 2,141, FNV c84bb7f7eaa0f230;
residual:   20,506, FNV 3cd0863a93c7602e;
next top: endpoint 588, 66 rows, FNV 18cf9a572cf9a5be. (9)
~~~

The partition remains disjoint and exhaustive in the frozen 22,647-row typed
universe. Equation `(9)` is a proof-graph ledger advance, not physical entry.

The theorem is finite-exact for the displayed pool, carrier, typed universe,
and endpoint layer. It proves no arbitrary-pair theorem, carrier termination,
physical owner/entry map, or LRC(14). **QED.**
