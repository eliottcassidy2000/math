---
id: THM-4321
title: "LRC(14) endpoint-589 q=50 response nineteen-channel obstruction"
status: >
  PROVED RELATIVE TO THM-4320 + FINITE-EXACT + INDEPENDENTLY AUDITED.
  After the fixed q=96 two-mask peel, the 19,134-body q=50 residual response
  hypergraph has cover number between 19 and 95. The lower certificate is a
  nineteen-body disjoint-responder packing; the upper certificate is an exact
  95-mask cover. The full 20,025-body q=50 family has bounds 19..96. The
  95-cover is also a two-for-one local minimum. This diagnostic is not the
  endpoint-589 closure, an optimum computation, a physical entry, or LRC(14).
source: root + endpoint592_capacity + endpoint590_structure_theory / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4320-lrc14-endpoint-589-direct-literal-closure
related:
  - THM-4318-lrc14-endpoint-590-exact-nine-response-size-preserving-exchange
  - THM-4319-lrc14-endpoint-590-fortyone-body-response-critical-core
artifact_root: 05-knowledge/results/lrc14_endpoint589_q50_response_channels_thm4321
artifact_manifest: 05-knowledge/results/lrc14_endpoint589_q50_response_channels_thm4321/SHA256SUMS
artifact_manifest_sha256: f67f749ccf380eaa05c901ae7e418a893dca19a881d1cf311630e6a93a1d2574
primary_scripts:
  - 04-computation/lrc14_endpoint589_q50_response_channels_thm4321/audit_capacity_residual_19_95.cpp
  - 04-computation/lrc14_endpoint589_q50_response_channels_thm4321/q50_cover_two_for_one_direct.py
  - 04-computation/lrc14_endpoint589_q50_response_channels_thm4321/verify_packet.py
audit: >
  PASS / ACCEPT for the frozen manifest. A no-import C++ replay rebuilds the
  wall geometry, streams all 20,160,075 rank-eight/rank-nine masks for the
  lower certificate, and directly checks all 95 upper-certificate masks.
  O2/O3 transcripts agree byte-for-byte. A separate direct Python audit checks
  every possible 2->1 exchange normally and under -O. The bracket, unlike the
  local-minimum sidecar, received an independent clean-room replay.
---

# THM-4321 -- LRC(14) endpoint-589 q=50 response nineteen-channel obstruction

**PROVED RELATIVE TO THM-4320 + FINITE-EXACT. LRC(14) REMAINS OPEN.**

## 1. Statement and representation map

Use THM-4320's 20,025 rank-nine carrier-failure bodies at `(50,589)`. The
fixed q=96 peel masks are

~~~text
0000b3a5, 0220932c.                                    (1)
~~~

The first is inactive at q=50. The second is active and covers exactly 891
q=50 failures. Let `U` be the remaining 19,134 bodies. For an active
rank-eight or rank-nine pool mask `R`, put

~~~text
E_R={B in U : R intersect B is empty}.                 (2)
~~~

If `tau(U)` denotes the minimum number of edges `(2)` covering `U`, then

~~~text
19 <= tau(U) <= 95.                                    (3)
~~~

For the unpeeled 20,025-body family, the corresponding bounds are `19..96`.
Under the architecture that keeps both q=96 peel masks, the total number of
endpoint additions lies in `21..97`. These are bounds, not an optimum.

The map is explicit:

~~~text
source:  labelled endpoint-589 carrier failures;
map:     active response mask -> its disjoint residual bodies;
target:  a finite set-cover hypergraph;
kept:    sufficiency of a response cover and packing duality;
lost:    q=96 coupling, exchange behavior, higher-rank responses,
         literal-mass magnitude, and physical owner/entry data;
sidecar: THM-4320's unprojected literal wall mass.       (4)
~~~

The lost literal coordinate is decisive: THM-4320 already closes endpoint
589, whereas `(3)` says that this sufficient-certificate representation is
expensive.

## 2. Nineteen-channel lower certificate

The frozen packing contains nineteen residual bodies:

~~~text
013c6401 071c6400 052c6401 07186401 031c5600
27046401 23165400 27146008 171a5000 15107401
35105401 07485401 31106601 0518f008 236c4001
055c4408 2504f400 23943001 14386401.                   (5)
~~~

The independent audit reconstructs the 2,383 low-rank q=50 wall classes and
streams all

~~~text
binom(30,8)+binom(30,9)=5,852,925+14,307,150           (6)
~~~

masks. Exactly 480,409 rank-eight and 4,112,383 rank-nine masks are active.
No active response is disjoint from two bodies in `(5)`: the maximum load is
one. Hence unit weight on the nineteen bodies is an integral dual certificate,
and every response cover of `U` has size at least nineteen.

There are 172 active rank-eight and 17,049 active rank-nine responders hitting
one of the packed bodies. They form nineteen mandatory channels. A complete
cover must choose at least one response from each channel; responders outside
the channels cannot discharge those obligations.

## 3. Ninety-five-mask upper certificate

The frozen artifact gives 95 distinct active rank-nine masks. Direct replay
against `U` covers all 19,134 bodies; its ordered mask FNV is
`fa44f9bfad76cfe7`. The smallest exact activity surplus among its masks is

~~~text
1,550,209,054,968 ticks at 0a624049.                    (7)
~~~

This proves the upper half of `(3)`. The active peel mask `0220932c` is not
one of the 95, and the 95 alone do not cover the 891 peeled bodies. Adjoining
it proves the full-family upper bound 96. Adjoining the two fixed q=96 masks
to the residual cover gives the architecture-specific upper bound 97.

The 95-mask cover is an exact two-for-one local minimum. For each of its
`binom(95,2)=4,465` pairs, remove the pair and union the bodies thereby left
uncovered. Any one-mask replacement must be a rank-eight/rank-nine subset of
the complement of this union. Across all pairs there are 40,468 such possible
masks, counted with pair multiplicity. Exact wall activity finds none active.
Thus no `2 -> 1` response exchange improves this incumbent. This does not rule
out a coordinated larger exchange or a smaller unrelated cover.

## 4. Inheritance, mechanism, and scope

The closest proved mechanism is THM-4318's exact-nine response exchange. The
canonical hostile object is THM-4320's q=50 failure family: its literal Haar
inequality is strictly positive even though the inherited carrier misses all
20,025 bodies. The corrected near miss is to distinguish response-cover cost
from literal difficulty. The least-used sidecar is the actual wall mass lost
by the response quotient.

The lower mechanism is not low average coverage but an exact obstruction:
nineteen labelled bodies have pairwise disjoint responder sets. The upper
mechanism is a displayed cover, not a heuristic claim about optimality. The
local audit explains why the cheapest pairwise repair stalls, while leaving
global optimization open.

This theorem is finite-exact for the fixed pool, endpoint pair, peel,
rank-eight/rank-nine response family, and displayed certificates. It proves no
new carrier exchange, endpoint closure beyond THM-4320, physical entry,
terminating descent, or LRC(14). **QED.**
