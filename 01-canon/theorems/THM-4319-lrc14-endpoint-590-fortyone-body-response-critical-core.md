---
id: THM-4319
title: "LRC(14) endpoint-590 forty-one-body response-critical core"
status: >
  PROVED RELATIVE TO THM-4318 + FINITE-EXACT + INDEPENDENT O2/O3 SEARCHES.
  A 41-body restriction of the endpoint-590 response hypergraph is
  9-cover-critical: its cover number is nine, while deleting any one body
  lowers the cover number to eight. Its incompatibility graph has chromatic
  number four, so the missing five units are genuinely higher-order response
  constraints. No new carrier exchange, physical entry, or LRC(14) follows.
source: root + endpoint590_structure_theory + endpoint592_structure / LRC14 continuation session, 2026-09-01
depends_on:
  - THM-4318-lrc14-endpoint-590-exact-nine-response-size-preserving-exchange
related:
  - THM-4313-lrc14-endpoint-592-fortythree-response-size-preserving-exchange
  - THM-4305-lrc14-endpoint-response-capacity-landscape-and-deletion-obstruction
artifact_root: 05-knowledge/results/lrc14_endpoint590_fortyone_body_response_critical_core_thm4319
artifact_manifest: 05-knowledge/results/lrc14_endpoint590_fortyone_body_response_critical_core_thm4319/SHA256SUMS
artifact_manifest_sha256: 3b8e55fcc95d93fdd9456ab80f7dd442f72fb9ada7aef9c04d4653a8f4105e51
primary_scripts:
  - 04-computation/lrc14_endpoint590_fortyone_body_response_critical_core_thm4319/endpoint590_obstruction41_no8.cpp
  - 04-computation/lrc14_endpoint590_fortyone_body_response_critical_core_thm4319/build_core41_ledgers.py
  - 04-computation/lrc14_endpoint590_fortyone_body_response_critical_core_thm4319/verify_endpoint590_core41_packet.py
audit: >
  PASS / ACCEPT, subject to the frozen SHA256SUMS. O2/O3 exact-search
  transcripts agree byte-for-byte. A second agent rebuilt the search and
  independently checked the 328-row deletion-witness ledger, direct
  nine-cover, dual support, saturation identity, clique censuses, proper
  four-coloring, and manifest closure. Normal and optimized Python verifier
  outputs are byte-identical and contain no assert dependency.
---

# THM-4319 -- LRC(14) endpoint-590 forty-one-body response-critical core

**PROVED RELATIVE TO THM-4318 + FINITE-EXACT. LRC(14) REMAINS OPEN.**

## 1. Critical core

Let `F=(f_0,...,f_99)` be THM-4318's ordered endpoint-590 failure ledger and
let `R` be its complete family of 14,368 distinct nonempty response signatures
realized by active rank-eight/rank-nine masks. For `X subseteq {0,...,99}`,
write `tau(X)` for the minimum number of signatures in `R` whose union
contains `X`.

Set

~~~text
U = {0,5,9,12,14,18,20,22,24,25,29,32,35,37,38,40,43,47,48,
     53,55,57,63,65,68,71,73,76,77,79,80,83,84,88,89,90,94,
     95,96,97,99}.
~~~

Then

~~~text
tau(U) = 9,
tau(U \ {u}) = 8 for every u in U.
~~~

Thus `U` is an inclusion-minimal nine-response obstruction: all 41 bodies are
essential to its cover number.

## 2. Exact proof

Restriction maps a full response `S` to `S intersect U`. It yields 2,041
distinct nonempty traces. Quotienting by inclusion leaves 395 maximal traces.
This quotient preserves `k`-cover existence because any discarded trace can
be replaced by a containing trace.

A deterministic search branches on every maximal trace containing a selected
uncovered pivot. Local dominated gains are removed only in favor of a
containing gain; the sum and integer-dual prunes are necessary conditions.
Independent O2/O3 builds agree on

~~~text
U at depth 8:
  nodes 1,968,373; dead states 1,698,303;
  sum prunes 530,603; dual prunes 1,009,466; UNSAT.

all 41 sets U \ {u} at depth 7, aggregated:
  nodes 278,730; dead states 270,531;
  sum prunes 52,383; dual prunes 198,392; 41 UNSAT.
~~~

A direct nine-response ledger covers `U`. A separate 328-row ledger supplies
eight distinct atlas representatives for every `U \ {u}`; each union is
checked to equal that exact deleted core on the core coordinates. These upper
witnesses and exhaustive lower searches prove both equalities. The numerical
optimizer used to discover the witnesses is not a proof dependency.

## 3. The missing unit beyond the dual

THM-4318's integer dual has a 21-point support `D` contained in `U`, total
weight 22, and response load at most three. The other 20 core points are
essential satellites. Among the 395 maximal traces, the dual-load distribution
is

~~~text
load 0: 13,  load 1: 88,  load 2: 177,  load 3: 117.
~~~

For a hypothetical eight-cover, let `l_j` be response `j`'s dual load, let

~~~text
Delta = sum_j (3-l_j),
Omega = sum_{i in D} w_i(c_i-1),
~~~

where `c_i` is the number of selected responses covering dual point `i`.
Double-counting dual weight gives

~~~text
sum_j l_j = 22+Omega,
Delta+Omega = 24-22 = 2.
~~~

Both terms are nonnegative integers. Hence at least six of the eight responses
would have dual load three and their weighted overlap on `D` would be at most
two. The exact search proves that the 20 satellites exclude every such
near-partition. This isolates the last unit that the dual bound alone misses.

## 4. Pairwise information loses the obstruction

Join two failure bodies when no realized response covers them together. This
is an intrinsic incompatibility graph, not a tournament: compatible pairs are
ties, and orienting them would add information not present in the response
relation.

Finite-exact censuses give

~~~text
full 100-point graph: 740 edges, clique number 4, 50 maximum cliques;
core 41-point graph:  175 edges, clique number 4,  7 maximum cliques.
~~~

An explicit proper four-coloring of the core, checked edge by edge, combines
with a four-clique to prove its chromatic number is exactly four. Thus the
complete pairwise-incompatibility abstraction permits four pairwise-compatible
classes, whereas the actual response hypergraph requires nine responses. The
five-unit gap is higher-order. No chromatic-number claim is made for the full
100-point graph.

## 5. Geometry and scope

Every core body avoids labels `1,2,4,5,8`; label 26 occurs in 40 of 41 bodies,
label 20 in 36, and label 12 in 31. These are exact frequencies, not asserted
symmetries: restricting to failure bodies forgets mask rank, multiplicity,
margin geometry, the other 59 failures, and every deletion-safety coordinate.

The theorem concerns only the fixed THM-4318 response hypergraph. It proves no
new carrier exchange, endpoint-589 statement, physical owner/entry map,
terminating descent, or LRC(14).
