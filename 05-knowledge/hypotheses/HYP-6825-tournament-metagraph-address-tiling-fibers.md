---
id: HYP-6825
title: Canonical merged-metagraph addresses and exact tiling fibers
status: IN PROGRESS — representation and executable audit reserved
source: codex-2026-07-14-S4
related:
  - HYP-2245
  - HYP-2989
  - HYP-3808
  - HYP-3809
  - HYP-3813
  - HYP-6815
  - THM-550
  - THM-646
---

# HYP-6825 — Canonical metagraph addresses and tiling fibers

This session is testing whether every tournament isomorphism-class node in the
merged tournament metagraph can be given an objective address, rooted at the
transitive class, and whether the tournament-tiling explorer's tilings can be
mapped exactly to and from the corresponding node fibers.

The candidate address is not raw graph distance.  It is a lexicographically
minimal colored path word, refined by recursive parent addresses, node
invariants, and a canonical adjacency-bit code.  Blue/black lines are retained
as structural letters rather than flattened into uncolored adjacency.

The exact relation under study is two-sorted:

~~~text
labelled tiling/tournament realization
  -> canonical tournament isomorphism class
  -> merged complement-orbit node,
~~~

with inverse maps returning fibers, not a fictitious unique tiling.

Known at reservation time:

- the explorer already enumerates tournaments and tiling matrices at small n;
- the merged metagraph has BLUE/BLACK/GREEN structural edges plus independent
  score, parent-fiber, H-step, parity, and multiplicity labels;
- boolean profiles cease to distinguish nodes at n=6;
- complement merging and vertex relabelling are distinct quotients; and
- past LRC work repeatedly shows that raw tournament class loses metric,
  endpoint-owner, scale, and threshold information.

Missing before promotion:

1. reconstruct the exact explorer and merged-edge schemas;
2. define and prove invariance of the canonical node address;
3. execute tiling -> class -> merged-node and node -> tiling-fiber maps;
4. quantify collisions under each controlled forgetting step;
5. state which extra LRC sidecars must decorate a node/fiber; and
6. run Tournament Analysis on candidate address coordinates and route
   carriers, with an explicit tie Hamiltonian path.

