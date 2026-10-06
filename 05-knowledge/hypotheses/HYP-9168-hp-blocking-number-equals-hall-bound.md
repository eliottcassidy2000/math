---
id: HYP-9168
title: "The HP-blocking number of a tournament equals its Hall-deficiency bound: the fewest arcs whose deletion leaves no Hamiltonian path is min over X (|X| >= 2) and Y (|Y| <= |X| - 2) of the number of arcs entering X from outside Y (or, dually, leaving X to outside Y)"
status: >
  OPEN. The inequality beta <= hall is PROVED: after the deletion, the
  vertices of X can have predecessors only in Y, but an HP needs |X| - 1
  distinct predecessors (THM-4524).
  FINITE-EXACT:
  - every class N <= 9 (the lane's exact subset tables for N <= 7 and
    branch and bound for N = 8, 9; the orchestrator's audit independently
    covers N <= 7);
  - the two-source/two-sink special case sigma(T) first fails at N = 7
    (beta = hall = 3 < sigma = 4);
  - the hall < sigma classes number 1 at N = 7, 13 at N = 8 and 180 at
    N = 9.
  Corollary if true: every regular tournament survives the deletion of any
  N - 2 arcs (FINITE-EXACT N <= 9).
  UPDATE 2026-10-06 (mac-mini chessboard-weave session, hp_blocking lane):
  - FINITE-EXACT N = 10: beta = hall for every one of the 9,733,056 classes
    (gentourng 10; C branch and bound with HP-guided branching, protection and a
    packing lower bound, validated against brute force on all classes N <= 7 and on
    5000 damaged digraphs; 2440 s). Hall histogram at N = 10:
    1:13727 2:203041 3:1161890 4:3007400 5:3434830 6:1664367 7:234468 8:13333;
    hall < sigma (two sources/sinks not optimal) in 3368 classes. Independent
    cross-check (orchestrator; different code and algorithm: Konig-form hall plus
    SAT with lazy Hamiltonian-path cuts): beta = hall on 111 sampled N = 10 classes,
    11 of them with hall = 8.
  - VERIFIED (random): uniform random labeled tournaments N = 11 (10000), 12 (10000),
    13 (1000), 14 (200); 300 almost-regular N = 10.
  - PROVED reformulation (Konig): hall(T) = min over A, B with |A| + |B| = N + 2 of the
    number of arcs from A to B = the fewest deletions killing every spanning
    1-path-cycle factor; so HYP-9168 says the cheapest way to kill all Hamiltonian paths
    of a tournament is to kill all path-plus-disjoint-cycles factors.
  - SCOPE (FINITE-EXACT counterexamples): the equality fails for oriented graphs that are
    not tournaments (a tournament minus one arc can have beta = 2 < hall = 3), so any
    proof must use semicompleteness.
  - Proof routes tested and failing as stated: Gutin's multipartite theorem (HP iff
    1-path-cycle factor) holds, but for some classes no minimum blocking set has
    multipartite type; the single-factor merge count (blocking all simple merges of a
    path with each cycle) falls short of hall at N = 6 (3 classes) and N = 7 (65
    classes, deficit up to 3).
  Notes: 05-knowledge/results/chessboard_weave_20261006.md section 7; scripts
  04-computation/experiments/chessboard_weave_20261006_hp_blocking_*.
source: collatz-procgen-20260922 session, selfie lane (2026-10-01), Conjecture D1; promoted with THM-4524
related:
  - 01-canon/theorems/THM-4524-selfie-tournaments-loop-gauge-arc-parity-and-the-odd-mallows-sloane-partner.md
  - 05-knowledge/results/procgen_selfie_20261001_selfie_tournaments.md
---

# HYP-9168 — the HP-blocking number equals the Hall bound

**Setting.** Let `beta(T)` be the fewest arcs of the tournament `T` whose deletion leaves a digraph with no Hamiltonian path. This answers the owner's question of how far a tournament can be shaved before it loses its last Hamiltonian path.

**The bound.** In an HP, every vertex but the first has a predecessor. Delete every arc that enters `X` from outside `Y`, where `|Y| <= |X| - 2`. Then the `|X| - 1` vertices of `X` that need predecessors can find them only in `Y`, which is too small, so no HP survives. The dual uses successors. Hence `beta <= hall`.

**The conjecture.** Equality always holds: the cheapest way to kill every HP is to starve a set of predecessors (or successors). If true, this is a König–Hall-type min-max theorem for Hamiltonian paths in tournaments.

**Evidence.**
- `beta = hall` for every class with `N <= 9`.
- The naive bound with `|X| = 2` (two sources or two sinks) fails first at `N = 7`, on the tournament `0 -> A => B -> 0` with `A`, `B` cyclic triangles. Deleting the three arcs of `A` starves `A` of predecessors.
