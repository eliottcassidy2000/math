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
