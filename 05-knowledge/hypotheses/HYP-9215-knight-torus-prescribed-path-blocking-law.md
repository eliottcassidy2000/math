---
id: HYP-9215
title: "Prescribed-path blocking law on the knight torus: for n >= 6 and every path P of at most 3 moves on the n x n knight torus G_n, the fewest moves (disjoint from P) whose deletion kills every closed tour through P equals the cheapest one-square obstruction lambda_P (7, 6, 6 for 1, 2, 3 moves), and every minimum blocking set is local (a starved square or an early closing)"
status: >
  OPEN. FINITE-EXACT (session reader, C enumeration + stored-cycle pools + CP-SAT, not re-run here):
  beta_P = lambda_P for every class of 1-move paths at n = 6 and every class of 2- and 3-move paths at n = 5, 6, 7, 8;
  minimum blocking sets are exactly the local ones for three sampled classes at n = 6 (43, 45, 84 sets);
  every path of <= 5 moves lies in a closed tour of G_n, n = 5..8. PROVED: beta_P <= lambda_P with the stated lambda values for n >= 6.
source: mac-mini-2026-10-06-oaimath2, 05-knowledge/results/oai2_openai_math_second_reading_20261006.md (section on openai/math #180, Barnette)
related:
  - HYP-9211 (the 0-move case: beta = 7 by stars only); THM-4552 (G_n, restricted edge connectivity 14)
  - THM-4550 (the weave law: the first non-local obstruction to extending paths on the 8 x 8 board)
  - openai/math #180 (Barnette's conjecture, unrefereed): Cor 7.1 / 7.2 are the cubic (k = 3) case of the same law
---

# HYP-9215 — the knight torus blocks prescribed paths only locally

**Definitions.**
* `G_n` is the `n × n` knight torus (8-regular for `n ≥ 5`).
* For a path `P`, `β_P` is the least number of edges disjoint from `P` whose deletion leaves no Hamiltonian cycle through `P`.
* `λ_P` is the least cost of a one-square obstruction. Either a square is starved to at most one usable move (moves into `P`'s interior do not count), or `P` is forced to close early at a common neighbour of its ends.

**Conjecture.** For `n ≥ 6` and every `P` with at most 3 moves:
* `β_P = λ_P`, with `λ_P = 7, 6, 6` for 1, 2, 3 moves;
* every minimum blocking set is local.

**Why this is the right shape.**
* In a cubic graph (`k = 3`), #180's corollaries say `β_∅ = k − 1 = 2` and `β_(P_4) = k − 2 = 1`. That is the same local law with `k = 3`.
* HYP-9211 is `β_∅ = k − 1 = 7` at `k = 8`.
* "Stars only" is not Barnette-like. In a 14-vertex Barnette graph with a nontrivial 3-edge cut, 15 of the 57 minimum blocking pairs are not stars. So the locality here is a consequence of the knight torus's restricted edge connectivity 14 (THM-4552 v), not of cubicity or planarity.

**Evidence.** See the status line. The enumerations covered `4.48·10^8` sets per class at `n = 6` and `8.47·10^9` at `n = 8`, and none blocked below `λ_P`.

**What would prove it.**
* A 2-factor step: an analogue of #180's Proposition 3.1, which produces states by planar density and flows.
* Plus a merging step that joins the 2-factor's cycles into one tour through `P`. This step is open; #180's disk identity needs planarity, and `genus(G_n) ≥ 1 + n²/2`.
