---
id: HYP-9215
title: "Prescribed-path blocking law on the knight torus: for n >= 6 and every path P of at most 4 moves on the n x n knight torus G_n, the fewest moves (disjoint from P) whose deletion kills every closed tour through P equals the cheapest one-square obstruction lambda_P (7, 6, 6 for 1, 2, 3 moves; 5 or 6 for 4 moves), and every minimum blocking set is local (a starved square or end, or an early closing)"
status: >
  OPEN. FINITE-EXACT (session reader: C enumeration of all (lambda-1)-sets, stored-cycle pools, CP-SAT; not re-run here):
  beta_P = lambda_P for the 1-move class at n = 5..8 (lambda = 6 at n = 5, 7 at n >= 6; up to 3.6e11 six-sets),
  every class of 2- and 3-move paths at n = 5..8, every class of 4-move paths at n = 6, 7, 8 (173, 184, 181 classes);
  every minimum blocking set is local for all 35 classes of 1- to 3-move paths at n = 6;
  every path of <= 6 moves lies in a closed tour of G_n, n = 5..8 (P7-Hamiltonicity). PROVED: beta_P <= lambda_P with
  the stated lambda values for n >= 6. Partial: 216 of the 1088 classes of 5-move paths at n = 6 have beta = lambda = 5.
  INDEPENDENT AUDIT (2026-10-06): well-posed and consistent with HYP-9211; spot check n = 6, P = (0,0)-(1,2): beta_P = lambda_P = 7
  (a CEGAR hitting-set proof over 412 tours shows no 6 moves block).
source: mac-mini-2026-10-06-oaimath2, 05-knowledge/results/oai2_openai_math_second_reading_20261006.md (section 6, openai/math #180, Barnette)
scripts: 04-computation/experiments/oai2_20261006_readers/barnette_cancellation/ (kp_enum.c, kp_certify.py, kp_minsets.py, kp_extend.py; kc_*.out, km_*.out)
related:
  - HYP-9211 (the 0-move case: beta = 7 by stars only); THM-4552 (G_n, restricted edge connectivity 14)
  - THM-4550 (the weave law: the first non-local obstruction to extending paths on the 8 x 8 board)
  - openai/math #180 (Barnette's conjecture, unrefereed): its Cor 7.1 / 7.2 (for Barnette graphs; 7.2 in the Pfaffian case) are the k = 3 case of the same law
---

# HYP-9215 — the knight torus blocks prescribed paths only locally

**Definitions.**
* `G_n` is the `n × n` knight torus (8-regular for `n ≥ 5`).
* For a path `P`, `β_P` is the least number of edges disjoint from `P` whose deletion leaves no Hamiltonian cycle through `P`.
* `λ_P` is the least cost of a one-square obstruction. Either a square (or an end of `P`) is starved of usable moves, or `P` is forced to close early at a common neighbour of its ends. A move is useless if it enters `P`'s interior or closes `P` early.

**Conjecture.** For `n ≥ 6` and every `P` with at most 4 moves:
* `β_P = λ_P`, where `λ_P = 7, 6, 6` for 1, 2, 3 moves and `λ_P ∈ {5, 6}` for 4 moves;
* every minimum blocking set is local.

**Why this is the right shape.**
* In a Barnette graph (cubic, bipartite, planar, 3-connected; `k = 3`), #180's corollaries (CONDITIONAL on #180) give `β_∅ = k − 1 = 2`, and in the Pfaffian case `β_(P_4) = k − 2 = 1`. That is the same local law with `k = 3`.
* HYP-9211 is `β_∅ = k − 1 = 7` at `k = 8`, the 0-move case. The 1-move case of this conjecture implies HYP-9211's lower half.
* "Stars only" is not Barnette-like. In the 14-vertex Barnette graph "cube with one vertex replaced by `Q_3 − v`", 15 of the 57 minimum blocking pairs are not stars (re-confirmed by the independent audit). The locality here plausibly reflects the knight torus's restricted edge connectivity 14 (THM-4552 v) rather than cubicity or planarity. This is a heuristic, not a proof.
* Every path of at most 6 moves lies in a closed tour (`n = 5..8`). This "P7-Hamiltonicity" is the degree-8 analogue of #180's P4 property.

**Contrast with the ordinary boards (FINITE-EXACT, reader).**
* On the 8 × 8 board, every non-extendable path of at most 5 moves is caught by iterated local forcing. The first genuinely global obstruction is the weave law of THM-4550 (9 moves inside rings 1–2).
* On the 10 × 10 board, every failure up to 4 moves is local.
* On the 6 × 6 board, non-local failures already appear at 4 moves (1 of 594 classes, e.g. `(0,2),(1,4),(3,5),(2,3),(3,1)`). These were not analysed.

**What would prove it.**
* A 2-factor step: an analogue of #180's Proposition 3.1, which produces states by planar density and flows.
* Plus a merging step that joins the cycles of the 2-factor into one tour through `P`. This step is open; #180's disk identity needs planarity, and `genus(G_n) ≥ 1 + n²/2`.
