---
id: THM-4474
title: "The strategy cube of 3n+-1: bounded-lookahead provability and residue divergence are decided by the extreme cycle densities of a finite parity graph (vs log_3 2), not by the drift; Collatz has the maximal window [0,1] at every level"
status: >
  PROVED + INDEPENDENTLY AUDITED (Theorems A, B, C, E and Proposition F);
  FINITE-EXACT (the classification of all 65,814 strategies at levels 1-5).
  A level-k sign strategy sigma assigns +-1 to the odd residues mod 2^k,
  and T_sigma(n) = n/2 (n even), (3n + sigma(n mod 2^k))/2 (n odd). This
  is Althofer's 3n+-1 game with the sign frozen into residues; Collatz is
  sigma = +, and 3n-1 is sigma = -. The parity graph G_sigma has nodes
  Z/2^k, with edges from s to the two lifts of T(s) mod 2^(k-1). A cycle
  with a odd nodes and length p is expanding iff 3^a > 2^p, i.e. its odd
  density exceeds c = log_3 2.
  (A) Bounded-lookahead provability. T_sigma has a bounded-lookahead
  descent certificate iff its exceptional set is empty at some finite
  level, iff every cycle of G_sigma is contracting (rho_max < c).
  (B) Residue divergence. A residue-class divergence certificate exists
  iff some closed class of G_sigma has every cycle expanding
  (rho_min > c). Positive drift on a trap is necessary but not sufficient.
  (C) Drift sandwich. On every closed class, rho_min <= pi(odd) <= rho_max.
  So the drift, the quantity every heuristic uses, lies strictly inside the
  window that decides provability.
  (E) Negation. Negation x -> -x is an isomorphism of parity graphs
  sigma <-> nu sigma, with (nu sigma)(m) = -sigma(-m), so classes (i) and
  (ii), the densities, the drifts and the exceptional counts are all
  nu-invariant.
  (F) Fixed points. Class (i) forces sigma(1) = + and sigma(-1) = -, so
  that both fixed points +-1 are contracting.
  Collatz has rho_min = 0 (the loop at 0) and rho_max = 1 (the loop at -1)
  at every level, hence no residue certificate of either kind.
  Class counts at levels 2..5: (i) 1, 1, 16, 1052 and (ii) 1, 1, 3, 32.
  UPDATE 2026-09-26: Collatz's flip distance to class (i) tends to 0 at
  the sharp rate 2^(-(1-h)k) (HYP-9138 PROVED by THM-4479: flip exactly the
  undecided residues).
  The OPEN fraction stabilizes near 0.435, which is the probability of no
  extra cycle in the k = infinity random-sign model.
source: collatz-procgen-20260922 session, strategy-cube lane (2026-09-25), extending the Kuratowski lane's strategy square; audited and promoted by the session orchestrator 2026-09-25
depends_on:
  - 01-canon/theorems/THM-4471-kawasaki-fixed-point-collatz-proof-refuted.md (Banach periodic points x_w = c_w/(2^p - 3^a))
related:
  - 01-canon/theorems/THM-4470-collatz-pairing-ladder-am-fair-and-defect-blind.md (the pair-0 obstruction)
  - 05-knowledge/results/procgen_kuratowski_20260925_tait_kempe_triples.md (the strategy square)
  - 05-knowledge/results/collatz_procgen_20260922_choice_ladder.md (exceptional set, dimension h(log_3 2))
  - 05-knowledge/hypotheses/HYP-9136-provable-pairing-price-tends-to-zero.md
script: 04-computation/experiments/procgen_cube_20260925_core.py
script_engine: 04-computation/experiments/procgen_cube_20260925_engine.c
script_audit: 04-computation/experiments/procgen_cube_20260925_orchestrator_check.py
output: 05-knowledge/results/procgen_cube_20260925.out
output_audit: 05-knowledge/results/procgen_cube_20260925_orchestrator_check.out
script_sha256: fb6a1de7b945fb4dc39c4579d6c7eb0aefccd462703a51082ef519f1f50ed31a
script_engine_sha256: a29f78e351acc71b86f613e4aea10f73a2d9c6159cedb5f71f04a7e2b8b4198b
script_audit_sha256: 581be1b9fd8746d4a5a502570b54ae13ac3264a01c5dc17478680319cb546bae
output_sha256: 0d51cec223dc8050332229ea8e9a09a394bba25a56eba1d45eea20b8afcf3b47
output_audit_sha256: 1ea820f5a1ecf4b4598c6f6a31244f375240045f0f0bfc685838b49184114426
hash_basis: raw LF bytes
audit: >
  The orchestrator read the proofs of Lemmas 1-2 (Terras bijection;
  Banach periodic points plus the cycle lemma), Theorems A, B, C and E, and
  Proposition F, and found them sound.
  Independent code (procgen_cube_20260925_orchestrator_check.py, written
  without reading the lane's scripts) implements the parity graph and
  Karp's maximum cycle mean. It confirms:
  * the class-(i) counts 1, 1, 16 at levels 2, 3, 4;
  * that every class-(i) strategy there sends each n < 3000 to {1,2};
  * Proposition F on all class-(i) strategies;
  * nu-invariance of the maximum density;
  * rho_max = 1 for Collatz and for 3n-1 at levels 2-6.
  The lane's full pipeline (C engine, exact Python reference, MaxSAT
  distances) was re-run. Its output is identical except for timing lines
  (peak memory 534 MB).
---

# THM-4474 -- the strategy cube: provability is decided by extreme cycle densities

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_cube_20260925_strategy_cube](../../05-knowledge/results/procgen_cube_20260925_strategy_cube.md).

## 1. Setting

* **Strategies.** A level-`k` sign strategy is a map `sigma : {1, 3, ..., 2^k - 1} -> {+1, -1}`, and
  `T_sigma(n) = n/2` for even `n`, `(3n + sigma(n mod 2^k))/2` for odd `n`.
* **The parity graph `G_sigma`.** Its nodes are `Z/2^k`. From `s` there are two edges, to the two lifts mod `2^k` of `T(s) mod 2^(k-1)`.
* **Cycle densities.** A cycle with `a` odd nodes among `p` has odd density `a/p`. It is **expanding** iff `3^a > 2^p`, i.e. `a/p > c = log_3 2`.
* `rho_max(X)` and `rho_min(X)` are the extreme cycle densities inside a node set `X`. The closed classes are the bottom strongly connected components.
* **Terras bijection.** For `m >= k`, `T` maps each class mod `2^m` affinely and bijectively onto a class mod `2^(m-1)`. So the classes mod `2^(k+L)` correspond to the `L`-paths of `G_sigma`.
* **Banach periodic points.** Each cycle `gamma` of `G_sigma` carries exactly one periodic point `x_gamma = c_gamma/(2^p - 3^a)` in `Z_2`. If `gamma` is expanding, a rotation of it lies in the exceptional set `Bad_inf`.

## 2. The theorems

**Theorem A.** The following are equivalent:
* (a) every cycle of `G_sigma` is contracting;
* (b) `Bad_L` is empty for some `L`;
* (c) every integer `n > n_0` descends below itself within `L` steps, for some `L` and `n_0`.

When they hold, `L_min <= 2^k (1 + log(3/2)/|mu|) + 1`, where `mu = rho_max log 3 - log 2`.

*Proof sketch.*
* (a) ⇒ (b): a walk decomposes into cycles plus a simple path of at most `2^k` nodes.
* (b) ⇒ (c): take `n_0` to be the largest of the finitely many descent thresholds.
* (c) ⇒ (a): an expanding cycle's periodic point gives, at every level, a class of `Bad_L` that contains arbitrarily large integers which rise for `L` steps. ∎

**Theorem B.** A residue-class divergence certificate exists (a nonempty closed node set in which every `M`-path multiplies by more than 1) iff some closed class has `rho_min > c`. Then every integer above an explicit threshold in the trap diverges.
* Positive drift on a trap is necessary, but it is not sufficient. At level 3, `(+,+,-,+)` has drift `+0.222` per step, yet its trap contains the contracting trivial cycle `{1,2}`. There are 1, 11 and 1005 such OPEN strategies at levels 3, 4 and 5.

**Theorem C (sandwich).** On every closed class `C`, `rho_min(C) <= pi_C(odd) <= rho_max(C)`.
* *Proof:* the stationary edge flow is a circulation, hence a positive combination of simple cycles, so `pi(odd)` is a mediant of cycle densities. ∎

**Theorem E (negation).** `s -> -s` is a parity-preserving isomorphism `G_sigma ≅ G_(nu sigma)`. So every residue-level datum is `nu`-invariant. The positive cycles of `nu sigma` are the negatives of the negative cycles of `sigma`.

**Proposition F.** Class (i) forces `sigma(1) = +1` and `sigma(-1) = -1`.
* With `sigma(1) = -1`, the point 1 is an expanding fixed point.
* With `sigma(-1) = +1`, the point `-1` is one: the Collatz loop at node `2^k - 1`.

This is the cube's version of THM-4470(5)'s pair-0 obstruction. A sign strategy can avoid both expanding fixed points, whereas no periodic pairing can.

## 3. What it says about Collatz

* **The window is maximal.** Collatz (`sigma = +`) has `rho_min = 0` (the loop at 0) and `rho_max = 1` (the loop at `-1`) at every level, so it lies in neither provable class. Its drift `pi_odd = 1/2` sits in the middle of the maximal window `[0, 1]`.
* **No residue argument can settle it.** By Theorem A, no bounded-lookahead argument settles Collatz at any modulus. By Theorem E, no residue-level datum distinguishes it from its sheet partner.
* **Distance to provability.** Collatz is adjacent only to class (iii) (a new cycle) at levels `>= 3`. The fewest flips reaching class (i) are `1, 2, 2, 4, 5, 9, 14, 23` for `k = 2..9`, a Haar fraction falling from `0.5` to `0.090`. Whether it tends to 0 is OPEN, the residue-periodic analogue of HYP-9136.
* **Microcosm and macrocosm.** In the `k = infinity` random-sign model, 43.5% of fields have no extra cycle; the extra-cycle counts are `{0: 43.5%, 1: 47.6%, 2: 8.5%, 3: 0.4%}`. The OPEN fraction at levels 5–12 matches it.
  * The macrocosm is the drift, the exceptional dimension and the density window. It converges to Collatz's values and is the same on both sheets.
  * The microcosm is the signs at the few smallest integers. It decides OPEN versus extra cycles: 97.8% of the level-5 cycle witnesses have minimum `<= 100`.

## 4. Classification (FINITE-EXACT)

| level | (i) | (ii) | (iii) | OPEN |
|---|---|---|---|---|
| 2 | 1 | 1 | 1 | 1 |
| 3 | 1 | 1 | 5 | 9 |
| 4 | 16 | 3 | 115 | 123 |
| 5 | 1,052 | 32 | 35,399 | 29,069 |

* Counts are totals, lifts included. A strategy may be in both (ii) and (iii).
* Class (iii) was searched with starts `<= 131,072` and cap `2^60`. OPEN strategies were re-searched to `2^20` with no new cycle.
