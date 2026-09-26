---
id: THM-4486
title: "The min-max cycle density of sign strategies is a mean-payoff game value; every sign strategy of every q n +- 1 has a cycle of density >= log_(q+1) 2 (the negative-integer adversary); for 5n+1 the min-max density is 1/2, 3/7, 5/12 and then exactly 2/5 = the density of the sporadic cycle (1,3,8,4,2) for every k >= 15; stationary-law (entropy) floors cannot exceed 1/3"
status: >
  PROVED + INDEPENDENTLY AUDITED; FINITE-EXACT; CITED (Ehrenfeucht-Mycielski
  positional determinacy, used only for the existence of matching
  certificates). Setting: the q n +- 1 strategy cube (THM-4474, THM-4481).
  rho*(q,k) is the least rho_max over level-k sign strategies. Class (i) at
  level k is nonempty iff rho*(q,k) < log_q 2.
  (G) rho*(q,k) is the value of a two-player mean-payoff game: Min fixes
  the signs, and Max fixes a lift at every residue pair. Lower bounds come
  with compact certificates (a lift strategy plus an integer potential,
  checked exactly); upper bounds come with Lemma P potentials.
  (S) rho*(q,k) depends only on q mod 2^(k-1), and rho*(-q,k) = rho*(q,k).
  (N) Theorem N: for every odd q >= 3, every k and every sign strategy,
  G_sigma has a cycle of density >= log_(q+1) 2. The adversary always
  takes the top lift, i.e. plays the negative integers -u; u grows by at
  most (q+1)/2 per odd step and halves per even step. It is tight at
  q = 2^j - 1 (the free cycle 1 -> 2^(j-1) -> ... -> 1). It gives
  rho*(3,k) = 1/2, the pin by 1 -> 2 -> 1. It beats THM-4481's entropy
  floor p0 = 0.2271 exactly for q <= 19, but settles no new q, since
  log_(q+1) 2 < log_q 2.
  (5) For q = 5, rho*(5,k) = 1/2 (k <= 6), 3/7 (7..10), 5/12 (11..14) and
  2/5 for every k >= 15. The lower bound at all levels uses the potential
  Phi(u) = u^2 (Phi(3) = 8) on the negative-integer graph. The value 2/5
  is the density of the sporadic 5x+1 cycle (1,3,8,4,2) (THM-4484), which
  no sign strategy can avoid. This is the first q >= 5 whose min-max
  density is known at every level; 2/5 lies in (log_6 2, log_5 2).
  (F) Stationary-law floors are capped at 1/3. The level-2 max-halving
  strategy has stationary odd frequency exactly 1/3 and rho_max 1/2, so no
  argument through uniform-chain stationary laws (entropy, merges, gains)
  can settle q <= 7. Certified strategies close that route for q = 9, 11.
  REFUTED: rho*(q,k) >= 3/7 for all q >= 5 (it is 2/5 for q = 5); the
  q-independence of rho* beyond k = 9.
  FINITE-EXACT: both certificates for q = 5 to k = 21, q = 7 to k = 22
  (value 14/37 = 0.3784 > log_7 2), and q = 9..21 to k = 18. Hence no
  provable sign strategy exists for 7n+-1 at k <= 22, for q = 9..21 at
  k <= 18, or for any q >= 7 at k <= 11.
  OPEN: whether lim rho*(7,k) < log_7 2 (proved floor 1/3; at k = 22 still
  0.022 above); a floor above p0 valid for all q.
  Collatz is untouched: q = 3 is pinned at 1/2 < log_3 2.
source: collatz-procgen-20260922 session, floor lane (2026-09-26); audited and promoted by the session orchestrator 2026-09-26
depends_on:
  - 01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md
  - 01-canon/theorems/THM-4481-entropy-merge-law-sign-strategies.md
related:
  - 01-canon/theorems/THM-4484-free-and-sporadic-cycles.md (the free and sporadic cycles that the adversary finds)
  - 05-knowledge/hypotheses/HYP-9141-5n-plus-1-haar-closure-concentration.md
note: 05-knowledge/results/procgen_floor_20260926_density_floor.md
scripts: 04-computation/experiments/procgen_floor_20260926_{run,lib}.py, procgen_floor_20260926_game.c
script_audit: 04-computation/experiments/procgen_floor_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_floor_20260926.out
output_audit: 05-knowledge/results/procgen_floor_20260926_orchestrator_check.out
output_sha256: 34374e14665bb374c0c27c22e91d867515abb780c15a3c466084d607d8de12d5
hash_basis: raw bytes
audit: >
  The orchestrator checked, line by line and found sound, Lemmas G1/G2
  (certificates), Lemma L (least fixed point), Theorem N (the top-lift
  adversary, closure of W = [H, N), the growth bound v <= (q+1)u/2),
  Corollary 5 (the u^2 potential with the correction at 3), and
  Proposition F.
  Independent code (procgen_floor_20260926_orchestrator_check.py, written
  without reading the lane's scripts) confirms:
  * rho*(5,k) = rho*(7,k) = 1/2 for k <= 4 by exhaustion;
  * Theorem N on 260 random strategies (q = 5..13, k = 3..7);
  * Corollary 5's edge inequalities on the negative-integer graph for all
    k <= 16;
  * the 5x+1 cycle's density 2/5;
  * Proposition F's stationary law (1/3, 1/6, 1/3, 1/6) for q = 3..23.
  The lane's pipeline was re-run (708 s, 407 MB, 41 checks). Its output
  is identical up to timing fields, including every value and certificate
  hash. Script hashes match the note.
---

# THM-4486 -- the min-max cycle density of sign strategies

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_floor_20260926_density_floor](../../05-knowledge/results/procgen_floor_20260926_density_floor.md).

## 1. A game

* **The game.** Provability of a sign strategy is a question about its densest cycle (THM-4474). The least possible densest-cycle density, `rho*(q,k)`, is the value of a mean-payoff game: Min chooses the signs, Max chooses which lift each orbit takes.
* **Certificates.** A lower bound comes with Max's lift strategy and a potential, both checkable exactly. This replaces SAT searches, which stopped at `k = 8`, and reaches `k = 22`.

## 2. The adversary that plays the negative integers

* **Negative integers are what no strategy can escape.** Max always takes the top lift, so the orbit runs through the negative integers `-u`. There the map is `u -> u/2` or `u -> (qu ∓ 1)/2` reduced, and `log u` is a potential. Every cycle then satisfies `2^p <= (q+1)^a`.
* **Collatz.** For `q = 3` this is the pin `rho* = 1/2` from the cycle `{1, 2}`. That is the identity `3 + 1 = 4` once more (THM-4484 §2).
* **For `q = 5`.** The squared potential `u^2`, corrected at `u = 3`, proves that every strategy at every level has a cycle of density `>= 2/5`.
  * The bound is attained by the sporadic 5x+1 cycle `(1, 3, 8, 4, 2)`, whose shape `(5,2)` has gap `2^5 - 5^2 = 7`.
  * In the Kuratowski–Tutte reading: the obstruction that fixes the provability floor of 5n+1 is its sporadic cycle, just as the free cycle `{1,2}` fixes 3n+1's.

## 3. What entropy cannot do

* **The cap.** The entropy law (THM-4481) bounds `rho_max` through stationary odd frequencies. Those can be as low as `1/3` while `rho_max = 1/2` (the max-halving strategy).
* **Why it matters.** The obstruction is a worst-cycle phenomenon, invisible to averages: the discrete (cycle) side cannot be reached from the continuous (stationary-measure) side beyond `1/3`.
* **What remains.** Whether 7n±1 ever becomes provable is open. It is not provable at `k <= 22`; the proved floor is `1/3` and the target `log_7 2 = 0.356`.
