---
id: THM-4508
title: "The structure of optimal qn±1 sign strategies. Parity-graph cycles are exactly the periodic points of the 2-adic map, so a finite rule is a level-independent object. Max-halving has density 1/2, attained on a full 2-shift S_inf of rational cycles in [-1/(q-4), 1/(q-4)]. A flip gains halvings only on the pattern (s,2),(s,2), which for q = 7 is the classes of ±1/3. An explicit 46-class rule is optimal for 5n±1 at every level >= 15, and an 18-class rule gives 2/5 for 7n±1. rho*(7,30) = 37/100 and rho*(7,31) >= 7/19, so 7n±1 is not provable at any level <= 31"
status: >
  PROVED + INDEPENDENTLY AUDITED (Lemma C, Lemma MH, Lemma F with
  Corollary F and the q = 5 analogue); FINITE-EXACT (the explicit rules,
  re-checked independently; the level-30/31 certificates, checked by the
  lane's separately written checker); EMPIRICAL (statistics of certified
  optima). Setting of THM-4486.
  (C) The closed walks of the parity graph G_sigma correspond one-to-one
  to the periodic points of T_sigma on Z_2 (rationals c/(2^p - q^a)),
  with the same parities. So rho_max(sigma) is the largest odd density of
  a periodic orbit, and a rule reading d bits has the same rho_max at
  every level k >= d.
  (MH) For odd q >= 5, max-halving (the sign making qx+s = 0 mod 4) has
  rho_max = 1/2 at every level, attained exactly on
  S_inf = {all accelerated valuations = 2}. S_inf is a full 2-shift with
  exactly 2^n points of period dividing n, all rationals in
  [-1/(q-4), 1/(q-4)]; the fixed points are ±1 (q = 5) and ±1/3 (q = 7).
  Every strategy with rho_max < 1/2 deviates from max-halving on every
  periodic orbit of S_inf.
  (F) Flipping the max-halving sign at x gives y = 2^(v_1-1) x_1 - s_1,
  whose valuation is an explicit function of the first three itinerary
  symbols. Relative to max-halving, a flip gains halvings exactly on the
  pattern (s,2),(s,2) (for q = 7, the classes 11, 21 mod 32 of ±1/3); it
  is neutral exactly on (s,4),(-s,.), where the orbit rejoins; it loses
  otherwise. At the fixed points the flip gives, for q = 5, the sporadic
  cycle (1,3,8,4,2), of density 2/5 = the q = 5 limit; for q = 7, a cycle
  of density 1/3.
  (Rules)
  - q = 5: max-halving flipped on 46 explicit classes (a 122-leaf tree of
    depth 15) has rho_max = 2/5 at every level k >= 15, an explicit
    level-independent optimal family. Depth 15 separates residue
    collisions: pairs of rationals on conflicting cycles, congruent mod
    2^14 but not mod 2^15. That every strategy needs level 15 is
    rho*(5,14) = 5/12.
  - q = 7: max-halving flipped on 18 explicit classes (depth 10) has
    rho_max = 2/5 at every k >= 10. Its residual obstruction is the (5,2)
    cycles of denominator 17.
  (Values) rho*(7,30) = 37/100 exactly, and 7/19 <= rho*(7,31) <= 37/100.
  So 7n±1 has no provable sign strategy at any level k <= 31
  (7^37 > 2^100, 7^7 > 2^19).
  OPEN: whether lim rho*(7,k) < log_7 2 (fits of 1/rho* land on both sides
  of log_2 7). Collatz is OPEN.
source: collatz-procgen-20260922 session, seven2 lane (2026-09-26), answering the owner's "keep going on the open problems, especially 7n±1"; audited and promoted by the session orchestrator 2026-09-26
depends_on:
  - 01-canon/theorems/THM-4486-min-max-cycle-density-game.md
  - 01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md
related:
  - 01-canon/theorems/THM-4484-free-and-sporadic-cycles.md (the sporadic cycle (1,3,8,4,2) of 5x+1)
  - 05-knowledge/results/procgen_seven_20260926_seven_n_plus_one_provability.md (levels 23-29)
  - 05-knowledge/results/procgen_floor_20260926_density_floor.md (the game, Lemmas G1, G2, L)
note: 05-knowledge/results/procgen_seven2_20260926_seven_structure.md
scripts: 04-computation/experiments/procgen_seven2_20260926_{rhomax,restrict,lean8,verify8}.c and procgen_seven2_20260926_{lib,run}.py
script_audit: 04-computation/experiments/procgen_seven2_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_seven2_20260926.out
output_sha256: a60b79a44a4e2189625d1532eb43b432f8916efc36a7199a388a005a39f349db
output_audit: 05-knowledge/results/procgen_seven2_20260926_orchestrator_check.out
output_audit_sha256: bf8cd6019e39a6b68c663ebf927a8ac7f1ccd8737b3cd8b1dd5ddbed701f7c9c
hash_basis: raw bytes
audit: >
  The orchestrator read the proofs of:
  - Lemma C (inverse branches are 2-adic contractions; induction on the
    residue depth);
  - Lemma MH (itinerary bijection; S_inf as the attractor of the
    contractions g_s(y) = (4y - s)/q);
  - Lemma F and Corollary F (case analysis on v_1, v_2, s_2, s_3),
  and found them sound.
  Independent code (procgen_seven2_20260926_orchestrator_check.py,
  written from the note's statements and THM-4486's certificate lemmas
  without reading the lane's scripts or checkers) confirms:
  - A separate numpy game solver (least fixed points of Min's and Max's
    operators on pairs) reproduces rho*(7,k) exactly at every k = 8..21:
    both least fixed points are finite at the published value. It also
    gives rho*(5,13) = rho*(5,14) = 5/12 and rho*(5,15) = rho*(5,16) = 2/5,
    and an exact Lemma G1 re-check of the q = 7, k = 16 potential.
  - The 46-class q = 5 rule (disjoint, negation-closed) has
    rho_max <= 2/5 at k = 15, 16, 17, and the periodic point 1 of density
    2/5.
  - The 18-class q = 7 rule has rho_max <= 2/5 at k = 10, 12, 14, 16, and
    the periodic point -9/17 of density 2/5.
  - Lemma MH: rho_max(MH) = 1/2 for q = 5, 7, 9, 11 at k = 8, 12, 16;
    exactly 2^n S_inf points of period dividing n (n <= 8, q = 5, 7).
  - Lemma C(b): 5555 closed walks of random strategies.
  - The Lemma F table and gain classification on 100000 random 2-adic
    integers each for q = 7 and q = 5.
  Not re-checked by the orchestrator: the level-30 and level-31
  certificates. They need potential arrays of 512 MiB and were deleted
  after checking. The lane's checker verify8, written separately from its
  engines, accepted all four. The level-30 10/27 upper certificate was
  also accepted by the seven lane's checker; the two 37/100 certificates
  and the level-31 7/19 certificate by verify8 only. The lane re-ran its
  full pipeline itself after a memory fix (28 checks, ALL CHECKS PASSED,
  1411 s, largest process 581 MiB). The orchestrator did not re-run it,
  because of memory pressure from two concurrent lanes on the 8 GB
  machine.
  Scope notes:
  - The explicit trees come from a greedy search and are not canonical.
    Leaf counts are upper bounds on description length.
  - The Fibonacci values at k = 19-21 are the passage of 1/rho* through
    phi^2. No renormalization is claimed.
---

# THM-4508 — the structure of optimal sign strategies, and 7n±1 at level 31

**PROVED + INDEPENDENTLY AUDITED; FINITE-EXACT.** Full note: [procgen_seven2_20260926_seven_structure](../../05-knowledge/results/procgen_seven2_20260926_seven_structure.md).

## 1. What a strategy is, once it is written down

THM-4486 made provability a game on the parity graph at level `k`. Lemma C removes the level:
- a closed walk of the parity graph is the same thing as a rational periodic point of the 2-adic map;
- so a sign rule that reads finitely many bits is one map, and one number `rho_max`.

"Optimal at every level" can therefore be proved by one finite object. For 5n±1 it is max-halving with 46 flipped classes: optimal from level 15 on, forever.

## 2. The skeleton and the only profitable move

**The skeleton.** Max-halving is the natural greedy rule. Its worst cycles are a full 2-shift `S_inf` of rational cycles, all of valuation 2, sitting in the small interval `[-1/(q-4), 1/(q-4)]`. Any better rule must break every one of them.

**The only profitable move.** The flip calculus says:
- a local sign flip pays off only on the pattern `(s,2),(s,2)`;
- it is free only when the orbit rejoins;
- everywhere else it loses.

**What the flips produce.**
- For 5n±1, the very first profitable flip, at the fixed point 1, already produces the optimum. It creates the sporadic cycle `1, 3, 8, 4, 2` of density `2/5`.
- For 7n±1, the first profitable flips (at `±1/3`) reach only `2/5`. The remaining improvements need flips whose descriptions grow fast (46 → 300 → 508 leaves for `2/5 → 15/38 → 7/18`).

## 3. 7n±1

`rho*(7,30) = 37/100` and `rho*(7,31) >= 7/19`, still above `log_7 2 = 0.3562`. So no bounded-lookahead proof exists for 7n±1 at any level up to 31. The mean valuation `1/rho*` has climbed to `2.703`; provability needs it above `log_2 7 = 2.807`. Whether it ever gets there is open.
