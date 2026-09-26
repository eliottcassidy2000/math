---
id: THM-4482
title: "Ranks are tensions: a sign strategy is bounded-lookahead provable iff it has a rank a log n + h with h periodic (equivalently merely bounded); the least rank defect is the maximum cycle mean; for Collatz it is log(3/2), from the loop at -1, so no periodic correction beats log n, and no finite bank of 2-adic valuation counters helps"
status: >
  PROVED + INDEPENDENTLY AUDITED; FINITE-EXACT. Setting: THM-4474's
  strategy cube. T_sigma(n) = n/2 or (3n + sigma(n mod 2^k))/2; G_sigma is
  the parity graph with weights w = log(3/2) on odd and -log 2 on even
  nodes.
  (R) The following are equivalent:
  (a) class (i);
  (b) a periodic rank R(n) = a log n + h(n mod 2^k) with
  R(Tn) <= R(n) - eps for all n >= n_1;
  (b') a rank a log n + h(n) with h merely bounded, of any period or none,
  and R(Tn) <= R(n) for n >= n_1;
  (c) an integer potential certificate at some F < log_3 2.
  Under (c) one may take a = 1, h = (log 3/r) psi, eps = |mu|/2 and
  n_1 = ceil(2/(3|mu|)).
  (D) min over periodic h of max over edges [w(s) + h(t) - h(s)] equals
  lambda* = rho_max log 3 - log 2. Bounded corrections do no better, and
  lambda* is the maximum mean weight of a normalized circulation (LP
  duality): a descent tension exists iff no positive circulation does.
  (C) For Collatz, lambda* = log(3/2) at every level, from the self-loop at
  -1, attained by h = 0. The Bernoulli-boundary obstruction (witnesses
  2^H - 1) is exactly this loop's integer shadow.
  (K) Adding any finite bank sum c_i v_2(n - beta_i) of 2-adic valuation
  counters (beta_i in Q_2 not positive integers) to any periodic correction
  still fails for Collatz, with the same least defect log(3/2). This
  contains the Kuratowski reframe's section 5 theorem.
  (M) The maximum-density cycles of THM-4479's sigma_k are exactly the
  periodic concatenations of density-F_k first-descent blocks. These are
  the words strictly above slope F_k, one per necklace. The upper
  Christoffel word is the unique balanced one. For k >= 4 there are
  non-Christoffel maximizers (1100 at k = 4, 11100 at k = 5..7), and the
  maximizing set has positive entropy.
  (S) A maximal cycle can be chosen with its nodes distinct mod 2^(k-1).
  Hence the largest rho_max of a class-(i) level-k strategy lies in
  [F_k, G_(2^(k-2))], where G_m is the best lower approximation of log_3 2
  with numerator <= m. The value is exactly 1/2, 1/2, 3/5, 5/8, 5/8 for
  k = 2..6. So "rho_max <= F_k for every provable strategy" is REFUTED for
  4 <= k <= 26.
  These are statements about modified maps. Collatz is OPEN.
source: collatz-procgen-20260922 session, tension lane (2026-09-26); audited and promoted by the session orchestrator 2026-09-26
depends_on:
  - 01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md
  - 01-canon/theorems/THM-4479-strategy-cube-distance-to-provability-sharp.md
related:
  - 05-knowledge/results/bernoulli_boundary_20260925.md (the bounded-correction obstruction, Corollary B)
  - 05-knowledge/results/kuratowski_reframe_20260925.md (section 5, contained in (K))
  - 01-canon/theorems/THM-4481-entropy-merge-law-sign-strategies.md
note: 05-knowledge/results/procgen_tension_20260926_ranks_christoffel_duality.md
scripts: 04-computation/experiments/procgen_tension_20260926_{lib,q1,q2,q3,q4,mip,run}.py (sha256 in the note, section 8; re-checked unchanged)
script_audit: 04-computation/experiments/procgen_tension_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_tension_20260926.out
output_audit: 05-knowledge/results/procgen_tension_20260926_orchestrator_check.out
output_sha256: 05de1d98418511a7736fde97ec8fa4f0632a62de173191c44afbd2472df030b5
hash_basis: raw bytes
audit: >
  The orchestrator checked the following line by line and found them
  sound:
  * Theorem R. The node-path Terras lemma for arbitrary strategies (a path
    of length j fixes n mod 2^(k+j)); the chained inequalities along
    gamma^m; the integer correction log(1 + sigma/(3n)) <= 1/(3n).
  * Theorem D (telescoping plus Lemma P; the circulation-tension
    alternative).
  * Corollary K (the 2-adic shadow of an expanding orbit of 1^(L-1)0 that
    avoids the bank returns every periodic and valuation term to its value
    after one period).
  * Theorem M (cycle lemma; best-lower-approximation argument; heredity).
  * Lemma S (siblings share predecessors, so a cycle through both splits).
  Independent code (procgen_tension_20260926_orchestrator_check.py, written
  without reading the lane's scripts) confirms:
  * Theorem D(ii) by LP against exact Karp on 316 strategies;
  * Collatz's defect log(3/2);
  * Corollary K on random 12-center banks (L = 5, 7, 9);
  * the non-Christoffel maximizers 1100 (k = 4) and 11100 (k = 5..7);
  * the refutation at level 4 (a provable strategy with rho_max = 3/5 >
    F_4, equal to the numerator bound).
  The lane's pipeline was re-run (183 s). Its output is identical up to
  timing and RSS lines, including the CP-SAT and flow-MIP decisions of the
  level-5/6 value sets.
---

# THM-4482 -- ranks are tensions

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_tension_20260926_ranks_christoffel_duality](../../05-knowledge/results/procgen_tension_20260926_ranks_christoffel_duality.md).

## 1. The rank theorem

A **rank** is a Lyapunov function for descent. On the strategy cube, Theorem R says:

> `sigma` is bounded-lookahead provable (class (i)) iff there is `R(n) = a log n + h(n mod 2^k)` with `R(T n) <= R(n) - eps` for all large `n`; and even an arbitrary *bounded* `h` cannot do better than a periodic one.

**The dictionary** (exact, via LP duality):

| object | meaning |
|---|---|
| `h` (rank correction) | a tension (potential difference) on `G_sigma` |
| `w + dh <= lambda` | a descent certificate with margin `-lambda` |
| expanding cycle | a positive circulation, the dual obstruction |
| `lambda* = rho_max log 3 - log 2` | the least possible worst-step defect (max-plus eigenvalue / maximum cycle mean) |
| Lemma P's integer `psi` | the rational certificate |

*Proof sketch.*
* **(c) ⇒ (b).** Rescale the potential. The integer correction `log(1 + sigma/(3n))` is absorbed above `2/(3|mu|)`.
* **(b') ⇒ (a).** Integers following an expanding cycle `gamma` `m` times form one class mod `2^(k+mp)`. On it `a m log(3^a/2^p) <= 2||h||` for all `m`, which is impossible.

## 2. What it says about Collatz

* **The defect.** Collatz's least rank defect is `log(3/2)` at every level. It comes from the self-loop at `-1`, where `-1 -> -1` is an expanding fixed point. So **no periodic or bounded correction improves on the plain logarithm**.
* **Three obstructions are one.** Three obstructions found by three sessions are the same fact:
  * the bounded-correction obstruction (Bernoulli boundary), whose witnesses `2^H - 1` shadow `-1`;
  * the finite-center-bank obstruction (reframe §5), now with a one-period proof for any finite bank;
  * bounded-lookahead non-provability (THM-4474, Proposition F: `sigma(-1) = +` forbids class (i)).
* **What a rank for Collatz would need.** The obstruction list is infinite, since every expanding necklace needs its own counter. So a Collatz rank must use unboundedly many 2-adic centers, adaptive centers, or height-dependent (archimedean) information. This is the "two places" boundary again: THM-4480 shows that height-aware edits are exponentially cheaper than periodic ones.

## 3. The maximizers and the refuted rigidity

* **(M) Maximizers of `sigma_k`.** The maximum-density cycles of THM-4479's `sigma_k` form a positive-entropy family of block concatenations. The Christoffel word is one of them, the unique balanced one, not the only one.
* **(S) Rigidity fails.** "`rho_max <= F_k` for every provable level-`k` strategy" is false. The truth is `F_k <= max rho_max <= G_(2^(k-2))`, which is exact for `k <= 6`. The value for `k = 7` (`17/27` or `29/46`) is open: CP-SAT did not decide it in 2 hours.
