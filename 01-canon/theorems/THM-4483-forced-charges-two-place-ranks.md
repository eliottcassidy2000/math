---
id: THM-4483
title: "What a Collatz rank must look like: a rank a log n + bounded h + a bank of 2-adic valuation counters must charge every point of the backward tree of every expanding cycle x with at least a chi(x); so no height-summable bank works for Collatz, nonnegative banks need infinite height moment in every open set of Z_2, and nonnegative rational-center banks with strict descent exist iff every orbit reaches 1"
status: >
  PROVED + INDEPENDENTLY AUDITED; FINITE-EXACT; EMPIRICAL (adaptive
  centers). Setting: T is the Collatz shortcut on Z_(2), the rationals with
  odd denominator. A rank is R(n) = a log n + h(n) + Phi(n), with a > 0,
  h bounded, and Phi(n) = sum_beta c_beta v_2(n - beta) a bank of 2-adic
  valuation counters at rational centers beta that are not positive
  integers (the 2-adic logarithmic potential of the charge
  sum c_beta delta_beta). "Height-summable" means
  sum |c_beta|(1 + log_2 H(beta)) < infinity. (S) means R(Tn) <= R(n) for
  n >= n_1; (L) is the bounded-lookahead version.
  (A) Forced charges. If R satisfies (S) or (L) with a height-summable
  bank, then c_z >= a chi(x) for every expanding periodic point x of T
  (chi = log(3^a/2^p)/p) and every z in its backward orbit. Under (S)
  the charges are non-increasing along forward orbits. The same holds
  for every sign strategy T_sigma.
  (A1) Hence no height-summable signed bank works for Collatz: the
  backward tree of -1 has 2^(j-1) points at depth j, each forcing charge
  >= a log(3/2). Explicit violating integers are given.
  (R+) In the strategy cube, a height-summable bank rank exists iff the
  strategy is class (i). This extends THM-4482: infinite banks add
  nothing to periodic corrections.
  (C) For nonnegative banks finite at one odd and one even integer,
  (S) or (L) forces finite total mass and infinite height moment in every
  nonempty open subset of Z_2.
  (D) A nonnegative bank at rational centers, finite on integers, with
  strict descent R(Tn) <= R(n) - eta for all n >= 2, exists iff every
  Collatz orbit reaches 1. The construction stores each integer's stopping
  time in the height of a private center m + 2^(D_m)/3. This is an
  equivalence, not progress on Collatz.
  (Q2) Adaptive-center ranks fail on explicit seam families:
  2^(V+1) - 2 for bounded period; 2^(Q+1)(2^V - 1) for bounded
  preperiod; hovering integers from near-critical words for unbounded
  period, where the counter is an excursion time.
  Every tested violation on real orbits (to 10^12) is a reset, never a
  shadow step (EMPIRICAL).
  OPEN: nonnegative period-only weight profiles in the thin band
  sum f(p)P(p) < infinity = sum p f(p)P(p); whether any sign strategy has
  finitely many (but some) expanding periodic points (class (i')); the
  future-peak horizon (every orbit point reaches the maximum of its
  remaining orbit within about 13 log_2 m steps on all tested data).
  Collatz is OPEN.
source: collatz-procgen-20260922 session, rank lane (2026-09-26); audited and promoted by the session orchestrator 2026-09-26
depends_on:
  - 01-canon/theorems/THM-4482-ranks-are-tensions-strategy-cube.md
  - 01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md
related:
  - 05-knowledge/results/kuratowski_reframe_20260925.md (section 5 finite centers; section 7 open resource)
  - 05-knowledge/results/bernoulli_boundary_20260925.md
  - 01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md (two places)
  - 01-canon/theorems/THM-4480-peak-discounted-provability-price.md
note: 05-knowledge/results/procgen_rank_20260926_two_place_lyapunov.md
scripts: 04-computation/experiments/procgen_rank_20260926_{lib,q1,cube,q2,q3,run}.py (sha256 in the note, section 7; re-checked unchanged)
script_audit: 04-computation/experiments/procgen_rank_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_rank_20260926.out
output_audit: 05-knowledge/results/procgen_rank_20260926_orchestrator_check.out
output_sha256: ab32013335a0348243729640aa13d3d334c2ec31bca5325755db758490552bf2
hash_basis: raw bytes
audit: >
  The orchestrator checked the following line by line and found them
  sound:
  * Lemma L: the local expansion Phi(y) = c_z D + Phi*(z) + o(1) at good
    integers, using the Liouville/product-formula bound
    v_2(n - beta) <= log_2((n+1) H(beta)).
  * Theorem A: shadows at exact depth D - i; dividing the step
    inequalities by D gives monotone charges; chaining m periods gives
    c >= a chi; the lookahead case via a pigeonhole on jump patterns.
  * Corollary A1 and Theorem R+.
  * Theorem C: (iii) the unconditional shadow lower bound for nonnegative
    banks; (ii) the density of the backward tree of -1.
  * Theorem D: the leak bound and the induction on D_m.
  Independent code (procgen_rank_20260926_orchestrator_check.py, written
  without reading the lane's scripts) confirms:
  * the 2^(j-1) growth of the backward tree of -1 (j <= 12);
  * the shadow violation of the bank {-1: 0.9 kappa} (exact increment
    0.1 kappa at D = 30, 60, 120);
  * the seam violation 1.2 kappa (D - 2) + O(1) when entering -1's shadow
    from the uncharged preimage -4;
  * the forced charge a chi on the cycle -5, -7, -10: 0.95 chi fails over
    one period, 1.05 chi passes.
  The lane's pipeline was re-run (93.5 s). Its output is identical up to
  timing lines.
---

# THM-4483 -- what a Collatz rank must look like

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_rank_20260926_two_place_lyapunov](../../05-knowledge/results/procgen_rank_20260926_two_place_lyapunov.md).

## 1. Two places

A rank for Collatz combines the two places of the rationals:
* the **real** place, through `a log n` (height);
* the **2-adic** place, through valuation counters `v_2(n - beta)`, the 2-adic logarithmic potential of a charge distribution on `Z_2`.

THM-4482 showed that bounded corrections cannot work, because of the loop at `-1`. This theorem shows what the 2-adic part must pay.

## 2. Forced charges

Let `x` be an expanding periodic point, with `chi(x) = log(3^a/2^p)/p > 0`.
* **Inside the shadow.** An integer `y` 2-adically close to `x` (exact depth `D`) follows the cycle, and its height grows by `chi` per step. Its counter at `x` drops by 1 per step. So the charge at `x` must be at least `a chi`.
* **At the seam.** An integer entering that shadow from a preimage `z` also enters with its counter at `z` at depth `D`. If `c_z < c_x` the rank jumps by `(c_x - c_z) D`.
* **Consequence.** So every point of the backward tree of every expanding cycle must carry charge `>= a chi`: the backward tree of `-1` alone needs infinite mass.

| obstruction | source | the same fact |
|---|---|---|
| bounded corrections fail | Bernoulli boundary | the loop at `-1` (THM-4482) |
| finite center banks fail | reframe §5, THM-4482 K | an expanding cycle outside the bank |
| bounded lookahead fails | THM-4474 | an expanding cycle in the parity graph |
| height-summable infinite banks fail | this theorem | the backward tree of `-1`, which is dense in `Z_2` |

## 3. What remains

* **Nonnegative banks.** They must have infinite height moment everywhere (C). Such banks with strict descent exist *iff* Collatz holds (D): the rank can store each integer's stopping time in the height of a private center. So valuation-bank ranks are exactly as hard as Collatz, and no weaker.
* **Adaptive centers.** Ranks that use the currently shadowed center fail at resets (seams), never inside shadows. The failures are the reset families of the reframe (`n_H`, `x_N`) and the hovering near-critical integers.
* **What a successful rank would need.** Unbounded 2-adic resolution tied to archimedean height, as in (D), i.e. information equivalent to stopping times. This is the precise sense in which "discrete (2-adic) plus continuous (real)" is necessary and, in this class, sufficient.
