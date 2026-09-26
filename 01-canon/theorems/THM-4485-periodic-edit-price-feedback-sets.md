---
id: THM-4485
title: "The exact price of periodic Collatz edits is the feedback-vertex number of the expanding cycles of the de Bruijn parity graph; Golomb-Mykkeltveit fails for density thresholds (ratio >= 1.35 on a positive-density set of k for log_3 2); 5n+1's periodic deletion price is ~1/k (k price -> 1), 3n+1's is 2^(-(1-h)k); delta_k >= FVS^odd >= FVS >= nu >= N"
status: >
  PROVED + INDEPENDENTLY AUDITED; FINITE-EXACT; CITED (Mykkeltveit 1972,
  re-derived). A periodic edit of level k sends the residues R mod 2^k to a
  constant v_0 and follows T_q (n/2 or (qn+1)/2) elsewhere. FVS_c(k) is the
  least number of nodes of B(2,k) (the parity graph of T_q) meeting every
  cycle of odd density > c = log_q 2.
  (1) Price. A periodic edit has bounded-lookahead descent iff R meets
  every expanding cycle of B(2,k); the "only if" direction holds for any
  edit values. When it does, every positive orbit is eventually periodic.
  So the minimal Haar density of a provable periodic edit is exactly
  FVS_c(k)/2^k.
  (2) Golomb's conjecture (Mykkeltveit 1972; proof re-derived via the sine
  weight): the least set meeting every cycle of B(2,k) has Z(k) nodes,
  where Z(k) is the number of necklaces. Hence FVS_c(k) <= Z(k) - 1 for
  every c.
  (3) DRIFT for deletions.
  q = 5: 1/(2k) <= price <= (Z(k)-1)/2^k <= 1/k + 2^(-k/2), and
  k price -> 1.
  q = 3: 2^(-(1-h)k)/(3k^2) <= price <= 2^(-(1-h)k).
  Positive drift turns the exponential rate into 1/k; both tend to 0.
  (4) Chain: delta_k (sign flips) >= FVS^odd >= FVS >= nu (cycle packing)
  >= N (expanding necklaces).
  (5) The Golomb analogue for density thresholds is REFUTED.
  At k = 3 and every c in [1/3, 1/2), the cycles (1), (01), (0011) are
  disjoint and expanding while N = 2. At k = 5, c = log_3 2, four disjoint
  expanding cycles against N = 2.
  In general FVS_c(k) >= ceil(Nprim_c(k+1)/2). For c = log_3 2 this beats
  N_k at every k >= 8 with ceil(c(k+1)) = ceil(ck), a set of density
  0.369, with ratio >= 1.35 from k = 10 (limit 1.3548). For c < 1/2 the
  ratio tends to 1.
  The Golomb analogue does hold for closed thresholds j/k with
  j in {0, 1, k-2, k-1, k}.
  FINITE-EXACT:
  FVS_(log_3 2)(k) = 1, 2, 2, 4, 5, 8, 12, 20, 37 (k = 2..10);
  FVS^odd_11 = 58; 95 <= FVS^odd_12 <= 102;
  FVS_(log_5 2)(k) = 2, 3, 4, 6, 9, 15, 27, 45, 84 (k = 2..10).
  New sign-flip bounds: delta_11 >= 58 (was 52) and delta_12 >= 95 (was
  70). nu < FVS at (q,k) = (3,7), with a certified fractional cover of
  weight 22/3.
  These are prices of modified maps. Collatz is OPEN.
  UPDATE 2026-09-26 (THM-4495, opus): for q = 3 the whole chain
  N <= nu <= FVS <= FVS^odd <= delta_k <= |Bad_k| is Theta(2^(hk) k^(-3/2))
  with explicit constants (0.26 and 545), so the price is
  Theta(2^(-(1-h)k) k^(-3/2)) rather than 2^(-(1-h)k)/(3k^2) .. 2^(-(1-h)k).
source: collatz-procgen-20260922 session, mykk lane (2026-09-26); audited and promoted by the session orchestrator 2026-09-26
depends_on:
  - 01-canon/theorems/THM-4474-strategy-cube-provability-by-cycle-densities.md
  - 01-canon/theorems/THM-4479-strategy-cube-distance-to-provability-sharp.md
related:
  - 01-canon/theorems/THM-4480-peak-discounted-provability-price.md (arbitrary edits, which can see height)
  - 01-canon/theorems/THM-4481-entropy-merge-law-sign-strategies.md (why flips cost more than deletions: a flip merges two orbits, a deletion merges all)
note: 05-knowledge/results/procgen_mykk_20260926_expanding_cycle_feedback.md
scripts: 04-computation/experiments/procgen_mykk_20260926_{lib,search,run,certs}.py (sha256 in the note; re-checked unchanged)
script_audit: 04-computation/experiments/procgen_mykk_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_mykk_20260926.out
output_audit: 05-knowledge/results/procgen_mykk_20260926_orchestrator_check.out
output_sha256: e5d870f1c942b36241437051b905cc4a98f1fc26d5e53d61b43e9adcf3aa9ea0
hash_basis: raw bytes
audit: >
  The orchestrator checked the proofs of Theorem 3 (price = FVS), Theorem 4
  (DRIFT bounds), Proposition 5 (chain) and Theorem 2(a)-(b)
  (counterexamples; each node lies on at most two (k+1)-cycles). The
  Mykkeltveit construction is cited and the lane re-derived it; the
  orchestrator did not re-derive it and relied on the computational
  Golomb check below.
  Independent code (procgen_mykk_20260926_orchestrator_check.py, an
  implicit hitting set with HiGHS written without reading the lane's
  scripts) confirms:
  * FVS_(log_3 2)(k) for k = 2..8;
  * FVS_(log_5 2)(k) for k = 2..7;
  * the least set meeting every cycle of B(2,k) has Z(k) nodes (k = 2..6);
  * the k = 3 counterexample;
  * under an optimal periodic edit at q = 3, k = 6, every 2 <= n <= 2*10^5
    reaches 1.
  The lane's full pipeline was re-run (1874 s, 457 MB, ALL CHECKS
  PASSED). Its output is identical up to timing lines and one
  timing-dependent split of step-function segments between RC2 and the
  HiGHS fallback (82/8 against 81/9). Every segment is verified either way.
---

# THM-4485 -- the price of periodic edits, and why Golomb–Mykkeltveit fails for thresholds

**PROVED + INDEPENDENTLY AUDITED.** Full note: [procgen_mykk_20260926_expanding_cycle_feedback](../../05-knowledge/results/procgen_mykk_20260926_expanding_cycle_feedback.md).

## 1. Three kinds of modification, one table

| modification | what it can see | 3n+1 price | 5n+1 price | theorem |
|---|---|---|---|---|
| arbitrary edits | height and residue | `rho_L exp(-Theta(L^(1/3)))` | `2^(-0.0139 L + o(L))` | THM-4480 |
| periodic deletions (`G = v_0` on `R`) | residue only | `2^(-(1-h)k)`, polynomially sharp | `~ 1/k` (`k price -> 1`) | this theorem |
| periodic sign flips | residue only; merges only two orbits | `2^(-(1-h)k + O(log k))` | between `~2/k` and `~3.4/k` on data; limit OPEN | THM-4479, THM-4481 |

* **Drift.** Positive drift changes the *rate* for periodic modifications from exponential to polynomial, but not the limit for deletions.
* **Height.** Arbitrary edits see height, so the positive-drift price falls back to exponential.
* **Flips versus deletions.** Flips cost at least deletions (the chain `delta >= FVS^odd >= FVS`). By THM-4481 a flip merges only two orbits, whereas a deletion merges all edited orbits into one point.

## 2. Golomb, Mykkeltveit, and the failure for thresholds

* **Golomb–Mykkeltveit.** The least feedback vertex set of `B(2,k)` equals the number of necklaces. It is a perfect packing–covering duality, proved by Mykkeltveit with a continuous sine weight that selects the discrete feedback set.
* **For expanding cycles only** (density above a threshold `c`), the duality fails:
  * cycles of length `k+1` pack more densely than necklaces of length `k`, so `FVS_c(k) >= Nprim_c(k+1)/2`;
  * for `c = log_3 2` this exceeds the necklace count by a factor tending to `1.3548` along the `k` where `ceil(c(k+1)) = ceil(ck)`.
* **What the necklace bound reflects.** It gives the right exponent (THM-4479), but it is not the exact covering number.
