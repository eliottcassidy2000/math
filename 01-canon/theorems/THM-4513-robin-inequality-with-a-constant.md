---
id: THM-4513
title: "The Robin inequality with a constant: for all m >= 2 and L >= 1, the survivors of the reflected barrier satisfy N_m(L) <= 1.0039 A_(m+15)(L). This is HYP-9142 with the fixed shift 15, replacing the polynomial factors (m+2)^17 and 4e(m+3)L of the two earlier proofs. With pairpeak Theorem C, the private pairing price is O(L) rho^peak_L. For 2 <= m <= 24 the shifts 1, 2, 3 hold with constants 1.1562, 1.0474, 1.0148 for every L"
status: >
  PROVED (hand proofs; computer-assisted finite parts: exact integer
  computation, plus floating-point evaluation of explicit closed forms with
  relative margins >= 1.5e-3 against float error < 1e-12, not interval
  arithmetic) + INDEPENDENTLY AUDITED (the lane's fresh-code pre-audit;
  the orchestrator's own exact counters and proof reading); FINITE-EXACT
  (Theorem R_small; exact verification ranges); CITED (Hoeffding 1963;
  Fortuin-Kasteleyn-Ginibre 1971; Robbins 1955; Kolmogorov and Ville
  maximal inequalities; Morse-Hedlund 1940).
  OPEN: HYP-9142 with shift 1 for all m (observed max N_m/A_(m+1) =
  1.0032546 at (m, L) = (7, 50)); the constant 1 is false for every shift.
  Setting (HYP-9142): N_m(L) counts the survivors of the reflected barrier
  (the zone site m sends both letters down); A_M(L) counts the undecided
  words with peak slope < 3^M. In the coordinate V_t = U_t - floor(ct),
  both are lazy nearest-neighbour walks in the same Sturmian environment.
  (Lemma 1) Hybrid telescoping: N_m(L) - A_M(L) = sum over zone visits of
  the gradient V_k(m-1) - V_k(m) of the hard-wall continuation count.
  (Prop 1) If that gradient is <= 0 for every remaining time k >= k1, then
  N_m <= Gamma* A_M with Gamma* a supremum over the last k1 steps.
  (Lemma 3 + Prop 2, eventual monotonicity) One bridge inequality implies
  the gradient condition for every k >= n1 + 1, by TP2 and
  Chapman-Kolmogorov. The bridge inequality holds for
  n1(m) = ceil((m+1)/(mu - 10^-3)), mu = log_3 2 - 1/2. It is proved by a
  leading-zero injection, FKG on e-subsets, a real-valued cycle lemma and
  Hoeffding; exact for 2 <= m <= 25 over all Sturmian factors, analytic for
  m >= 26.
  (Lemma 2 + Prop 3, the last k1 steps) Robin <= half-line, and the
  hard-wall/half-line ratio is monotone in the start. So
  Gamma* <= 1/(1-q), with q = P(top | survive) <= 3.8517e-3 for shift 15,
  by four regimes.
  Consequence: pi_L <= (1.0039 3^15 (27L+18) + 4L 2^-L) rho^peak_L =
  O(L) rho^peak_L (constant about 3.9e8), improving O(L^18) (robin note)
  and O(L^3) (THM-4488).
  Collatz is not addressed: these are finite-horizon counting statements.
source: collatz-procgen-20260922 session, robin2 lane (2026-09-26), attacking the constant form of HYP-9142 (Conjecture R); audited and promoted by the session orchestrator 2026-09-26
depends_on:
  - 05-knowledge/hypotheses/HYP-9142-robin-inequality-barrier-price.md (definitions of N_m, A_M)
  - 05-knowledge/results/procgen_pairpeak_20260926_pairing_peak_price.md (Theorems B, C: the consequence for pi_L)
related:
  - 01-canon/theorems/THM-4488-private-pairing-peak-price-rational-bridge.md (codex crossroads223: the polynomial-factor proof, O(L^3))
  - 01-canon/theorems/THM-4480-peak-discounted-provability-price.md (the peak discount)
  - 05-knowledge/hypotheses/HYP-9140-pairing-price-is-peak-discounted.md
note: 05-knowledge/results/procgen_robin2_20260926_constant_robin.md
scripts: 04-computation/experiments/procgen_robin2_20260926_{lib,run}.py
script_audit: 04-computation/experiments/procgen_robin2_20260926_orchestrator_check.py
output: 05-knowledge/results/procgen_robin2_20260926.out
output_sha256: 02ab8aeada4b518de8defba7bfb253a9a1033955451a4400ccdcfe8fc58bfbfe
output_audit: 05-knowledge/results/procgen_robin2_20260926_orchestrator_check.out
output_audit_sha256: d05ec443e3b291273e566539c76f0b6a8b3837b289a70512b98b1c89cd326c40
hash_basis: raw bytes
audit: >
  Pre-audit (the lane's fresh agent, its own code, not importing the lane's
  library):
  - every lemma, proposition and the assembly rated SOUND, with no
    counterexample;
  - Proposition 1 checked exactly for seven (m, K) pairs, including shift 1;
  - Lemma 4 checked on 320 large random bridges;
  - regime D checked on exact instances;
  - the exact EM check extended to m = 26..45, 50, 60, 80;
  - all constants reproduced; one wording fix, applied.
  Orchestrator:
  - Read Lemma 1 (the hybrid telescoping identity: the zone step
    contributes 2V(m-1) against V(m-1) + V(m)), Proposition 1, Lemma 2
    (Robin below the half-line; the ratio monotonicity by induction) and
    Lemma 3 (TP2 by tail swapping at the first meeting;
    Chapman-Kolmogorov) and found them sound.
  - Read the structure of Propositions 2 and 3.
  - Independent exact counters (the orchestrator's own N_m and A_M
    implementations from the robin audit; the robin2 lane's code was not
    read) confirm:
    - N_1 = 0;
    - 10000 N_m(L) <= 10039 A_(m+15)(L) for 2 <= m <= 14, L <= 400;
    - the ratio profile (c0 = 1: 1.0032546 at (7,50); c0 = 2: 1.0000889 at
      (9,58); c0 = 3: 1 + 2.0e-6 at (12,77));
    - Theorem R_small's inequalities for 2 <= m <= 24, L <= 250-400.
  - The lane's runner was re-run (186.9 s, 28 checks, ALL CHECKS PASSED).
    It is identical to the committed .out up to timing fields.
  Scope notes:
  - The analytic constants are floating-point evaluations with large
    margins, not interval arithmetic.
  - Shift 1 is proved only for m <= 24.
---

# THM-4513 — the Robin inequality with a constant

**PROVED (computer-assisted finite parts) + INDEPENDENTLY AUDITED.** Full note: [procgen_robin2_20260926_constant_robin](../../05-knowledge/results/procgen_robin2_20260926_constant_robin.md).

## 1. What changed

HYP-9142 compares two ways of keeping the letter walk in a strip:
- a reflecting barrier at level `m`, which counts the orbits a periodic pairing strategy fails to catch;
- a hard wall slightly higher, which counts the undecided words below a given peak.

Two sessions had proved the comparison up to polynomial factors. It now holds with a constant factor: `N_m(L) <= 1.0039·A_(m+15)(L)` for all `m, L`. Through pairpeak Theorem C, this makes the private pairing price `O(L)·rho^peak_L`. So the price of a periodic pairing is the peak-discounted undecided mass times a linear factor.

## 2. How

In integer coordinates both walks live in the same Sturmian environment. The reflecting walk differs from the hard-wall walk only at one site. An exact telescoping identity charges the difference to the gradient of the hard-wall count at that site. A bridge inequality, proved by TP2, FKG, a cycle lemma and Hoeffding, makes the gradient nonpositive after the natural descent time `(m+1)/mu`. The remaining short window costs at most `1/(1 - q)`, where `q < 0.4%` is the chance of touching the far wall.

## 3. What stays open

Shift 1 for all `m` (the data say `1.0033`), and any version with constant exactly 1, which is false.
