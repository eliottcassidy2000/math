---
id: THM-4534
title: "The edge multiset dimension of the hypercube grows like exp(Theta(d^(1/3))) (Allikvere's Open Problem 2): 0.8146 <= liminf ln edim_m(Q_d)/d^(1/3) <= limsup <= 2.0526, where the upper constant C* = (3 sqrt2 ln2)^(2/3) is exactly the reach of the forest-lemma union bound. A closed-form estimate proves existence for every d >= 17, and explicit sets cover 6 <= d <= 16, so finiteness for all d >= 6 needs no union-bound computation (Open Problem 3, partly); density-1/2 sets overshoot by exp(d ln 2 - O(d^(1/3))) (Open Problem 4); 11 <= edim_m(Q_7) <= 19"
status: >
  PROVED + INDEPENDENTLY AUDITED (proofs read; numerical components re-checked):
  - Theorem A, the upper bound: a random set of density q with ln(q 2^d) >= (C* + eps) d^(1/3) resolves Q_d with
    probability -> 1. The proof uses a Fourier atom bound (Lemma A), Chebyshev star forests (Lemma S), path/zigzag forests
    and the forest lemma for every q;
  - Theorem B, the lower bound c* = (3 ln2/(2 sqrt2))^(2/3), a sharper evaluation of THM-4525's entropy inequality L5;
  - Theorem C (C* is the exact reach of the forest-lemma union bound);
  - Proposition K (weighted spanning-tree count of the parallel-pair level graph; see the correction below);
  - OP4 asymptotics.
  PROVED (lane-level numerics): Theorem D, the closed-form bound U_d < 1 for every d >= 17 (interval evaluations at
  d = 17, 18, 19); the Fourier-Hoelder multi-forest lemma and the certified U_10 <= 0.2649.
  VERIFIED + AUDITED: explicit resolving sets for every 6 <= d <= 16, new Q_13 <= 105, Q_14 <= 125, Q_15 <= 135,
  Q_16 <= 171.
  FINITE-EXACT (lane only):
  - no resolving set of Q_7 of size <= 10 (three complete normal forms), so edim_m(Q_7) >= 11;
  - k = 11 is unfinished except the last case;
  - interval-certified sparse bounds for 11 <= d <= 64 (e.g. M_11 = 361, M_64 = 31808).
  OPEN:
  - whether lim ln edim_m(Q_d)/d^(1/3) exists, and its value in [0.8146, 2.0526] (HYP-9169, updated);
  - OP3 in the strict form, one analytic estimate valid from d = 11 (the best forest bound at d = 11 has margin only
    about 6).
  CORRECTION (orchestrator audit): Proposition K's formula is right, but its corollary "the level graph is connected iff h
  is odd" fails at h = n. par(n) gives the matching {p, n-p}, which is disconnected even for odd n. Correct form: connected
  iff h is odd and h < n.
source: collatz-procgen-20260922 session, edim2 lane (2026-10-01; resumed after a reboot and a network drop), the owner's request to settle the remaining open problems of arXiv:2608.09983; audited and promoted by the session orchestrator 2026-10-01
depends_on:
  - 01-canon/theorems/THM-4525-edge-multiset-dimension-of-q6-is-15.md (OP1; lemmas L1-L5; explicit sets d <= 12)
  - external: J. Allikvere, arXiv:2608.09983v1 (forest lemma, Lemma 11; pair types, Lemma 13 and Prop. 14)
related:
  - 05-knowledge/hypotheses/HYP-9169-edge-multiset-dimension-of-hypercubes-grows-like-exp-cube-root.md (upgraded: Theta form proved)
note: 05-knowledge/results/procgen_edim2_20261001_growth_and_uniform_bounds.md
scripts: 04-computation/experiments/procgen_edim2_20261001_{run,lib}.py, procgen_edim2_20261001_{anneal,q7search}.c
script_audit: 04-computation/experiments/procgen_edim2_20261001_orchestrator_check.py
output: 05-knowledge/results/procgen_edim2_20261001.out (a --quick run, 71 checks)
output_sha256: f00aba94da0a4b74d71e2f5b90584d59d25adb56d7b2521bbecbef520d3f056b
output_audit: 05-knowledge/results/procgen_edim2_20261001_orchestrator_check.out
output_audit_sha256: cc542ad74c776fa343e8eefcbb1d308a134d8e17426ac6e9776375783eec800b
hash_basis: raw bytes
audit: >
  The orchestrator read and found sound:
  - Lemma A (Fourier inversion with |1 - q + q e^(i theta)|^N <= exp(-2x sin^2(theta/2))) and Lemma A2 (bathtub);
  - Lemma S (the conditional law given A = a is hypergeometric; Chebyshev puts at least 3/4 of the mass strictly closer
    to the centre);
  - the exponent estimates E1-E3 (level sums of lambda - 2x^2/n give (sqrt2/3) lambda^(3/2) sqrt n);
  - the case analysis of Theorem A (bulk, extreme and par(n) pairs, with (sqrt2/3) C*^(3/2) = 2 ln 2 balancing the
    d^2 4^d pairs);
  - Theorem C (Lemma A2 plus N_e <= 2^(d+1) b_n(a));
  - Theorem B (re-derived by hand: integrating log2 mu_r over the heavy levels gives (2 sqrt2/3) sqrt d lambda^(3/2)/ln 2 = d,
    so c* = (3 ln2/(2 sqrt2))^(2/3)).
  Theorem D's structure (tangent-line concavity reducing to four endpoint values, plus elementary bounds for d >= 20) was
  read; its d = 17..19 interval evaluations are lane-level.
  Independent code (procgen_edim2_20261001_orchestrator_check.py; only the certificate data was taken from the lane):
  - all 11 explicit sets resolve Q_6..Q_16 (up to 524288 edges);
  - Lemmas A and A2 on 300 random (n, n', q);
  - Proposition K's formula exact for every par(h) with d <= 8, which found the connectivity correction above;
  - the constants C*, c*, C*/c* = 4^(2/3) and (sqrt2/3) C*^(3/2) = 2 ln 2, and the entropy-integral identity.
  The lane's runner was re-run in --quick mode (165 s, 71 checks, ALL CHECKS PASSED); identical up to the timing line. The
  full mode (about 40 minutes: the Q_7 k = 10 search, sparse bounds for d >= 15) was not repeated.
---

# THM-4534 — exp(Θ(d^(1/3))) growth of the edge multiset dimension

**PROVED + INDEPENDENTLY AUDITED** (with the lane-level numerical parts listed above). Full note: [procgen_edim2_20261001_growth_and_uniform_bounds](../../05-knowledge/results/procgen_edim2_20261001_growth_and_uniform_bounds.md).

## The four open problems of arXiv:2608.09983

1. **`edim_m(Q_6)`.** It equals 15 (THM-4525).
2. **Growth rate: SOLVED in Θ-form.** `ln edim_m(Q_d)` lies between `0.8146 d^(1/3)` and `2.0526 d^(1/3)` for large `d`.
   - The upper constant is exactly the limit of the random-landmark union bound.
   - The two constants differ by the factor `4^(2/3)`.
   - Whether `ln edim_m(Q_d)/d^(1/3)` converges is open.
3. **A uniform estimate: PARTLY.**
   - A single closed-form inequality proves existence for every `d >= 17`.
   - Explicit resolving sets cover `6 <= d <= 16`.
   - Together these prove finiteness for all `d >= 6` without the paper's computer-checked union bounds for `11 <= d <= 50`.
   - One analytic estimate starting at `d = 11` is still open. The available bound at `d = 11` has a margin of only about 6.
4. **Density 1/2 versus sparse: ANSWERED.**
   - Density-1/2 sets are larger than optimal by a factor `exp(d ln 2 - O(d^(1/3)))`.
   - Sparse random sets are optimal up to the constant `4^(2/3)` in the exponent.

**Also: `11 <= edim_m(Q_7) <= 19`.**

**EXTERNAL CERTIFICATE AUDIT 2026-10-01 (see THM-4525's `external_audit`).** The paper's archived rational union bounds `U_d < 1` for `11 <= d <= 50` were checked to be fractions `< 1`; their derivation was not redone. THM-4534 does not use them: Theorem D covers every `d >= 17` and the explicit sets cover `6 <= d <= 16`.
