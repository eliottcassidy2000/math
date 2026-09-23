---
id: THM-4468
title: "Golden-zero super-blocks beat the golden constant for AMM 12592: C* <= 159/100 < 1 + log_5(phi^2)"
status: >
  PROVED (Lemma R re-derived; fold identities exact; interval-arithmetic
  contour certificate at 60-bit precision with outward rounding for every
  N >= 4096, uniform in N by monotonicity; exact finite certificates for
  N = 16, 64, 256, 1024; classical input CITED: Robbins' Stirling bounds) +
  INDEPENDENTLY AUDITED. There is an exactly fair, deterministic,
  complement-symmetric fair-coin extractor with T(L) <= ceil(159 L/100) for
  every L >= 16 and T(L) <= 2L for L < 16. Hence the uniform constant
  satisfies C* <= 159/100 < C_* = 1 + log_5(phi^2) = 1.59799. The golden
  constant, optimal for separately balanced dyadic blocks (THM-3009; Long
  2026-09-20), is not the uniform optimum. With THM-4467 (extended to
  gamma = 0.377): 1.377 <= C* <= 1.59.
source: collatz-procgen-20260923 AMM lane (HYP-9128), audited and promoted by the session orchestrator 2026-09-23
depends_on:
  - THM-2966-spine-normal-form-for-critical-run-fair-extractors
related:
  - THM-3009-archimedean-floor-for-balanced-block-extractors
  - THM-4467-uniform-polya-capacity-gap-amm12592
  - THM-3027-capacity-threshold-is-log-sqrt5-phi
script: 04-computation/experiments/amm12592_procgen_20260923_hyp9128_contours.py
script_finite: 04-computation/experiments/amm12592_procgen_20260923_hyp9128_finite.py
output: 05-knowledge/results/amm12592_procgen_20260923_hyp9128.out
script_sha256: 8ef7d42def100d0d4b95983a7c8864a96706fa95747fea57cc3ad5a5ba8e744c
script_finite_sha256: 206f5f7e3ed90ad307c6fc65935a13cf4ba3e5118de20dc8e908ef5fb12a439f
output_sha256: 8edace740e77533687d28004eae776473d437d3612f8fb671c6653892091afde
hash_basis: raw LF bytes
audit: >
  The orchestrator re-derived the super-block reduction (Lemma S: homogenize
  to degree BN; p<->q invariance iff P is palindromic; Phi_block = w^N S(w)).
  Independent code
  (04-computation/experiments/amm12592_procgen_20260923_hyp9128_orchestrator_check.py)
  rebuilt Phi_block(p) for the explicit N = 16 and N = 64 super-blocks
  directly from the fair-coin definition. It confirms, as exact polynomial
  identities, realizability (box and parity), deadlines <= ceil(1.59 L),
  Phi(p) = Phi(1-p), and Phi = w^N S_N(w). The orchestrator re-ran the
  interval contour certificate and reproduced every number: 152 cells, worst
  rate -0.01464, bottom margin 33.6, top 2415.6, middle 7e-5. The analytic
  lemmas (F1, K, L0) were read and checked for structure, not re-derived line
  by line.
---

# THM-4468 -- super-blocks with a golden zero beat the golden constant

**PROVED + INDEPENDENTLY AUDITED** (finite part re-verified from the definition; certificate re-run).
Full proof: [amm12592_procgen_20260923_hyp9128_proof](../../05-knowledge/results/amm12592_procgen_20260923_hyp9128_proof.md).

## 1. Inheritance

* **Closest proved mechanism.** THM-3009 (and Long's §11) show that
  *separately* balanced dyadic blocks cannot go below
  `C_* = 1 + log_5(phi^2)`.
* **Corrected near miss.** The working expectation `C* = C_*`
  (06-writeups remark 6; THM-3009 §4).
* **Sidecar.** THM-4467's Theorem B shows the lacunary 0-spine has natural
  boundary `|w| = 1`, and that any sub-golden extractor must be analytic at
  the golden point `w = -1`.

## 2. Construction and statement

Fix `c = 159/100`, `B = 4`, `N = 16*4^k`, `m = N/16`, and put

```text
S_N(w) = (1+w)^m (1 + w^m + ... + w^(15m)).
```

`S_N` has a zero of order `m` at `w = -1`. The super-block covers the levels
`L in [N, 4N)`, with deadlines `t_i = min(4N, ceil(c(N+i)))`. Lemma R, a lattice
rounding for palindromic targets that is re-derived rather than cited, produces
integer signed counts `e_(L,z)` with

```text
|e| <= C(R_L, z),    e == C(R_L, z) (mod 2),    Phi_block = w^N S_N(w).
```

So each block is fair (`w = pq` is `p<->q` symmetric), and the blocks tile
`[16, infinity)`. Long's ratio-2 blocks cover `L < 16`.

**Theorem.** The resulting extractor is exactly fair, and `T(L) <= ceil(159L/100)`
for `L >= 16`. Hence `C* <= 159/100`.

## 3. Proof structure

* **N >= 4096.** Cauchy estimates on explicit circles, through Long's fold
  identities (Lemmas F1, K, L0). Interval arithmetic certifies negative rates
  in all 152 level cells and in the level-0 bottom, top and middle regimes.
  The bounds have the form `poly(N) exp(-delta N)`, so they hold for every
  `N >= 4096`.
* **16 <= N <= 2048.** Exact margin certificates. Explicit integral blocks are
  built and verified for `N <= 256`.
* **Assembly.** Lemma S plus realizability turns each block into a balanced
  super-block.

## 4. Scope and open questions

* The constant `159/100` is limited only by a crude bottom-regime majorant.
  The family's measure-level threshold is about `1.578`.
* The lower bound is `C* >= 1.377` (THM-4467 at `gamma = 0.377`).
* **Open.**
  * Whether `C* < 3/2` (HYP-9129). Realizable handoff states stall near
    `1.567`.
  * The exact value of `C*`.
