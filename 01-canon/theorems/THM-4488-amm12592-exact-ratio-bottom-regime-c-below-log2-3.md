---
id: THM-4488
title: "AMM 12592: C* <= 197/125 = 1.576 < log_2 3, by replacing the crude bottom-regime majorant of THM-4468 with the exact binomial ratio; the log_2 3 coincidence of the recurrent-numbers atlas is a proof artifact"
status: >
  PROVED modulo the inherited machinery of THM-4468 (construction, Lemma R,
  Long's fold identities, the interval contour certificates; PROVED +
  INDEPENDENTLY AUDITED there) + FINITE-EXACT (finite certificates 32 <= N
  <= 2048 and explicit integral super-blocks N <= 256 at c = 197/125). One new
  step: in THM-4468's level-0 bottom regime the ratio
  C(A_0 + r - 1, r)/C(R_0, r), A_0 = K + 2m, was bounded by the largest
  factor to the r-th power, which needs c > 203/128 = 1.5859; the exact
  running product prod_(j<r)(A_0 + j)/(R_0 - j) is largest at r = 1 (its
  factors increase in j and decrease in N; the convexity of
  I(tau) = int_0^tau log((a+s)/(b-s)) ds with I(3/128) < 0 controls the
  range r > 3N_A/128), giving theta = A_0/R_0 = 0.952946 at N_A = 4096 in
  place of the crude 1.03643, and a bottom margin of 110.967 >= 5.
  Every other inequality of THM-4468 is the original interval certificate
  re-run at c = 197/125 (all rates negative, worst -0.00134; each bound
  poly(N) exp(-delta N), hence uniform in N >= N_A). Consequence: the
  uniform constant of AMM 12592 satisfies 1.377 <= C* <= 1.576 < log_2 3 =
  1.58496, so log_2 3 has no role in this problem: the earlier window
  [1.377, 1.59] contained it only because 159/100 was stuck at the crude
  threshold 203/128. Quick-mode probes: the certificate still passes at
  c = 197/125 = 1.576 and fails at c = 63/40 = 1.575 (no circle with a
  negative level rate), matching the family's measure-level threshold
  1.578 (THM-4468 section 10) up to the circle discretisation. HYP-9129
  (C* < 3/2?) is untouched; the lower bound remains THM-4467's 1.377.
source: collatz-exponent-atlas-20260926 session (opus), 2026-09-26; testing the "log_2 3 in the AMM window" lead of the recurrent-numbers atlas. Mechanism: THM-4468's certificate with one crude majorant replaced by the exact binomial ratio (the crossroads/procgen proposal H5 asked for a saddle-point bound; the exact ratio suffices).
depends_on:
  - 01-canon/theorems/THM-4468-golden-zero-superblocks-beat-golden-amm12592.md (construction, Lemma R, fold identities, certificate machinery, assembly)
  - 01-canon/theorems/THM-2966-spine-normal-form-for-critical-run-fair-extractors.md (the model and C* = 1 + gamma*)
related:
  - 05-knowledge/results/amm12592_opus_20260926_exactratio_c158.md (full note: the lemma, the runs, the probes)
  - 05-knowledge/results/amm12592_procgen_20260923_hyp9128_proof.md (THM-4468's proof note; sections 5.1-5.2 name the crude majorant as the sole limit)
  - 05-knowledge/results/constants_atlas_20260926_recurrent_numbers.md (the lead this settles)
  - 05-knowledge/hypotheses/HYP-9129-amm12592-below-three-halves.md (the open question the window bears on)
  - 01-canon/theorems/THM-4467-uniform-polya-capacity-gap-amm12592.md (the lower bound)
script: 04-computation/experiments/amm12592_opus_20260926_exactratio_contours.py
script_finite: 04-computation/experiments/amm12592_opus_20260926_exactratio_finite.py
output: 05-knowledge/results/amm12592_opus_20260926_exactratio_c1576.out
output_finite: 05-knowledge/results/amm12592_opus_20260926_exactratio_finite_c1576.out
script_sha256: 38923a4054af122210cb2441eb519e607f88ada021664dc06e87a4fa63c8c9ca
script_finite_sha256: b3d39ef160a37e2367bff5f32d1efbe694be27346db3bca1e55e826a83a3db7f
output_sha256: f97fe3e059baf07e3e0fdc0cb3c1b80c82c4cc7e4986811b5a45da3d6a83192c
output_finite_sha256: 7dcc4f0e51cf1c1689a47598ea8cff826f7e0d23bd767a1ca86fea611f6e48c7
hash_basis: raw LF bytes
audit: >
  Self-audited: the exact-ratio lemma re-derived (monotonicity of each
  factor in j and N, left Riemann sum of an increasing integrand, convexity
  of I with I(0) = 0, interval certificate of I(3/128) < 0); THM-4468's
  scripts changed only in the constant and in the theta step (documented in
  the docstrings); outputs retained. Independent audit not yet performed.
---

# THM-4488 -- AMM 12592 below log_2 3

**PROVED (modulo THM-4468's machinery).** Full note:
[amm12592_opus_20260926_exactratio_c158](../../05-knowledge/results/amm12592_opus_20260926_exactratio_c158.md).

## 1. Statement

There is an exactly fair, deterministic, complement-symmetric extractor for
an unknown-bias coin with pathwise deadline `T(L) <= ceil(197/125 * L)` for
every critical value `L >= 16` (and `T(L) <= 2L` for `L < 16`). Hence

```text
1.377  <=  C*  <=  197/125 = 1.576  <  log_2 3 = 1.58496.
```

## 2. The one new step

THM-4468's bottom regime (level 0, `r < r_1 = 3N/128`) uses Lemma L0(b),
`|[x^r] A| <= C(A_0 + r - 1, r)`, and needs `|e_(0,r)| <= (theta^r + eps_D) C(R_0, r)`
with `theta < 1`. The original certificate took
`theta = (A_0 + r_1 - 1)/R_0`, the largest factor of the ratio, which is
`< 1` only for `c > 203/128`. The exact ratio
`P_r = prod_(j<r) (A_0 + j)/(R_0 - j)` has increasing factors, so it first
decreases and then increases in `r`; with `a = 2 - c + 1/8 >= A_0/N`,
`b = c - 1 - 1/N_A <= R_0/N` and `I(tau) = int_0^tau log((a+s)/(b-s)) ds`,
one has `log P_r <= N I(r/N)`, `I` convex, `I(0) = 0`, and `I(3/128) < 0`
(certified), so `P_r <= max( max_(r <= 3N_A/128) P_r(N_A), exp(3N_A I(3/128)/128) )`
for every `N >= N_A` and every `r < 3N/128`. At `c = 197/125` this maximum is
`P_1 = A_0/R_0 = 0.952946`, and the bottom margin is `110.967`.

## 3. What is re-run, what is inherited

* Re-run at the new `c`: the 152 level cells (I1), the tail terms `eps_D`
  and `eps_T`, the 118 middle cells (I4), the finite certificates for
  `16 <= N <= 2048`, and the explicit Lemma-R roundings for `N <= 256`.
  At `N = 16` the margin certificate gives `4.87 < 5` but the explicit
  integral block exists and is verified, as in THM-4468.
* Inherited unchanged: the construction, Lemma R, Lemmas F1/K/L0, the
  assembly of blocks into an extractor.

## 4. Non-consequences

The lower bound `C* >= 1.377` is untouched; whether `C* < 3/2` (HYP-9129)
is open. The constant is limited by the level rates of this family
(`kappa = 1/16`), not by the bottom regime any more; `kappa` near `0.09`
is expected to reach about `1.574` (THM-4468 section 10).
