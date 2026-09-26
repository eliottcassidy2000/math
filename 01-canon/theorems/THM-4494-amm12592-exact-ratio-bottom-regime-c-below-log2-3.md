---
id: THM-4494
title: "AMM 12592: C* <= 197/125 = 1.576 < log_2 3, by replacing the crude bottom-regime majorant of THM-4468 with the exact binomial ratio; the log_2 3 coincidence of the recurrent-numbers atlas is a proof artifact"
status: >
  PROVED modulo the inherited machinery of THM-4468 (construction, Lemma R,
  Long's fold identities, the interval contour certificates; PROVED +
  INDEPENDENTLY AUDITED there) + FINITE-EXACT (finite certificates 32 <= N
  <= 2048 and explicit integral super-blocks N <= 256 at c = 197/125). One new
  step: in THM-4468's level-0 bottom regime the ratio
  C(A_0 + r - 1, r)/C(R_0, r), A_0 = K + 2m, was bounded by the largest
  factor to the r-th power, which needs c > 203/128 = 1.5859; the exact
  running product prod_(j<r)(A_0 + j)/(R_0 - j) is at most the majorant
  product prod_(j<r)(aN + j)/(bN - j), a = 2 - c + 1/8 >= A_0/N,
  b = c - 1 - 1/N_A <= R_0/N, whose factors increase in j and decrease in
  N, so it is largest at r = 1 (the convexity of
  I(tau) = int_0^tau log((a+s)/(b-s)) ds with I(3/128) < 0 controls the
  range r > 3N_A/128), giving theta = a/b = 0.953529 in place of the crude
  1.03643, and a bottom margin of 109.592 >= 5. (The first version used the
  exact factors at N_A and claimed they decrease in N; the independent
  audit showed they do not, since ceil(cN) - cN varies with N; repaired,
  conclusion intact, MISTAKES 2026-09-26.)
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
  INDEPENDENTLY AUDITED (2026-09-26): one repairable error found and
  repaired (see audit field); verdict after repair: conclusion stands.
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
script_sha256: b19bbee7ffdea5e252b51ef8382c8143edc7e9285b4f347143685e16b28a6bdf
script_finite_sha256: b3d39ef160a37e2367bff5f32d1efbe694be27346db3bca1e55e826a83a3db7f
output_sha256: 9fed56abe544df90dcf76ef482fadc153d82ad444c9025a48b7905df1c4f9df7
output_finite_sha256: 6765386c179f0c0387581cfe81a3b8e7610495d6a24d8664e126831c761d74c4
hash_basis: raw LF bytes
audit: >
  Self-audited: the exact-ratio lemma re-derived (monotonicity of each
  factor in j, left Riemann sum of an increasing integrand, convexity of I
  with I(0) = 0, interval certificate of I(3/128) < 0); THM-4468's scripts
  changed only in the constant and in the theta step (documented in the
  docstrings).
  Independent adversarial audit (auditor subagent, 2026-09-26;
  04-computation/experiments/amm12592_opus_20260926_exactratio_audit.py ->
  05-knowledge/results/amm12592_opus_20260926_exactratio_audit.out, exact
  rationals + mpmath 50 dps, nothing imported from the audited scripts).
  CONFIRMED: factor increasing in j; left-Riemann direction log P_r <= N I(r/N);
  I convex, I(0) = 0, chord bound; I(3/128) = -1.38315e-4 < 0 at 197/125
  (primitive = quadrature); chord range r > r_A; (I2) used in THM-4468's form
  with theta^r -> P_r; scripts differ from THM-4468's only in docstring/CLI/theta
  step; finite N <= 256 and quick contour re-runs reproduce the outputs; 63/40
  fails as stated; N = 16, 64 blocks at 197/125 re-verified from the fair-coin
  definition; assembly and N = 16 handling as in THM-4468.
  ERROR (repairable, conclusion intact): lemma claims 1/4 -- the EXACT factor
  (A_0+j)/(R_0-j) is not decreasing in N along N = 16*4^k (the residue of
  ceil(cN) varies): at 197/125, P_1(16384) = 8994/9437 = 0.953057 >
  P_1(4096) = 0.952946 = the first script's theta, and P_1(N) -> a/(c-1) =
  0.953125. Only the majorant factors (aN+j)/(bN-j) are monotone in N.
  Repair applied (2026-09-26, same session): theta = max_r G_r(N_A) = a/b =
  0.953529 at 197/125 (chord 0.5675 for r > r_A), bottom margin 109.592 >= 5;
  at 158/100 theta = 0.940051, margin 142.360; both full certificates re-run,
  outputs and hashes replaced. Provenance gap (header output hashes matched
  no committed file) fixed by recomputing after the re-run. Cosmetic: cell
  counts are 154/116; I(3/128) < 0 down to c = 1.57434, not 1.576.
  Verdict after repair: SOUND; C* <= 197/125 < log_2 3 stands.
---

# THM-4494 -- AMM 12592 below log_2 3

**PROVED (modulo THM-4468's machinery) + INDEPENDENTLY AUDITED (one repairable error found and repaired).** Full note:
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
`P_r = prod_(j<r) (A_0 + j)/(R_0 - j)` is at most the majorant product
`G_r(N) = prod_(j<r) (aN + j)/(bN - j)` with `a = 2 - c + 1/8 >= A_0/N`,
`b = c - 1 - 1/N_A <= R_0/N` (`N >= N_A`); each majorant factor increases in
`j` and decreases in `N`, so `G_r(N) <= G_r(N_A)`, which first decreases and
then increases in `r`. With `I(tau) = int_0^tau log((a+s)/(b-s)) ds` one has
`log G_r(N) <= N I(r/N)`, `I` convex, `I(0) = 0`, and `I(3/128) < 0`
(certified), so `P_r(N) <= max( max_(r <= 3N_A/128) G_r(N_A), exp(3N_A I(3/128)/128) )`
for every `N >= N_A` and every `r < 3N/128`. At `c = 197/125` this maximum is
`G_1 = a/b = 0.953529`, and the bottom margin is `109.592`. (The exact factors
`(A_0 + j)/(R_0 - j)` are *not* monotone in `N`: `ceil(cN) - cN` varies
along `N = 16 * 4^k`, and `P_1(16384) = 0.953057 > P_1(4096) = 0.952946`;
the first version of this step used them, and the independent audit
caught it.)

## 3. What is re-run, what is inherited

* Re-run at the new `c`: the 154 level cells (I1), the tail terms `eps_D`
  and `eps_T`, the 116 middle cells (I4), the finite certificates for
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
