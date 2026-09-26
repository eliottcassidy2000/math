# AMM 12592: C* ≤ 197/125 = 1.576 < log₂3, by the exact binomial ratio in THM-4468's bottom regime; the log₂3 lead of the recurrent-numbers atlas is dead

**Status: PROVED modulo the inherited machinery of THM-4468 (construction,
Lemma R, Long's fold identities F1/K/L0, interval contour certificates,
which are PROVED + INDEPENDENTLY AUDITED there) + FINITE-EXACT (the
finite certificates and explicit integral super-blocks re-run at
`c = 158/100`) + the one new step, an elementary replacement of a crude
majorant, proved below. Consequence: the uniform constant of AMM 12592
satisfies `1.377 <= C* <= 1.576 < log_2 3 = 1.58496`, so the coincidence
"`log_2 3` lies in the window `[1.377, 1.59]`" recorded in the
[atlas](constants_atlas_20260926_recurrent_numbers.md) is a proof artifact:
THM-4468's constant `159/100` was limited by a bound that needs
`c > 203/128 = 1.5859`, which happens to sit `0.001` above `log_2 3`.
Session `collatz-exponent-atlas-20260926` (opus), 2026-09-26.**

Scripts (derived from THM-4468's, changes documented in their docstrings):
`04-computation/experiments/amm12592_opus_20260926_exactratio_contours.py`
(analytic part, `c` a parameter, exact-ratio bottom regime),
`amm12592_opus_20260926_exactratio_finite.py` (finite part, only the
constant changed). Outputs: `amm12592_opus_20260926_exactratio_c158.out`,
`amm12592_opus_20260926_exactratio_finite_c158.out`.

## 0. What THM-4468 proves and what limited its constant

THM-4468 builds, for `N = 16 * 4^k`, a complement-symmetric super-block on
the levels `L in [N, 4N)` with deadlines `t_i = min(4N, ceil(c(N+i)))` from
the handoff state `S_N(w) = (1+w)^m (1 + w^m + ... + w^(15m))`, `m = N/16`,
and proves that the integer counts it needs exist (Lemma R) once four
families of inequalities hold: the level masses (I1), the level-0 bottom
(I2), top (I3) and middle (I4) margins. For `N >= N_A = 4096` these are
certified by interval arithmetic on explicit circles; for `16 <= N <= 2048`
exactly. The proof note states (section 5.2) that the constant `159/100` is
limited by (I2) alone: the bottom regime bounds the exact coefficient of the
auxiliary series `A(x) = S_N(x/(1+x)^2)(1+x)^(-K)`, `K = 2N - ceil(cN)`, by
Lemma L0(b), `|[x^r] A| <= C(A_0 + r - 1, r)` with `A_0 = K + 2m`, and then
compares with `C(R_0, r)`, `R_0 = ceil(cN) - N - 1`, through

```text
C(A_0 + r - 1, r) / C(R_0, r) <= ((A_0 + r - 1)/R_0)^r =: theta_1^r,
```

the largest factor raised to the `r`-th power. With `r < r_1 = 3N/128` this
needs `theta_1 < 1`, i.e. `2 - c + 1/8 + 3/128 < c - 1 - 3/128`, i.e.
`c > 203/128 = 1.5859`. Every other certified rate is comfortably negative
at `c = 1.59`, and the family's measure-level threshold is `1.578`.

## 1. The replacement: the exact ratio is largest at `r = 1`

**Lemma (exact ratio, majorant form).** Put
`P_r = C(A_0 + r - 1, r)/C(R_0, r) = prod_(j<r) (A_0 + j)/(R_0 - j)`.
Let `a = 2 - c + 1/8` and `b = c - 1 - 1/N_A` (rationals), so that
`A_0 <= aN` and `R_0 >= bN` for `N >= N_A`, and put
`G_r(N) = prod_(j<r) (aN + j)/(bN - j)` and
`I(tau) = int_0^tau log((a+s)/(b-s)) ds`. Then:

1. `P_r(N) <= G_r(N)` (each exact factor is at most the majorant factor);
   the majorant factor `(aN + j)/(bN - j)` is increasing in `j` and
   decreasing in `N` (`d/dN = -j(a+b)/(bN-j)^2 <= 0`), so `G_r(N) <= G_r(N_A)`
   for `N >= N_A`;
2. `log G_r(N) <= N I(r/N)` (left Riemann sum of the increasing integrand);
3. `I` is convex with `I(0) = 0`, so if `I(tau_1) < 0` for `tau_1 = 3/128`
   then `I(tau) <= (tau/tau_1) I(tau_1) < 0` on `(0, tau_1]`;
4. hence for `N >= N_A`: `P_r(N) <= G_r(N_A)` for `r <= r_A = 3N_A/128`,
   and `P_r(N) <= exp(r_A I(tau_1)/tau_1)` for `r_A < r < 3N/128`.

So `max_(1 <= r < r_1) P_r <= theta := max( max_(r <= r_A) G_r(N_A), exp(r_A I(tau_1)/tau_1) )`,
computed exactly (rationals for the products, interval arithmetic for `I`);
the inner maximum is `G_1 = a/b`.

*Correction after the independent audit (2026-09-26).* The first version
of this lemma bounded `P_r(N)` by the **exact** products `P_r(N_A)`,
claiming that the exact factor `(A_0 + j)/(R_0 - j)` decreases in `N`. It
does not: `A_0` and `R_0` depend on `ceil(cN)`, whose residue
`ceil(cN) - cN` varies along `N = 16 * 4^k` (`0.704` at `4096`, `0.816` at
`16384`, `0.264` at `65536`), and at `c = 197/125` one has
`P_1(16384) = 8994/9437 = 0.953057 > P_1(4096) = 2248/2359 = 0.952946`,
with `P_1(N) -> a/(c-1) = 0.953125`. Only the majorant factors are
monotone. The repair costs `0.0006` in `theta` (`0.952946 -> 0.953529` at
`197/125`, `0.939789 -> 0.940051` at `158/100`) and about `1.4` in the
bottom margin; every conclusion survives. Logged in `01-canon/MISTAKES.md`.

At `c = 158/100`, `N_A = 4096`: `[exact ratio] N_A = 4096: K = 1720, A0 = K + 2m = 2232, R_0 = 2375, r_A = 96; a = 109/200 >= A0/N_A, b = 59367/102400 <= R_0/N_A; max_r<=r_A majorant product G_r(N_A) = 0.940051 (G_1 = a/b = 0.940051; the exact P_1(N_A) = 0.939789 is not uniform in N); G_(r_A) = 1.389e-01; chord bound for r > r_A: 1.448e-01`; `[level 0] r1 = 3N/128;  qbar = 0.04211, pbar = 0.04041, g1 = 0.922447 (log g1 = -0.08073), theta1 majorant-ratio(N_A) = 0.94005 (r<=r_A value 0.94005, r=r1 value 1.448e-01; I(3/128) in [-0.000472,-0.000472]; crude theta1 would be 1.02178), sigma0 in [0.579756, 0.580000]`. The bottom margin is then
`((c-1)N_A - 1)(1 - theta - eps_D)`, which the certificate requires to be
at least `5`.

*Why this is the whole change.* (I2) is the only place where `theta_1`
enters; (I1), (I3), (I4) and the two tail terms `eps_D`, `eps_T` are the
original interval certificates, re-run at the new `c`. The finite part
(`N <= 2048`) is THM-4468's script with the constant changed: it computes
the level-0 packet exactly and needs no majorant.

## 2. Results

**`c = 158/100`** (full resolution, `N_A = 4096`):

```text
=> worst certified rate -0.00637;  max_i M_i <= exp(-22.76) = 1.307e-10 at N = N_A (bound decreasing in N since every rate < 0)
[V2D] D-contour (-0.3009,1.8625): rate <= -0.02999, extra <= 1.458, eps_D <= exp(-120.74)
[V2T] top-contour (-0.3009,1.8625): rate <= -0.02999, extra <= 2.887, eps_T <= exp(-119.32)
[bottom] margin >= ((c-1)N_A - 1)(1 - theta1 - eps_D) = 142.360 (need >= 5; increasing in N)
=> sqrt(2 R0) * sum <= exp(-6.84) = 1.070e-03 (need <= 1/2, then margin >= binom(R0,r)/2 >= 5)
All analytic inequalities certified at N = N_A = 4096; each bound is of the form poly(N) exp(-delta N) with delta >= 0.0064 (levels, tails) and the middle-regime exponent < 0, hence holds for all N >= N_A.
```

Finite part: `7` margin certificates PASS (`N = 32..2048`); explicit integral super-blocks verified `5` times (`N <= 256`);
`N=   16: max T/L = 31/19 = 1.63158; structure True; h=24; max_(i>=1) M_i = 2.973e-01; level-0 min margin 4.87 (r=1), max |e|/C = 0.4591, |e_(0,R0)| = 1.20e-02; MARGIN CERTIFICATE FAIL  [0s]`
`N= 2048: max T/L = 3270/2069 = 1.58047; structure True; h=3137; max_(i>=1) M_i = 1.382e-11; level-0 min margin 455.00 (r=1), max |e|/C = 0.6167, |e_(0,R0)| = 3.99e-32; MARGIN CERTIFICATE PASS  [31s]`

**`c = 197/125`** (full resolution, `N_A = 4096`):

```text
=> worst certified rate -0.00134;  max_i M_i <= exp(-2.18) = 1.125e-01 at N = N_A (bound decreasing in N since every rate < 0)
[V2D] D-contour (-0.1997,1.7568): rate <= -0.02802, extra <= 0.122, eps_D <= exp(-114.08)
[V2T] top-contour (-0.1997,1.7568): rate <= -0.02802, extra <= 2.334, eps_T <= exp(-111.86)
[bottom] margin >= ((c-1)N_A - 1)(1 - theta1 - eps_D) = 109.592 (need >= 5; increasing in N)
=> sqrt(2 R0) * sum <= exp(-5.72) = 3.281e-03 (need <= 1/2, then margin >= binom(R0,r)/2 >= 5)
All analytic inequalities certified at N = N_A = 4096; each bound is of the form poly(N) exp(-delta N) with delta >= 0.0013 (levels, tails) and the middle-regime exponent < 0, hence holds for all N >= N_A.
```

Finite part: `7` margin certificates PASS (`N = 32..2048`); explicit integral super-blocks verified `5` times (`N <= 256`);
`N=   16: max T/L = 13/8 = 1.62500; structure True; h=24; max_(i>=1) M_i = 2.973e-01; level-0 min margin 4.87 (r=1), max |e|/C = 0.4591, |e_(0,R0)| = 1.20e-02; MARGIN CERTIFICATE FAIL  [0s]`
`N= 2048: max T/L = 3257/2066 = 1.57648; structure True; h=3150; max_(i>=1) M_i = 1.025e-08; level-0 min margin 439.00 (r=1), max |e|/C = 0.6277, |e_(0,R0)| = 1.46e-30; MARGIN CERTIFICATE PASS  [30s]`


**Theorem.** THM-4468's construction with `c = 197/125` (and likewise with `c = 79/50`) is an exactly fair,
deterministic, complement-symmetric extractor with `T(L) <= ceil(197L/125)`
for every `L >= 16` (and `T(L) <= 2L` for `L < 16`, Long's ratio-2 blocks).
Hence `C* <= 197/125 = 1.576 < log_2 3`.

*Proof.* As in THM-4468 sections 3–4: the analytic certificate covers every
`N = 16 * 4^k >= 4096` (each bound is `poly(N) exp(-delta N)` with all rates
negative, monotone in `N`); the exact margin certificates cover
`32 <= N <= 2048`; at `N = 16` the margin certificate falls short of the
convenient threshold `5` (margin `4.87`) but the independent Lemma-R
rounding builds and verifies an explicit integral super-block, as it does
for every `N <= 256`; the blocks tile `[16, infinity)` and Long's blocks
cover `L < 16`. ∎

## 3. Where this certificate stops

| `c` | quick-mode certificate | detail |
|---|---|---|
| `789/500 = 1.578` | certified (quick mode) | `=> worst certified rate -0.00362;  max_i M_i <= exp(-11.48) = 1.030e-05 at N = N_A (bound decreasing in N since every ra` |
| `197/125 = 1.576` | certified (quick mode) | `=> worst certified rate -0.00125;  max_i M_i <= exp(-1.75) = 1.743e-01 at N = N_A (bound decreasing in N since every rat` |
| `63/40 = 1.575` | FAILS: no circle with a negative level rate (`assert target < 0` in the circle search) | `AssertionError` |

The exact-ratio step itself is valid as long as `I(3/128) < 0`, which holds
for `c` down to `1.57434` (root by bisection at `N_A = 4096`, audit script);
the failure of the certificate at `63/40` comes from the level rates (I1),
not from the lemma. Below `1.574` the middle-regime pairing choice
`r_1 = 3N/128` would also have to change. The family's measure-level threshold
is `1.578` (THM-4468 section 10), and the crossroads/procgen H5 proposal (a
saddle-point bound in the bottom regime) targets `1.58` for this family and
`1.575` for `kappa = 0.09`; the exact ratio reaches the first of these with
no analysis at all.

## 4. What it means for the atlas

* The window `[1.377, 1.59]` contained `log_2 3` only because the upper
  bound was stuck at a majorant threshold `203/128 = 1.5859`. With
  `C* <= 1.576` the constant `log_2 3` has no role in AMM 12592; the atlas
  entry moves from "lead" to "trap".
* The proved window is now `[1.377, 1.576]`; HYP-9129 (`C* < 3/2`?) is
  untouched, and the lower-bound side (Pólya capacity, reach `1.3775`)
  remains the frontier.

## 5. Independent audit (2026-09-26)

Auditor subagent, blind re-derivation with exact rationals and mpmath at
50 digits: `04-computation/experiments/amm12592_opus_20260926_exactratio_audit.py`
-> `05-knowledge/results/amm12592_opus_20260926_exactratio_audit.out`.
Confirmed: the binomial-ratio identity; factors increasing in `j`; the
left-Riemann direction; convexity of `I` and the chord; `I(3/128) < 0` at
both `c` (primitive = quadrature); the `r > r_A` range; that (I2) is used
in THM-4468's form with `theta^r -> P_r`; that the scripts differ from
THM-4468's only in docstring, CLI and the `theta` step; the finite re-runs
(`N <= 256`, explicit blocks rebuilt from the fair-coin definition at
`N = 16, 64`); the quick contour re-run; the `63/40` failure; the assembly.
Found: the monotonicity-in-`N` error above (repaired: majorant factors), a
provenance error (the theorem header's output hashes matched no committed
file; recomputed after the re-run), and two cosmetic slips (cell counts
`154`/`116`, and the root of `I(3/128) = 0` at `1.57434`). Verdict after
repair: the conclusion `C* <= 197/125 < log_2 3` stands with bottom margin
`109.592` (`c = 197/125`) and `142.360` (`c = 158/100`).
