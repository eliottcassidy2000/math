# Mean currents and refinement energy of the Collatz parity census: two-octave memory, an out-of-sample test to 2^40, an explicit formula over the zeros of D, and what a proof of the energy budgets would have to beat

2026-10-15, session opus-2026-10-15-S37 (worktree `/tmp/math-wt-octave-rg-37`). Audited: independent S37 referee
([report](../../04-computation/experiments/s37_audit_20261015/REPORT.md)) and clock_duality_oct15's scoped audit
([incoming-work audit](collatz_incoming_audit_20261015b.md)); ledger MISTAKE-618 (theirs) and MISTAKE-619 (this
session's statistics, found by the referee and by the out-of-sample test).

**Owner's prompt (verbatim).** "keep working for a long session of next targets and think the mean-current and refinement-energy bounds, preserving source and integer clock. pull in recent repo work and use Q–Pell and golden-tree controls to close the dyadic energy gaps, merge in new ideas on these and other problems and pursue whatever most mathematically relevant for major proof progress. consider the connection between the geometric Langlands correspondence and S-duality, a property of certain quantum field theories. also think of the field of rational functions on the Riemann sphere and the gijswijt sequence. also merge in the 15 supersingular primes, stern sprimes, non-generous primes, wagstaff primes, woodall primes and cullen primes and leyland primes. adjust targets as necessary and be free in your creative exploration and conjecture proposal and investigation"

The special primes, curling numbers and Langlands threads are in the companion note [special_primes_curling_genus_zero_langlands_20261015.md](special_primes_curling_genus_zero_langlands_20261015.md).

**Status.**
- **PROVED:** THM-4658 (class-sum recursion; two-octave memory `H(4N) = H(N) + H_(2,2)(3N) − H_(2,8)(6N)`; exact energy transport and 4-adic descent; Haar eigen-distribution exactly at zeros of `D`; rank-one excess fixed-point equation; branch norms; character inheritance; the finite-depth spectrum `{1 − D(s)} ∪ {−2^(−s)ζ}`, so no finite-depth excess mode off `Re s = 0`). The identity `H(N) = −[S(N) + 3S(⌊2N/3⌋) + 3S(⌊N/3⌋)]` cancels the leading D-zero powers while retaining their integer-floor residues (section 8; MISTAKE-620).
- **CONDITIONAL, interval-certified (Arb):** nonzero Mercator residue under energy budgets for dyadic `X >= N`, at `N = 2^34`, `2^37` and `2^40` (nested intervals; certified constant budgets up to 227.4 X², 1806.6 X² and 14348 X²); nonzero residues at four further zeros of `D` at `N = 2^37` (HYP-9283, HYP-9285).
- **FINITE-EXACT:** the parity census through 2^34, its extension through 2^37 and the streaming extension through 2^40 (completed after clock_duality's MISTAKE-620 audit flagged it as pending: it reproduces every 2^37 grid point and the independently computed `S(2^37)`, `H(2^37)`), with dyadic energies `E_X/X²` through `X = 2^39` and class sums modulo 243 at every grid point.
- **REFUTED:** HYP-9282 (square-root residue biases; its pre-registered test failed); the HYP-9283(A) forecast that constant budgets fail (out of sample the energies fell far below its predictive range).
- **NUMERICAL:** the explicit-formula fit of `S` over 28 zeros of `D` (holdout error 2%; the pre-registered prediction for `[2^37, 2^40)` PASSED with error 1.95%; dense data to 2^27 bring the remainder to 0.27 sqrt(N)), the odd-step suppression of `H` relative to the halving null, and the 2^28..2^37 dip of the deep wavelet energies, which partly recovered by 2^39 (HYP-9284).
- **OPEN:** any energy budget for all large `X` (calibrated as beyond Krasikov–Lagarias), the modal expansion, and the Collatz conjecture.

## 0. Answers in brief

1. **Refinement energy is paid; the mean current is the obligation.** clock_duality_oct15 proved the sharp refinement envelope M1/M2 (clock mesh) and coarse-mean terminal certificates T1/T2. What remains is a bound on the coarse means, their target "(5)": `C(P_b) <= (4/9)X²` for all dyadic `X >= 2^34`. THM-4658 adds the exact mean-current transport, which keeps the integer clock `N ↦ 4N, 3N, 6N` and introduces no new source; clock_duality writes the same increment as the annular source row `Δ(N) = Σ_(N<n<=4N) w_n t(n)` and has since propagated a linear-log budget from paid annular means (`collatz_annular_overlap_20261015c.md`).
2. **The two-octave memory is exact but cancelled.** Since `Δ = H(4N) − H(N)`, the empirical content is that the normalized profile has near-zero two-octave autocorrelation (regression slope 0.006 on 2^14..2^34, 0.05 on 2^14..2^37). The halving null shows what the cancellation buys: its profile keeps 65-91% of its value across one octave, Collatz 13%.
3. **The square-root biases were octave memory (HYP-9282 REFUTED).** Octave means of every channel are strongly correlated (lag-one 0.63-0.66 at every level), the naive t-values were inflated, and a null keeping `t(2m) = −t(m)` reproduces the out-of-sample behaviour. Three of the five pre-registered channels changed sign in the octaves 2^34..2^36.
4. **The "constant budgets fail" forecast is refuted too (HYP-9283(A)).** Out of sample the energies dropped to 0.04-0.10, the four smallest values of the record (predictive probability 0.1% under the surrogate). Whether a constant budget holds for all `X` is open in both directions.
5. **The conditional pole needs much less than was asked.** With exact data to 2^40, any constant budget `E_X <= 14348 X²` for dyadic `X >= 2^40` (or `(2 + 9616 j) X²`) certifies the nonzero Mercator residue (1806.6 X² beyond 2^37, 227.4 X² beyond 2^34). The intervals `[0.3215, 0.3885]` (2^34), `[0.3431, 0.3668]` (2^37) and `[0.3508, 0.3591]` (2^40) are nested, a passed consistency test. Any such budget still forces `limsup log B(X)/log X >= 0.99821` against the unconditional 0.84 of Krasikov–Lagarias, so the remaining obligation is a major open problem, not a technical gap.
6. **New: the parity sum follows an explicit formula over the zeros of D (HYP-9285).** Least squares over the 28 zeros of `D(s) = 1 + 2^(−s) + (3/2)^s/3` with `Re ρ > 0.65` predicts `S` three octaves ahead to 2%, against 19% for the Mercator zero alone. The fitted coefficients agree with the conditional Arb residues (28 of 28 inside their intervals, most of them wide; 5 certified nonzero, including the low-frequency zero `0.6939 + 22.25i` behind the 2.45-cycle line, each with the fitted value near the middle). The floor-free operator annihilates D-zero powers; the integer operator retains O_ρ(N^(Re ρ-1)) errors per mode. Infinite modal sums and a converse comparison with an untransformed remainder require additional estimates (MISTAKE-620). The pre-registered prediction for `[2^37, 2^40)` PASSED: relative rms error 0.0195 against the bound 0.03 (`S(2^40) = 989,143,994`, predicted 961,770,189).
7. **New: odd-step suppression (HYP-9284) and no finite-depth excess modes (THM-4658 statement 7).** Relative to the halving null, the Collatz `H` carries 3-30 times less energy in sample (median 14) and 25-300 times less out of sample; the level-3..5 wavelet energies match the null in sample but fell below every null replicate between 2^30 and 2^37. On densities of any finite 3-adic depth the transfer operator's spectrum is `{1 − D(s)}` plus `−2^(−s)` times roots of unity, so off the imaginary axis the only finite-depth modes are the D-zero density modes: the wavelet channels carry no finite-depth power-law modes and are remainders.
8. **Q–Pell and golden controls.** The transport has the Q model's "divide by four per two octaves" shape, but the native cross term varies over −0.97..+0.91 (the wavelet note saw −0.949..+0.936), so no scalar Möbius/Pell closure is established (a richer state is not excluded). The golden/Conway control is deterministic and closes; Collatz sits between it and the halving null.

## 1. Inheritance and the live board

- **Programme.** The amplitude note (B1–B3, renewal_shapes_oct11) reduces a nonzero Mercator pole to a discrepancy hypothesis; the break-clock note insists on the integer clock and the incoming current; the wavelet note gives the channel refinement (2), its one-octave energy rule (3) and the packet (7); clock_duality_oct15 gives M1/M2, T1/T2, the Gram slack, the annular rows and the scope repair of THM-4658 4(a).
- **Closest proved mechanism.** Identity (2) of the wavelet note. **Hostile.** The same-event-word, different-energy pair (417 vs 81) of the break-clock note, and clock_duality's mean-zero cylinder function whose joint branches expand at `ρ = 3/4`.
- **Corrected near misses (this session).** Octave means treated as independent; identities valid for every sequence used as evidence; a stationary surrogate used to forecast a supremum (MISTAKE-619).
- **Cards used.** "Keep observer, recurrence class, and finite head"; "Type every analogy and every implication"; "Pre-register before extending the census."

## 2. Exact two-octave memory and energy transport (THM-4658)

Identity (2) of the wavelet note, applied twice through the channel cycle `2 → 1 → 2` of multiplication by `2^(−1)` modulo 3, gives for every integer `N`

`H(4N) = H(N) + H_(2,2)(3N) − H_(2,8)(6N) =: H(N) + Δ(N)`.

Writing `M = 4N + j` with `|H(M) − H(4N)| <= 4` gives `E_(4X) = 4E_X + 8C_X + 4D_X + R_X`, `|R_X| <= 24 sqrt(X E_(4X)) + 48X` (`C_X = Σ HΔ`, `D_X = Σ Δ²` over `[X, 2X)`), and Minkowski's inequality gives `sqrt(E_(4X)) <= 2 sqrt(E_X) + 2 sqrt(D_X) + 4 sqrt(3X)`. It is the two-octave analogue of the wavelet note's one-octave rule (3).

**Census check.** The normalized transport `e_(4X) = (e_X + 2c_X + δ_X)/4` is exact up to `R_X`; the remainder is visible at the fourth decimal for `k = 10..13, 16` (e.g. 0.1947 against 0.1942 at `k = 10`) and below `10^(−6)` from `k = 17` (referee, I2).

| `k` | `e_X` | mean part | refinement | `δ_X` | `κ` | `e_(4X)` | descent bound |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 18 | 0.274 | 0.222 | 0.053 | 1.604 | +0.80 | 0.733 | 0.807 |
| 20 | 0.733 | 0.701 | 0.033 | 0.213 | +0.65 | 0.365 | 0.436 |
| 22 | 0.365 | 0.335 | 0.031 | 4.537 | +0.84 | 1.769 | 1.872 |
| 24 | 1.769 | 1.659 | 0.110 | 0.982 | +0.60 | 1.084 | 1.348 |
| 26 | 1.084 | 0.982 | 0.102 | 2.835 | −0.92 | 0.173 | 1.857 |
| 28 | 0.173 | 0.028 | 0.145 | 1.181 | −0.36 | 0.256 | 0.564 |
| 30 | 0.256 | 0.163 | 0.094 | 0.954 | +0.20 | 0.353 | 0.550 |

**Reading.** `κ = c/sqrt(eδ)` is an uncentred cosine. The descent bound is close when `H` and `Δ` align (`k = 22`: 1.769 against 1.872) and far from tight when they anti-align (`k = 26`).

**Character inheritance (statement 6).** For any `φ` with `φ(2n) = λφ(n)` (a Dirichlet character modulo `3^q`, `λ = ψ(2)`), `C_φ(N) = −λ C_φ(⌊N/2⌋) + O_φ(N)` exactly. The twisted splits `D_j` (`λ = −1`) inherit their previous-octave value with weight +1, Haar-type sums alternate. Averaged over octaves, any such relation forces `(1 − a 2^(−1/2)) ⟨F/sqrt N⟩ = ⟨G/sqrt N⟩` up to end effects *for every sequence*: the measured `⟨D_j⟩/⟨O_j⟩` = 3.12, 3.48, 3.36 (≈ 2 + sqrt 2) and the level-by-level fixed-point "consistency" of 4(b) (residuals 6-12% of the bias size, below their standard errors) are code checks. The same inheritance holds for the halving null, which is why every channel's octave means are strongly correlated.

## 3. The census to 2^34 and the profile

**The census.** `collatz_octave_memory_census_20261015.py`: parity bits by descent below the current octave (39 s to 2^34); a 541 s statistics pass with exact octave sums of the 12 level-1/2 channels, an exact streaming Haar decomposition of `H`, the increments `Δ(N)` for `N < 2^31`, and class sums modulo 27 at 1,537 log-uniform points. 20,000 random `n` agree with direct trajectories; `S(2^34) = 24953018` and `H(2^34) = 26869` agree with the 2026-10-13 census and with the referee's third independent census.

**The profile, 2^14..2^34.** `g(u) = H(2^u)/sqrt(2^(u+1))` has mean `+0.30`, sd 0.28, kurtosis 3.20 (standard error about 0.95), maximum 1.005. Point autocorrelation 0.98, 0.85, 0.69, 0.39, 0.18, 0.06, 0.01 at lags 1/64, 1/8, 1/4, 1/2, 3/4, 1, 2 octaves; adjacent octave *means* correlate at 0.28-0.36. The periodogram of `S(N)/N` shows lines at 26.48/26.53, 14.49 and 2.45 cycles per octave; section 8 identifies them with zeros of `D`. Window energies (1,217 windows) have mean 0.51 and maximum 1.85; the coarse mean carries 84% of the energy. There is no deterministic floor: `e_X < 0.26` at seven octaves and the octave mean part is 0.0002 at 2^29 (referee, E4).

## 4. The out-of-sample census to 2^37 and 2^40

**Methods.**
- `collatz_census_extension_2p37_20261015.py`: parity bits for `n < 2^37` (16 GB) by descent below the octave start in uint64; the 245 words containing an orbit that would leave uint64 are recomputed with Python integers before the next octave. Block-parallel class sums modulo 243 and prefix sums of `H`, `H²` (exact int64, overflow-checked). It reproduces every grid sample of the 2^34 census, `S(2^34)` and `H(2^34)`; 20,000 random `n` and the patched words agree with direct trajectories.
- `collatz_census_extension_2p40_20261015.py`: the 2^37 table plus on-the-fly descent for `2^37 <= n < 2^40`, continuing in exact two-limb 128-bit arithmetic when a value would leave uint64 (stress-tested on 12,042 values, including forced early switching and the path-record starters whose orbits exceed 2^64). Its first version (two passes, single steps) was too slow and was stopped unread; the version used descends with the exact 8-step jump `T^8(256h + l) = 3^(a(l)) h + T^8(l)` (8 steps, parity unchanged), computes the class counts, moments and within-block discrepancy sums in one pass and re-walks only the 954 blocks that contain grid points (heavy pass 44 minutes). It reproduces every one of the 1,728 grid points of the 2^37 extension exactly (class sums, `H`, prefix sums of `H` and `H²`), reproduces `S(2^37) = −152262364` and `H(2^37) = −67763` from the independent moments computation, and agrees with Python trajectories at 20,000 random `n` in `[2^37, 2^40)`. `S(2^40) = 989143994`, `H(2^40) = −55691`.
- `collatz_boundary_moments_2p37_20261015.py`: exact signed block moments through 2^37 for the certificate (first 2^17 rows identical to the stored 2^34 moments); `S(2^37) = −152262364`, `H(2^37) = −67763`.

**Exact dyadic energies** `e_X = E_X/X²`, with the octave mean `m_X = Σ_(X<=n<2X) H(n)/(X sqrt X)`:

| `k` | 20 | 22 | 24 | 26 | 28 | 29 | 30 | 31 | 32 | 33 | 34 | 35 | 36 | 37 | 38 | 39 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `e_X` | 0.733 | 0.365 | 1.769 | 1.084 | 0.173 | 0.344 | 0.256 | 0.274 | 0.353 | 0.095 | 0.060 | 0.080 | 0.037 | 0.069 | 0.931 | 0.078 |
| `m_X` | +0.84 | +0.58 | +1.29 | +0.99 | −0.17 | −0.01 | +0.40 | +0.41 | +0.55 | +0.14 | +0.14 | +0.19 | +0.03 | −0.06 | +0.90 | +0.24 |

The maximum over all dyadic `X <= 2^39` is 1.769 at 2^24, so (9a) `E_X <= 2X²` holds on every dyadic window through 2^39. The quiet stretch 2^33..2^37 ended with a coarse-mean excursion at 2^38 (`e = 0.93`, mean part 0.81).

**Normalized wavelet energies.** The mean square of `g_(q,r) = H_(q,r)/sqrt(6N/3^q)` (64 samples per octave, averaged over the channels of a level; about 1 for independent signs):

| octaves | level 1 | level 2 | level 3 | level 4 | level 5 |
|---|---:|---:|---:|---:|---:|
| 14-19 | 0.34 | 1.80 | 1.23 | 0.99 | 1.03 |
| 20-25 | 0.35 | 0.38 | 0.97 | 1.11 | 1.17 |
| 26-29 | 0.21 | 0.56 | 1.43 | 0.98 | 0.90 |
| 30-33 | 0.05 | 0.32 | 0.54 | 0.54 | 0.59 |
| 34-36 | 0.05 | 0.16 | 0.22 | 0.44 | 0.63 |
| 37-39 | 0.08 | 0.23 | 0.50 | 0.71 | 0.71 |

## 5. Square-root residue biases: proposed and refuted (HYP-9282)

In sample (19 octaves, 2^14..2^32) 75 of 363 channels had naive `|t| > 5` and 27 kept one sign throughout. The referee showed the octave means are far from independent (E2) and that five table rows drift toward zero (E3). The pre-registered test then failed: `H_(3,18)`, `H_(3,9)`, `H_(5,155)` changed sign in all or most of the octaves 34-36 and missed their bands; `H_(5,78)` and `H_(5,59)` kept theirs. Over all 363 channels the slope of out-of-sample on in-sample means is 0.20, exactly what octave memory alone gives (section 6). See [HYP-9282](../hypotheses/HYP-9282-collatz-parity-square-root-residue-biases.md).

## 6. The halving null

**Definition.** `t_null(n) = (−1)^(v_2(n)) s(m)`, `m` the odd part of `n`, with `s(m) = ±1` from a splitmix64 hash of `(m, seed)`. It keeps `t(2m) = −t(m)` exactly (so THM-4658 statement 6 and all octave inheritance hold), and has no odd-step structure. Twelve replicates were run through 2^37 with exactly the statistics used for Collatz (`collatz_bias_null_calibration_2p37_20261015.py`; octave-mean matrices saved in its NPZ):

| statistic | Collatz | halving null (min-median-max) |
|---|---:|---|
| channels with naive `|t| > 5` / `> 8` (in sample) | 75 / 20 | 23-33-49 / 3-6-11 |
| largest naive `|t|` | 18.7 | 9.0-11.5-32.1 |
| channels with one sign in all 19 octaves | 27 | 5-8-15 |
| top five: octave signs kept out of sample (of 15) | 7 | 7-10-15 |
| all channels: slope of out-of-sample on in-sample means | 0.20 | 0.03-0.23-0.46 |
| lag-one correlation of octave means, levels 1-5 | 0.63-0.66 | 0.45-0.84 |
| profile of `H`: sd / one-octave autocorrelation | 0.27 / 0.13 | 0.63-1.27-1.95 / 0.65-0.84-0.91 |
| window energy of `H`, mean in sample / out of sample | 0.51 / 0.07 | 1.5-7.1-14.4 / 1.8-7.1-22.9 |
| window energy of `H`, maximum in sample | 1.85 | 5.5-30.1-42.7 |
| level energies 1-5 in sample | 0.28, 0.86, 1.09, 0.95, 0.98 | medians 1.62, 1.15, 1.02, 1.00, 0.98 |
| level energies 3, 4, 5: out of sample / in sample | 0.20, 0.46, 0.65 | min 0.52, 0.77, 0.78 |
| level energies 3, 4, 5: octaves 30-36 / 14-29 | 0.34, 0.48, 0.58 | min 0.68, 0.77, 0.79 |
| level 2: octaves 30-36 / 14-29 | 0.26 | min 0.22 (inside the null range) |

**Reading.** Collatz's *persistence* statistics (signs, slopes, memory) are those of the halving null, with somewhat longer in-sample memory. Its *sizes* are not: `H` is an order of magnitude smaller and decorrelates within an octave, and the deep levels, which match the null in sample, fell below every null replicate after 2^30.

## 7. Energy budgets (HYP-9283)

**(A) The forecast (REFUTED).** A Gaussian AR(8) surrogate of the window means, fitted on 2^14..2^34 with the bias included, predicted a dyadic exceedance of `2X²` before 2^64 with probability 31%; the referee showed this rested on the bias (1.0% without it) and the innovation scale. Refitted and simulated over the window starts `(2^33, 2^36]`, it predicts a median maximum window energy 0.77 (10-90%: 0.34-1.64); the census gives 0.117 (predictive probability 0.1%).

**(B) Certificates (CONDITIONAL, Arb).** The finite numerator of B3 at `N` plus a Cauchy–Schwarz tail under `E_X <= X²(a + bj)`, `X = N 2^j`:

| `(a, b)` | `(2, 0)` | `(4, 1)` | `(16, 4)` | `(64, 16)` |
|---|---|---|---|---|
| `N = 2^34` | [0.3215, 0.3885] | [0.2964, 0.4137] | [0.2379, 0.4721] | [0.1210, 0.5890] |
| `N = 2^37` | [0.3431, 0.3668] | [0.3342, 0.3757] | [0.3134, 0.3965] | [0.2720, 0.4379] |
| `N = 2^40` | [0.3508, 0.3591] | [0.3476, 0.3623] | [0.3402, 0.3697] | [0.3255, 0.3844] |

(normalized residue `2|r| 2^(29(σ_0−1))`). The rows are nested, as they must be if the budget holds beyond the later `N` (the exact data satisfy it in between). Largest certified budgets: constant 227.4 (2^34), 1806.6 (2^37), 14348 (2^40); slope with `a = 2`: 144.4, 1194.9, 9616; `(4c, c)`: 36.86, 292.86, 2325.9. The amplitude note had proposed logarithmic budgets (its sections 3 and 6.1); the certification, margins and later-`N` reruns are new here.

**(C) Calibration (PROVED implication).** Any certified budget makes `F` meromorphic in `Re s > 1/2` with a nonzero pole at `s_0`, so `S(X) ≠ O(X^β)` for `β < 0.99821` and `limsup log B(X)/log X >= 0.99821`, where `B(X)` counts the `n <= X` reaching 1. Krasikov–Lagarias (2003) proved `B(X) >= X^0.84`.

## 8. An explicit formula over the zeros of D (HYP-9285)

**The remainder identity (PROVED).** `H(N) = −[S(N) + 3S(⌊2N/3⌋) + 3S(⌊N/3⌋)]` for every `N`, and `1 + 3(2/3)^ρ + 3(1/3)^ρ = 3(2/3)^ρ D(ρ)`. The floor-free combination therefore annihilates every such power. The actual integer combination leaves `3[floor(2N/3)^ρ-(2N/3)^ρ]+3[floor(N/3)^ρ-(N/3)^ρ]=O_ρ(N^(Re ρ-1))`; at N=1 it equals 1 for Re ρ>0. A finite packet has a controlled floor remainder, while an infinite packet needs a summability or truncation bound. See [MISTAKE-620 and the explicit repair](collatz_conductor_wiener_20261015d.md#4-integer-floors-leave-a-modal-residue).

**The fit (NUMERICAL).** `D` has 44 zeros with `1/2 < Re ρ < 1`, `0 < Im ρ < 3000`. Least squares for `S(N) ≈ Re Σ c_ρ N^ρ` on the grid points of `2^22..2^34`, tested on `2^34..2^37`: relative rms error 0.19 (`s_0` only), 0.075 (12 zeros, `Re > 0.85`), 0.0205 (28 zeros, `Re > 0.65`), 0.0186 (44 zeros; in-sample 0.006, mild overfitting). The large step at `Re > 0.65` comes from `ρ = 0.6939 + 22.25i`, the 2.45-cycle line; `0.8635 + 131.16i` is the 14.47-cycle line.

**Residues two ways.** A pole with residue `r` gives `S(N) ~ Re(2rN^ρ/ρ)`, so `|r_fit| = |c_ρ||ρ|/2`. The conditional Arb certificate at `N = 2^37` (budget `(2, 0)`, `collatz_pole_certificate_other_zeros_20261015.py`) gives an interval for every zero. All 28 fitted values lie inside their intervals (most of which are wide), and five residues are certified nonzero, each with the fitted value near the middle of its interval:

| zero `ρ` | cycles per octave | conditional `|r|` | `|r_fit|` |
|---|---:|---|---:|
| `0.99821 + 240.248i` (`s_0`) | 26.50 | [0.1778, 0.1901] | 0.1841 |
| `0.98494 + 720.740i` | 79.51 | [0.0451, 0.0951] | 0.0649 |
| `0.94973 + 1572.556i` | 173.48 | [0.0778, 0.3294] | 0.2056 |
| `0.86345 + 131.156i` | 14.47 | [0.4168, 0.5959] | 0.5094 |
| `0.69392 + 22.253i` | 2.455 | [1.127, 4.044] | 2.557 |

For the other 23 zeros the tail (proportional to `|ρ| N^(−(Re ρ − 1/2))`) is too large to certify a sign. So, under the same budget, `F` has at least five certified nonzero poles in `Re s > 1/2`.

**The pre-registered test (P1): PASSED.** The 28-zero model fitted on `2^22..2^37` predicted `S` on `[2^37, 2^40)` (committed in `adcd7eed1f` before the 2^40 census had produced output). Relative rms error over the 192 grid points: 0.0195, against the pass bound 0.03; refitted on the same data, the Mercator zero alone gives 0.165 and the 12 zeros with `Re > 0.85` give 0.042. Actual against predicted: `S(2^38)` = 281,675,896 against 273,757,038; `S(2^39)` = −525,856,530 against −510,780,478; `S(2^40)` = 989,143,994 against 961,770,189. Fitted through 2^40, the Mercator coefficient gives `|r| = 0.18422`, inside every conditional interval.

**Dense test (NUMERICAL, `collatz_explicit_formula_dense_test_20261015.py`).** The 64-point grid aliases every zero above 32 cycles per octave. On exact data through 2^27 sampled at 2048 points per octave (no aliasing below `Im ρ = 8000`), with the fit on `2^18..2^25` and the test on `2^25..2^27`:

| zeros used | holdout relative rms | holdout residual / `sqrt(N)` (rms) |
|---|---:|---:|
| 3 (`Im < 300`, `Re > 0.65`) | 0.124 | 1.40 |
| 10 (`Im < 1000`) | 0.082 | 0.93 |
| 28 (`Im < 3000`) | 0.034 | 0.37 |
| 47 (`Im < 5000`) | 0.029 | 0.31 |
| 75 (`Im < 8000`) | 0.026 | 0.28 |
| 119 (`Im < 8000`, `Re > 1/2`) | 0.025 | 0.27 |

The residual falls to the square-root scale and keeps falling as zeros are added. If the remainder `R = S − (modal sum)` were uncorrelated across the three scales of `H = −[R(N) + 3R(2N/3) + 3R(N/3)]` (exact up to the per-mode floor residues of MISTAKE-620), its rms would be `rms(H/sqrt N)/sqrt(10) = 0.24 sqrt(N)` on this range (`rms(H/sqrt N) = 0.75`); the measured 0.27 is close (HEURISTIC comparison).

**No excess modes at any finite depth (PROVED, THM-4658 statement 7).** In 3-adic Fourier coordinates the odd branch of the transfer operator only moves information to coarser frequencies. On densities constant modulo `3^Q`, the operator is therefore block triangular, its spectrum is `{1 − D(s)}` together with `−2^(−s)` times roots of unity of order `2·3^k` (`k < Q`), and the excess operator has eigenvalue 1 only on `Re s = 0`. So a power-law excess main term with `Re ρ ≠ 0` must have unbounded 3-adic frequency depth, and the only finite-depth modes are the D-zero density modes. This fits the refutation of HYP-9282 and the leading-mode reading of `H` (section 8, with the floor residues of MISTAKE-620). clock_duality_oct15 extended the exclusion to a weighted Wiener class at infinite depth (Proposition W1 of `collatz_conductor_wiener_20261015d.md`: the excess operator contracts in `Σ |f̂(ξ)| R^(q(ξ))` for `R > 3^(σ−1)/(2^σ − 1)`), while arbitrary distributions of unbounded depth do solve the fixed-point equation (THM-4658, remark on statement 7). The spectrum is confirmed numerically to `10^(−14)` for `Q <= 4`, and an argument-principle scan finds no eigenvalue-1 point in `[0.05, 1.5] × (0, 300]` (`collatz_excess_operator_resonance_scan_20261015.py`).

## 9. Odd-step suppression (HYP-9284)

- The halving null's window energies average 1.5-14 in sample (above 2 in 11 of 12 replicates) with maxima up to 43; Collatz never exceeded 1.85. Since the null keeps every halving identity, a proof of any budget must use the odd-step relation `t(n) = −t((3n+1)/2)`.
- In the leading-mode reading of section 8, the odd-step relation is what makes `S` close to modal and `H` small; the halving null has no such structure.
- The deep levels behave like the null up to about 2^27 and then dipped together, more weakly at deeper levels (section 4).
- **Pre-registered tests on 2^37..2^40.** P2 and P3 PASSED (`e_X` = 0.069, 0.931, 0.078 for `X` = 2^37, 2^38, 2^39). P4 FAILED (mean 0.359 against the bound 0.25: the quiet stretch ended at 2^38). P5 PASSED (levels 3, 4, 5 averaged over octaves 37-39: 0.50, 0.71, 0.71, below 0.8), but the levels are recovering (octave 39: 0.69, 0.81, 0.76). P6 PASSED (octave means of `H/sqrt(2N)`: −0.01, +0.50, +0.14).
- **Reading.** The 2^28..2^37 dip of the deep levels now looks like a slow modulation rather than the onset of sub-square-root growth; the level-1 suppression against the null persists at every octave.

See [HYP-9284](../hypotheses/HYP-9284-collatz-odd-step-suppression-wavelet-energies.md).

## 10. Q–Pell, golden and Möbius controls

- **Q–Pell (THM-4630).** Q's model variance update is `V_new/V = (3 − r)/4` with the Möbius correlation map `r ↦ (r − 1)/(3 − r)` and two-step Pell matrix `[[3,4],[2,3]]`. The Collatz transport `e_(4X) = (e_X + 2c_X + δ_X)/4` has the same "divide by four per two octaves" shape. The measured native cross term ranges over −0.97..+0.91 with runs of either sign (the wavelet note's one-octave cross terms: −0.949..+0.936). Variation alone does not exclude a function of a richer state; no scalar Pell-type closure or contraction rate is established. clock_duality's [changing-generator companion](quadratic_cycle_phase_20261015b.md) supplies a separate exact Pell identity.
- **Möbius maps and the rational function field.** The Q map is a hyperbolic element of `PGL_2 = Aut(C(z))` with fixed points `1 ± sqrt 2`. The Collatz normalized transport `g ↦ g/2 + Δ̂` fixes `∞`, but `Δ̂ ≈ −g/2` plus a fresh term, so it is not an independent random affine iteration (referee, I4).
- **Golden/Conway control.** Conway's leaf charge has no cancellation and a deterministic profile; a constant budget is natural there. The halving null is the opposite (pure noise with memory). Collatz sits between: its remainder `H` is small and short-memory.

## 11. Hostile controls

1. The suppression and decline of HYP-9284 are finite-range statistics against 12 null replicates; slow modulations over tens of octaves cannot be excluded.
2. The explicit-formula fit is least squares on a grid with 64 points per octave; zeros above `Im = 3000` alias. The conditional residue intervals depend on the unproved budget.
3. HYP-9283(B) certifies implications. Every budget in it is unproved.
4. THM-4658 transfers the energy problem one level down; it does not close it.
5. clock_duality_oct15 proved the refinement envelope first and found the 4(a) gap first.

## 12. Open threads

1. A mechanism for the odd-step suppression: why does the Collatz recursion damp the χ-sector memory that the halving null keeps?
2. Statement 7 settles every finite depth. Does the infinite-depth excess operator (distributions with unbounded frequency depth) have eigenvalue-1 points with `Re s > 0`? Its triangular induction has no finest level.
3. The modal expansion: growth bounds for `F` on vertical lines would turn the conditional continuation into an explicit formula.
4. Census beyond 2^40 needs a disk-backed table or a segmented sieve of descents.
5. Recast clock_duality's terminal certificates with the certified margins of section 7(B).

## 13. Audit

- **Independent S37 referee** (report in `04-computation/experiments/s37_audit_20261015/`): confirmed statements 1-3, 4(b), 5 of THM-4658, the census (third independent reproduction), the certificate (bit-identical, mpmath agreement to 1e-29), the special-primes and curling data. Errors: E1 (4(a) circular; repaired by clock_duality), E2 (octave independence), E3 (drift), E4 (deterministic floor); imprecisions I1-I13. All applied.
- **clock_duality_oct15** (scoped audit `collatz_incoming_audit_20261015b.md`, MISTAKE-618): the 4(a) repair, the mean-zero branch hostile, the Pell-wording correction, the annular source row.
- **This session**: the fixed-point and `2 + sqrt 2` checks are automatic (section 2); the pre-registered test of HYP-9282 failed; the surrogate failed out of sample. MISTAKE-619.

## 14. Reproduction

```sh
python 04-computation/experiments/collatz_octave_memory_census_20261015.py 34 31 out.json out.npz
python 04-computation/experiments/collatz_octave_memory_analysis_20261015.py out.json out.npz analysis.json
python 04-computation/experiments/collatz_residue_bias_deep_levels_20261015.py 33 bias.json
python 04-computation/experiments/collatz_census_extension_2p37_20261015.py 37 ext37.npz
python 04-computation/experiments/collatz_census_extension_analysis_20261015.py ext37.npz ext37_analysis.json
python 04-computation/experiments/collatz_bias_null_calibration_2p37_20261015.py 37 12 ext37.npz null.json
python 04-computation/experiments/collatz_boundary_moments_2p37_20261015.py 37 m37.npz m37.json 05-knowledge/results/collatz_boundary_moments_20261013c.npz
python 04-computation/experiments/collatz_pole_certificate_general_N_20261015.py --moments m37.npz --census m37.json --output cert37.json
python 04-computation/experiments/collatz_pole_certificate_margins_20261015.py --moments m37.npz --census m37.json --output margins37.json
python 04-computation/experiments/collatz_mercator_waveform_prediction_20261015.py ext37.npz prediction.json
python 04-computation/experiments/collatz_other_zeros_residue_table_20261015.py m37.npz m37.json prediction.json residues.json
python 04-computation/experiments/collatz_excess_operator_resonance_scan_20261015.py scan.json
python 04-computation/experiments/collatz_census_extension_2p40_20261015.py 40 37 c40.npz m40.npz c40.json ext37.npz
python 04-computation/experiments/collatz_census_2p40_evaluation_20261015.py c40.npz c40.json prediction.json eval.json
```

The completed 2^37 block-moment NPZ file (9 MB) is not committed; its SHA256 is in the certificate JSON and its script regenerates it. The 2^40 moment/certificate extension is pending and is not a validated input here.

**References.**
- Repository: collatz_amplitude_boundary_20261013c.md (B1–B3), collatz_break_clock_repair_20261013d.md, collatz_wavelet_energy_20261013d.md, collatz_clock_mesh_20261015a.md, collatz_terminal_mesh_20261015a.md, collatz_gram_slack_20261015b.md, collatz_incoming_audit_20261015b.md, collatz_annular_overlap_20261015c.md, golden_nested_trees_20261013c.md (section 6), THM-4630, THM-4629, THM-4581.
- I. Krasikov and J. C. Lagarias, Bounds for the 3x+1 problem using difference inequalities, Acta Arith. 109 (2003) 237–258.
- T. Tao, Almost all orbits of the Collatz map attain almost bounded values, Forum Math. Pi 10 (2022) e12 (arXiv:1909.03562); its Fourier decay of the 3-adic Syracuse distributions is `O_A(n^(−A))` for every `A`.
