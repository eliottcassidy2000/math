# Mean currents, refinement energy and square-root biases of the Collatz parity census: two-octave memory, log budgets for the Mercator pole, and what a proof would have to beat

2026-10-15, session opus-2026-10-15-S37 (worktree `/tmp/math-wt-octave-rg-37`). Pre-audit draft.

**Owner's prompt (verbatim).** "keep working for a long session of next targets and think the mean-current and refinement-energy bounds, preserving source and integer clock. pull in recent repo work and use Q–Pell and golden-tree controls to close the dyadic energy gaps, merge in new ideas on these and other problems and pursue whatever most mathematically relevant for major proof progress. consider the connection between the geometric Langlands correspondence and S-duality, a property of certain quantum field theories. also think of the field of rational functions on the Riemann sphere and the gijswijt sequence. also merge in the 15 supersingular primes, stern sprimes, non-generous primes, wagstaff primes, woodall primes and cullen primes and leyland primes. adjust targets as necessary and be free in your creative exploration and conjecture proposal and investigation"

The special primes, curling numbers and Langlands threads are in the companion note [special_primes_curling_genus_zero_langlands_20261015.md](special_primes_curling_genus_zero_langlands_20261015.md).

**Status.**
- **PROVED (THM-4658):**
  - the class-sum recursion in measure form;
  - two-octave memory `H(4N) = H(N) + H_(2,2)(3N) − H_(2,8)(6N)`;
  - the exact energy transport and the 4-adic energy descent;
  - power-law main terms are eigen-distributions of a 3-adic transfer operator, whose branches are balanced in `L²(Z_3)` exactly on `Re s = 1/2` and `Re s = 1`;
  - the calibration implication of HYP-9283(C).
- **CONDITIONAL, interval-certified (HYP-9283(B)):** the nonzero Mercator residue under logarithmic energy budgets up to `E_X <= X²(64 + 16 log_2(X/2^34))`.
- **FINITE-EXACT:**
  - an independent census of `t(n) = (−1)^τ(n)` through `2^34`, reproducing `S(2^34) = 24953018` and `H(2^34) = 26869` exactly;
  - octave means, energies, the Haar refinement spectrum and level-2 increments;
  - class sums modulo 243 through `2^33`.
- **NUMERICAL:**
  - decorrelation and cancelled memory of the normalized profile;
  - square-root residue biases (HYP-9282);
  - a Gaussian surrogate predicting that constant budgets fail (HYP-9283(A), HEURISTIC).
- **OPEN:** any energy budget for all large `X` (calibrated as beyond Krasikov–Lagarias), the bias mechanism, and the Collatz conjecture.

## 0. Answers in brief

1. **Refinement energy is paid; the mean current is the obligation.** The concurrent session `clock_duality_oct15` (credited throughout) proved a sharp interval envelope that pays every scale below a square-root mesh (M1/M2), together with source-preserving terminal certificates (T1/T2). What remains is a bound on the coarse means, their target "(5)": `C(P_b) <= (4/9)X²` for all dyadic `X >= 2^34`.
   - THM-4658 adds the exact mean-current transport: the octave profile of `H` at `4X` equals the profile at `X` plus a level-2 increment. A bound on that increment's energy bounds the top energy with the same constant (4-adic descent). This keeps the integer clock `N ↦ 4N, 3N, 6N` exactly and introduces no new source.
2. **The memory is real but cancelled.** In the census the level-2 increment almost exactly cancels the inherited value: the regression slope of `H(4N)/sqrt(4N)` on `H(N)/sqrt(N)` is 0.006. The normalized profile behaves like a short-memory stationary process in `log N`, with Gaussian-like marginal and lag-one-octave correlation 0.06.
3. **Constant budgets are the wrong target.** A stationary Gaussian-like process exceeds every constant eventually. A surrogate fitted to the census puts the probability that some dyadic `E_X` exceeds `2X²` before `2^64` at about 31%, and before `2^200` at about 88% (HYP-9283(A), HEURISTIC). So hypothesis (9a), and the constant coarse target (5), are expected to fail eventually.
4. **But the conclusion survives logarithmic budgets.** With the same finite numerator, the nonzero Mercator residue is interval-certified under `E_X <= X²(a + b log_2(X/2^34))` for `(a, b)` up to `(64, 16)`, giving `[0.1210, 0.5890]` (HYP-9283(B)). The dyadic-energy "gap" is therefore a gap in the formulation of the target. The repaired target, a log-growing budget, is consistent with the data by a wide margin.
5. **How hard is the repaired target?** Any such budget forces `limsup log B(X)/log X >= 0.99821` for the count `B(X)` of integers reaching 1 (HYP-9283(C)). The best unconditional exponent is 0.84 (Krasikov–Lagarias). The packet route of the wavelet note would need exponential Fourier decay of the 3-adic carry distribution, where Tao (2019) has polynomial decay. This is honest calibration, not a proof.
6. **New structure: square-root residue biases.** The census shows persistent signed main terms of order `sqrt(N)` on many residue classes modulo `3^q`, `q <= 5` (HYP-9282).
   - Examples: `H_(3,18)` and `H_(3,9)` have opposite signs in all 19 octaves, with `t = −18.7` and `+14.5`.
   - `H = H_(1,2)` itself averages `+0.30` in normalized units.
   - By THM-4658 statement 4, such main terms would be a fixed point of the transfer operator at its balanced line `Re s = 1/2`.
7. **Q–Pell and golden controls.**
   - The energy transport `e_(4X) = (e_X + 2c_X + δ_X)/4` has exactly the algebraic shape of the Q covariance model's variance update `V_new = V(3 − r)/4`.
   - The measured cross term varies widely (`κ` swings between −0.97 and +0.91), so these data do not establish the same scalar Möbius/Pell closure. They do not exclude a richer deterministic state or a different closure; the observed cancellation must be retained.
   - The golden/Conway control is the opposite regime: a deterministic profile, so bounded and decaying.

## 1. Inheritance and the live board

- **Programme.**
  - The amplitude note (B1–B3, renewal_shapes_oct11) reduces a nonzero Mercator pole to a discrepancy hypothesis.
  - The break-clock repair note insists on the integer clock and the incoming current, giving the (5a)–(5c) block compiler and mean/innovation split.
  - The wavelet note gives the exact channel refinement (2), the Gram rule and the packet (7).
  - The clock-mesh note (M1/M2) gives the sharp refinement envelope.
  - The terminal-mesh note (T1/T2) gives coarse-mean certificates without grounding every source.
- **Closest proved mechanism.** Identity (2) of the wavelet note.
- **Hostile.** The same-event-word, different-energy pair (417 vs 81) of the break-clock note.
- **Corrected near miss.** An expected energy treated as a uniform bound.
- **Least-used sidecar.** The four-step (two-octave) composite of identity (2), which closes the top channel on itself modulo level 2.
- **Concept board.**
  - ports and incoming current;
  - integer lengths;
  - coarse means;
  - refinement innovations;
  - the 3-adic transfer operator;
  - the cancelled two-octave memory;
  - residue biases.
- **Cards used.** "Keep observer, recurrence class, and finite head"; "Retain local incidence through counting"; "Type every analogy and every implication."

## 2. Exact two-octave memory and energy transport (THM-4658)

Identity (2) of the wavelet note, applied twice through the channel cycle `2 → 1 → 2` of multiplication by `2^(−1)` modulo 3, gives for every integer `N`

`H(4N) = H(N) + H_(2,2)(3N) − H_(2,8)(6N) =: H(N) + Δ(N)`.

Writing `M = 4N + j` with `|H(M) − H(4N)| <= 4` gives the exact transport

`E_(4X) = 4E_X + 8C_X + 4D_X + R_X`, `|R_X| <= 24 sqrt(X E_(4X)) + 48X`,

with `C_X = Σ H(N)Δ(N)` and `D_X = Σ Δ(N)²` over `[X, 2X)`. Minkowski's inequality gives the descent

`sqrt(E_(4X)) <= 2 sqrt(E_X) + 2 sqrt(D_X) + 4 sqrt(3X)`.

So the top-channel energy is controlled, with the same constant in the limit, by the energy of a level-2 quantity. This is a source-preserving rule in the sense of the break-clock note: it never resets the current and never reorders sources.

**Census check.** The normalized transport `e_(4X) = (e_X + 2c_X + δ_X)/4` holds to four decimals on every octave `2^10..2^30`. For example, at `X = 2^22`: `e_(4X) = 1.7688`, prediction `1.7688`.

| `k` | `e_X` | mean part | refinement | `δ_X` | `κ` | `e_(4X)` | descent bound |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 18 | 0.274 | 0.222 | 0.053 | 1.604 | +0.80 | 0.733 | 0.807 |
| 20 | 0.733 | 0.701 | 0.033 | 0.213 | +0.65 | 0.365 | 0.436 |
| 22 | 0.365 | 0.335 | 0.031 | 4.537 | +0.84 | 1.769 | 1.872 |
| 24 | 1.769 | 1.659 | 0.110 | 0.982 | +0.60 | 1.084 | 1.348 |
| 26 | 1.084 | 0.982 | 0.102 | 2.835 | −0.92 | 0.173 | 1.857 |
| 28 | 0.173 | 0.028 | 0.145 | 1.181 | −0.36 | 0.256 | 0.564 |
| 30 | 0.256 | 0.163 | 0.094 | 0.954 | +0.20 | 0.353 | 0.550 |

**Reading.** The descent bound is attained when `H` and `Δ` align (`κ ≈ +1`, as at `k = 22`) and is far from tight when they anti-align (`κ ≈ −1`, as at `k = 26`). The large energies at `2^24` and `2^26` are an aligned increment arriving on an already large profile. The next two octaves show the cancellation.

## 3. The census and the cancelled memory

**The census.** `collatz_octave_memory_census_20261015.py` computes parity bits by descent below the current octave, in parallel, in 39 seconds through `2^34`. Its uint64 guard aborts rather than certifies. A statistics pass of 541 seconds then records:
- exact octave sums and sums of squares for the 12 level-1/2 channels;
- a streaming exact Haar decomposition of `H` on every octave;
- the increments `Δ(N) = H(4N) − H(N)` for `N < 2^31`;
- samples of all class sums modulo 27 at 1,537 log-uniform points.

**Checks.** 20,000 random `n` were checked against direct trajectories, with no mismatch. The Haar identity holds to `1e−16`. `S(2^34)` and `H(2^34)` agree with the repository's independent census.

**The profile.**
- `g(u) = H(2^u)/sqrt(2^(u+1))` on `u ∈ [14, 34]` has mean `+0.30`, standard deviation 0.28, kurtosis 3.20 and maximum 1.005.
- Its autocorrelation is 0.98, 0.85, 0.69, 0.39, 0.18, 0.06, 0.01 at lags 1/64, 1/8, 1/4, 1/2, 3/4, 1 and 2 octaves.
- Despite the exact memory identity, the regression of `g(u + 2)` on `g(u)` has slope 0.006 and correlation 0.007.

So the level-2 increment almost exactly cancels the memory: `Δ(N) ≈ −H(N)` plus a fresh term. This is the quantitative content of "the mean current must be paid at every scale": nothing carries over to help.

**Spectrum and energy split.**
- The periodogram of `g` has 92% of its power below one cycle per octave and no line.
- `S(N)/N` shows the Mercator line at 26.48/26.53 cycles per octave, the base line at 14.49, and a strong line at 2.45 cycles per octave (positive controls).
- Window energies (1,217 windows at 64 phases per octave) have mean 0.51 and maximum 1.85, at `log_2 X = 23.89`, slightly above the dyadic maximum 1.77. The coarse mean carries 84% of the energy and refinement 16%, consistent with the clock-mesh note's split at `2^24` and `2^26`.

## 4. Square-root residue biases (HYP-9282)

The class sums modulo `3^q`, `q <= 5`, sampled through `2^33`, show persistent signed main terms. In units of the independent-sign scale `sqrt(6N/3^q)`, many channels keep one sign in all 19 octaves `2^14..2^32`: `H_(3,18)` (`−2.44`, `t = −18.7`), `H_(5,155)` (`−1.42`, `−14.6`), `H_(5,78)` (`+1.55`, `+14.6`), `H_(3,9)` (`+2.03`, `+14.5`), `H_(5,59)` (`+0.92`, `+14.1`).
- 75 of the 363 channels have `|t| > 5`; 20 have `|t| > 8`.
- The top channel `H` averages `+0.30`, with 17 of 19 octaves positive.
- The splits of the multiples of `3^j` by the next ternary digit give `−0.37, −0.62, +0.70, −0.02` times `sqrt(N)` for `j = 0..3`.
- Relative drifts within ±2.5% per octave are consistent with exponent `1/2 ± 0.04`.

**Interpretation (CONJECTURAL).** The white-sign model gives zero for every channel, and the Mercator modes cancel from every wavelet channel exactly. A nonzero `b_(q,r) sqrt(N)` is therefore a deterministic feature of the Collatz parity at exactly the square-root scale. By THM-4658 statement 4(b), its coefficient distribution `β` on `Z_3` would be a fixed point of the rank-one modified transfer operator `P̃_(1/2)`.

THM-4658 statement 5 shows that `Re s = 1/2` is precisely the line on which the two branches of the 3-adic transfer operator have equal `L²(Z_3)` norms and squared norms summing to 1. That is the variance-neutral line of the white model. The biases sit exactly on the critical line of the operator that propagates them.

**For the energy programme.** The bias contributes a deterministic floor to every octave energy, about `3·0.30² ≈ 0.26` of the mean `e ≈ 0.5`. It is harmless for the pole: a real `sqrt(x)` main term only adds a pole at `s = 1/2`.

## 5. Constant versus logarithmic budgets (HYP-9283)

**(A) Constant budgets fail eventually (HEURISTIC).** The profile in section 3 is a short-memory, Gaussian-like stationary process in `log N`, so its supremum is infinite. A Gaussian AR(8) surrogate (quarter-octave steps, fitted on `2^14..2^34`, bias included, refinement energies resampled) gives:
- P(some dyadic `E_X > 2X²`) from `2^34` to `2^40, 2^50, 2^64, 2^100, 2^200`: 5.7%, 16.4%, 30.7%, 56.6%, 87.7%;
- the same at all window phases: 9.4%, 27.6%, 46.6%, 75.9%, 97.4%.

Hypothesis (9a), the pointwise `1.5 sqrt(x)` bound, and the constant coarse target `(4/9)X²` of the clock-mesh and terminal-mesh notes are therefore expected to fail at some large dyadic `X`. They probably hold well beyond `2^34`, but not for all `X`.

**(B) Logarithmic budgets suffice (CONDITIONAL, interval-certified).** `collatz_pole_certificate_logenergy_20261015.py` reuses B3's exact block moments, root disk, quadrature and `R`-tail, and replaces only the future `E`-tail by a Cauchy–Schwarz sum over dyadic blocks with budget `X²(a + b j)` at `X = 2^(34+j)`. The resulting intervals for `2|r| 2^(29(σ_0−1))`:
- `(2, 0)`: `[0.3215, 0.3885]`, the published constant budget;
- `(4, 1)`: `[0.2964, 0.4137]`;
- `(16, 4)`: `[0.2379, 0.4721]`;
- `(64, 16)`: `[0.1210, 0.5890]`.

All are nonzero. In clock-mesh language, a coarse budget `C(P_b) <= X²(15 + 4j)` plus the proved fine allowance lies inside the `(64, 16)` row.

**(C) Calibration (PROVED implication).** Any of these budgets makes `E` holomorphic in `Re s > 1/2` and `F` meromorphic with a certified nonzero pole at `s_0`. Hence `S(X) ≠ O(X^β)` for `β < 0.99821`, and so `limsup log B(X)/log X >= 0.99821`, where `B(X)` counts the `n <= X` reaching 1. Krasikov–Lagarias (2003) proved `B(X) >= X^0.84` for large `X`. A proof of any energy budget, even logarithmic, would therefore improve the known population exponent along a sequence of scales.

## 6. Q–Pell, golden and Möbius controls

- **Q–Pell (THM-4630).**
  - Q's model variance update is `V_new/V = (3 − r)/4` with a Möbius correlation map `r ↦ (r − 1)/(3 − r)` and a two-step Pell matrix `[[3,4],[2,3]]`.
  - The Collatz energy transport `e_(4X) = (e_X + 2c_X + δ_X)/4` has the same "divide by four per two steps" shape.
  - Q's closure comes from prescribed source correlations. The measured Collatz cross term `κ_X` ranges from −0.97 to +0.91 with runs of either sign. Variation alone does not prove that it is not a function of a richer previous state.
  - No scalar Pell-type closure or contraction rate is established here. A proposed closure must retain the observed cross term and be tested against it.
  - This is an obstruction to the analogy without extra data, not a theorem excluding all Pell-based constructions. The [changing-generator companion](quadratic_cycle_phase_20261015b.md) supplies a separate exact Pell identity without identifying it with native Collatz transport.
- **Möbius maps and the rational function field.** The Q map is an element of `PGL_2 = Aut(C(z))`: a hyperbolic automorphism of the Riemann sphere, with fixed points `1 ± sqrt(2)` and multiplier `(3 − 2 sqrt(2))` at the attractor. The Collatz normalized transport `g ↦ g/2 + Δ̂` is an affine Möbius map fixing `∞`, but with a random translation part. Its stationary law is that of a random affine iteration, which is unbounded when the increments are Gaussian; this is the content of (A). See the companion note for the rational-function-field thread.
- **Golden/Conway control** (golden nested-tree note; `gold_heavy_child_20261014j.md`). Conway's leaf charge has no cancellation, and its profile is deterministic: the integer limit is closed. A constant budget is natural there because nothing random enters. The Collatz census sits in the opposite regime.

## 7. Hostile controls

1. The biases of HYP-9282 are finite-range statistics. A slowly varying fluctuation over 19 octaves cannot be excluded, and the exponent is fixed only to `1/2 ± 0.04`.
2. HYP-9283(A) is a surrogate model, not a theorem; its probabilities depend on the AR fit and on stationarity.
3. HYP-9283(B) certifies an implication. Every budget in it is unproved.
4. THM-4658 transfers the energy problem one level down; it does not close it.
5. The concurrent clock-duality session proved the refinement envelope first. This note's contributions are the mean-current transport, the data, the reformulated target and the calibration.

## 8. Open threads

1. A mechanism for the square-root biases, perhaps an explicit `P̃_(1/2)`-fixed distribution built from the root forcing `t(1) = 1`. Test against levels `q <= 5`.
2. Census octaves `2^35..2^40`, to test the pre-registered bias signs and the surrogate's 6% exceedance prediction.
3. Recast the terminal-controller certificates of the clock-duality notes with logarithmic budgets. The finite certificates already found (`2^10..2^20`) are far inside any such budget.
4. A transfer inequality for the joint vector of level-1 and level-2 octave means that captures the cancellation `Δ ≈ −H`.

## 9. Audit

The original independent session audit is pending. A concurrent Codex proof audit
repaired statement 4(a) of THM-4658: a scalar main term does not alone justify
Haar-leading residue-class main terms; `H(N)=o(N^σ)` is now explicit. It also
narrowed the Q--Pell nonexistence wording above and supplied a mean-zero
finite-cylinder hostile to confusing branch squared norms with a norm of the
sum. These scoped checks are recorded in
[the incoming-work audit](collatz_incoming_audit_20261015b.md); they are not
an audit of the entire numerical census or surrogate model.

## 10. Reproduction

```sh
python 04-computation/experiments/collatz_octave_memory_census_20261015.py 34 31 out.json out.npz
python 04-computation/experiments/collatz_octave_memory_analysis_20261015.py out.json out.npz analysis.json
python 04-computation/experiments/collatz_residue_bias_deep_levels_20261015.py 33 bias.json
python 04-computation/experiments/collatz_pole_certificate_logenergy_20261015.py --moments 05-knowledge/results/collatz_boundary_moments_20261013c.npz --census 05-knowledge/results/collatz_discrepancy_20261013c.json --output cert.json
```

Outputs in `05-knowledge/results/`: `collatz_octave_memory_census_20261015.{json,npz,out}`, `collatz_octave_memory_analysis_20261015.{json,out}`, `collatz_residue_bias_deep_levels_20261015.json`, `collatz_pole_certificate_logenergy_20261015.json`.

**References.**
- Repository: collatz_amplitude_boundary_20261013c.md (B1–B3), collatz_break_clock_repair_20261013d.md, collatz_wavelet_energy_20261013d.md, collatz_clock_mesh_20261015a.md, collatz_terminal_mesh_20261015a.md, golden_nested_trees_20261013c.md (section 6), THM-4630, THM-4629, THM-4581.
- I. Krasikov and J. C. Lagarias, Bounds for the 3x+1 problem using difference inequalities, Acta Arith. 109 (2003) 237–258.
- T. Tao, Almost all orbits of the Collatz map attain almost bounded values, Forum Math. Pi 10 (2022) e12 (arXiv:1909.03562).
