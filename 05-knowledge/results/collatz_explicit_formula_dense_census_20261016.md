# The explicit formula of the Collatz parity sum on a dense census: the Mercator family of zeros made uniform, aliasing at resonant grid densities, a second pre-registered test, and the conditional certificate at the 1054/665 zero

2026-10-16, session opus-2026-10-16-S38 (worktree `/tmp/math-wt-explicit-modal-38`). Continues
[the S37 note](collatz_octave_memory_bias_energy_budgets_20261015.md) (HYP-9283, HYP-9284, HYP-9285). Audit: independent
adversarial referee, [report](../../04-computation/experiments/s38_audit_20261016/REPORT.md): the mathematics of THM-4659
and THM-4660 confirmed; credit, aliasing scope, residue numbers and census reliability corrected below (MISTAKE-623).

**Prior art (credited).** The family equation of THM-4659, the parity rule `q` odd / `p` even, the leading constant `C`
and the zero `ρ*` itself (the 332.5-cycle line, decay `5.4756e-7` per octave) are in Codex's
[resonance-pair note](collatz_signed_resonance_pair_mass_20261013b.md) (equation (9) and its root table) and in THM-4629
statement 3; HYP-9230 and HYP-9260 statement 4 carry the same constant. The identity `F D = G` and the conditional
continuation of THM-4660 are Codex's P4 and renewal_shapes_oct11's Lemma B1 / Theorem B3; the boundary residue
estimator used in section 3 is renewal_shapes_oct11's formula (4). clock_duality_oct15's W1/W2 use the same function
(W1's `γ_1(σ)` is `f(σ)`, W2 uses `f'(1) = −κ/2`). THM-4632 is the line-side counterpart of THM-4659 statement 4.

**Owner's prompt (verbatim).** "keep going with the next targets for another long session"

**Targets taken from S37's open obligations.** (i) the modal expansion of `S(N)` (HYP-9285): denser data, more zeros,
a second pre-registered out-of-sample test; (ii) the residues at many zeros, two ways; (iii) the growth of the boundary
numerator on vertical lines; (iv) whether the conditional pole certificate can reach zeros closer to `Re s = 1` than the
Mercator zero.

**Status.**
- **PROVED:** THM-4659, the uniform version of the prior-art family statement. The zeros of `D(s) = 1 + 2^(−s) + (3/2)^s/3`
  near `Re s = 1` are exactly the parity-matched near-returns `2^p ≈ 3^q`, `p` even, `q` odd, uniformly in `q`;
  `1 − Re ρ = Cε² + O(ε⁴)` with `C = π² ln²2/ln³(4/3) = 199.1647`; no zeros with `1 <= Re s <= 2`; approximate aliasing
  at resonant grid densities. Unconditional lemmas of THM-4660 on `D` (uniformly simple zeros, bounded zero density,
  a uniform lower bound off the zeros).
- **CONDITIONAL (complete proof):** THM-4660. If `H` has dyadic energy `O(X^(2+ε))` for every `ε > 0` (a square-root
  hypothesis weaker than every budget of HYP-9283), then `F(s) = Σ t(n) n^(−s)` is meromorphic in `Re s > 1/2` with
  simple poles at zeros of `D` (the continuation is the prior-art one, under the weaker hypothesis), and the second
  Riesz mean of `S` has the explicit formula (new)
  `Σ_(n<=x) t(n)(1 − n/x)² = Σ_(Re ρ > 1/2+δ) 2r_ρ x^ρ/(ρ(ρ+1)(ρ+2)) + O(x^(1/2+δ))`.
- **FINITE-EXACT:** a dense census (1024 log-uniform points per octave) through 2^40, reproducing the S37 census at all
  1,920 shared grid points (a replication by the same code) and, below 2^34, the referee's independent census at all
  24,489 grid points; an irrational-density census (1000 sqrt 2 per octave) through 2^38. The resumed extension through
  2^42 is a single run: see section 5 for its integrity checks.
- **NUMERICAL:** the zeros of `D` with `1/2 < Re < 1`, `0 < Im < 20000` (295; floating-point argument principle,
  independently re-counted by the referee).
- **NUMERICAL:** residues at all 295 zeros from the exact 2^40 block moments; agreement with the dense least-squares
  fit for the dominant low modes; the holdout calibration of the explicit formula on the dense grid (0.17% with the
  69 zeros below the Nyquist height); the pre-registered test on `[2^40, 2^42)`: five of six tests PASSED (M2 predicts
  two new octaves to 0.115%), the energy bound `e_X <= 1` FAILED at `2^41` (section 5); the second-stage test on
  `[2^42, 2^43)` with replication is PENDING.
- **CONDITIONAL (Arb):** the residue at `ρ* = 0.99999921 + 3014.0202i` (the `1054/665` member, prior art: Codex's
  332.5 line) is `0.00993 + 0.00113i` before the future tail, with radius 0.0261 under the budget `E_X <= 2X²` for
  `X >= 2^40`. The radius scales like `N^(−1/2)` (0.0130 at 2^42, 0.0092 at 2^43), so the certificate cannot succeed
  before 2^43, and at 2^43 the expected margin is about 8%.

## 1. The zeros of D in the critical strip up to height 20000

`collatz_D_zeros_and_boundary_residues_20261016.py` finds 295 zeros with `1/2 < Re ρ < 1` and `0 < Im ρ < 20000`
(vectorised Newton plus a dense local search where needed); the argument principle on the forty rectangles
`[1/2, 1] × [500k, 500(k+1)]` confirms the count in every rectangle (floating point: NUMERICAL; the referee's lattice
enumeration gives the same count). About 14.75 zeros per unit 1000 of height.

**Three structures.**
1. **The Mercator family near `Re s = 1` (THM-4659; prior art as credited above).** 23 zeros with `Re > 0.98` below height 20000, all at
   frequencies within 0.012 of half-integers (26.5036, 79.5105, ..., 332.5001, ... cycles per octave). They are the
   parity-matched near-returns: `84/53` gives `s_0 = 0.99820922 + 240.248i`, `1054/665` gives
   `ρ* = 0.99999921 + 3014.0202i`, its odd multiples give `0.99999289 + 9042.06i` and `0.99998025 + 15070.10i`, and the
   other members sit around `2 Im ρ*` and its neighbours, at near-returns such as `(2·1054 ± 84, 2·665 ± 53)` (e.g.
   `2533.5`, `3494.5`, `5787.8`, `6268.3`). `(1 − Re ρ)/ε²` is 197.3 at `84/53` and 199.1638 at `1054/665` against
   the limit `C = 199.1647` (approached from below, as the `ε⁴` coefficient `+206120.8` requires). The next matched convergent `50508/31867` predicts a zero `0.999999978121 + 144432.720754i`.
   The zero `0.693918 + 22.2532i` behind the 2.45-cycle line of `S(N)/N` is the analytic continuation of the family to
   the convergent `8/5`.
2. **A near-1/2 family.** Zeros at `Re = 0.500806, 0.5012, 0.5004, 0.5016, 0.50006, 0.5019` with heights 1527.9,
   4500.1, 7555.9, 10528.2, 13584.0, 16556.2 (spacings alternating 2972.25 and 3055.79, mean 3014.02). Near them
   `|D(1/2 + it)|` drops to `2.6·10^(−5)` (at 13584.0): the multiplier `1/(3(2/3)^(1/2)|D|)` that converts the
   square-root-scale content of `H` into the remainder of `S` reaches 15,645 there. This is a concrete reason why a bound
   on `H` does not by itself bound the remainder of the explicit formula (MISTAKE-620's "converse" caveat).
3. **Almost periodicity.** Shifting the zeros by `6028.04i = 2πi·665/ln 2` lands within 0.00061 (median) of another
   zero: `2^(−6028.04i) = 1` and `3^(6028.04i) = e^(2πiε)` with `ε = 6.3·10^(−5)`, so `D(s + 6028.04i) ≈ D(s)`.

## 2. The dense census and its reliability

`collatz_census_dense_resumable_20261016.py` reuses the exact kernels of the S37 streaming census (16 GB parity table
below 2^37, exact 8-step jumps, two-limb fallback) on a grid of `PER` log-uniform points per octave and can resume from
an exact state file. Validation: a split run (2^28, then resumed to 2^30) reproduces a direct 2^30 run in every output
(grid values, prefix sums, block moments); every shared grid point agrees with the S37 64-point census.

**A corrupted run.** The first dense run to 2^40 ended in an internal inconsistency: one final class count was off by
exactly `2^56` and the chained `H` by 8, two bit-level errors, with no software difference to the second run. The script
now checks, for every block, that the class counts, the register sum of `t` and the sub-block moments agree, and that
the chained `H` equals `3C(2 mod 3) − S` from the final counts; the rerun passed both with 0 inconsistent blocks,
gave `S(2^40) = 989143994`, `H(2^40) = −55691` (as in S37) and agrees exactly with the S37 census at all 1,920 shared grid
points. This replicates the S37 2^40 census on `[2^37, 2^40)` (same code, independent run), which had no second run
before. The referee's fault injection (report E6) shows what the invariants cover: single-bit faults in the class
counts, the register total, `Σt` per sub-block and the `H` increment are detected; faults in the weighted moments
`Σ(n − lo)t` (which feed the residues and certificates) and in the energy sums are not, and neither the parity table
nor the grid pass is covered. For `[0, 2^40)` this is closed by the replication and the referee's independent census
below 2^34; for later single runs it is addressed in section 5.

## 3. Residues two ways

The boundary estimate `r_hat(ρ) = −(3/2)^ρ numer(ρ)/D'(ρ)` (S37 certificate formula, complex128, block moments through
2^40, without the tail) and the dense least-squares fit (69 zeros, `[2^22, 2^40)`, 1024 points per octave) agree for the
dominant low modes: `s_0` 0.18398 / 0.18427, `0.8635 + 131.16i` 0.5065 / 0.5081, `0.9849 + 720.74i` 0.0701 / 0.0648 (aliased with 8561.6),
`0.9497 + 1572.56i` 0.2036 / 0.2025, `0.6939 + 22.25i` 2.570 / 2.579, `0.8925 + 611.62i` 0.1361 / 0.1357. They disagree
for the Mercator members above height 2500 (`ρ*`: 0.00999 / 0.00712; `3494.5`: 0.0128 / 0.0326; `3975.0`: 0.0135 /
0.0592; `4455.5`: 0.0124 / 0.151). THM-4659 statement 5 explains why: 1024 is a resonant density (`1024 log2 3 =
1623.0016`), so the high members alias approximately onto the low ones (`6268.3 → 3014.0`, `5787.8 → 3494.5`,
`5307.3 → 3975.0`, `4826.8 → 4455.4`, common detuning 0.00386 cycles per octave), and the fitted coefficients of those
zeros mix several modes. A grid with a non-resonant density resolves them: the referee's synthetic fits recover the
residues at 1000 points per octave as well as at `1000 sqrt 2`. (The drop of the holdout error from 0.0118 to 0.00167
on adding the near-Nyquist zeros is absorption of aliased power, not resolution of the high modes: report I5.)

**The test (FINITE-EXACT census, NUMERICAL fit).** A second census on the grid `round(2^(10 + j/(1000 sqrt 2)))`
(39,208 points to `2^38`; `S(2^38) = 281675896`, `H(2^38) = 130565`; it agrees with the 1024-per-octave census at
all 2,429 shared points) has Nyquist height 6409.7. The fit with the 94 zeros below height 6400 on `[2^24, 2^38)`
([collatz_irrational_grid_residue_fit_20261016.json](collatz_irrational_grid_residue_fit_20261016.json))
recovers the boundary residues of the whole Mercator family:

| zero (height) | fit, irrational grid | boundary (2^40 moments) | fit, 1024-per-octave grid |
|---:|---:|---:|---:|
| 240.248 | 0.18398 | 0.18398 | 0.18427 |
| 720.740 | 0.07009 | 0.07011 | 0.0648 |
| 2533.525 | 0.03701 | 0.03703 | |
| 3014.020 (`rho*`) | 0.01019 | 0.00999 | 0.00712 |
| 3494.515 | 0.01410 | 0.01284 | 0.0326 |
| 5307.300 | 0.05783 | 0.05736 | |
| 5787.793 | 0.02476 | 0.02514 | |
| 6268.288 | 0.00433 | 0.00532 | |

So the disagreement on the integer grid was aliasing, as THM-4659 statement 5 predicts. With the aliasing removed,
the least-squares residues agree with the boundary formula to 0.1% at the dominant members and to 2% at `rho*`. The
fit's holdout (fit to `2^36`, test `[2^36, 2^38)`) is 1.09% and its residual is `1.99 sqrt(N)`.

## 4. The vertical-line behaviour of the boundary numerator

The median of `|numer(σ + it)|` over the height bands `[10, 10³)`, `[10³, 4·10³)`, `[4·10³, 10⁴)`, `[10⁴, 2·10⁴)` is
267, 64, 49, 37 at `σ = 0.55`; 3.9, 0.93, 0.71, 0.55 at `σ = 0.7`; 0.055, 0.014, 0.010, 0.008 at `σ = 0.85`; 8e-4, 2e-4,
1e-4, 1e-4 at `σ = 1`. The finite numerator decays slowly with height. This is not evidence about the analytic function
it approximates: the certified tail bound grows like `|s| N^(−(σ−1/2))`, which at `σ = 0.85` for `N = 2^40` first
exceeds these values near `t ≈ 528` on the grid, and does so at 96% of the grid points with `t >= 200`.

## 5. The pre-registered test on [2^40, 2^42)

[The pre-registration](collatz_census_2p42_preregistration_20261016.md) was committed in `dc233fcf45` (11:43:32) before
the resumed run started (11:43:48; relaunched 14:20:09 after the session restart, with the same grid; errata section
added after the audit, criteria unchanged). The run gave `S(2^42) = 3573799150`, `H(2^42) = 1268051`, 0 inconsistent
blocks ([census](collatz_census_dense_1024_2p42_20261016.json)). It is a single run, so its published outputs were
checked against bounds that any significant corruption would violate: Lipschitz bounds on all 2048 grid values,
monotone prefix sums of `H²`, bounded increments of the prefix sums of `H`, moment bounds and the moment total against
`S` ([integrity](collatz_census_dense_2p42_integrity_20261016.out); all pass). Results
([evaluation](collatz_census_2p42_evaluation_20261016.json)):

| test | criterion | value | verdict |
|---|---|---:|---|
| Q1 | M1 (S37's 28-zero model) relative rms <= 0.03 | 0.0205 | PASSED |
| Q2 | M2 (69 zeros below Nyquist) relative rms <= 0.005 | 0.00115 | PASSED |
| Q3 | `e_X <= 2` at `X = 2^40, 2^41` | 0.139, 1.142 | PASSED |
| Q4 | `e_X <= 1` at both | 1.142 at `2^41` | FAILED |
| Q5 | mean window energy of `H` over starts in `[2^40, 2^41)` below 1 | 0.817 (max 1.25) | PASSED |
| Q6 | M2 residual / `sqrt(N)` rms <= 3.0 | 1.73 | PASSED |

M2 predicted two new octaves of `S(N)` to 0.115% relative rms, better than its own holdout (0.167%), with a residual
of `1.73 sqrt(N)`. The S37 model M1, fitted on a coarse grid with 28 zeros, held at 2.05%. The energy budget
`e_X <= 1` failed at `2^41` (as P4 failed at `2^38` in S37), while `e_X <= 2` held: the dyadic energy of `H` wanders
between 0.04 and 1.8 times `X²` with no visible drift (HYP-9283). The aliasing caveat applies to individual
coefficients, not to the prediction: 1024 is a resonant density, and the near-Nyquist zeros absorb aliased power in
and out of sample alike (report I5).

**Second stage.** [The second pre-registration](collatz_census_2p43_preregistration_20261016.md) (pushed 17:53:38)
precedes Run C (launched 17:54:00): the checked census of `[2^40, 2^43)` with checksummed accumulators. Its tests are
R1 (M2 refitted to `2^42` predicts `[2^42, 2^43)` within 0.0035), R2 (`e_X <= 2` at `2^42`), R3 (the residue estimate
at `ρ*` within 1% of 0.00999), R4 (the Arb certificate at `ρ*` with `N = 2^43`) and R5 (exact replication of Run B).
PENDING.

## 6. Toward a certified pole at ρ*

The boundary estimate of the residue at `ρ*` (renewal_shapes_oct11's formula (4), block moments) is stable:
`0.009924 + 0.001131i` (modulus 0.009988) at `N = 2^42` ([certificate](collatz_pole_certificate_rhostar_2p42_20261016.json)),
within `[0.009969, 0.009995]` over `N = 2^32, ..., 2^42`. The conditional Arb radius under the budget `(2, 0)` scales like
`N^(−1/2)`: 0.0261 at `2^40`, 0.0130 at `2^42` (interval `[−0.00304, 0.0230]`, not certified, as expected), 0.0092
expected at `2^43`. So `N = 2^43` is the first size at which the certificate can succeed, with a margin of about 8%.
That is test R4 of the second stage, and it needs exact moments through `2^43`, which is why Run C recomputes
`[2^40, 2^42)`. For comparison, the Mercator pole `s_0` stays certified at `2^42` (interval `[0.3529, 0.3571]` for
`2|r| 2^(29(σ−1))`, nested inside the `2^40` interval: a code check, not a test, by MISTAKE-621).

## 7. The conditional explicit formula (THM-4660)

The identity `F(s)D(s) = G(s) = −(3/2)^s(F_H(s)/3 + E(s))` (even/odd split of the Dirichlet series, the odd-step
relation `t(n) = −t((3n+1)/2)` and `[m ≡ 2 (3)] = (w_m + 1)/3`) puts all the arithmetic into `F_H(s) = s∫H(x)x^(−s−1)dx`
and an absolutely convergent correction `E`. Under the square-root hypothesis (W) `G` is analytic in `Re s > 1/2` with
`|G(s)| = O(|s|)`. Three unconditional facts about the exponential polynomial `D` then make the contour argument work:
its zeros are uniformly simple (`|D'(ρ)| >= c_0`; numerically `min |D'| = 0.1438` over the 295 zeros), they have bounded
density (at most one per unit height in the strip, numerically; Jensen gives fewer than 12 per rectangle of height 1),
and `|D|` is bounded below at fixed distance from the zeros (compactness of the torus of phases and Hurwitz). The
second Riesz mean then has an absolutely convergent modal expansion with error `O(x^(1/2+δ))`. The unsmoothed sum needs
much more (report E4): vertical-line control `∫ |G(σ' + it)| dt/|t| = O(T^ε)`, horizontal bounds of order `T^(1/2+δ)` at
`T ≍ x`, and convergence of `Σ r_ρ x^ρ/ρ`, which is not absolute if `|r_ρ| ≍ 1`. That is the remaining content of HYP-9285,
together with the square-root size of the remainder seen numerically.

## 8. Open threads

1. A certified nonzero residue at `ρ*` (the `1054/665` member) would give `limsup log B(X)/log X >= 0.99999921` under
   a budget beyond the census; PENDING at `N = 2^43`.
2. Non-vanishing of `G` at infinitely many Mercator members would give `limsup log B(X)/log X = 1` under (W).
3. The unsmoothed explicit formula (HYP-9285): vertical-line and horizontal bounds for `G` as in section 7.
4. Done: the irrational-density grid resolves the Mercator family (section 3).
