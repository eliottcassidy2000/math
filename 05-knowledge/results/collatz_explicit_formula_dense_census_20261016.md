# The explicit formula of the Collatz parity sum on a dense census: the Mercator family of zeros, aliasing on integer grids, a second pre-registered test, and the near-critical zero 1054/665

2026-10-16, session opus-2026-10-16-S38 (worktree `/tmp/math-wt-explicit-modal-38`). Continues
[the S37 note](collatz_octave_memory_bias_energy_budgets_20261015.md) (HYP-9283, HYP-9284, HYP-9285). Audit: PENDING.

**Owner's prompt (verbatim).** "keep going with the next targets for another long session"

**Targets taken from S37's open obligations.** (i) the modal expansion of `S(N)` (HYP-9285): denser data, more zeros,
a second pre-registered out-of-sample test; (ii) the residues at many zeros, two ways; (iii) the growth of the boundary
numerator on vertical lines; (iv) whether the conditional pole certificate can reach zeros closer to `Re s = 1` than the
Mercator zero.

**Status.**
- **PROVED:** THM-4659 (the zeros of `D(s) = 1 + 2^(−s) + (3/2)^s/3` near `Re s = 1` are exactly the parity-matched
  near-returns `2^p ≈ 3^q`, `p` even, `q` odd; `1 − Re ρ = Cε² + O(ε⁴)`, `C = π² ln²2/ln³(4/3) = 199.16`; aliasing of the
  family on integer-density log grids). Unconditional lemmas of THM-4660 on `D` (uniformly simple zeros, bounded zero
  density, a uniform lower bound off the zeros).
- **CONDITIONAL (complete proof):** THM-4660. If `H` has dyadic energy `O(X^(2+ε))` for every `ε > 0` (a square-root
  hypothesis weaker than every budget of HYP-9283), then `F(s) = Σ t(n) n^(−s)` is meromorphic in `Re s > 1/2` with
  simple poles at zeros of `D`, and the second Riesz mean of `S` has the explicit formula
  `Σ_(n<=x) t(n)(1 − n/x)² = Σ_(Re ρ > 1/2+δ) 2r_ρ x^ρ/(ρ(ρ+1)(ρ+2)) + O(x^(1/2+δ))`.
- **FINITE-EXACT:** a dense census (1024 log-uniform points per octave) through 2^40, reproducing the S37 census at all
  1,920 shared grid points, and its resumed extension through 2^42 (PENDING); the zeros of `D` with `1/2 < Re < 1`,
  `0 < Im < 20000` (295, argument principle).
- **NUMERICAL:** residues at all 295 zeros from the exact 2^40 block moments; agreement with the dense least-squares
  fit for the dominant low modes; the holdout calibration of the explicit formula on the dense grid (0.17% with the
  69 zeros below the Nyquist height); the pre-registered test on `[2^40, 2^42)` (PENDING).
- **CONDITIONAL (Arb):** the residue at `ρ* = 0.99999921 + 3014.0202i` (the `1054/665` member) is `0.00993 + 0.00113i`
  before the future tail, with radius 0.026 under the budget `E_X <= 2X²` for `X >= 2^40`: not yet certified. PENDING at
  2^42 and 2^43.

## 1. The zeros of D in the critical strip up to height 20000

`collatz_D_zeros_and_boundary_residues_20261016.py` finds 295 zeros with `1/2 < Re ρ < 1` and `0 < Im ρ < 20000`
(vectorised Newton plus a dense local search where needed); the argument principle on the forty rectangles
`[1/2, 1] × [500k, 500(k+1)]` confirms the count in every rectangle. About 14.75 zeros per unit 1000 of height.

**Three structures.**
1. **The Mercator family near `Re s = 1` (THM-4659).** 23 zeros with `Re > 0.98` below height 20000, all at
   frequencies within 0.012 of half-integers (26.5036, 79.5105, ..., 332.5001, ... cycles per octave). They are the
   parity-matched near-returns: `84/53` gives `s_0 = 0.99820922 + 240.248i`, `1054/665` gives
   `ρ* = 0.99999921 + 3014.0202i`, its odd multiples give `0.99999289 + 9042.06i` and `0.99998025 + 15070.10i`, and the
   sidebands at `±480.5` and `±240.25` (e.g. `2533.5`, `3494.5`, `5787.8`, `6268.3`) are the near-returns
   `(p ± 2·84, q ± 2·53)`-type combinations. `(1 − Re ρ)/ε²` is 197.3 at `84/53` and 199.17 at `1054/665` against
   the limit `C = 199.16`. The next matched convergent `50508/31867` predicts a zero `0.999999978121 + 144432.720754i`.
   The zero `0.693918 + 22.2532i` behind the 2.45-cycle line of `S(N)/N` is the analytic continuation of the family to
   the convergent `8/5`.
2. **A near-1/2 family.** Zeros at `Re = 0.5006, 0.5012, 0.5004, 0.5016, 0.50006, 0.5019` with heights 1527.9,
   4500.1, 7555.9, 10528.2, 13584.0, 16556.2 (spacings alternating 2972.25 and 3055.79, mean 3014.02). Near them
   `|D(1/2 + it)|` drops to `2.6·10^(−5)` (at 13584.0): the multiplier `1/(3(2/3)^(1/2)|D|)` that converts the
   square-root-scale content of `H` into the remainder of `S` reaches 15,645 there. This is a concrete reason why a bound
   on `H` does not by itself bound the remainder of the explicit formula (MISTAKE-620's "converse" caveat).
3. **Almost periodicity.** Shifting the zeros by `6028.04i = 2πi·665/ln 2` lands within 0.0012 (median) of another
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
points. This also verifies the S37 2^40 census on `[2^37, 2^40)`, which had no independent check before.

## 3. Residues two ways

The boundary estimate `r_hat(ρ) = −(3/2)^ρ numer(ρ)/D'(ρ)` (S37 certificate formula, complex128, block moments through
2^40, without the tail) and the dense least-squares fit (69 zeros, `[2^22, 2^40)`, 1024 points per octave) agree for the
dominant low modes: `s_0` 0.18398 / 0.18427, `0.8635 + 131.16i` 0.5065 / 0.5081, `0.9849 + 720.74i` 0.0701 / 0.0648,
`0.9497 + 1572.56i` 0.2036 / 0.2025, `0.6939 + 22.25i` 2.570 / 2.579, `0.8925 + 611.62i` 0.1361 / 0.1357. They disagree
for the Mercator members above height 2500 (`ρ*`: 0.00999 / 0.00712; `3494.5`: 0.0128 / 0.0326; `3975.0`: 0.0135 /
0.0592; `4455.5`: 0.0124 / 0.151). THM-4659 statement 5 explains why: with 1024 points per octave the high members alias
exactly onto the low ones (`6268.3 → 3014.0`, `5787.8 → 3494.5`, `5307.3 → 3975.0`, `4826.8 → 4455.4`), so the fitted
coefficients of those zeros mix several modes. Fits that resolve the family need a grid with a non-integer number of
points per octave.

## 4. The vertical-line behaviour of the boundary numerator

The median of `|numer(σ + it)|` over the height bands `[10, 10³)`, `[10³, 4·10³)`, `[4·10³, 10⁴)`, `[10⁴, 2·10⁴)` is
267, 64, 49, 37 at `σ = 0.55`; 3.9, 0.93, 0.71, 0.55 at `σ = 0.7`; 0.055, 0.014, 0.010, 0.008 at `σ = 0.85`; 8e-4, 2e-4,
1e-4, 1e-4 at `σ = 1`. The finite numerator decays slowly with height. This is not evidence about the analytic function
it approximates: the certified tail bound grows like `|s| N^(−(σ−1/2))`, which exceeds these values beyond `t ≈ 200` at
`σ = 0.85` for `N = 2^40`.

## 5. The pre-registered test on [2^40, 2^42)

PENDING (see [the pre-registration](collatz_census_2p42_preregistration_20261016.md), committed in `dc233fcf45`
before the resumed run started).

## 6. Toward a certified pole at ρ*

PENDING.

## 7. The conditional explicit formula (THM-4660)

The identity `F(s)D(s) = G(s) = −(3/2)^s(F_H(s)/3 + E(s))` (even/odd split of the Dirichlet series, the odd-step
relation `t(n) = −t((3n+1)/2)` and `[m ≡ 2 (3)] = (w_m + 1)/3`) puts all the arithmetic into `F_H(s) = s∫H(x)x^(−s−1)dx`
and an absolutely convergent correction `E`. Under the square-root hypothesis (W) `G` is analytic in `Re s > 1/2` with
`|G(s)| = O(|s|)`. Three unconditional facts about the exponential polynomial `D` then make the contour argument work:
its zeros are uniformly simple (`|D'(ρ)| >= c_0`; numerically `min |D'| = 0.1438` over the 295 zeros), they have bounded
density (at most one per unit height in the strip, numerically; Jensen gives fewer than 12 per rectangle of height 1),
and `|D|` is bounded below at fixed distance from the zeros (compactness of the torus of phases and Hurwitz). The
second Riesz mean then has an absolutely convergent modal expansion with error `O(x^(1/2+δ))`. The unsmoothed sum needs
`|G(s)| = o(|s|)` on average along horizontal lines; that is the remaining content of HYP-9285, together with the
square-root size of the remainder seen numerically.

## 8. Open threads

1. A certified nonzero residue at `ρ*` (the `1054/665` member) would give `limsup log B(X)/log X >= 0.99999921` under
   a budget beyond the census; PENDING at `N = 2^43`.
2. Non-vanishing of `G` at infinitely many Mercator members would give `limsup log B(X)/log X = 1` under (W).
3. The unsmoothed explicit formula (HYP-9285): a bound for `G(σ + it)/(σ + it)` on average along horizontal lines.
4. Residues of the high Mercator members from a grid with an irrational number of points per octave (in progress).
