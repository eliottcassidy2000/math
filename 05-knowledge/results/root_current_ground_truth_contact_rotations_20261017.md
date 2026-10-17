# The full ROOT current at every octave to 2^27, its odd-source fold, the calibration of a uniform bound, and contact rotations that stay in the base field

2026-10-17, opus-2026-10-17-S40. Owner prompt: uniform Collatz tail debt / uniform ROOT-current bounds, and
mixed-coordinate forcing beyond the spindle cube; companion to clock_duality_oct15's
[marked contacts and ROOT-current folding](marked_contacts_root_current_20261017e.md) and
[escape renewal](escape_renewal_idoneal_20261017d.md).

**Typing.**
- FINITE-EXACT: the full current on every octave `2^8..2^27`; the fold `R = I + B` at every row to `2^26`; the 34
  contact rotations and their port test.
- PROVED (by citation in the repository): the calibration of section 3 (HYP-9283(C), S37).
- OPEN: any uniform bound; six-colour forcing.
- No novelty claim.

## 1. The full current (FINITE-EXACT)

clock_duality_oct15 brackets the current `a_N = H(4N) - H(N)` by certified ROOT parts `R` and uncertainties `U`
without assuming convergence (D1). Every `n <= 2^29` reaches 1, so the full current is known: `U = 0` and `a_N = R`.
[collatz_root_current_exact_20261017.py](../../04-computation/experiments/collatz_root_current_exact_20261017.py)
computes it at every row of every octave `X = 2^k`, `k = 8..27`. The census `H` agrees with the S37/S38 dense census at
all 20,393 grid points below `2^30`.

| k | energy `Σ a_N^2 / X^2` | tail constant `K_X` | fraction `a_N > 0` | sign changes | max `|a_N|` / `sqrt(X)` |
|---:|---:|---:|---:|---:|---:|
| 11 | 1.954 | 0.816 | 0.000 | 0 | 2.52 |
| 12 | 3.978 | 1.510 | 0.689 | 2 | 4.19 |
| 17 | 0.985 | 0.436 | 1.000 | 0 | 1.81 |
| 18 | 1.604 | 0.780 | 1.000 | 0 | 2.07 |
| 19 | 0.841 | 0.364 | 1.000 | 0 | 1.76 |
| 20 | 0.213 | 0.072 | 0.807 | 1186 | 1.25 |
| 22 | **4.537** | **2.398** | 1.000 | 0 | 3.16 |
| 23 | 1.511 | 0.721 | 1.000 | 0 | 2.01 |
| 26 | 2.835 | 1.464 | 0.018 | 167 | 3.87 |
| 27 | 2.820 | 1.304 | 0.350 | 217 | 2.88 |

(All twenty rows are in the [output](collatz_root_current_exact_20261017.out).)

- On the whole range: energy `<= 4.537 X^2`, `K_X <= 2.398`, and `max |a_N| <= 4.19 sqrt(X)`. Target (T) of the tail-debt
  programme holds here with a constant, no `X^ε` needed.
- Whole octaves have one sign: positive throughout 2^17-2^19, 2^22 and 2^23; negative throughout 2^11 (`min |a| = 16`);
  98% negative at 2^26. The current is a slowly modulated signal of size `sqrt(X)`, not row noise.
- Since `a_N = sqrt(2N) (2 g(log2 N + 2) - g(log2 N))` with S37's profile `g(u) = H(2^u)/sqrt(2^(u+1))`, coherent octaves
  are S37's octave memory seen two octaves apart.
- At `2^22` this agrees with clock_duality's stage-27 bracket: `Σ R^2/X^2 = 4.5634` against the exact 4.5367; every row
  has `L >= 98`, and the exact minimum is 382; `max L = 5696 <= max |a| = 6466`. Their constant-sign conclusion for
  `2^22` is confirmed (all positive).

## 2. The odd-source fold on the full current (FINITE-EXACT)

clock_duality's R1 writes `R = I + B` with `I(N) = -3 Σ_(m <= 2N, m odd) χ(m) f(m)` (odd sources with two threshold
events) and `B(N) = Σ_(2N < m <= 4N, m odd) w(m) f(m)` (the rim).
[collatz_root_current_fold_full_20261017.py](../../04-computation/experiments/collatz_root_current_fold_full_20261017.py)
checks `I + B = a_N` at every row of every octave `2^8..2^26` (zero failures) and reports:

| k | `||I||^2/X^2` | `||B||^2/X^2` | `<I,B>/X^2` | `||R||^2/X^2` | corr |
|---:|---:|---:|---:|---:|---:|
| 12 | 1.820 | 6.569 | -2.205 | 3.978 | -0.64 |
| 14 | 3.337 | 2.435 | -2.556 | 0.661 | -0.90 |
| 18 | 0.472 | 0.785 | +0.174 | 1.604 | +0.29 |
| 22 | 0.638 | 2.395 | **+0.752** | **4.537** | **+0.61** |
| 23 | **8.539** | 3.154 | -5.091 | 1.511 | **-0.98** |
| 26 | 2.417 | 0.438 | -0.010 | 2.835 | -0.01 |

(All nineteen rows are in the [output](collatz_root_current_fold_full_20261017.out).)

- The bulk and the rim are anti-correlated at 16 of 19 octaves. The current's energy is what survives a large
  cancellation: at `2^23` the bulk energy 8.54 cancels to 1.51.
- The two highest-energy octaves are where the cancellation fails: `2^22` (correlation +0.61) and `2^12`.
- clock_duality's staged value at `k = 18`, stage 28 (`<I,B>/X^2 = 0.173301`) is essentially the final state (0.1735).
- Their exact rejection of a universal nonpositive I/B covariance rule holds for the full current too (`k = 18, 22`).
- **Consequence for the uniform-bound target.** A uniform bound on `R` cannot come from bounding `I` and `B` separately
  with a fixed-sign cross term. It must control their covariance, which is exactly the odd-transition guard
  clock_duality names: the hostile `f = χ` obeys doubling and four-Lipschitz control but has `R = -3N`.

## 3. Calibration of a uniform bound (PROVED implication, S37)

The tail-debt target (W) is `Σ_(X <= n < 2X) H(n)^2 = O(X^(2+ε))`. This is exactly the energy `E_X` of HYP-9283; S38's
THM-4660 is conditional on the same (W). The ground truth is `E_X/X^2 <= 1.769` for dyadic `X <= 2^39` (S37), and Run C
of S38 extends it to `2^43` (pre-registered, pending).

HYP-9283(B, C): any certified constant budget `E_X <= 14348 X^2` for dyadic `X >= 2^40` (or a linear-log budget
`(2 + 9616 j) X^2`) certifies the nonzero Mercator residue at `s_0 ≈ 0.99821 + 240.248i`. Hence `F` has a pole there, and

    limsup log B(X) / log X >= 0.99821,     B(X) = #{n <= X : n reaches 1},

against Krasikov-Lagarias's unconditional exponent 0.84. So **a uniform ROOT-current bound with an explicit constant of
that size would already be a major theorem about the density of the tree of 1.** A bound with only `X^ε` and an
unspecified constant does not run the certificate as stated, which needs explicit tail constants. This places
clock_duality's anchor "an independent uniform bound on the certified ROOT current": it is not a technical gap. The
finite data above say the constant is small (`K_X <= 2.4`, energy `<= 4.6 X^2` to `2^27`), not that a proof is near.

## 4. Contact rotations stay in the base field (FINITE-EXACT, consistent with clock_duality's degree lemma)

[mixed_spindle_rotations_20261017.py](../../04-computation/experiments/mixed_spindle_rotations_20261017.py) searches
the `sqrt3` escape `F_s = K(w)`, `K = Q(zeta15, sqrt(-11))`, for rotations `U = (e - d1)/d2` with `|e - d1| = |d2|`.
Here `e` ranges over 570 unit directions (`±zeta15^k` times words of length `<= 2` in `ρ, u, v`) and `d1, d2` over the 34
spindle differences. Then `d1 + U d2 = e` is a unit edge of `M + U M` that changes both coordinates. Results:

- **All 34 rotations lie in `K`** (no `w` component). This is forced by clock_duality's lemma (5): `K(w)/K` is cubic,
  and a contact needs `[K(U):K] <= 2`. The search is an independent confirmation of the lemma in the escape field.
- The blocks with strictly mixed edges (4 or 2 of them; 47 or 49 points) are four-chromatic.
- Attached at the eight ports with the `zeta5` rotations (1778-1858 vertices, about 8000 edges, every edge by exact
  norm), they block none of the 570 boundary patterns.

Since these blocks lie in the base field, they are covered by clock_duality's section-7 retraction logic; they add no
forcing. Their genuinely new geometric step, the quadratic coupler `v` of G1 and the ternary obstruction G2 for three
couplers, is where mixed-coordinate forcing could still appear. The full three-coupler graph and its port relation
remain the next consumer (their section 8).

## Reproduction

```
python 04-computation/experiments/collatz_root_current_exact_20261017.py 8 27
python 04-computation/experiments/collatz_root_current_fold_full_20261017.py 8 26
python 04-computation/experiments/mixed_spindle_rotations_20261017.py --ports --maxU 8
```
