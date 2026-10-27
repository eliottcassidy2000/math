# The full ROOT current at every octave to 2^27, its odd-source fold, the calibration of a uniform bound, and contact rotations that stay in the base field

2026-10-17, opus-2026-10-17-S40. Owner prompt: uniform Collatz tail debt / uniform ROOT-current bounds, and
mixed-coordinate forcing beyond the spindle cube; companion to clock_duality_oct15's
[marked contacts and ROOT-current folding](marked_contacts_root_current_20261017e.md) and
[escape renewal](escape_renewal_idoneal_20261017d.md).

**Typing.**
- FINITE-EXACT: the full current on every octave `2^8..2^27`; the fold `R = I + B` at every row to `2^26`; the 34
  contact rotations and their port test.
- Independently audited 2026-10-17 (s40_audit_20261017/REPORT.md): every number for `k = 8..24` and the fold at every
  row reproduced; two prose errors corrected (MISTAKE-626).
- PROVED (by citation in the repository): the calibration of section 3 (HYP-9283(C), S37).
- OPEN: any uniform bound; six-colour forcing.
- No novelty claim.

## 1. The full current (FINITE-EXACT)

clock_duality_oct15 brackets the current `a_N = H(4N) - H(N)` by certified ROOT parts `R` and uncertainties `U`
without assuming convergence (D1). Every `n <= 2^29` reaches 1, so the full current is known: `U = 0` and `a_N = R`.
[collatz_root_current_exact_20261017.py](../../04-computation/experiments/collatz_root_current_exact_20261017.py)
computes it at every row of every octave `X = 2^k`, `k = 8..27`. The census `H` agrees with the S37/S38 dense census at
all 20,393 grid points below `2^30`. Energies and tail constants are exact integer numerators over `X^2` (int64
accumulation, with an overflow bound asserted before summing; the numerators are in the JSON). The decimals below are
roundings. The first version used float64 reductions, which clock_duality rightly declined to treat as integer
receipts.

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

- The bulk and the rim are anti-correlated at 17 of 19 octaves, positive only at 2^18 and 2^22 (corrected after the
  audit; the first version and commit 309af1b454 said 16 of 19). The current's energy is what survives a large
  cancellation: at `2^23` the bulk energy 8.54 cancels to 1.51.
- The highest-energy octave, `2^22`, is where the cancellation fails (correlation +0.61). At `2^12` bulk and rim still
  anti-correlate (-0.64); there the energy comes from a large rim (`||B||^2 = 6.57 X^2`).
- clock_duality's staged cross term at `k = 18`, stage 28 (`<I,B>/X^2 = 0.173301`) is close to the final 0.1735. The
  staged energy 1.634206 differs from the final 1.6042, because the stages have different terminal coverage (as
  clock_duality notes in root_braid_logistic_20261017f.md).
- Their exact rejection of a universal nonpositive I/B covariance rule holds for the full current too (`k = 18, 22`).
- **Consequence for the uniform-bound target.** A uniform bound on `R` cannot come from bounding `I` and `B` separately
  with a fixed-sign cross term. It must control their covariance, which is exactly the odd-transition guard
  clock_duality names: the hostile `f = χ` obeys doubling and four-Lipschitz control but has `R = -3N`.

## 3. Calibration of a uniform bound (PROVED implication, S37)

The tail-debt target (W) is `Σ_(X <= n < 2X) H(n)^2 = O(X^(2+ε))`. This is exactly the energy `E_X` of HYP-9283; S38's
THM-4660 is conditional on the same (W). The ground truth is `E_X/X^2 <= 1.769` for dyadic `4 <= X <= 2^39` (S37; `X = 2`
gives 3.25), S38's first stage adds 0.139 at `2^40` and 1.142 at `2^41`, and Run C of S38 extends the census to `2^43`
(pre-registered, pending).

HYP-9283(B, C): any certified constant budget `E_X <= 14348 X^2` for dyadic `X >= 2^40` (or a linear-log budget
`(2 + 9616 j) X^2`) certifies the nonzero Mercator residue at `s_0 ≈ 0.99821 + 240.248i`. Hence `F` has a pole there, and

    limsup log B(X) / log X >= 0.99821,     B(X) = #{n <= X : n reaches 1},

against Krasikov-Lagarias's unconditional exponent 0.84. So **a uniform ROOT-current bound with an explicit constant of
that size would already be a major theorem about the density of the tree of 1.** A bound with only `X^ε` and an
unspecified constant does not run the certificate as stated, which needs explicit tail constants. This places
clock_duality's anchor "an independent uniform bound on the certified ROOT current": it is not a technical gap. The
finite data above say the constant is small (`K_X <= 2.4`, energy `<= 4.6 X^2` to `2^27`), not that a proof is near.

**From the current to `E_X`.** A budget for the current reaches `E_X` through THM-4658:
`sqrt(E_4X) <= 2 sqrt(E_X) + 2 sqrt(D_X) + 4 sqrt(3X)`, with `D_X = Σ_(X <= N < 2X) a_N^2`. So `D_X <= C X^2` for
`X >= 2^40` gives `E_X <= (max(e_(2^40), e_(2^41), C) + o(1)) X^2`. A constant tail constant `K_X` gives only
`D_X = O(X^2 log X)` (graded_tail_carrier (T3)), i.e. a linear-log budget, which is one of the certified shapes.

**A floor for every certificate.** clock_duality's `K_X` (graded_tail_carrier (T1)) is computed on certificate
envelopes `C(N) >= |a_N|`, so the exact values above are floors for any certificate. At `2^22` the floor is 2.398,
against the certified stage-27 value 3.928 (escape_renewal_idoneal_20261017d.md section 3).

## 4. Contact rotations stay in the base field (FINITE-EXACT, consistent with clock_duality's degree lemma)

[mixed_spindle_rotations_20261017.py](../../04-computation/experiments/mixed_spindle_rotations_20261017.py) searches
the `sqrt3` escape `F_s = K(w)`, `K = Q(zeta15, sqrt(-11))`, for rotations `U = (e - d1)/d2` with `|e - d1| = |d2|`.
Here `e` ranges over 570 unit directions (`±zeta15^k` times words of length `<= 2` in `ρ, u, v`) and `d1, d2` over the 34
spindle differences. Then `d1 + U d2 = e` is a unit edge of `M + U M` that changes both coordinates. Results:

- **All 34 rotations lie in `K`** (no `w` component). This is forced by clock_duality's lemma (5): `K(w)/K` is cubic,
  and a contact needs `[K(U):K] <= 2`. The search script does not test it; it was checked separately in S40, and the
  audit's audit_rotations_field.py confirms 34/34 with zero `w`-coordinates, `|U| = 1`, and six roots of unity.
  Together these confirm the lemma independently in the escape field.
- The non-root-of-unity blocks with mixed edges (4 or 2 of them; 47 or 49 points) are four-chromatic. The six
  root-of-unity rotations also give 3-4 mixed edges on 35-37 points.
- Attached at the eight ports with the `zeta5` rotations (1778-1858 vertices, about 8000 edges, every edge by exact
  norm), they block none of the 570 boundary patterns.

Since these blocks lie in the base field, they are covered by clock_duality's section-7 retraction logic; they add no
forcing. Their genuinely new geometric step, the quadratic coupler `v` of G1 and the ternary obstruction G2 for three
couplers, is where mixed-coordinate forcing could still appear. The full three-coupler graph and its port relation
remain the next consumer (their section 8).

**The eleventh-root escape.** clock_duality's other escape, `Q(zeta55, w)`, keeps its fixed 19-places, so its
residue colouring survives. It also cannot host this construction. The Moser spindle needs `sqrt(-3)`, which
`Q(zeta55)` lacks: it contains `sqrt(-11)` and `sqrt5`, but not `sqrt(-3)`. And the same degree lemma forces any contact
rotation in `Q(zeta55)(w)` back into `Q(zeta55)`. So the sqrt3 escape is the only escape where spindle-based mixed
contacts can be tried, and there they stay in the base field.

## Reproduction

```
python 04-computation/experiments/collatz_root_current_exact_20261017.py 8 27
python 04-computation/experiments/collatz_root_current_fold_full_20261017.py 8 26
python 04-computation/experiments/mixed_spindle_rotations_20261017.py --ports --maxU 8
```
