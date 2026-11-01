# Verification receipt

Source snapshot: `2a18212bfe0e8f4e82444ca2b96d3d708cb520e0`.
This was also the HEAD of `/tmp/session135-source` when checked. No source
checkout or shared repository file was edited by this lane.

## Reproduced after the final degree-nine bound change

| Check | Coverage | Result |
|---|---|---|
| `verify.py` | 59,049 full exponent vectors 0..2; 9,072 onehot states data 0..5; integers 1..30,000; 300 seeded vectors 0..8; source pins and rule deltas | 583,291 explicit gates; normal/-O byte-identical |
| `verify_entry_bound.py` | 2,760 inputs, 2,744 with positive control; literal positive phases and every claimed reset letter | 330,773 gates; normal/-O byte-identical |
| `density_bounds.py` | exact onehot/T=0 masses, inherited six-cycle point density, complementary upper bounds | normal/-O byte-identical |
| `density_approximation.py --cap 1` | all 1,024 binary exponent vectors, common-denominator weights and omitted-mass identity | normal/-O byte-identical |

The largest producer decision prefix in these controls is 1,566 actual steps.
The largest explicit phase/reset entry trace in the separate entry suite is
1,040 steps. Neither observed maximum is asserted as a universal bound.
All proof bounds and quantifiers are in `THEORY.md`.

The parent session's independent auditor separately confirmed the phase-budget
proof, onehot predicate, degree-nine bound, raw-entry bound, density existence,
degree-nineteen precision algorithm, exact lower bounds, and independent capped
state/cylinder tests. Its separately owned audit packet should accompany this
lane in the session's integration; it is not fabricated or duplicated here.

## Numerical scope

The finite integer range 1..30,000 contains 17,486 standard-orbit basin inputs,
12,478 two-cycle basin inputs, and 36 six-cycle basin inputs. These are exact
finite counts, **not** asymptotic estimates. The exact lower bounds in
`DENSITY.normal.json` instead follow from infinite-domain formulas.

`SHA256SUMS` pins every retained file except itself. Verification commands are
listed in `THEORY.md`; all scripts use only the Python standard library.
