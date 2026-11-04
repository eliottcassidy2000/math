# Residual probe: Paley-7 with C3 children (EXPLORATORY, FINITE-EXACT single instance)

This probe was NOT pre-registered. It tests the HYP-9456 residue.

## Weighted model

The modulo-4 weighted model (`../wmodel.py`) was applied to a one-child test.
- **Instance.** Paley-7 quotient. One C3 child is named; every other child is contracted.
- **Finding.** With the C3 kept in its orientation and the quotient conversed, the configuration is **blind to every
  query named in one or two C3 children**. The reversed-C3 target is visible.
- **Scope across quotients.**
  - Paley-7 has no "cyclic ordered pair-module" witness at any position. So the S188 cyclic-child quartic does not
    apply.
  - 199 of 280 positions in random zero-parity prime 7-quotients also lack one.

## Exact integer test

`paley_c3.py` computes this with `runword.cpp`.
- **Method.** A run-word substitution DP treats the six other children as C3 modules with exact ordered-cover profile
  (q1,q2,q3) = (3,6,6). It was validated 12/12 against direct named Hamilton counts (`rw.py`).
- **Setup.**
  - T = Paley7[C3^7].
  - The two targets are U = Paley7^op[C3^7] (C3 unchanged) and Paley7^op[C3^op at the named child, C3 elsewhere].
- **Data.** Exact coefficients over all (1+9+36)^2 = 2116 downset queries. Variables: 6 coarse block variables and
  the 3 named vertices.

## Results

| target | nonzero exact coefficient differences | 2-adic valuations present |
|---|---:|---|
| unchanged C3 | 870 | {2, 3, 4, 6} (minimum 2) |
| reversed C3 | 636 | {1, 2, 3, 4, 6} |

- H(T) = H(U) = 123,156,134,945,235.
- **The unchanged case is not a collision.** It is separated at the integer level, already on a coarse coefficient:
  [alpha_1 gamma_2] differs by -18,257,826,464,676 = -4 x odd, so the residue is 4 mod 8. Modulo 4 it is blind,
  which agrees with the weighted model.

## Reading

The residual frontier HYP-9456 needs the modulo-8 layer, for which the S188 coarse_mod8_affine five-type law is the
template. No further modulo-4 idea is needed.

## Reproduction

    g++ -O2 -std=c++17 runword.cpp -o runword.exe
    python rw.py
    python paley_c3.py

`rw.py` validates against `fcoef_named.cpp`, which must be compiled to `fcoef.exe` beside `lib.py`.
