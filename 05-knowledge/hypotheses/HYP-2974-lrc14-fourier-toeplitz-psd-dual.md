---
id: HYP-2974
title: LRC14 Fourier-Toeplitz PSD dual
status: STUB / claimed phase-sensitive proof-interface route; computation pending
source: codex-2026-06-24-S157
script: 04-computation/lrc14_fourier_toeplitz_psd_dual_codex_s157.py
result: 05-knowledge/results/lrc14_fourier_toeplitz_psd_dual_codex_s157.out
related:
  - HYP-2973
  - HYP-2972
  - HYP-2971
  - HYP-2970
  - HYP-2969
  - HYP-2965
  - HYP-2963
  - HYP-2901
  - THM-406
  - OPEN-Q-108
---

# HYP-2974: LRC14 Fourier-Toeplitz PSD Dual

This slot is reserved for the Fourier version of the multiplicity-dual route.
For a 13-speed row `S`, let

```text
C_S(t) = #{v in S : ||v t|| < 1/14}.
```

A strict LRC14 counterexample would force the open danger arcs to cover the
circle, hence

```text
F_S(t) = C_S(t) - 1 >= 0
```

almost everywhere.  Therefore every finite Toeplitz moment matrix built from
the Fourier coefficients of `F_S` must be positive semidefinite:

```text
T_d(S) = [hat F_S(i-j)]_{0<=i,j<=d} >= 0.
```

If any low-degree `T_d(S)` has a negative eigenvalue, the associated
trigonometric square `|P(e^{2 pi i t})|^2` gives a dual certificate:

```text
integral F_S(t) |P(e^{2 pi i t})|^2 dt < 0,
```

so `F_S` is not nonnegative and `S` has a strict safe interval.

The planned computation will test whether the hard HYP-2963/LRC14 rows fail
this PSD necessary condition at low harmonic degree, compare the failures with
HYP-2973's danger-count duals, HYP-2971's scalar moment barriers, and
HYP-2972's twist witnesses, and record which rows remain PSD-blind.

This route is phase-sensitive in a way HYP-2973 deliberately is not.  HYP-2973
keeps only the distribution of the count random variable `C_S`; HYP-2974 keeps
the Fourier locations of that count function through the Toeplitz matrices.

No computational evidence is claimed by this stub yet.
