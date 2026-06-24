# LRC14 Fourier-Toeplitz Dual Scout

## Claim

HYP-2974 tries a non-packet proof route, complementary to HYP-2970's
endpoint-credit winding-cycle dual, HYP-2971's multiplicity-moment barriers,
HYP-2972's twist ladders, and HYP-2973's danger-count moment duals.

If the danger arcs cover the circle, then

```text
F_S(t) = sum_v 1_{||v t|| < 1/14} - 1
```

is nonnegative.  Therefore every Toeplitz moment matrix
`(hat F_S(i-j))` must be positive semidefinite.  A negative eigenvalue is a
Fourier/Farkas certificate of a safe open interval.

## Default Run

```text
04-computation/lrc14_fourier_toeplitz_dual_scout_codex_s156.py
05-knowledge/results/lrc14_fourier_toeplitz_dual_scout_codex_s156.out
```

Summary:

```text
rows audited                 52
positive-open rows           50
boundary-only rows           2
dual-certified rows          48
PSD-through-degree-90 rows    4
```

The four PSD-through-degree-90 rows are AP, GW, K33 near `12->36`, and
`P10+GW`.  So the harmonic lens sees most positive-open covering/migration rows
while naturally handing K33/petal exceptions back to the labelled stack.

## New Target

Prove bounded-degree Toeplitz negativity after removing AP/GW boundary atoms
and named K33/petal exits.  Unit-apex residue harmonics dominate the observed
negative eigenvectors.
