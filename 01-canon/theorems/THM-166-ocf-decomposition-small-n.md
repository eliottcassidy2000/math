# THM-166: OCF Decomposition for Small n

**Status:** VERIFIED (n=3-6 exhaustive, n=7 regular)
**Session:** kind-pasteur-2026-03-13-S61

## Statement

### The Truncated OCF

For n-vertex tournaments with n <= 8:

    H(T) = 1 + 2*alpha_1(T) + 4*alpha_2(T)

where:
- alpha_1 = total number of directed odd cycles (all lengths 3, 5, 7, ...)
- alpha_2 = number of vertex-disjoint pairs of directed odd cycles
- Higher alpha_k = 0 because 3+3+3 = 9 > 8 (no three disjoint 3-cycles fit)

### The Additive Decomposition

    H(T) = 1 + 2*c3_dir + 2*c5_dir + 2*c7_dir + ... + 4*alpha_2

This splits into:
- **Score-determined**: c3_dir is determined by score sequence (Rao's formula)
- **5-cycle contribution**: 2*c5_dir (at n=5, this is the Hamiltonian cycle count)
- **7-cycle contribution**: 2*c7_dir (at n=7, this is the Hamiltonian cycle count)
- **Disjointness contribution**: 4*alpha_2

### The Coefficient-2 Principle

Any directed odd cycle that uses more than n/2 vertices conflicts with EVERY other odd cycle (since any 3-cycle needs 3 vertices, and with > n/2 already used, the remaining < n/2 vertices can't support a 3-cycle on their own).

Consequence: such a cycle contributes exactly 2^1 = 2 to I(Omega, 2).

At n=5: ALL cycles use >= 3 = ceil(5/2) vertices, and alpha_2 = 0 always.
At n=6: 5-cycles use 5 > 3 = 6/2 vertices, so each contributes 2.
At n=7: 7-cycles use 7 > 3.5 vertices (and 5-cycles use 5 > 3.5).

More precisely: a k-cycle uses k vertices. It's disjoint from no cycle of length 3 when k > n-3 (i.e., fewer than 3 vertices remain).

## Verified Cases

### n=5 (exhaustive, 1024 tournaments)

    H = 1 + 2*(c3_dir + c5_dir)

- alpha_2 = 0 always (3+3=6 > 5)
- c3_dir = score-determined
- c5_dir in {0, 1, 2, 3} varies within score class (1,2,2,2,3)
- H = 1 + 2*alpha_1, all values odd

### n=6 (exhaustive, 32768 tournaments)

    H = 1 + 2*(c3_dir + c5_dir) + 4*alpha_2

- alpha_2 >= 0 (3+3=6=n, disjoint 3-3 pairs possible)
- Within each (c3_dir, alpha_2) subclass, H = const + 2*c5_dir:
  - (c3=4, a2=0): H = 9 + 2*c5_dir
  - (c3=5, a2=1): H = 15 + 2*c5_dir
  - (c3=6, a2=0): H = 13 + 2*c5_dir
  - (c3=6, a2=1): H = 17 + 2*c5_dir
  - (c3=6, a2=2): H = 21 + 2*c5_dir
  - (c3=7, a2=1): H = 19 + 2*c5_dir
  - (c3=8, a2=1): H = 21 + 2*c5_dir
- All residuals = 1 + 2*c3_dir + 4*alpha_2

### n=7 regular (2640 tournaments)

    H = 1 + 2*(14 + c5_dir + c7_dir) + 4*disj_33

With the rigid constraint c5_dir + 2*disj_33 = 56:
- H = 1 + 28 + 2*c5_dir + 2*c7_dir + 4*disj_33
    = 29 + 2*(c5_dir + 2*disj_33) + 2*c7_dir
    = 29 + 112 + 2*c7_dir = 141 + 2*c7_dir

| H | c3_dir | c5_dir | c7_dir | disj_33 | nc | c5+2*disj |
|---|--------|--------|--------|---------|-----|-----------|
| 171 | 14 | 36 | 15 | 10 | 65 | 56 |
| 175 | 14 | 28 | 17 | 14 | 59 | 56 |
| 189 | 14 | 42 | 24 | 7 | 80 | 56 |

## Connections

- **THM-164** (Fourier-Cycle Bridge): H_2 captures c3 (score), H_4 captures c5 (non-score). This theorem gives the OCF explanation: H_2 = score-dependent part of 1 + 2*c3_dir, and H_4 + H_6 = 2*c5_dir + 4*alpha_2 + 2*c7_dir.

- **THM-027** (BIBD Maximization): Paley/BIBD minimizes alpha_2 but maximizes c5_dir and c7_dir. The net effect is maximum H because the coefficient on cycles (2) outweighs the loss from fewer disjoint pairs (4 per pair, but fewer pairs).

## Verification Scripts

- ocf_directed_fix.py: OCF with directed cycles, verified n=3-6
- cycle_disjointness_constraint.py: Full decomposition at n=5,6
- n7_regular_cycle_spectrum.py: Regular n=7 analysis
