---
id: THM-550
title: Half-tiling dimension half(n)=floor((n-1)^2/4) (odd k^2 / even k(k-1)); phi-self-converse tilings = 2^half(n); grid-sym fraction = 2^((floor((n-1)/2)-C(n-1,2))/2); folding identity full(n)=2*half(n)-d
status: PROVED (orbit counting on THM-549; verified exhaustively n=3..7 and to n=40 numerically)
source: kind-pasteur-2026-06-20
depends_on:
  - THM-549   # rho is a pure reflection with d=floor((n-1)/2) fixed cells
related:
  - THM-551   # recursions
  - the grid-sym fraction line in CLAUDE.md (corrected kind-pasteur-S20ex)
external: A002620 (quarter-squares); Burnside / orbit counting.
---

# THM-550 — The half-tiling dimension and the grid-symmetric fraction

## Statement

Let `m = C(n-1,2)` (full-tiling dimension) and `d = floor((n-1)/2)` (THM-549 fixed
cells).  `rho` splits the tiles into `d` fixed cells and `(m-d)/2` transposed pairs.
The number of `rho`-orbits — the dimension of a fundamental domain, the **half-tiling** —
is
```
    half(n) = d + (m-d)/2 = (m+d)/2 = floor((n-1)^2/4),
```
the quarter-square sequence `A002620(n-1) = 0,1,2,4,6,9,12,16,20,25,30,...`
(`n = 2,3,4,...`), with the parity split
```
    odd  n = 2k+1 :  half = k^2          (a "gnomon square")
    even n = 2k   :  half = k(k-1)       (a pronic / oblong).
```

The `phi`-self-converse (grid-symmetric) tilings number exactly
```
    #{t : t = t . rho} = 2^{half(n)},
```
and the grid-symmetric **fraction** is
```
    2^{half(n) - m} = 2^{(d - m)/2} = 2^{(floor((n-1)/2) - C(n-1,2))/2},
```
giving exponents `0, -1, -2, -4, -6, -9, -12, ...` for `n = 3,4,5,6,7,8,9` —
**exactly the canon grid-sym fraction** (CLAUDE.md, corrected by kind-pasteur-S20ex from
the earlier wrong `2^{-(n-2)}`).

**Folding identity.**  The full and half dimensions are tied by
```
    full(n) = C(n-1,2) = 2*half(n) - d,      d = floor((n-1)/2):
```
the full staircase is two copies of the half-tiling glued along the `d` diagonal cells
(which are shared, hence counted once).

## Proof

Orbit count: `#orbits = (#fixed) + (#pairs) = d + (m-d)/2 = (m+d)/2`.  Evaluate:
- `n = 2k+1`: `m = C(2k,2) = k(2k-1)`, `d = k`, so `(m+d)/2 = (2k^2)/2 = k^2 = floor((2k)^2/4)`.
- `n = 2k`:   `m = C(2k-1,2) = (2k-1)(k-1)`, `d = k-1`, so `(m+d)/2 = (k-1)*2k/2 = k(k-1) = floor((2k-1)^2/4)`.

Grid-symmetric tilings are exactly those constant on each `rho`-orbit (THM-549
corollary): one free bit per orbit gives `2^{#orbits} = 2^{half(n)}`.  The fraction is
`2^{half(n)}/2^m = 2^{half-m} = 2^{(d-m)/2}`.  Folding: `2*half - d = (m+d) - d = m`. **QED.**

## Connections

- **Burnside.**  `2^{half(n)}` is the Burnside fixed-count `Fix(phi)` for the reversal
  anti-automorphism `phi` acting on the `2^m` tilings.  It is the tiling-model analogue
  of the canon Burnside law "all-odd cycle types contribute `2^{orbit-pairs}`".  It
  enters the self-complementary count `SC_n` as the `phi`-term of the anti-automorphism
  Burnside sum (see HYP-2686).
- **Two-sheeted structure.**  `2^{half}` of `2^m` is a fraction `2^{-(m-half)} =
  2^{-(m-d)/2}` — the codimension is the number of off-diagonal `rho`-pairs.
- **A halving made explicit.**  This answers, at the labeled-tiling level, the
  "intrinsic halving structure of SC tournaments" flagged as open question #4 in
  07-reflections/two-models-staircase-recursion.md.

## Verification

`04-computation/half_tiling_framework_kps.py` blocks [A],[C]: half(n) table for n=2..12
matches `floor((n-1)^2/4)` and the odd/even `k^2 / k(k-1)` split; exhaustive grid-sym
counts for n=3..7 are 2,4,16,64,512 = `2^{half}`, equal to the `phi`-self-converse
counts, with log2-fractions 0,-1,-2,-4,-6 matching canon.
