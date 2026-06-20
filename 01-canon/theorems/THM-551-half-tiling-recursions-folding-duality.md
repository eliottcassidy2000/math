---
id: THM-551
title: Half-tiling recursions (even half(n)=2A-C; odd half(n)=2A-2E+G) and the full<->half folding duality -- the user's A+B-C and A+B-C+D-E-F+G are the anti-diagonal fold of the full A+B+C-D-E-F+G (three overlapping (n-1)-triangles)
status: PROVED (exact finite-difference identities; verified to n=40)
source: kind-pasteur-2026-06-20
depends_on:
  - THM-549   # the fold rho
  - THM-550   # half(n), folding identity full=2half-d
related:
  - two-models-staircase-recursion.md   # full Mode-B: overlap+wiring+apex
  - everything-is-the-triangle.md        # Mode A / Mode B
external: finite-difference calculus; quarter-square A002620.
---

# THM-551 — Half-tiling recursions and the fold of the full-tiling recursion

## Statement

Write `full(n)=C(n-1,2)` and `half(n)=floor((n-1)^2/4)`.

**(a) Full-tiling recursion (the canon "A+B+C-D-E-F+G").**
```
    full(n) = 3*full(n-1) - 3*full(n-2) + full(n-3).
```
This is the inclusion-exclusion of THREE overlapping `(n-1)`-staircases
`A,B,C`, with pairwise overlaps `D,E,F` of size `full(n-2)` and triple overlap `G` of
size `full(n-3)`.  Equivalently it is `Delta^3 = 0` for the degree-2 polynomial
`C(n-1,2)`.

**(b) Half-tiling recursions (the user's framework).**
```
    even n :  half(n) = 2*half(n-1) - half(n-2)                       [ A + B - C ]
    odd  n :  half(n) = 2*half(n-1) - 2*half(n-3) + half(n-4)         [ A+B-C+D-E-F+G ]
```
with `A,B ~ size n-1`, `C,D ~ size n-2`, `E,F ~ size n-3`, `G ~ size n-4`; in the odd
sum the `-C + D` pair cancels (`C = D = half(n-2)`), leaving `2A - 2E + G`.
Region layout of the half-region (a triangle): corners `A, D, B`; edges `A+D-E`,
`B+D-F`, `A+B-C`; center `A+B+G-E-F`.  Both identities are exact for `floor((n-1)^2/4)`.

**(c) Folding duality.**  Under the anti-diagonal fold of THM-549/550
(`full = 2*half - d`, `d = floor((n-1)/2)`), the 3-fold-symmetric full decomposition
descends to the 2+1 asymmetric half decomposition: **the symmetry axis (anti-diagonal)
absorbs one of the three corners**, shrinking that corner from size `n-1` (corner `C`
in the full) to size `n-2` (corner `D` in the half).  The second difference of
`half(n)` is `0` for even `n` and `1` for odd `n`; the `+1` is precisely the new
diagonal (fixed) cell created when `n` grows by 2 (so `d` grows by 1).  Thus:
```
    odd half-recursion = even-style 3-term  +  the apex/center diagonal cell.
```

## Proof

**(a)** `C(n-1,2) = (n-1)(n-2)/2` is quadratic, so its 3rd finite difference vanishes:
`full(n) - 3 full(n-1) + 3 full(n-2) - full(n-3) = 0`.  The geometric reading: a
side-`s` triangle of unit cells is covered by three side-`s` copies overlapping in a
central inverted triangle; inclusion-exclusion over the three copies and their
pairwise/triple intersections gives `3*full(n-1) - 3*full(n-2) + full(n-3)`.

**(b)** Direct algebra on `q(n) := floor((n-1)^2/4)`.  For even `n`, `q(n)=(n^2-2n)/4`;
`2q(n-1)-q(n-2) = 2*(n-2)^2/4 - (n^2-6n+8)/4 = (n^2-2n)/4 = q(n)`.  For odd `n`,
`q(n)=(n-1)^2/4`; `2q(n-1)-2q(n-3)+q(n-4) = [2(n-2)^2 - 2(n-4)^2 + (n-5)^2]/4
= (n^2-2n+1)/4 = q(n)`.  The `-C+D` cancellation is by construction (`C=D=q(n-2)`).
Verified to `n=40` for both, and the literal 7-term/3-term signed sums to `n=29/28`,
in `04-computation/half_tiling_framework_kps.py` block [D].

**(c)** From THM-550, `full = 2*half - d` with `d` increasing by 1 each time `n`
increases by 2.  The second difference `q(n)-2q(n-1)+q(n-2)` equals `0,1,0,1,...`
(alternating, `1` at odd `n`); the unit jump at odd `n` is the extra fixed cell.  The
full Mode-B decomposition (two-models-staircase-recursion.md) is
`delta_{n-2} = delta_{n-4}(overlap) + bottom-wiring + top-wiring + apex`; folding the
two legs onto each other by `rho` identifies bottom-wiring with top-wiring (they are
`rho`-images: `(x,1) <-> (n, n+1-x)`) and fixes the apex `(n,1)` only when it is on the
anti-diagonal, i.e. when `n+1 = n+1` — always — so the apex `(n,1)` is a fixed cell iff
`n+1 = a+b = n+1` ✓; it is the outermost gnomon corner. **QED.**

## The gnomon picture (odd n)

For odd `n = 2k+1`, `half = k^2 = 1 + 3 + 5 + ... + (2k-1)`: the half-region is built
from `k` L-shaped gnomons.  Passing `n -> n+2` (odd to odd) adds the gnomon `2k+1`,
which is exactly `2*half(n+1) - half(n) ...` — the Mode-B (n -> n-2) slow time scale.
This is the half-tiling shadow of the Cayley-Dickson / Mode-B descent
(everything-is-the-triangle.md): odd steps add an odd-area shell.

## Why this matters

The user asked to connect the half framework to the full `A+B+C-D-E-F+G`.  The answer:
**they are the same inclusion-exclusion, before and after folding across the converse
symmetry axis.**  The full triangle is 3-fold (three equal corners); folding it across
its axis of symmetry merges two corners' worth of cells pairwise and leaves the diagonal
corner at half size, turning `(+A+B+C, -D-E-F, +G)` into `(+A+B-C, +D-E-F, +G)` with the
asymmetric corner sizes `(n-1, n-1, n-2)`.  The odd/even dichotomy (extra center cell at
odd `n`) is the geometric source of the second-difference `0/1` alternation and mirrors
the deep even-odd tournament split (THM-016/017).
