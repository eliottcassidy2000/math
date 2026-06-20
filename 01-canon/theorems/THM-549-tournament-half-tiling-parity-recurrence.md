---
id: THM-549
title: Tournament half-tilings are the mirror-folded staircase fundamental domain and obey parity-split recurrences
status: PROVED; verified by `half_tiling_recursion_codex_s56.py`
source: codex-2026-06-20-S56
depends_on:
  - THM-280
  - THM-442
related:
  - THM-513
  - HYP-2685
  - HYP-2660
---

# THM-549 - Half-Tiling Parity Recurrence

## Setup

Fix the tournament tiling model with base Hamiltonian path

```text
n -> n-1 -> ... -> 1.
```

The free cells are pairs `(x,y)` with `1 <= y < x <= n` and `x-y >= 2`.
Their number is

```text
F_n = binom(n-1, 2).
```

By THM-280, the grid reflection

```text
R_n(x,y) = (n+1-y, n+1-x)
```

induces tournament complement up to relabeling.  Define the half-tiling carrier
as the fixed line plus one side of this reflection:

```text
HCells_n = {(x,y): x-y>=2 and x+y>=n+1}.
```

## Statement

Let

```text
h_n = |HCells_n|.
```

Then

```text
h_n = (F_n + floor((n-1)/2))/2 = floor((n-1)^2/4).
```

Equivalently, beginning with tournament size `2`,

```text
0, 1, 2, 4, 6, 9, 12, 16, 20, 25, 30, ...
```

The parity-split recurrences are exact.

For even `n`:

```text
h_n = A + B - C = 2h_{n-1} - h_{n-2}.
```

For odd `n`:

```text
h_n = A + B - C + D - E - F + G
    = h_{n-1}+h_{n-1}-h_{n-2}+h_{n-2}-h_{n-3}-h_{n-3}+h_{n-4}.
```

Here `A,B` are half-tiling carriers of size `n-1`, `C,D` of size `n-2`,
`E,F` of size `n-3`, and `G` of size `n-4`.  The two `n-2` terms cancel in
cardinality but not in geometry: they occupy different overlap slots.

Thus the user's half-tiling recurrences are the mirror-folded parity refinement
of THM-442's full-staircase inclusion-exclusion recursion.

## Proof

The reflection fixed cells satisfy

```text
(x,y) = (n+1-y,n+1-x),
```

so

```text
x+y=n+1,  y<x,  x-y>=2.
```

There are exactly `floor((n-1)/2)` such cells.  Every nonfixed cell belongs to
a two-element reflection orbit.  Keeping the fixed cells and one member of each
two-cell orbit gives

```text
h_n = fixed_n + (F_n-fixed_n)/2 = (F_n+fixed_n)/2.
```

Substituting `F_n=binom(n-1,2)` and `fixed_n=floor((n-1)/2)` gives
`floor((n-1)^2/4)`.

For even `n=2k`:

```text
h_n = k(k-1),       h_{n-1} = (k-1)^2,       h_{n-2} = (k-1)(k-2),
```

so

```text
2h_{n-1}-h_{n-2}
  = 2(k-1)^2 - (k-1)(k-2)
  = k(k-1)
  = h_n.
```

For odd `n=2k+1`:

```text
h_n     = k^2,
h_{n-1} = k(k-1),
h_{n-3} = (k-1)(k-2),
h_{n-4} = (k-2)^2.
```

Since the `-h_{n-2}+h_{n-2}` terms cancel in cardinality,

```text
2h_{n-1}-2h_{n-3}+h_{n-4}
  = 2k(k-1)-2(k-1)(k-2)+(k-2)^2
  = k^2
  = h_n.
```

This proves both recurrences.

## Verification

`04-computation/half_tiling_recursion_codex_s56.py` verifies the orbit count,
the sequence, the full THM-442 third-difference count, and the two half-tiling
recurrences for `n<=15`.  The stored output is

```text
05-knowledge/results/half_tiling_recursion_codex_s56.out
```

## Scope

This theorem is a cell-count/fundamental-domain theorem.  It does not imply
that Hamiltonian-path count, OCF, maximum `H`, or LRC sector measures satisfy
the same additive recursion.  THM-442 already records that `H` is not
cell-affine.  The half-tiling carrier is instead an address quotient: it tells
which mirror-side and fixed-line data should be retained before passing to
cycle-space or proof-obligation invariants.
