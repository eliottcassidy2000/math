---
id: THM-557
title: LRC single-block decorrelated extremality and diagonal-freeze error
status: PROVED; finite exact coherent-block quotient plus elementary BV error
source: codex-2026-06-20-S61
depends_on: []
related:
  - HYP-2694
  - HYP-2684
  - HYP-2675
  - HYP-2695
  - THM-546
  - THM-556
---

# THM-557 - LRC Single-Block Decorrelated Extremality

## Statement

Fix the HYP-2694 coherent-block quotient for LRC(14).  The speed `0` is an
anchored singleton cluster.  The `m=k-1` nonzero speeds are partitioned into
far coherent consecutive blocks, and all far blocks share the slow coordinate
`x` but have independent carrier phases.

For `m=7,8,9,10,11`, the exact shared-`x` decorrelated cover is maximized by
the one-part partition `[m]`.  The exact single-block values and cap margins
are:

```text
k=8,  m=7:  D_m=283/1470,    cap-D_m=1111/5880
k=9,  m=8:  D_m=629/2058,    cap-D_m=111019/588588
k=10, m=9:  D_m=16969/41160, cap-D_m=102803/535080
k=11, m=10: D_m=30551/61740, cap-D_m=184957/802620
k=12, m=11: D_m=71111/123480, cap-D_m=34729/123480
```

The closest split is always `[m-1,1]`, with exact split gaps:

```text
m=7:  1111/10290
m=8:  374/5145
m=9:  6561/96040
m=10: 42661/864360
m=11: 9047/172872
```

For the actual shifted single-block row

```text
E_M = {0} union {M,M+1,...,M+m-1},
```

the finite-scale error from the decorrelated single-block value satisfies

```text
|p0(E_M)-D_m| <= 7*binom(m,2)/M.
```

Consequently the shifted single-block branch is rigorously below `cap_k` once

```text
M > 7*binom(m,2)/(cap_k-D_m).
```

The corresponding sufficient cutoffs for `k=8..12` are:

```text
k=8:  779
k=9:  1040
k=10: 1312
k=11: 1367
k=12: 1369
```

## Proof

For the finite decorrelated quotient, enumerate all integer partitions of
`m=k-1`, attach the anchored cluster `(0,)`, and compute the shared-`x`
carrier integral exactly.  For each fixed `x`, a far block contributes the law
of its covered inner-sector set as its carrier phase varies uniformly on
`[0,1)`.  Convolving these laws and integrating over the exact `x`-breakpoints
gives a rational value for every partition.  The finite exact enumeration in
`04-computation/lrc14_single_block_extremality_margin_codex_s61.py` checks all
`15,22,30,42,56` partitions for `m=7..11` and verifies that `[m]` is the unique
maximum in each case.

For the finite-scale error, write

```text
H(x,phi) = 1{ {phi + d x : 0<=d<m} hits all six inner sectors }.
```

Then

```text
p0(E_M) = integral_0^1 H(x, Mx) dx,
D_m     = integral_0^1 integral_0^1 H(x, phi) dphi dx.
```

Split `[0,1)` into cells `[a/M,(a+1)/M)`, set `x=(a+t)/M`, and use that
`phi=Mx mod 1=t`.  On each cell, compare `H((a+t)/M,t)` with the cell average
`integral_0^1 H((a+s)/M,t) ds`.  For fixed carrier phase `t`, the function
`x -> H(x,t)` changes only when one of the `m` block points crosses a sector
wall.  The point with offset `d` has at most `7d` such crossings, so the total
variation in `x` is at most

```text
J_m = sum_{d=1}^{m-1} 7d = 7*binom(m,2).
```

Summing the cell variations and multiplying by the cell width gives
`|p0(E_M)-D_m| <= J_m/M`.

## Scope

This theorem proves the exact coherent-block quotient and the single shifted
block's diagonal-freeze error.  It does not prove the full HYP-2694
multi-carrier discrepancy bound for arbitrary split cluster shapes.  That
remaining task is the joint decorrelation step in HYP-2684/HYP-2694.

## Computational Check

The stored output is
`05-knowledge/results/lrc14_single_block_extremality_margin_codex_s61.out`.
