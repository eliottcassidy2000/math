---
id: HYP-2849
title: LRC origin-decoy divisor-depth ledger -- bounded Farey-neighborhood subtraction collapses to denominator obligations plus the p=1 decoy
status: PROOF REDUCTION; exact interval verification through q=7, sampled through q=12, deepest-branch scans through q=10
source: codex-2026-06-22-S89
related:
  - HYP-2847
  - HYP-2844
  - HYP-2846
  - HYP-2848
  - HYP-2842
  - THM-527
---

# HYP-2849 -- origin decoy and divisor-depth ledger

## Claim

For the bounded lower-threshold model in HYP-2847, the proper Farey
neighborhood part is not really an interval-union problem.  At span
`S=2q-1`, every proper neighborhood around a reduced center `a/b`,
`2 <= b < q`, is touched by a `G_P` hole for speed `p` only when `b | p`.
Therefore the proper residual is exactly a divisor-depth ledger:

```text
sum_{2 <= b < q} phi(b) * 2 * max(0, r_b - 1/(2q p_b)),
where r_b = (q-b)/(b q (2q-1)),
and p_b = min {p in P : b | p}, or infinity if no such p exists.
```

The origin neighborhood is not just an extra interval.  It is a decoy slot:
if `p=1` is absent, origin mass survives; if `p=1` is present, the origin is
killed, but one scarce `P` slot has been spent and the divisor-depth cover of
proper denominators is weaker.

## Proof Skeleton

Two elementary separation inequalities drive the reduction.

First, proper neighborhoods are pairwise disjoint.  Distinct reduced Farey
centers `a/b` and `c/d` satisfy `|a/b-c/d| >= 1/(bd)`.  Their bounded-span
radii are

```text
r_b = (q-b)/(b q S),   r_d = (q-d)/(d q S),   S=2q-1.
```

The inequality `r_b+r_d < 1/(bd)` follows after clearing denominators from

```text
d(q-b) + b(q-d) < q(2q-1),
```

with `2 <= b,d < q`; the bilinear maximum occurs on the boundary of the
rectangle and is below `q(2q-1)`.

Second, a non-divisor hole cannot touch a proper `a/b` neighborhood.  If
`b ∤ p`, then every `j/p` has

```text
|a/b - j/p| >= 1/(bp).
```

The hole radius is `1/(2qp)`, so

```text
1/(bp) - r_b - 1/(2qp)
  = [q(2q-1) - p(q-b) - b(2q-1)/2] / [b q (2q-1) p]
  >= 1/(2q(2q-1))
```

because `p <= 2q-1`.  Thus only selected multiples of `b` matter, and among
those multiples only the smallest selected multiple matters.

This gives the exact ledger above.  The origin interval has radius
`R_0=(q-1)/(q(2q-1))`.  Since `R_0 < 1/(2q)`, the `p=1` hole kills it
completely.  If `p=1` is absent, the origin interval becomes the fallback
positive mass.

## Computed Evidence

Script: `04-computation/lrc_origin_decoy_depth_ledger_codex_s89.py`

Transcript: `05-knowledge/results/lrc_origin_decoy_depth_ledger_codex_s89.out`

The script verifies the divisor-depth reduction by exact interval subtraction:

| q | status |
|---:|---|
| 3..7 | exhaustive exact verification, `0` mismatches |
| 8..12 | deterministic sampled exact verification, `0` mismatches |

The observed non-divisor separation margin is exactly `1/(2q(2q-1))` in every
row checked.

Deepest bounded dense branch scan (`k=q+1`, `|P|=q-2`):

| q | N | minimum | worst P | live proper contributor |
|---:|---:|---:|---|---|
| 3 | 6 | `1/15` | `(1,)` | `b=2` |
| 4 | 8 | `5/168` | `(1,6)` | `b=2` |
| 5 | 10 | `1/60` | `(1,3,4)` | `b=2` |
| 6 | 12 | `5/264` | `(1,3,4,5)` | `b=2` |
| 7 | 14 | `2/273` | `(1,2,3,4,5)` | `b=6` |
| 8 | 16 | `1/90` | `(1,2,3,4,5,7)` | `b=6` |
| 9 | 18 | `1/204` | `(1,2,3,5,6,7,8)` | `b=4` |
| 10 | 20 | `1/342` | `(1,2,3,4,5,7,8,18)` | `b=6` |

In every scanned deepest row, the global minimizer includes `p=1`.  The
origin interval is killed, but the lost `P` slot leaves a small proper
denominator obligation alive.  Without `p=1`, the proper residual can be zero,
but the origin tail is much larger; for example at `q=7`, no-`p=1` minimum is
`11/182`, while the global minimum is `2/273`.

## LRC14 Readout

For the HYP-2847 stress row `q=7`, `k=8`, `|P|=5`, this gives a compact proof
route:

1. If `1 ∉ P`, the origin interval contributes at least `11/182`.
2. If `1 ∈ P`, the origin interval is killed but the divisor-depth ledger over
   denominators `2,3,4,5,6` has minimum `2/273`, attained by
   `P=(1,2,3,4,5)` with only denominator `b=6` live.
3. Hence the origin+proper bounded certificate has residual at least `2/273`.

This is weaker than the HYP-2844/HYP-2846/HYP-2848 witness-floor and
q-widening route, but it is much easier to prove directly in the bounded
stress row: interval geometry reduces to two separation inequalities plus a
finite divisor-depth optimization.

## Tournament Analysis And Assumption Challenge

Pairwise observable: exact residual after subtracting `G_P` holes.

Switch/gauge: proof compression wins when it preserves exact residual.

Hamiltonian path:

```text
divisor_depth_ledger > origin_decoy_slot_tradeoff > proper_interval_geometry > low_denominator_resonant_route > raw_interval_union > exact_farey_centers > raw_runner_vertices
```

Challenged assumption: the useful vertices are not runners, arcs, or exact
Farey centers.  They are denominator obligations and decoy slots.  This
quotient preserves the exact proper residual and the origin/no-origin
dichotomy; it discards irrelevant center identities and most raw interval
geometry.
