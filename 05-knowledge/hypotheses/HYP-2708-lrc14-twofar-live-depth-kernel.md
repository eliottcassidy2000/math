---
id: HYP-2708
title: LRC14 two-far deviation has only four live survival depths
status: OPEN; live-depth kernel lemma proved, deviation bound still open
source: codex-2026-06-20-S67
tangent: T941
depends_on:
  - HYP-2701
  - HYP-2704
  - HYP-2705
  - THM-548
  - THM-556
related:
  - HYP-2706
  - HYP-2684
  - HYP-2675
  - OPEN-Q-108
---

# HYP-2708 - Two-Far Live-Depth Kernel

## Claim

The HYP-2701 two-far survival deviation is lower-dimensional than the full
seven-depth missed-sector profile suggests.  For the survival currency

```text
C = p1+p2+p3+p4-4p6
```

and two far hits, only before-depths

```text
t = 1, 2, 5, 6
```

can contribute to the actual-minus-decorrelated-boundary deviation.  Depths
`3` and `4` are exactly silent.

Therefore the sharp two-far proof obligation is not a seven-depth discrepancy
bound.  It is a signed four-depth kernel bound on the residual deletion
automaton.

## Exact Scout

Script:

```text
04-computation/lrc14_twofar_live_depth_kernel_codex_s67.py
```

Output:

```text
05-knowledge/results/lrc14_twofar_live_depth_kernel_codex_s67.out
```

The script decomposes exact wall atoms by:

```text
t = number of inner sectors missed by the bounded core,
h = number of distinct missed sectors hit by the two far colors.
```

It verifies the exact checksum

```text
actual C(B union {u,v}) - C_boundary(B,2)
  = sum_atoms mass * ( C(t-h) - K2(t) )
```

where `K2(t)` is the two-iid-hit death-chain boundary coefficient.

## Proved Kernel Lemma

For two iid hits into seven colors, let `K2(t)` be the expected value of `C`
after two hits starting from missed-depth `t`.  The exact coefficients
`C(t-h)-K2(t)` for `h=0,1,2` are:

```text
t=0: K2=0      coeffs  0, 0, 0        silent
t=1: K2=36/49  coeffs  13/49, -36/49, -36/49
t=2: K2=47/49  coeffs  2/49, 2/49, -47/49
t=3: K2=1      coeffs  0, 0, 0        silent
t=4: K2=1      coeffs  0, 0, 0        silent
t=5: K2=45/49  coeffs  -45/49, 4/49, 4/49
t=6: K2=26/49  coeffs  -222/49, -26/49, 23/49
```

Proof is immediate from the deterministic currency values:

```text
C(0)=0, C(1)=C(2)=C(3)=C(4)=1, C(5)=0, C(6)=-4.
```

If `t=3` or `t=4`, then after at most two hits the possible depths stay inside
`{1,2,3,4}`, where `C=1`, so actual and boundary both contribute `1`
regardless of the hit kernel.  Those depths cannot spend boundary margin.

## Risk-Bank Evidence

The S67 finite bank re-audits the tight S64 two-far rows.  Examples:

```text
k=8 floor-fail row (0,3,6,9,12,14,15,18):
  slack = -107/8820
  actual-boundary deviation = -103403/432180
  live_mass = 569/882
  deviation by depth:
    t=1: -4484/108045
    t=2: -56489/432180
    t=5: -2153/72030
    t=6: -803/21609

k=9 leader (0,4,6,8,10,12,14,15,16):
  slack = 29/980
  actual-boundary deviation = -6395/28812
  live_mass = 2077/2940
  deviation by depth:
    t=1: -46217/288120
    t=2: -53/13720
    t=5: -477/38416
    t=6: -1739/38416

k=10 leader (0,2,4,6,8,10,12,14,15,16):
  slack = 29/588
  actual-boundary deviation = -15763/144060
  live_mass = 17/35
  deviation by depth:
    t=1: -811/14406
    t=2: 219/48020
    t=5: -477/38416
    t=6: -1739/38416
```

The tight rows spend the boundary through the live depths, not through the
bulk middle depths.  In the compact S64 risk bank, aggregate live-depth
deviation was concentrated as:

```text
t=1: -156504763/191023560
t=2: -2607955393/4202518320
t=5: -216232799/1050629580
t=6: -1054201/2938824
```

The mass values in that aggregate are not normalized probabilities because
they are summed across several different rows; they are a risk-bank ledger.

## Proof Consequence

The remaining two-far lemma can be stated as a four-depth inequality:

```text
sum_{t in {1,2,5,6}} signed_kernel_debt_t(B; u,v)
  <= boundary_margin(B,k).
```

Depths `3` and `4` should be removed from the analytic target entirely.  This
matters because those middle depths are the bulk of the comfortable
decorrelated mass; bounding them wastes constants and hides the actual
obstruction.

## Post-S23 Live-Depth Ladder

After pulling the KPS S23 OPEN-Q-108 localization, the global wide `p0`
residual is now mostly a bounded-base single-far finite window, while this
HYP-2708 route remains the two-far survival branch.  The addendum script

```text
04-computation/lrc14_survival_live_depth_ladder_codex_s67.py
```

with output

```text
05-knowledge/results/lrc14_survival_live_depth_ladder_codex_s67.out
```

puts those branches in the same exact depth-kernel language:

```text
direct p0, one far:       live depths = {1}
survival C, one far:      live depths = {1,5,6}
survival C, two far:      live depths = {1,2,5,6}
survival C, three far:    live depths = {1,2,3,5,6}
survival C, four+ far:    live depths = {1,2,3,4,5,6}
```

Thus S23's single-far finite window is exactly a `t=1` closure problem for
direct `p0`, whereas the one-far survival gate adds only high-tail debt
`t=5,6`.  HYP-2708 is the first branch where shallow `t=2` debt becomes live.
For `r=3`, only the `t=4` middle plateau remains silent; after four far hits
every positive missed depth is live.

This refines the proof order:

1. **Single-far p0 window:** use S23's finite window plus a sharp transfer
   bound on the pure `t=1` closure kernel.
2. **One-far survival gate:** add only `t=5,6` tail terms.
3. **Two-far true-wide survival:** prove the HYP-2708 four-depth inequality.
4. **Three-plus far:** use aggregate death-chain margin/synchronization, not a
   middle-depth silence shortcut.

## Candidate Bound Shape

For each bounded core state atom, let `R` be its missed-sector set and let
`H_{u,v}(R)` be the actual two-color hit law restricted to `R`.  Let
`I_t` be the iid two-hit law on a set of size `t`.  The needed object is:

```text
Debt_t = <C(t-h)-K2(t), H_{u,v}(R)-I_t>
```

averaged over all core atoms with `|R|=t`, and only for `t in {1,2,5,6}`.

Expected proof split:

1. **Nonresonant far pairs:** use the signed Abel/BV or projective phase
   angle route from HYP-2705/HYP-2684 to show the live-depth hit law is close
   enough to iid.
2. **Low relation-distance pairs:** route to a finite atlas keyed by far-pair
   difference, mod-7 phase word, and bounded-core live-depth state word.
3. **k=8 dividend:** allow the known floor failures and compare to the true
   cap dividend rather than the floor.

## Tournament Analysis

Vertices are row/proof obligations, not runners.  Pairwise observable:
smaller cap slack, then more negative live-depth deviation.  Switch/gauge:
keep the live-depth kernel before scalar row invariants.

The S67 risk-bank tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3cycles=0
```

Hamiltonian risk path begins:

```text
1. k=8 (0,3,6,9,12,14,15,18)
2. k=8 (0,3,6,7,9,12,15,18)
3. k=8 (0,2,3,6,9,12,15,18)
4. k=9 (0,4,6,8,10,12,14,15,16)
```

This confirms that the live-depth quotient preserves the current cap risk
ordering.

## Status

This does not prove LRC(14).  It proves a useful exact reduction of the
two-far deviation target:

```text
seven depths -> four live depths.
```

The next proof should bound these four live-depth signed kernel debts against
the already-positive HYP-2701 boundary margins.
