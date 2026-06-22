---
id: HYP-2864
title: LRC14 perturbed sheet discrepancy bridge -- exact sheet spacing becomes gcd-discrepancy, with bounded-quotient fallback to THM-405/finite atlases
status: PROVED bridge lemma for exactly dilated bases plus parked runners; open as full perturbed-cluster closure
source: codex-2026-06-22-S92
related:
  - HYP-2858
  - HYP-+2863
  - HYP-2853
  - HYP-2861
  - HYP-2862
  - THM-405
  - THM-562
  - THM-527
  - OPEN-Q-108
---

# HYP-2864 -- Perturbed Sheet Discrepancy Bridge

## Claim

The exact sheet spacing in the pure-dilation proof is the zero-loss case of a
finite discrepancy lemma.  Let `bE` be an exactly dilated base and add parked
runners `W_1,...,W_h`.  Write

```text
tau = (n+u)/b,    n=0,...,b-1,
r_l = W_l mod b, g_l = gcd(r_l,b).
```

For each fixed `u` in the base-good set `G_E`, the base `bE` is safe on every
sheet.  The runner `W_l` can spoil at most

```text
b/7 + g_l
```

of the `b` sheets.  Hence at least

```text
b*(1 - h/7) - sum_l g_l
```

sheets survive.  Integrating over `u in G_E`,

```text
L(bE union {W_l}) >= meas(G_E) * (1 - h/7 - (sum_l gcd(W_l,b))/b).
```

This is a discrete finite-`V` substitute for an Erdos-Turan/Koksma step.  It
turns exact equal spacing into a controlled gcd-discrepancy term.

## Proof

For one parked runner `W`, write `r = W mod b` and `g = gcd(r,b)`.  The sheet
map

```text
n -> frac(r n / b + theta)
```

visits the `b/g` equally spaced points with multiplicity `g`.  Any danger
interval for LRC14 has length `1/7`.  On an equally spaced `N`-grid an interval
of length `ell` contains at most `N ell + 1` grid points.  Multiplying by `g`
gives

```text
#danger <= g*((b/g)/7 + 1) = b/7 + g.
```

Union bound over the `h` parked runners gives the surviving-sheet count above.
No limiting argument is used.

## Bounded-Quotient Fallback

If the sheet floor is not positive, then

```text
sum_l g_l >= b*(7-h)/7.
```

Thus at least one parked runner has

```text
g_l >= b*(7-h)/(7h),
```

so its reduced sheet quotient satisfies

```text
b/g_l <= 7h/(7-h).
```

For `h=1,...,6`, the fallback quotient bounds are:

| h | quotient bound |
|---:|---:|
| 1 | `7/6` |
| 2 | `14/5` |
| 3 | `21/4` |
| 4 | `28/3` |
| 5 | `35/2` |
| 6 | `42` |

So the failure of the discrepancy branch is not amorphous.  It produces a
small quotient/resonance channel.  That channel should route to a bounded
residue atlas; if the actual speed ratio is at most `13`, THM-405 closes it
immediately, and otherwise it belongs to the existing bounded-ratio/finitary
LRC<=13 handoff stack (THM-562-style finite resonance atlas plus discrepancy
tail).

## Scope

This proves a bridge lemma for exactly dilated bases with parked runners.  It
does not by itself close arbitrary perturbed clusters where every runner
leaves the `bE` lattice.  The intended next split is:

1. normalize a coordinated-growth row to `bE` plus a bounded number of parked
   residue runners;
2. use the sheet-discrepancy floor when the residues have small gcds with `b`;
3. route large-gcd failures to a bounded quotient/resonance atlas;
4. route truly bounded speed ratio to THM-405.

This is the concrete version of the user's proposed path: exact spacing is the
lossless case; perturbed spacing is controlled by a discrepancy term; the
`V <= Cb` residue cases are bounded-ratio or finite-atlas work rather than a
new asymptotic obstruction.

## Computation

Script:

```text
04-computation/lrc14_perturbed_sheet_discrepancy_codex_s92.py
```

Output:

```text
05-knowledge/results/lrc14_perturbed_sheet_discrepancy_codex_s92.out
```

The script verifies the one-runner bound on `2000` exact random rows with zero
failures.  A tight sampled row has

```text
b=217, r=528, gcd=1, count=31, bound=32.
```

It also prints the bounded-quotient fallback table above and sample multi-runner
floors.

## Tournament Analysis And Assumption Challenge

Pairwise observable: surviving sheet fraction after parked-runner damage.

Switch/gauge: a proof carrier wins if it converts finite sheet counts into a
positive witness without invoking an asymptotic limit.

Hamiltonian path:

```text
sheet_discrepancy > bounded_quotient_gate > THM405_first_window > finite_residue_atlas > raw_Erdos_Turan_bound > runner_vertices
```

Challenged assumption: the Node-1 finite-`V` step must be treated as a
continuous equidistribution theorem.  In this quotient, the finite sheet count
is exact up to the arithmetic term `gcd(r,b)`, and failure of the bound exposes
a bounded quotient rather than an analytic wall.
