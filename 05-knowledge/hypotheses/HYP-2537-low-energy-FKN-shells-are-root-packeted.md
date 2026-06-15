# HYP-2537 - Low-energy FKN shells are root-packeted, not Hamming-uniform

**Status:** OPEN synthesis, with exact shell-1 and shell-2 anchors.
**Source:** codex-2026-06-15-S11, extending T823 / THM-511 / THM-513.
**Companions:** THM-513, THM-511, THM-284, THM-301, T823.
**Computation:** `04-computation/fkn_root_packet_shells_codex_s11.py`;
stored output
`05-knowledge/results/fkn_root_packet_shells_codex_s11.out`.

## Statement

The fixed-base-path tiling cube has `m=C(n-1,2)` free coordinates, but those
coordinates are not featureless Boolean directions. They are the non-simple
positive interval roots of `A_{n-1}`. The emerging conjectural principle is:

```text
near the transitive ground state, the meaningful local object is not
Hamming weight alone but the interval-root packet carried by the active tiles.
```

More concretely, for any low-shell tiling `S` with active root set
`R(S)={t_1,...,t_k}`, define the **packet graph**

```text
P(S):
  vertices  = active tiles t_i, weighted by gap g_i = x_i-y_i-1,
  edges     = tagged by same-end / cross-end / disjoint incidence,
  higher faces = relay / overlap patterns among triples and larger packets.
```

**HYP-2537.** For fixed low degree and low shell radius, the restriction of the
natural tournament invariants (`score`, `c3`, low-degree Fourier coefficients
of `H`, OCF packet counts) is governed by `P(S)` up to bounded-depth packet
corrections, not by raw Hamming weight alone.

Equivalently: the local FKN neighborhood should be read as a **packet algebra of
interval roots**.

## Exact Anchors Already Landed

THM-513 gives the first two shells exactly.

1. **Shell 1.** Every free tile is its own iso class. The one-flip shell has
   exactly `m=C(n-1,2)` classes, one per interval root. Their easy invariants
   collapse only by root gap:

   ```text
   score defect = e_y - e_x,
   c3           = x-y-1,
   H            = 2^(x-y-1) + 1.
   ```

   So `H` and `c3` see only root height, while score remembers the root address.

2. **Shell 2.** The first nonlinear cyclic statistic is already packeted:

   ```text
   same-end  -> c3 = g1+g2-1,
   cross-end -> c3 = g1+g2+1,
   disjoint  -> c3 = g1+g2.
   ```

   And `H`'s quadratic defect obeys the same incidence split:
   same-end negative, cross-end positive, disjoint nonnegative in the stored
   data through `n=8`.

These two exact shells already show why the naive count `2^m` is too coarse.
Even at radius `2`, "which two tiles were flipped?" matters through the packet
incidence graph.

## Why This Is the Right FKN Extension

T823 and THM-511 translate FKN into the tournament project as:

```text
level-1 / dictator / score content  <-> near-transitive ranking content,
level-2+ / packet / cycle content   <-> genuine cyclic content.
```

HYP-2537 sharpens the **dictator side** of that statement.

In the tiling chamber, the dictators are not uniform coordinates. They are
interval roots of varying height and address. The one-flip shell is therefore a
triangular atlas of possible single upsets, not a homogeneous Boolean sphere.
The first obstruction to Hamming-uniformity is shell 1 itself; the first packet
interaction is shell 2.

## Recursive Form

The `n -> n+1` growth matches the root-layer recursion.

- Old packet data is preserved for old tiles (THM-299).
- The new vertex contributes the new top row of roots `(n+1,y)`, `1 <= y <= n-1`.
- So the local shell geometry grows by **adding one new root row** to the
  triangular carrier.

This gives a recursive explanation for why the local FKN neighborhood is
triangular/interval-structured rather than just Boolean.

## Prediction

The usable metrics near the transitive state should be:

1. `energy =` number of active tiles.
2. `height mass = sum_i g_i`.
3. `packet incidence counts =` number of same-end, cross-end, and disjoint
   edges in `P(S)`.
4. `relay / overlap faces =` higher packet data controlling the first genuine
   shell-3 and shell-4 corrections.

I expect low-degree Fourier coefficients and low-cycle statistics on shell `k`
to stabilize as explicit functions of these packet metrics.

## Tournament Analysis

Rejected vertex set: raw Hamming coordinates. That quotient destroys the
interval-root address, which THM-513 shows is already visible at shell 1.

Chosen vertex set: active interval roots / packet graph.

Pairwise observable:

```text
incidence(t_i, t_j) in {same-end, cross-end, disjoint}.
```

This is not itself a tournament relation, which is why the default
tournament-analysis wrapper is not the clean tool here. The preserved structure
is interval incidence, not pairwise dominance.

## Next Moves

1. Prove the exact disjoint-pair quadratic `H` law beyond the sign pattern.
2. Identify the first shell-3 packet correction and compare it to the
   Möbius/Walsh `A+B+C-D-E-F+G` derivative.
3. Test whether low-degree Fourier data on shell `k` is determined by the
   packet graph `P(S)` plus bounded higher-face data.
