---
id: HYP-2989
title: LRC14 Haar-product discrepancy and tournament-tiling square
status: PROOF-INTERFACE / synthesis scout and quotient guardrail, not a proof
source: codex-2026-06-24-S165
artifacts:
  - 04-computation/lrc14_haar_product_discrepancy_tiling_codex_s165.py
  - 05-knowledge/results/lrc14_haar_product_discrepancy_tiling_codex_s165.out
related:
  - HYP-2988
  - HYP-2987
  - HYP-2986
  - HYP-2985
  - HYP-2984
  - HYP-2981
  - HYP-2978
  - HYP-2963
  - HYP-2951
  - HYP-2949
  - HYP-2595
  - HYP-2594
  - HYP-2736
  - HYP-2745
  - OPEN-Q-108
---

# HYP-2989: Haar Product Discrepancy And Tournament Tiling

This pass pulls the recent S162/S163/S164 packet work back into the older
colored-discrepancy thread.  The key observation is elementary but structural:

```text
h_I(x) h_J(y) = h_{I x J}(x,y)
```

On the two children of each dyadic interval this is the sign matrix

```text
[[+1, -1],
 [-1, +1]].
```

That same matrix is also:

```text
the elementary 2-by-2 fixed-margin switch,
the mixed Haar discrepancy atom,
the tope wall-crossing square,
and the oriented square tile hidden in the tournament tiling model.
```

So the Haar product rule is not just an analogy to tournament tiling.  It is
the minimal square on which a row/column quotient can look identical while the
mixed orientation sign flips.

## Computed Artifact

Script:

```text
04-computation/lrc14_haar_product_discrepancy_tiling_codex_s165.py
```

Stored output:

```text
05-knowledge/results/lrc14_haar_product_discrepancy_tiling_codex_s165.out
```

The script records:

```text
h_x=(1,-1), h_y=(1,-1)
h_x tensor h_y = [[1,-1],[-1,1]]
fixed-margin switch = [[1,-1],[-1,1]]
```

The minimal quotient collision is:

```text
diagonal_packet      [[1,0],[0,1]] rows=(1,1), cols=(1,1), Hxy= 2
anti_diagonal_packet [[0,1],[1,0]] rows=(1,1), cols=(1,1), Hxy=-2
```

Thus fixed row and column margins are not theorem-safe unless the mixed Haar
sign or an equivalent switch/cocircuit label is retained.  The wall-crossing
move is exactly:

```text
anti_diagonal + [[1,-1],[-1,1]] = diagonal
Hxy jump = 4
```

Exhaustive 2-by-2 stress through total mass `6` finds that every table with a
nontrivial fixed-margin mate has a mate changing `Hxy`.  This is the smallest
possible warning version of HYP-2978's quotient rule.

## Synthesis Across Recent Agents

HYP-2594 counted continuous interval components `K` and proved the crude bound

```text
actual_count >= V*Sigma - K.
```

HYP-2595 then showed that only color-compatible resonances survive in the
Fourier/color expansion.  HYP-2989 explains the geometry of that improvement:
the raw boundary components are row/column shadows, while the surviving
defect lives in mixed Haar product coefficients.

HYP-2986's tope/cocircuit pass is the same square in cyclic endpoint language.
An open all-safe tope is the positive-cell side of the switch; AP/GW boundary
cocircuits are zero-dimensional wall atoms; a genuine residual would have to
survive the mixed-product switch without exposing an open witness or a boundary
owner law.

HYP-2987's certificate handoff atlas supplies the gluing rule: a handoff may
forget margins, endpoint lists, or raw component counts only if it preserves
the mixed product packet or reconstructs it by Fejer, Ramanujan, endpoint,
Kaczynski, or state-lift labels.

HYP-2985's smoothing dispatcher supplies the analytic version: a smoothing
kernel chosen on a packet family must preserve the mixed Haar product mode or
route its loss as a boundary/exceptional approach class.

HYP-2981's Fejer interval certificates are high-order relatives of the same
mixed product atom.  A Fejer quadratic form is a positive trigonometric
packet that detects a signed mixed defect after lower quotients have been
fixed.

## Theorem-Facing Target

The next useful lemma is a Haar-product replacement for the raw component
count in HYP-2594:

```text
Colored Haar-product discrepancy lemma.
After fixing the row/column packet fibers in the colored CRT placement layer,
the deficit V*Sigma - actual_count is controlled by the number and strength of
color-compatible mixed Haar switches, not by the raw continuous component
count K.
```

A plausible bound remains the HYP-2595 shape:

```text
Delta(P,E,V) <= C * (k + c_GP)
```

where `k=|E|` and `c_GP` is the number of small-speed safe components.  HYP-2989
adds the mechanism: `k+c_GP` should count independent mixed-product switch
families after the row/column margins are fixed, while `K` counts every
micro-boundary component before cancellation.

If formalized, this connects directly to the HYP-2987 zipper theorem arrows:

```text
O1 source-kernel exclusion:
  reduce arbitrary rows to colored packet fibers.

O3 family compression:
  compress many endpoint rectangles into a bounded set of mixed Haar switches.

O4 admissible smoothing:
  allow Fejer/Kaczynski smoothing only after the mixed product mode is retained.

O6 F7 definition:
  any residual sector must be a named mixed-product harmonic sector, not an
  anonymous row/column-margin failure.
```

## Tournament Analysis

Vertices are proof carriers, not runners:

```text
labelled_haar_square_packet
fixed_margin_switch_cocircuit
mixed_haar_discrepancy
colored_resonance_congruence
tope_cocircuit_wall_label
row_column_margin_shadow
raw_component_count_K
```

Pairwise observable:

```text
retention of LRC predicate,
row/column margins,
mixed Haar sign,
fixed-margin switch,
color resonance,
Fejer/dual compatibility,
tope/cocircuit labels,
quotient guardrail.
```

Stored fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
SCC_sizes=[1,1,1,1,1,1,1]
hamiltonian_paths=1
canonical_path=
  labelled_haar_square_packet
  > fixed_margin_switch_cocircuit
  > mixed_haar_discrepancy
  > colored_resonance_congruence
  > tope_cocircuit_wall_label
  > row_column_margin_shadow
  > raw_component_count_K
```

## Assumption Challenge

Alternate vertex sets considered: runners, colors, dyadic rectangles,
row/column margin fibers, 2-by-2 switches, endpoint topes, boundary
cocircuits, Fourier modes, Fejer atoms, and proof obligations.  The chosen
vertices are proof carriers because the question is which quotient is allowed
to forget the mixed product sign.

The quotient preserves the LRC placement predicate only when it retains the
mixed Haar product packet or routes its loss through a labelled certificate.
It destroys runner identity, raw endpoint order, and continuous component
counts.  That destruction is acceptable only after the mixed switch packet is
reconstructed, annihilated by a dual certificate, or parked in a named residual
sector.

## Bottom Line

The 2D Haar product rule creates the same structure as the tournament tiling
model because both are controlled by one elementary square:

```text
observer row choice x observer column choice -> mixed orientation sign.
```

The fixed Hamiltonian path chooses the observer axes.  The Haar product tile
records the sign that the row/column quotient would forget.  This is exactly
the guardrail the LRC14 proof stack has been converging toward.
