---
id: HYP-2992
title: LRC14 Haar-product tile discrepancy lemma
status: PROOF-INTERFACE / exact finite product-table scout; not a proof of LRC14
source: codex-2026-06-24-S165
artifacts:
  - 04-computation/lrc14_haar_product_tile_discrepancy_codex_s165.py
  - 05-knowledge/results/lrc14_haar_product_tile_discrepancy_codex_s165.out
  - 07-reflections/lrc14-haar-product-tile-discrepancy-codex-s165.md
related:
  - HYP-2989
  - HYP-2988
  - HYP-2987
  - HYP-2986
  - HYP-2985
  - HYP-2984
  - HYP-2981
  - HYP-2979
  - HYP-2975
  - HYP-2948
  - HYP-2951
  - HYP-2963
  - HYP-2908
  - THM-572
  - OPEN-Q-108
---

# HYP-2992: LRC14 Haar-Product Tile Discrepancy Lemma

**Status:** PROOF-INTERFACE / exact finite product-table scout; not a proof of LRC14.

## Claim

The two-dimensional Haar product rule supplies the local algebra that the LRC14
tournament-tiling model has been missing.  On a labelled LRC14 packet fiber,
rectangular Haar coefficients over an endpoint/scale or wall/owner grid should
fall into the same finite interaction classes as fixed-Hamiltonian-path
tournament tiles:

```text
orthogonal zero,
same tile / boundary atom,
one-coordinate owner strip,
nested refinement,
cross-coordinate handoff.
```

The proposed proof target is a Haar-tile discrepancy lemma:

> every primitive non-AP/GW zero-open residual must expose a nonzero signed
> two-dimensional Haar tile coefficient in an owner-strip, cross-handoff, or
> nested-refinement class; if all such coefficients vanish, the only remaining
> local address is a same-tile boundary atom, forcing the AP/Goddyn-Wong
> boundary skeleton or a new THM-572/F7 state-lift atom.

## Evidence

Script:
`04-computation/lrc14_haar_product_tile_discrepancy_codex_s165.py`

Output:
`05-knowledge/results/lrc14_haar_product_tile_discrepancy_codex_s165.out`

The script enumerates all ordered products of dyadic Haar rectangles through
depth `3`.  With `225` rectangles and `50625` ordered pairs, the product classes
are:

```text
same_tile_indicator       225   0.00444
vertical_owner_strip     1020   0.02015
horizontal_owner_strip   1020   0.02015
cross_handoff            2312   0.04567
nested_refinement        2312   0.04567
orthogonal_zero         43736   0.86392
```

Every nonzero non-atom class has perfectly balanced signs:

```text
owner strips:       -510 / +510
cross handoff:     -1156 / +1156
nested refinement: -1156 / +1156
```

This is exactly the discrepancy-theory shape: most interactions are
orthogonal, while the nonzero signal is sparse, signed, local, and typed.

## Product Rule

For unnormalised one-dimensional Haar functions,

```text
h_I h_J = 0           if I and J are disjoint
h_I h_I = 1_I         if I = J
h_I h_J = +/- h_J     if J is contained in one child of I
h_I h_J = +/- h_I     if I is contained in one child of J
```

In two dimensions,

```text
h_{I x J} h_{I' x J'} = (h_I h_{I'})(h_J h_{J'}).
```

The coordinatewise rule is the whole point: equality in one coordinate and
nesting in the other creates an owner strip; nesting in both creates a
recursive state-lift packet; opposite nesting creates a zipper handoff.

## Tournament Analysis

Vertices are Haar-product interaction classes, not runners.

Pairwise observable:
labelled information retained by a Haar-product quotient.

Switch/gauge:
orient toward the carrier that preserves more of the LRC proof packet:
bulk locality, endpoint owner, boundary atom, zipper crossing, state-lift
refinement, and quotient guardrail.

Tie Hamiltonian path:

```text
same_tile_indicator
-> vertical_owner_strip
-> horizontal_owner_strip
-> cross_handoff
-> nested_refinement
-> orthogonal_zero
```

Fingerprint:

```text
score_histogram={0:1,1:1,2:1,3:1,4:1,5:1}
directed_3cycles=0
SCC_sizes=[1,1,1,1,1,1]
Hamiltonian_path_count=1
```

## Relationship to Recent LRC14 Work

HYP-2986 gives the endpoint tope/cocircuit wall language.  HYP-2987's handoff
atlas says the proof is a zipper of local certificates.  HYP-2984/HYP-2985 say
kernel and smoothing changes are admissible only when packet labels survive or
named defects are emitted.  HYP-2948/HYP-2951 give the 1D Haar/Baire boundary
split.

HYP-2992 is the proposed 2D local algebra underneath those routes.  HYP-2989's
independent Haar-square scout is the minimal `2 x 2` fixed-margin switch
version of the same phenomenon, HYP-2991 records the later local zipper
cocycle `zeta`, and HYP-2988 records the exposure-poset router; this hypothesis
records the larger dyadic rectangle interaction atlas.
The one
dimensional circle event is not discarded; it becomes one coordinate of a
rectangular packet grid.  The other coordinate is the owner/period/scale/fiber
clock that prevents the quotient from forgetting the reason a local certificate
works.

## Assumption Challenge

This session does not use runners as tournament vertices.  Plausible vertices
are Haar rectangles, endpoint-owner pairs, wall-crossing events, exact-period
packets, Fejer atom banks, K33/state-lift obligations, or proof carriers.

The quotient preserves typed local discrepancy and wall-packet interaction.
It destroys raw runner magnitude and hidden C27/Goddyn-Wong transfer data unless
those labels are attached.  This is why `same_tile_indicator` is not enough:
AP and GW can be Haar-mass-zero boundary atoms with different hidden legality.

## Next Steps

1. Build the actual LRC14 packet grid for HYP-2963 rows: one coordinate for
   endpoint or tope wall, the other for exact-period / packet-family / Fejer
   scale.
2. Compute Haar tile coefficients for AP, GW, K33, petal, splice, and covering
   rows and verify that known certificates land in the predicted classes.
3. Prove a finite vanishing lemma: if all owner-strip, cross-handoff, and
   nested-refinement coefficients vanish on a primitive zero-open packet, then
   the endpoint-owner skeleton is AP/GW or a new THM-572/F7 state lift.
4. Compare the Haar product table with the fixed-path tournament staircase
   tile cube to formalize the dictionary between rectangular coefficients and
   tile flips.
