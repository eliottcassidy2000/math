# LRC14 Taut Bridge Graph Curvature

Agent: codex-2026-06-24-S155  
Related: HYP-2975, HYP-2974, HYP-2971, HYP-2970, HYP-2969, HYP-2965, HYP-2951, HYP-2949, HYP-2908, THM-572, OPEN-Q-108

## New Angle

The current packet route already knows how to measure positive Haar-open safe
mass.  HYP-2970, which landed concurrently, adds an endpoint-credit winding
cycle dual for open danger covers.  The new S155 angle is smaller and more
local: inspect what remains exactly at the boundary where the open-cover graph
degenerates.

For a row `S`, the danger arcs are open intervals

```text
I(s,m) = ((m - 1/14)/s, (m + 1/14)/s).
```

S155 records two objects before scalarization:

```text
positive bridge:  uncovered interval between an end-owner and a start-owner
taut vertex:      isolated safe endpoint, covered on both sides, safe at the point
```

The quotient preserves the LRC predicate

```text
positive open witness / boundary-only equality / fully covered residual
```

and destroys row identity away from labelled endpoint owners.  That is why it
must remain downstream of exact `M`, qdiv, Farey branch, HYP-2970 endpoint
credit, HYP-2974 Toeplitz PSD, and C27/K33 labels.

Post-rebase HYP-2977 is the complementary global quotient.  It keeps the
safe-set Fourier shadow and high-frequency relation tails while discarding
endpoint ownership.  HYP-2975 keeps the endpoint-owner current while discarding
spectral energy.  A real proof engine should make those two shadows commute.

## Evidence

The named hard rows split cleanly.

```text
AP:          safe_mu=0, taut=6, positive bridges=0
GW 12->24:  safe_mu=0, taut=6, positive bridges=0
12->36:     safe_mu=1/1260, positive bridges=2
10->20:     safe_mu=1/980,  positive bridges=2
13->26:     safe_mu=1/182,  positive bridges=2
P10+K33:    safe_mu=4/2205, positive bridges=4
12->84:     safe_mu=563/105105, positive bridges=8
```

The AP/GW taut vertices are exactly the denominator-14 unit points:

```text
1/14:  1 -> 13
3/14:  5 -> 9
5/14:  3 -> 11
9/14:  11 -> 3
11/14: 9 -> 5
13/14: 13 -> 1
```

Every endpoint-owner pair sums to `0 mod 14`, and the total owner-current is
zero.  This is the taut-boundary version of the HYP-2970 zero-credit condition.

The bounded bank is also clean.  With one-swaps through added value `160` and
two-swaps through added value `36`, the script scans `21645` primitive rows:

```text
positive-open bridge rows: 21644
zero-curvature taut rows: 1
zero row: GW 12->24
smallest positive safe mass: 1/1260 at 12->36
```

## What This Does Not Prove

This is not a standalone LRC14 proof.  A positive bridge is just the exact
regular-open witness already recognized by the Haar/Baire route.  The useful
new thing is the boundary language for rows with no positive bridge.

The falsifier is now sharper:

```text
primitive qdiv>=14 row
no positive bridge
no AP/GW zero-sum taut current
no bounded-degree Toeplitz-negative exit
no C27/K33 state-lift debt
```

The S155 bank found no such row, but that is evidence, not a theorem.

## Tournament Analysis

The useful tournament vertices are not runners.  I considered:

```text
runners
raw endpoints
endpoint-owner labels
positive safe intervals
isolated taut points
boundary-current states
missed-depth sectors
HYP-2970 endpoint-credit cycles
HYP-2974 Toeplitz dual sections
C27/K33 proof obligations
```

Chosen layer:

```text
vertices = boundary/proof carriers
```

Pairwise observable:

```text
retention of open-witness status,
boundary equality,
endpoint ownership,
dual-certifiability,
state-lift address,
scalar-decoy resistance.
```

Switch/gauge:

```text
componentwise majority of the six retention scores,
with ties along
exact_M -> open bridge -> taut current -> endpoint credit.
```

Fingerprint:

```text
score_hist={0:1,2:2,3:1,4:2,6:1,7:1}
directed_3cycles=3
SCC sizes=[5,1,1,1]
leading order:
  endpoint_credit_winding
  taut_vertex_current
  Toeplitz_PSD_dual
  C27_K33_state_labels
  positive_open_bridge
  exact_M_Farey_scale
  multiplicity_moment_dual
  raw_runner_vertices
```

The nontrivial SCC is a warning: endpoint credit, taut current, harmonic PSD,
C27/K33 labels, and open bridges should travel as a packet.  Picking only one
of them recreates the scalarization problem that earlier routes ran into.

## Next Move

Implement the actual HYP-2970 potential solver and run it side by side with
S155's taut-current extractor and HYP-2974's negative-eigenvector localization.
The most useful comparison row is `12->36`: it is positive-open, K33-labelled,
Toeplitz-negative only at a higher degree, and has the smallest positive bridge
in the S155 bank.
