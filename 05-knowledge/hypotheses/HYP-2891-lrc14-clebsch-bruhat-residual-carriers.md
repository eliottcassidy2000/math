---
id: HYP-2891
status: STRUCTURAL SIGNAL / proof-target; exact finite carrier audit
source: codex-2026-06-22-S105
tags: [lrc14, clebsch, bruhat, truncated-octahedron, unital, residual-masks, jensen, tournament-analysis]
related:
  - HYP-2890
  - HYP-2889
  - HYP-2885
  - HYP-2887
  - HYP-2702
  - HYP-2698
  - OPEN-Q-108
results:
  - 04-computation/lrc14_clebsch_bruhat_unital_carriers_codex_s105.py
  - 05-knowledge/results/lrc14_clebsch_bruhat_unital_carriers_codex_s105.out
---

# HYP-2891: Clebsch/Bruhat carriers for the signed LRC14 residual

The user suggested three finite-geometry graph prompts: unital, Clebsch graph,
and truncated octahedral graph.  The productive reading is not that any of them
directly proves LRC14.  They supply two exact quotient carriers for the current
HYP-2889/HYP-2890 state.

```text
Clebsch graph              -> tangent residual-mask covariance quotient
truncated octahedral graph -> Bruhat S4 local compression/swap quotient
unital/design language     -> tangent/secant incidence discipline
```

## Clebsch as tangent residual quotient

S105 builds the Clebsch graph as the folded 5-cube.  It verifies:

```text
vertices = 16
edges    = 40
degree   = 5
strongly regular parameters (16,5,0,2)
```

Closed neighborhoods give an exact pair-balanced design:

```text
16 blocks, each of size 6;
each point lies in 6 blocks;
each pair lies in exactly 2 blocks.
```

This matches the design signal from the Clebsch page and is the useful LRC
object: a finite covariance carrier.

For LRC residual masks, choose a tangent/anchor sector among the six inner
sectors, delete it, and fold the remaining five Boolean coordinates by
complement.  The 64 residual masks map to 16 Clebsch classes, four masks per
class.  S105 finds:

```text
classes = 16
preimage size = 4 for every class
mixed-depth classes = 16
```

Every Clebsch class mixes missed-depth values.  Therefore Clebsch is **not** a
scalar missed-depth quotient.  It is a signed pair-covariance quotient.  This
fits the S31l/HYP-2890 message: keep the labelled signed functional; do not
collapse to `q_t` or scalar additive energy too early.

## Bruhat S4 as local compression quotient

S105 builds the truncated-octahedral graph as the Bruhat/Cayley graph on `S4`
with adjacent transposition generators.  It verifies:

```text
vertices = 24
edges    = 36
degree   = 3
commutation squares = 6
braid hexagons      = 8
automorphism group size = 48
edge orbit sizes    = 12, 24
```

The two edge orbits are exactly the warning needed for HYP-2889:

```text
outer swaps s1/s3: one orbit of size 24
middle swaps s2:   one orbit of size 12
```

Local compression cannot be treated as an unlabeled monotone path.  The middle
swap is structurally different.  This mirrors HYP-2889's exact compression
counterexamples: profile-improving moves can lower `p0` unless sector/Fourier
labels are retained.

## Unital-style lesson

The relevant unital idea is the tangent/secant discipline of a 2-design: pair
incidence must be controlled before scalarizing.  LRC transfer:

```text
choose a tangent sector / anchor;
use Clebsch closed neighborhoods as the pair-balanced covariance design;
then analyze Bruhat adjacent-swap faces before collapsing to scalar energy.
```

This is weaker than a literal unital, but stronger than an analogy.  It gives a
finite carrier that records exactly the kind of labelled pair incidence scalar
additive energy forgets.

## Proof target

Refine the HYP-2890 residual-leak inequality

```text
R_sf(E)-R_sf(AP) <= Gamma_sf(A*(AP)-A*(E))
```

by decomposing `R_sf` into:

1. tangent Clebsch closed-neighborhood components;
2. Bruhat `S4` compression faces.

Candidate theorem:

```text
Clebsch square/commuting-face components are nonpositive by pair design balance;
remaining braid-hexagon components are finite AP/Freiman low-depth packets;
high-depth packets route to HYP-2636/HYP-2887 signed cancellation.
```

The concrete next computation is to take the HYP-2890 worst residual-leak rows,
track compression moves in Bruhat face coordinates, and test whether the
near-violation mass sits on braid hexagons rather than commutation squares.

## Tournament Analysis

Vertices are proof carriers:

```text
AP_Jensen_labelled_functional,
Clebsch_tangent_residual_design,
Bruhat_S4_compression_carrier,
unital_tangent_secant_incidence,
scalar_additive_energy.
```

The observable is `(preserves signed labels, finite exactness, compression
relevance, scalar-loss warning)`.  The transitive path is:

```text
AP_Jensen_labelled_functional
  > Clebsch_tangent_residual_design
  > Bruhat_S4_compression_carrier
  > unital_tangent_secant_incidence
  > scalar_additive_energy.
```

Assumption challenged: residual-mask vertices, compression edges, or additive
energy can be unlabeled.  The exact carriers say the opposite: Clebsch preserves
pair covariance while destroying missed depth, and Bruhat preserves swap labels
precisely where local compression monotonicity fails.
