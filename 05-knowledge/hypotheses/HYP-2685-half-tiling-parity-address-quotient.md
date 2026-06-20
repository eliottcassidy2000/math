---
id: HYP-2685
title: Half-tiling parity address quotient for tournament recursion and proof ledgers
status: OPEN; THM-550 proves the parity carrier count, structural transfer remains exploratory
source: codex-2026-06-20-S56
depends_on:
  - THM-550
  - THM-549
  - THM-280
  - THM-442
  - THM-513
related:
  - THM-551
  - HYP-2689
  - HYP-2660
  - HYP-2683
  - HYP-2684
  - HYP-2681
  - HYP-2682
  - OPEN-Q-108
---

# HYP-2685 - Half-Tiling Parity Address Quotient

## Claim Being Tested

The half-tiling model is a useful proof quotient because it keeps exactly the
data that survives mirror/complement folding:

```text
fixed-line tiles + one side of the mirror line.
```

THM-550 proves the carrier count and the parity recurrences.  The hypothesis is
that this quotient should be used as an address layer before scalarizing
tournament or LRC-style proof states.

The lesson is the same as HYP-2660 and HYP-2683: choose the closure/gauge first,
then ask scalar questions.

## Connections

### Full tiling recursion

THM-442 identifies the full tiling recursion

```text
A+B+C-D-E-F+G
```

as the third finite difference of the quadratic full-cell count.  THM-550 says
the mirror-folded half carrier splits by parity:

```text
even n:  A+B-C
odd n:   A+B-C+D-E-F+G.
```

The odd half-tiling still has three visible corner slots `A,D,B`, but the two
`n-2` terms have different geometric roles despite equal cardinality.  This is
the first place the half model adds information beyond the full cell count.

### OCF and Hamiltonian paths

The additive recursion is not a recursion for `H(T)`.  THM-442 already warns
that Hamiltonian-path count and OCF live in cycle space, not cut/cell-affine
space.  The half-tiling quotient should therefore be paired with THM-513's
root-packet view:

```text
half-domain cell -> interval-root atom -> packet incidence -> cycle-space effect.
```

This suggests a new OCF proof route: compute Walsh/OCF contributions by mirror
orbits of interval roots, treating fixed-line roots separately.  The fixed line
is the likely parity-root analogue of the tournament/even-graph bijection in
HYP-2660.

### Tile-address refinement

THM-551 turns the quotient into a local coordinate system.  A tile `(a,b)` has
address

```text
(beta,tau)=(a,a+b-1),
```

where `beta` is the full birth strip and `tau` is the half-tiling fixed-line
crossing.  This refines the address-quotient claim: the half carrier is not
only a global Burnside count, but a layer clock for individual interval roots.
HYP-2689 asks whether complement-even invariant computations can update by the
`tau=n` crossing layer before importing OCF/Hamiltonian-path packet data.

### Relation to `A+B+C-D-E-F+G` in LRC work

HYP-2681 used the seven packets

```text
A,B,C; D,E,F; G
```

for three-far LRC corrections, where the user's recursion
`A+B+C-D-E-F+G` is a pair-tax shadow, not the actual residual.  The half-tiling
recursion gives a cleaner geometric model for why these signed seven-term
expressions keep appearing: they are inclusion-exclusion packets over a folded
carrier, and the sign pattern only becomes meaningful after the address slots
are fixed.

This does not repair the corrected HYP-2675 target by itself.  Incoming KPS S19
now shows that `Q(k-1)` is the decorrelated limit, not a finite wide bound; the
finite target is still cap-level via joint ET-Koksma.  The half-tiling analogy
should be used for routing and address repair, not for claiming an additive
LRC inequality.

## Proposed Proof Uses

1. **Mirror fundamental-domain OCF.**  Rewrite low-degree Walsh/OCF coefficients
   by reflection orbits of interval roots.  Fixed-line cells should be treated
   like forced parity-root edges in the tournament/even-graph bijection.
2. **Half-shell FKN packets.**  Refine THM-513's one-flip and two-flip root
   packets by which side of the mirror line they occupy and whether they lie on
   the fixed line.
3. **Self-complementary class census.**  Use half-tiling fixed-line data to
   separate genuinely self-complementary classes from classes paired only after
   relabeling.
4. **LRC address transfer.**  Use the half-tiling distinction "fixed line versus
   one side" as an analogy for HYP-2683's state-mass/private-sector address:
   keep boundary/fixed data before scalar plateau or decorrelation estimates.
5. **Even/odd recursion warning.**  Any recursive proof that ignores tournament
   parity loses corner geometry; even and odd tournament sizes have different
   half-carrier laws.

## Tournament Analysis

Vertices are proof lenses, not tournament vertices:

```text
odd_half_three_corner
mirror_orbit_fold
even_half_two_corner
full_IE_third_difference
root_packet_shell
OCF_cycle_space
LRC_address_analogy
```

Pairwise observable:

```text
(number of preserved predicates, overlap with opponent, declaration order).
```

Switch/gauge: larger observable wins.  The stored script gives a transitive
tournament with one Hamiltonian path:

```text
odd_half_three_corner
> mirror_orbit_fold
> even_half_two_corner
> full_IE_third_difference
> root_packet_shell
> OCF_cycle_space
> LRC_address_analogy.
```

This ranking is not a theorem about mathematical importance.  It is a proof
workflow reminder: establish the odd half-carrier geometry before importing
cycle-space or LRC analogies.

## Assumption Challenge

Do not assume tournament recursion vertices are runners, arcs, or even all
free tiles.  Alternate vertices considered:

- reflection orbits of free cells;
- fixed-line cells;
- half-tiling corner slots;
- interval roots;
- OCF odd-cycle packets;
- even-graph parity roots;
- LRC sector-state address slots;
- proof obligations.

The half-tiling quotient preserves complement/mirror structure and fixed-path
cell geometry.  It destroys the discarded side's independent bit data, so it is
not a full tournament encoding.  Its use is as a fundamental-domain address,
not as a replacement for full tiling enumeration.

## Status

THM-550 proves the count and recurrence layer.  The OCF/even-graph/LRC transfers
are open and should be tested by exact small-`n` scripts before being promoted
to canon.
