---
id: HYP-2948
title: LRC14 Borel-Baire-Haar any-angle carrier
status: PROOF-INTERFACE / category-measure boundary carrier, not a proof
source: codex-2026-06-24-S145
related:
  - HYP-2946
  - HYP-2947
  - HYP-2945
  - HYP-2944
  - HYP-2943
  - HYP-2942
  - HYP-2940
  - HYP-2937
  - HYP-2932
  - HYP-2920
  - HYP-2908
  - THM-572
  - OPEN-Q-108
---

# HYP-2948: LRC14 Borel-Baire-Haar Any-Angle Carrier

The prompt's Borel/Baire/Haar and any-angle path-planning ideas sharpen a
useful distinction at the LRC14 threshold.

For a finite runner set `S`, the strict safe set

```text
{t in T : ||v t|| > 1/14 for all v in S}
```

is an open Borel/Baire subset of the circle.  The threshold-safe set with
`>= 1/14` is closed Borel/Baire.  Haar's theorem supplies the canonical
normalized translation-invariant measure on the circle, and on any finite
orbit-closure subtorus cut out by integer relations.  In the compact metric
settings used here, the finite arc events are both Borel and Baire; the
important split is measure versus category/interior.

## Computation

The script

```text
04-computation/lrc14_borel_baire_haar_anyangle_codex_s145.py
```

stores output at

```text
05-knowledge/results/lrc14_borel_baire_haar_anyangle_codex_s145.out
```

At threshold `delta=1/14`, the exact rational audit gives:

```text
AP                         strict Haar 0       closed pts 6
GW 12->24                  strict Haar 0       closed pts 6
near 12->36                strict Haar 1/1260  components 2
petal 10->20               strict Haar 1/980   components 2
petal 13->26               strict Haar 1/182   components 2
two-swap 10,12->20,24      strict Haar 1/980   components 2
two-swap 10,12->20,36      strict Haar 4/2205  components 4
```

Thus AP and Goddyn-Wong are boundary-only at the threshold: they have no strict
Haar mass and no Baire-open interval, but they have closed threshold support at
six denominator-14 unit times.  The near-miss and petal rows already have
positive strict Haar mass and hence nonempty Baire interior.

## Haar-Baire Wave*

The proposed sixth any-angle proof carrier is:

```text
Haar-Baire Wave*
```

It propagates Borel interval fronts on the circle or on a relation-lattice
subtorus, labelled by

```text
(strict Haar mass, Baire interior, closed boundary support).
```

The path-planning dictionary is:

```text
line of sight  -> no unsafe arc blocks a witness interval
taut path      -> every heading change is a cover-arc boundary event
wavefront      -> exact denominator combs and wall crossings
Haar label     -> invariant mass on the orbit-closure torus
Baire label    -> open/nonmeager interval versus meager boundary support
```

This imports the useful parts of Field D*, Theta*, AP Theta*, Lazy Theta*,
Block A*, ANYA, and CWave without pretending a grid path algorithm is itself a
number-theoretic proof.

## Tournament Analysis

Tournament vertices are proof carriers, not runners:

```text
Haar-Baire Wave*
ANYA interval-taut
CWave primitive-front
AP Theta* angle-prop
Theta* line-of-sight
Field D* interpolation
Lazy Theta*
Block A*
```

The pairwise observable is the criterion vector

```text
exactness, interval_nodes, continuous_geometry, dynamic_update,
lrc_transfer, haar_baire_label, anti_scalar_guard.
```

The switch/gauge orients `A -> B` when `A` wins more coordinates than `B`;
ties follow the displayed Hamiltonian path.  The result is transitive:

```text
Haar-Baire Wave*
> ANYA interval-taut
> CWave primitive-front
> AP Theta* angle-prop
> Field D* interpolation
> Theta* line-of-sight
> Lazy Theta*
> Block A*.
```

Fingerprint:

```text
score_histogram = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
directed_3cycles = 0
hamiltonian_paths = 1
```

## Assumption Challenge

For this LRC route, tournament vertices need not be runners or arcs.  Alternate
vertex sets considered:

```text
runners
unsafe arcs
strict-safe interval components
closed boundary points
wall-crossing events
subtorus relation cells
proof obligations
```

The useful quotient for S145 uses interval-front/proof-obligation vertices.  It
preserves the predicate "strict open witness versus boundary-only threshold
support" and the measured Haar mass of the witness region.  It destroys runner
ownership, individual active-pair arithmetic, and most C27 shell labels unless
they are reattached as packet metadata.  The challenged assumption is that a
lonely-runner tournament should start from runner vertices; here the important
frontier is the cover boundary itself.

## Proof Targets

1. **Boundary-support lemma.**  Classify threshold-safe but strict-Haar-zero
   rows by exact denominator boundary points.  The current reduction predicts
   that AP and GW are the only q=14 survivors.
2. **Haar subtorus lift.**  For relation-lattice residuals, push unsafe arc
   covers to the orbit closure and integrate by normalized Haar measure there,
   not by an independent-block product measure.
3. **Baire-open strictness.**  Show every non-AP/GW reduced residual has a
   nonempty open safe interval.  In finite arc arrangements, this is the
   category version of positive strict Haar mass.
4. **Any-angle finite frontier.**  Use ANYA/CWave style interval fronts as
   nodes: expand only wall-crossing intervals where a cover arc becomes tight,
   with lazy exact-`M` checks.
5. **State-lift bridge.**  If Haar-Baire Wave* cannot produce an open
   interval, the remaining boundary-only packet should retain enough labels to
   feed HYP-2908/THM-572.

This is not a proof of LRC14.  It is a carrier split that says the AP/GW tight
locus is a closed-boundary phenomenon while the known low-gap perturbations are
already open/Haar-positive.

## Relation To HYP-2947

Incoming HYP-2947 claims the broader measurable rank recombination namespace:
exact `M`/Farey branch, C27 shell transfer, q=3 unital chart, affine depth,
K33/Kuratowski packet, measurable address tax, and HYP-2908/THM-572 state lift.
HYP-2948 is the boundary-front companion: it supplies the exact small-row
Haar/Baire readout saying when the measurable route has open interior and when
it must switch to boundary-owner labels.
