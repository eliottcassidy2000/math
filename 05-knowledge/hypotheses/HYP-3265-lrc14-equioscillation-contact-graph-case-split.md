---
id: HYP-3265
title: LRC14 equioscillation contact graph and case split
status: SYNTHESIS / exact contact-graph scout; not an LRC14 proof
source: kind-pasteur-2026-06-27-S255 plus codex-2026-06-28 contact-graph addendum
tangent: T1355
technique: LTI-355
tournament_technique: LTT-255
script: 04-computation/lrc14_equioscillation_contact_graph_codex_20260628.py
result: 05-knowledge/results/lrc14_equioscillation_contact_graph_codex_20260628.out
reflection: 07-reflections/lrc14-equioscillation-contact-graph-case-split-codex-20260628.md
related:
  - HYP-3300
  - HYP-3260
  - HYP-3259
  - HYP-3258
  - HYP-3257
  - HYP-3256
  - HYP-3255
  - HYP-3253
  - HYP-3252
  - HYP-3251
  - HYP-3250
  - HYP-3248
  - HYP-3246
  - HYP-3245
  - HYP-3243
  - HYP-3242
  - HYP-3241
  - HYP-3240
  - HYP-3238
  - HYP-3237
  - HYP-3236
  - HYP-3218
  - HYP-3214
  - HYP-3132
  - HYP-2928
  - HYP-2909
  - THM-523
  - THM-530
  - OPEN-Q-108
---

# HYP-3265: LRC14 Equioscillation Contact Graph And Case Split

## Claim

The AP safety function observation is not just an extremal-theory analogy.
It is a small exact contact graph that organizes the LRC14 case split.

For a 13-speed row `S`, write

```text
f_S(t) = min_{s in S} ||s t||,     M(S) = max_t f_S(t).
```

For the AP `{1,...,13}`, and also for the Goddyn-Wong tight row
`{1,...,11,13,24}`, the closed witnesses at threshold `1/14` are exactly

```text
t = a/14,  a in (Z/14)* = {1,3,5,9,11,13}.
```

They form three antipodal pairs

```text
(1,13), (3,11), (5,9)
```

under `t -> 1-t`.  At a unit `a`, the binding runner residues are

```text
s == +/- a^{-1} mod 14.
```

Thus the six contacts quotient to three complement-pair binders:

```text
{1,13} -> {1,13}
{3,11} -> {5,9}
{5,9}  -> {3,11}
```

This is the graph-theoretic form of the equioscillation picture: a
six-vertex unit contact set, a three-vertex complement-pair binder set, and
a `2-to-1` contact map that becomes a matching after antipodal quotient.

## Structural Case Split

The unit contacts give an elementary proof fragment.

At a unit `a`, multiplication by `a` permutes the residues modulo `14`.
Therefore

```text
f_S(a/14) < 1/14
iff some s*a == 0 mod 14
iff some s == 0 mod 14.
```

So:

```text
14-free row:
  no speed divisible by 14
  -> every unit a has f_S(a/14) >= 1/14
  -> M(S) >= 1/14.

tight row:
  the six unit contacts are global maxima, with no higher off-unit peak.

covering row:
  a speed divisible by 14 kills all six unit contacts
  -> the witness must move off the units into the covering/floor route.
```

This is THM-523 at `q=14` reread through the Chebyshev contact graph.  The
non-covering branch is the unit-contact branch; the covering branch is exactly
where the contact graph is destroyed and a separate sidecar, such as the
multi-far decorrelation floor or promoted `Phi_{14d}` witness, is required.

## Topology / Geometry / Graph Reading

The AP picture has three simultaneous meanings.

```text
Topology:
  the danger-cover complement has six zero-measure boundary holes;
  after antipodal quotient these are three indexed holes.

Geometry:
  the safety profile is a max-min saddle with six active contacts;
  Kolmogorov/Chebyshev local minimality should be a convex-hull condition
  on the active contact gradients.

Graph:
  the contact graph is K_{1,2} + K_{1,2} + K_{1,2};
  after quotient by t -> 1-t it is a perfect matching between three
  unit-pair vertices and three binder-pair vertices.
```

This suggests a sharper proof target than "show the AP is special":

```text
classify rows whose unit-contact graph survives globally.
```

If the graph survives and no off-unit peak exceeds `1/14`, the row is in the
tight AP/Goddyn-Wong core.  If the graph is killed by a `0 mod 14` speed, the
proof must name the off-unit witness route.  If the graph survives but an
off-unit peak opens, the row is already strict.  That converts the visual
picture into a finite chamber classifier.

## Assumption Challenge

The tournament vertices should not be runners here.  The useful vertices are:

```text
unit contacts, antipodal contact pairs, complement-pair binders,
danger-cover boundary holes, Kolmogorov active-gradient obligations,
covering kill-switches, bulk/core gluing obligations, and proof carriers.
```

The quotient preserves the LRC predicate only as:

```text
unit-contact witness OR known tight equality core OR off-unit covering floor.
```

It destroys individual endpoint owners and the exact off-unit profile unless
those are retained as sidecars.

## Tournament Analysis

The scout uses proof carriers as vertices.

Pairwise observable: retained equioscillation proof payload minus destroyed
sidecars.

Switch/gauge: orient toward the carrier with larger weighted retained payload;
ties use fewer destroyed coordinates and then the more finite/formal carrier.

Exact fingerprint:

```text
vertices = 8
directed_3cycles = 0
singleton SCCs
Hamiltonian path:
  unit_contact_matching
  -> antipodal_three_pair_quotient
  -> d7_index_degree_packet
  -> danger_nerve_boundary_holes
  -> covering_kill_switch
  -> kolmogorov_convex_hull
  -> bulk_equidistribution_glue
  -> raw_safety_scalar
```

The graph route is therefore not a scalar ranking.  It says the contact
matching is currently the least lossy local carrier for the user's six-touch
observation, while the topology, index, Kolmogorov, covering, and bulk routes
are the sidecars needed to turn the local carrier into a proof.

## Incoming Mainline Integration

Rebasing over HYP-3250/HYP-3251/HYP-3252 sharpens the role of this contact
graph.  HYP-3250's finite-and-isolated tight-locus evidence is exactly the
chamber classifier target: AP/Goddyn-Wong/dilations should be the equality
cells, and everything else should carry a uniform positive margin.  HYP-3251
and HYP-3252 supply the guardrail: the Borsuk-Ulam/Gauss-sum index is an
ambient description of the AP equioscillation saddle, and the floor remains
the `S`-dependent proof content.  The contact graph should therefore be glued
to the covering floor, decorrelation margin, or finite chamber discharge
whenever the six unit contacts are killed or stop being global.

The later mac-mini S81/S82 tightening entries strengthen the same split.  Unit
witnesses for AP/Goddyn-Wong/dilations are now treated as the rigorous
construction branch; bounded single-swap/two-swap evidence supports the margin
branch; and S82 reframes resonant multiples of `14` as an on-grid core hit
whose witness moves to off-grid bulk.  HYP-3256 then separates layers:
Q(sqrt(-7)) organizes residues and the three binding pairs, while the
tight-versus-loose census is magnitude-level.  Thus the contact graph is the
finite equioscillation system to solve on the equality side, but non-global or
killed contacts still need an `S`-dependent floor, off-grid bulk, or finite
discharge.

The later mac-mini S83 hidden-`C_3` note refines the graph symmetry: the three
binding complement-pair slots are one orbit under multiplication by `3` modulo
`14`, i.e. the `(Z/14)^*/{+/-1}` / real-cyclotomic Galois orbit.  A structural
proof should therefore try to prove one binding-pair local lemma and transport
it equivariantly across the three-slot contact graph.

HYP-3259 adds the geometric reading: the tight locus is a real manifold whose
unit/binding coordinates are infinitesimally rigid while covering coordinates
flex.  That is exactly the contact-graph split in differential form.
HYP-3260 supplies the complementary nullspace warning: the same contact graph
has only rank `3` over unit residues, so nonunit residue, height, and covering
data must remain explicit sidecars.
HYP-3300 packages that warning as an observability/Morse requirement: the
contact graph can be a boundary column, not a scalar certificate that forgets
the blind layer.

## Next Pull

Formalize the contact-graph chamber theorem:

```text
primitive row
-> surviving unit-contact graph with global maximum 1/14 (AP/GW core)
   or surviving unit-contact graph with higher off-unit peak (strict open)
   or killed unit contacts with promoted Phi_{14d}/covering-floor witness
   or finite chamber discharge / named residual debt.
```

This is the contact-graph refinement of HYP-3243's topology/geometry/graph
atlas and HYP-3246's Chebyshev frame.
