---
id: HYP-2887
status: COMPUTATIONAL SIGNAL / proof-target
source: codex-2026-06-22-S103
tags: [lrc14, coimage, packet-graph, octahedron, line-graph, k4, hodge, reciprocal-tail, wall-ledger, abel-summation]
related:
  - HYP-2885
  - HYP-2886
  - HYP-2884
  - HYP-2883
  - HYP-2636
  - HYP-2633
  - HYP-2632
  - HYP-2828
  - HYP-2605
  - OPEN-Q-108
---

# HYP-2887: the repeated-packet lift is an octahedral current-realizability problem

The S101/S102 packet graph should not be forced into an arbitrary tournament
orientation.  Its hidden finite realizability structure is more rigid:

```text
nonzero packet graph = K6 minus the affine-zero perfect matching
                     = octahedron graph
                     = line graph L(K4).
```

Thus the six residue packet vertices are naturally the six edges of a
tetrahedron.  The affine-zero pairs are opposite `K4` edges.

One convenient address is:

```text
0 -> 12
2 -> 34
3 -> 13
6 -> 24
4 -> 14
5 -> 23
```

Then `edge_weight(a,b) != 0` exactly when the corresponding `K4` edges share a
vertex.  The missing affine lane

```text
(0,2), (3,6), (4,5)
```

is exactly the opposite-edge matching.

## Computation

Script:

- `04-computation/lrc14_octahedral_current_realizability_codex_s103.py`
- output: `05-knowledge/results/lrc14_octahedral_current_realizability_codex_s103.out`

The script verifies the `L(K4)` incidence and enumerates all `3^6 = 729`
integer layer cochains

```text
h : V(L(K4)) -> {0,1,2}
```

for the fixed core `(1,8,15,22)`.  For each gauge it computes the actual
support-six reciprocal lift through `H=10`, the local divergence

```text
div_H(v) = loop_H(v) + sum_{u~v} edge_H(u,v),
```

and several realizability invariants: layer sum, diameter, zero/opposite-edge
stretch, nonzero-edge stretch, weighted graph Laplacian, low-height wall
incidence, and divergence sign profile.

Only `126` unique lifted packet supports are needed for all `729` gauges.

## Findings

The `L(K4)` identification has no incidence mismatches.  The octahedron has
eight triangular faces: four `K4` vertex-stars and four `K4` faces.  Its graph
cycle rank is

```text
12 - 6 + 1 = 7.
```

Equivalently, the eight triangular face curls have one global dependence.  This
is a striking match with the seven-sector `14 = 2*7` quotient: the repeated
packet graph has a seven-dimensional face-curl module before analytic lifting.

The best scanned layer gauge at `H=10` is:

```text
layer word 210210
L1 divergence     0.00225361
max divergence    0.00072254
net divergence   -0.00204067
sign profile      1 positive / 5 negative vertices
wall incidence    1863
wall max          69
```

For comparison, the constant gauges are:

```text
000000: L1 divergence 0.0219283, wall max 86
111111: L1 divergence 0.00754806, wall max 68
222222: L1 divergence 0.00954186, wall max 64
```

So the improvement is not just "raise every lift."  A nonconstant realizable
height cochain cuts the `H=10` divergence by about a factor `9.7` against the
start-aligned gauge and by about a factor `3.35` against the raised-pair gauge.

The rough correlations across all `729` gauges are also suggestive:

```text
wall_max vs L1 divergence               +0.467
wall_incidence vs L1 divergence          +0.416
negative_vertices vs L1 divergence       -0.442
wall_max within fixed layer-sum          +0.451
```

The sign-profile correlation is the "opposite bounded signs" signal: gauges
that create mixed local divergence signs are easier to cancel than gauges with
all-positive leakage.

## Proof target

The next lemma should be stated as an octahedral Hodge estimate, not as an
abstract graph or tournament claim.

For a realizable layer cochain `h`, let `J_h` be the lifted reciprocal packet
current on the octahedron `L(K4)`.  After deleting the finite low-height wall
ledger, control:

```text
1. vertex divergence       delta J_h,
2. triangular face curls   d J_h,
3. additive-frequency tails from HYP-2636.
```

Because the octahedron is a sphere, there is no independent harmonic
one-current:

```text
C^1(L(K4)) = gradients + face-curls
```

up to the usual finite-dimensional gauge choice.  Therefore, after coherent
face curls are charged to low-height walls, the remaining spread current should
have no hidden global cycle obstruction.  It should be Abel-summable in the
HYP-2636 additive channels.

The theorem-shaped target is:

```text
realizable layer cochain on L(K4)
  + coefficient-height <=2 wall deletion
  + face-curl ledger on the 8 octahedral triangles
  + HYP-2636 channelwise Abel/Cauchy bound
  => support-six reciprocal divergence <= LRC(14) margin.
```

This refines HYP-2884.  HYP-2884 found the lifted divergence defect; HYP-2887
identifies the finite topology that a realizable defect must obey.

## Relation To Tournament Analysis

Tournament Analysis remains useful as a discipline, but the object here is not
a tournament on runners or residues.  It is a realizable current on a fixed
octahedral carrier.  The analogue of directed cyclic excess is the face-curl
around octahedral triangles.

This also clarifies the KPS winding-scar result: non-lonely phases have
cyclic/no-sink excess.  In the repeated-packet reciprocal tail, that cyclic
excess should appear as coherent octahedral face curl.  Coherent face curl is
finite wall structure; incoherent face curl is spread current and should be
handled by Abel summation.

## Assumption challenge

Candidate vertices considered: runners, residues, finite packets, packet
edges, octahedral faces, `K4` edges, integer lift layers, low-height walls,
additive-frequency channels, and proof obligations.

Chosen quotient: residue packets as the six edges of `K4`, with integer lifts
as height cochains on the octahedron `L(K4)`.

This preserves the finite local packet balance, the affine-zero matching, the
integer-lift realizability constraint, and the cycle-rank-seven face-curl
module.  It destroys raw runner labels and arbitrary orientation data.

The challenged assumption is that the missing LRC structure should be another
tournament.  The packet graph says otherwise: the repeated-residue obstruction
is a line-graph current with a spherical Hodge decomposition.
