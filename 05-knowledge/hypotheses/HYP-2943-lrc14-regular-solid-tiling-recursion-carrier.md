---
id: HYP-2943
title: LRC14 regular-solid and Euclidean-tiling recursion carrier
status: PROOF-INTERFACE / annular 28-point 1/14 carrier plus curvature/norm guardrails; not a proof
source: codex-2026-06-24-S141
related:
  - HYP-2942
  - HYP-2941
  - HYP-2940
  - HYP-2939
  - HYP-2938
  - HYP-2937
  - HYP-2936
  - HYP-2935
  - HYP-2932
  - HYP-2894
  - HYP-2892
  - HYP-2891
  - HYP-2887
  - HYP-2908
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_platonic_tiling_recursion_codex_s141.py
  - 05-knowledge/results/lrc14_platonic_tiling_recursion_codex_s141.out
---

# HYP-2943: The 14-gonal prism/antiprism is the cyclic companion to the q=3 unital

S141 tests the user's Platonic/Archimedean/Johnson and
square/triangle/hexagonal tiling prompt as an LRC14 carrier audit.  The
synthesis of the two parallel S141 attempts is:

```text
regular maps: curvature boundary and dual skeletons
Euclidean tilings: zero-curvature recursion/norm labels
Platonic/Archimedean solids: finite uniform positive-curvature seeds
Johnson solids: finite nonuniform residual atlases
14-gonal prism/antiprism: exact LRC14 annular carrier
```

The last item is the new high-leverage object.  The norm and solid data are
guardrails telling us which geometric labels are legitimate; the annulus is the
candidate proof carrier to test next.

## Computation

The script
`04-computation/lrc14_platonic_tiling_recursion_codex_s141.py` stores output at
`05-knowledge/results/lrc14_platonic_tiling_recursion_codex_s141.out`.

It uses normalized vertex defect

```text
kappa = 1 - sum_p (p-2)/(2p)
```

for regular-polygon vertex configurations.  For a uniform spherical
configuration,

```text
V * kappa = 2.
```

This keeps vertex count, curvature denominator, and face counts as exact
Fractions.

## Curvature Boundary

The regular map notation `{p,q}` separates the regimes:

```text
(p-2)(q-2) < 4: spherical finite maps, the Platonic solids.
(p-2)(q-2) = 4: Euclidean zero-curvature tilings.
(p-2)(q-2) > 4: hyperbolic maps.
```

The five Platonic solids are the positive-curvature finite regular maps:

```text
tetrahedron  {3,3}, V,E,F=(4,6,4),    self-dual
cube         {4,3}, V,E,F=(8,12,6),   dual octahedron
octahedron   {3,4}, V,E,F=(6,12,8),   dual cube
dodecahedron {5,3}, V,E,F=(20,30,12), dual icosahedron
icosahedron  {3,5}, V,E,F=(12,30,20), dual dodecahedron
```

The plane regular tilings sit exactly on the zero-curvature wall:

```text
triangular tiling {3,6}, dual hexagonal tiling
square tiling     {4,4}, self-dual
hexagonal tiling  {6,3}, dual triangular tiling
```

Proof-use rule: do not collapse finite positive-curvature skeletons and
repeatable Euclidean recursion packets into one scalar "regular geometry"
object.

## Euclidean Recursion Labels

The prompt's recursion counts are real, but they live in different carriers.

Square tiling:

```text
Gaussian axis scaling m+0i.
indices m^2: 4,9,16,25,36,49,...
```

This is a self-dual square-grid carrier and is useful only if the two grid axes
remain labelled.

Triangular tiling:

```text
Eisenstein scaling m+0*omega.
general indices m^2.
dyadic spine emphasized by prompt: 4,16,64,256,...
```

The clean LRC-style carrier is the dyadic triangular spine, because it retains
orientation-sector labels.

Triangle-hexagon bridge:

```text
local index 6.
a regular hexagon dissects into 6 equilateral triangles.
the triangular tiling has 6 triangles around a vertex.
```

This belongs to support-six packet language, compatible with the
HYP-2887 octahedral/current carrier, but it is not a replacement for exact
`q`, C27, or K33 labels.

Hexagonal tiling:

```text
Eisenstein norm N(a+b*omega)=a^2-a*b+b^2.
N(3+omega)=7.
norm-7 powers: 7,49,343,2401,...
centered hex ring counts: 7,19,37,61,...
```

Guardrail: the prompt's `7,49,...` chain is a norm-index self-recursion, not
the same as ordinary centered hexagonal patch counts.  The two carriers agree
at `7` and then diverge.

The `7` is not automatically the forbidden tournament value in THM-572.  It
becomes relevant only if a labelled LRC residual constructs the right
tournament state-lift packet.

## Platonic, Archimedean, Johnson

Platonic solids are finite positive-curvature seeds:

```text
tetrahedron     3^3       kappa=1/2   V=4
octahedron      3^4       kappa=1/3   V=6
cube            4^3       kappa=1/4   V=8
icosahedron     3^5       kappa=1/6   V=12
dodecahedron    5^3       kappa=1/10  V=20
```

The 13 Archimedean solids retain one vertex-figure word.  The script records
their exact vertex defects; the deficit range is:

```text
1/60 .. 1/6.
```

So Archimedean solids are good analogues for labelled local quotients: one
local word is preserved globally.

Johnson solids supply the opposite guardrail:

```text
count = 92
regular polygon faces
convex
mixed vertex figures
not Platonic/Archimedean/prism/antiprism
```

For LRC14 this is the finite residual-atlas role.  After uniform carriers fail,
one should expect Johnson-like bounded catalogues of residual atoms, analogous
to AP/GW/petal/K33 frontier tables, not a single symmetric quotient.

## Prism / Antiprism Annular Families

The decisive computation is the infinite Archimedean prism/antiprism family:

```text
n-gonal prism:      vertex configuration (4,4,n)
n-gonal antiprism:  vertex configuration (3,3,3,n)
```

Both satisfy:

```text
kappa = 1/n
V = 2n.
```

Therefore:

```text
n=7:
  V=14, kappa=1/7

n=14:
  V=28, kappa=1/14
```

This is the live bridge.  HYP-2942 gave a 28-point q=3 unital that preserves
pair incidence but not cyclic order.  The 14-gonal prism/antiprism gives a
28-vertex object that preserves cyclic annular order and has exactly the LRC14
threshold curvature `1/14`.

```text
q=3 unital:
  28 points, pair-unique incidence, detects repeated H12 pairs

14-gonal prism/antiprism:
  28 vertices, two 14-cycles, per-vertex curvature 1/14, cyclic annular order
```

So the next test is not "which object is correct?" but:

```text
combine unital pair incidence with annular cyclic order.
```

## Tournament Analysis

S141 uses proof carriers as vertices, not raw solids.  The observable is:

```text
exact-LRC retention,
C27/branch retention,
pair-incidence retention,
cyclic/defect retention,
anti-scalar guard.
```

The carrier tournament is transitive:

```text
exact M/Farey branch
> C27 shell transfer
> q=3 unital pair incidence
> 14-prism/antiprism annulus
> Euclidean tiling norm recursions
> triangle/hex support-six bridge
> square Gaussian self recursion
> Platonic/Archimedean uniform atlas
> hex norm-7 recursion
> Johnson finite defect atlas
> raw solid numerology
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1}
c3=0
hp=1
```

## POKE Proof Target

Build a 28-vertex annular label model with two 14-cycles.  Attach

```text
AP, GW, H1..H13, D1..D13
```

and compare it to the HYP-2942 unital chart:

```text
GW  H12->D3
K33 H12->D9
```

Questions:

```text
Does the H12 conflict become a twist?
Does it become a diameter?
Does it force a two-chart obstruction?
Does annular cyclic order plus q=3 pair incidence produce the missing
state-lift packet for HYP-2908/THM-572?
```

This is not a proof of LRC14.  It is a sharper, typed POKE target: after exact
`M`/Farey and C27/unital labels are attached, use the 14-annulus to test cyclic
order constraints, and use Euclidean norm labels only as residual packet
classifiers.
