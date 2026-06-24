---
id: HYP-2943
title: LRC14 regular-solid and Euclidean tiling recursion carrier
status: PROOF-INTERFACE / annular 28-point 1/14 carrier identified; not a proof
source: codex-2026-06-24-S141
related:
  - HYP-2942
  - HYP-2941
  - HYP-2940
  - HYP-2938
  - HYP-2894
  - HYP-2892
  - HYP-2908
  - THM-572
results:
  - 04-computation/lrc14_platonic_tiling_recursion_codex_s141.py
  - 05-knowledge/results/lrc14_platonic_tiling_recursion_codex_s141.out
---

# HYP-2943: The 14-gonal prism/antiprism is the cyclic companion to the q=3 unital

S141 tests the user's prompt:

```text
Platonic, Archimedean, and Johnson solids;
Euclidean square, triangular, and hexagonal tilings;
square self-recursion 4,9,16,25;
triangle/hex recursion by 6;
hex self-recursion 7,49;
triangle self-recursion 4,16.
```

The conclusion is a split carrier atlas:

```text
Euclidean tilings: zero-defect recursion laws
Platonic/Archimedean solids: uniform positive-curvature seeds
Johnson solids: finite nonuniform defect packets
14-gonal prism/antiprism: exact LRC14 annular carrier
```

The last item is the new high-leverage result.

## Computation

The script
`04-computation/lrc14_platonic_tiling_recursion_codex_s141.py` stores output at
`05-knowledge/results/lrc14_platonic_tiling_recursion_codex_s141.out`.

It uses normalized vertex defect

```text
kappa = 1 - sum_p (p-2)/(2p)
```

for a regular-polygon vertex configuration.  For a uniform spherical
configuration, total curvature gives:

```text
V * kappa = 2.
```

This makes the vertex count, curvature denominator, and face counts exact
Fractions rather than visual guesses.

## Euclidean Regular Tilings

The three regular Euclidean tilings are exactly the zero-defect configurations:

```text
triangular tiling:  3^6     kappa = 0
square tiling:      4^4     kappa = 0
hexagonal tiling:   6^3     kappa = 0
```

The prompt's recursion counts split into labelled carriers:

```text
square self n^2:
  4, 9, 16, 25, 36, 49, ...

triangle self 4^k:
  4, 16, 64, 256, ...

triangle -> hex exchange 6*n^2:
  6, 24, 54, 96, 150, ...

hex centered patches 1 + 3r(r+1):
  7, 19, 37, 61, 91, 127, ...

hex self hexaflake 7^k:
  7, 49, 343, 2401, ...
```

Interpretation:

```text
square = self-dual product/grid recursion
triangle = dyadic self-subdivision
triangle/hex = dual exchange by six triangles per hexagon
hex = central cell plus six petals, either lattice patch or fractal packet
```

The `7` in the hex carrier is not by itself the forbidden tournament value
from THM-572.  It becomes relevant only if an LRC residual constructs the
right tournament state-lift packet.

## Platonic and Archimedean Solids

Platonic solids are uniform positive-curvature seeds:

```text
tetrahedron     3^3       kappa=1/2   V=4
octahedron      3^4       kappa=1/3   V=6
cube            4^3       kappa=1/4   V=8
icosahedron     3^5       kappa=1/6   V=12
dodecahedron    5^3       kappa=1/10  V=20
```

The finite Archimedean seeds add mixed regular faces.  They are useful
structural analogies, especially the old truncated-octahedron and
square/hex/triangular carriers, but none of the finite exceptional
Archimedean solids is the direct LRC14 scale.

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

This is the live bridge.  The q=3 unital from HYP-2942 also has `28` points,
but it preserves pair incidence rather than cyclic order.  The 14-gonal
prism/antiprism gives a complementary `28`-point annular object at exactly the
LRC14 threshold scale.

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

## Johnson-Solid Guardrail

Johnson solids are regular-faced but nonuniform.  The script checks small
diagnostic examples:

```text
square pyramid J1:
  apex      (3,3,3,3)   contribution 1/3
  base x4   (3,3,4)     contribution 5/3
  total = 2

pentagonal pyramid J2:
  apex      (3,3,3,3,3) contribution 1/6
  base x5   (3,3,5)     contribution 11/6
  total = 2
```

This is the right lesson: Johnson solids belong to finite residual-atlas work.
They are analogous to finite exceptional LRC packets after the uniform and
annular carriers fail, not a single global recursion law.

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
> Euclidean tiling zero-defect laws
> triangle/hex dual exchange 6
> square self n^2 recursion
> Platonic/Archimedean uniform atlas
> hex self 7^k recursion
> Johnson finite defect atlas
> raw solid numerology
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1}
c3=0
hp=1
```

## Proof-Route Consequence

HYP-2942 gave the pair-incidence side:

```text
q=3 unital detects that GW H12->D3 and K33 H12->D9 cannot globally share
the same raw H12 pair.
```

HYP-2943 adds the cyclic-annular side:

```text
the 14-gonal prism/antiprism has 28 vertices and local curvature 1/14,
so it is a candidate carrier for two 14-cycles, half-step twists, and
apex-clock order.
```

The next concrete experiment is:

```text
build a 28-vertex annular label model with two 14-cycles,
attach AP/GW/H/D labels,
and compare whether the H12/GW/K33 conflict seen by the unital becomes
a twist, a diameter, or a two-chart obstruction.
```

If successful, the pair-incidence unital and cyclic annulus could become two
orthogonal coordinates of the same C27/LRC14 state-lift packet.
