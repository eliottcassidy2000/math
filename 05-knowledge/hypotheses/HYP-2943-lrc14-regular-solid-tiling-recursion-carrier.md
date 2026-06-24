---
id: HYP-2943
title: LRC14 regular-solid and Euclidean-tiling recursion carriers
status: PROOF-INTERFACE / curvature and recursion guardrail; not a proof
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
  - HYP-2891
  - HYP-2887
  - HYP-2908
  - THM-572
  - OPEN-Q-108
results:
  - 04-computation/lrc14_platonic_tiling_recursion_codex_s141.py
  - 05-knowledge/results/lrc14_platonic_tiling_recursion_codex_s141.out
---

# HYP-2943: Regular-solid and Euclidean-tiling recursion carriers

S141 imports the prompt's Platonic, Archimedean, Johnson, square, triangular,
and hexagonal tiling ideas as POKE carriers for LRC14.  The point is to keep
the binary relations and recursion indices labelled before using the geometry.

The diagnostic script
`04-computation/lrc14_platonic_tiling_recursion_codex_s141.py` stores output
at `05-knowledge/results/lrc14_platonic_tiling_recursion_codex_s141.out`.

## Curvature boundary

The regular map notation `{p,q}` separates three regimes:

```text
(p-2)(q-2) < 4: spherical finite maps, the Platonic solids.
(p-2)(q-2) = 4: Euclidean zero-curvature tilings.
(p-2)(q-2) > 4: hyperbolic maps.
```

For the five Platonic solids the script records:

```text
tetrahedron  {3,3}, V,E,F=(4,6,4),   self-dual
cube         {4,3}, V,E,F=(8,12,6),  dual octahedron
octahedron   {3,4}, V,E,F=(6,12,8),  dual cube
dodecahedron {5,3}, V,E,F=(20,30,12), dual icosahedron
icosahedron  {3,5}, V,E,F=(12,30,20), dual dodecahedron
```

The plane regular tilings sit exactly on the zero-curvature wall:

```text
triangular tiling {3,6}, dual hexagonal tiling
square tiling     {4,4}, self-dual
hexagonal tiling  {6,3}, dual triangular tiling
```

Proof-use rule: Platonic solids are finite positive-curvature dual skeletons.
The Euclidean tilings are the recursion wall.  They should not be collapsed
into one "regular geometry" scalar.

## Euclidean recursion labels

The prompt's recursion counts are real, but they live in different carriers.

Square tiling:

```text
Gaussian axis scaling m+0i.
indices m^2 for m=2..5: 4,9,16,25.
```

This is a self-dual square-grid carrier.  It is useful only if the two grid
axes remain labelled, matching the general repo warning that square counts
alone are unsafe.

Triangular tiling:

```text
Eisenstein scaling m+0*omega.
general indices m^2 for m=2..5: 4,9,16,25.
binary chain emphasized by prompt: 4,16,64,256.
```

The clean LRC carrier is the dyadic triangular spine, because it retains
orientation-sector labels.  It is closer to the affine-depth packet in
HYP-2941 than to a raw face-count identity.

Triangle-hexagon bridge:

```text
local index 6.
a regular hexagon dissects into 6 equilateral triangles.
the triangular tiling has 6 triangles around a vertex.
```

This is the correct reading of "triangles and hexagons recurse with each other
in terms of 6."  It belongs to support-six packet language, compatible with
the HYP-2887 octahedral/current carrier, but it is not a replacement for exact
`q`, C27, or K33 labels.

Hexagonal tiling:

```text
Eisenstein norm N(a+b*omega)=a^2-a*b+b^2.
N(3+omega)=7.
norm-7 powers: 7,49,343,2401.
centered hex ring counts: 7,19,37,61.
```

This is the main guardrail.  The prompt's `7,49,...` chain is a norm-index
self-recursion, not the same thing as centered hexagonal patch counts.  The two
carriers agree at `7` and then diverge.

## Archimedean and Johnson roles

The 13 Archimedean solids retain one vertex-figure word.  The script records
their vertex figures and deficits, with range:

```text
deficit range = 1/60 .. 1/6.
```

They are useful analogues for labelled local quotients: one local word is
preserved globally.

Johnson solids supply the opposite guardrail:

```text
count = 92
regular polygon faces
convex
not Platonic/Archimedean/prism/antiprism
mixed vertex figures
```

For LRC14 this is the finite-residual-atlas role.  After uniform carriers fail,
one should expect a finite Johnson-like catalogue of bounded residual atoms,
not a single symmetric quotient.  This matches the AP/GW/petal/K33 frontier
tables and the S140 rule that the q=3 unital is branch-local rather than
global.

## Tournament Analysis

Tournament vertices are proof carriers:

```text
exact_M_Farey_branch
C27_unital_block_lift
Euclidean_tiling_recursion_indices
Archimedean_vertex_figure_words
Platonic_dual_curvature_skeleton
Johnson_residual_atlas
raw_polyhedron_visual_metaphor
```

Pair observable:

```text
LRC faithfulness,
label retention,
recursion exactness,
finite certifiability,
anti-scalar guard.
```

The role order is transitive:

```text
exact_M_Farey_branch
> C27_unital_block_lift
> Euclidean_tiling_recursion_indices
> Archimedean_vertex_figure_words
> Platonic_dual_curvature_skeleton
> Johnson_residual_atlas
> raw_polyhedron_visual_metaphor.
```

Fingerprint:

```text
c3=0
hp=1
```

## POKE proof target

The actionable target is a recursion-label test for future low-gap residuals:

```text
After exact M/Farey and C27/unital-block labels are attached, assign each
residual packet a Euclidean recursion type:
  square self-dual axis packet,
  triangular dyadic sector packet,
  triangle-hex support-six bridge,
  hex norm-7 self packet,
  or Johnson-like finite residual atom.
```

Then check whether unit-visible C27 defects always land in the triangular or
support-six discharge lanes, while nonunit K33 defects land in hex norm-7 or
Johnson-like residual lanes that can feed HYP-2908/THM-572.

This is not a proof of LRC14.  It is a structured POKE atlas: do not use
Platonic/Archimedean/Johnson/tiling language unless the preserved local word,
dual pair, or recursion index is named.
