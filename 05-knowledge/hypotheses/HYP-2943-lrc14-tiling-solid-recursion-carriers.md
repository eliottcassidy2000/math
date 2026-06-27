---
id: HYP-2943
title: LRC14 tiling and solid recursion carriers split by curvature and lattice norm
status: COMPUTATIONAL SCOUT / proof-interface carrier split; not a proof
source: codex-2026-06-24
tags: [lrc14, tilings, platonic, archimedean, johnson, eisenstein, gaussian, c27]
related:
  - HYP-2942
  - HYP-2939
  - HYP-2937
  - HYP-2936
  - HYP-2932
  - HYP-2920
  - HYP-2908
  - THM-554
  - THM-555
results:
  - 04-computation/lrc14_tiling_solid_recursion_codex.py
  - 05-knowledge/results/lrc14_tiling_solid_recursion_codex.out
---

# HYP-2943: tiling and solid recursion carriers for LRC14

This scout separates three carrier types that tend to get blended in the LRC14
analogy stack:

```text
local Euclidean incidence,
lattice self-similarity index,
finite spherical curvature defect.
```

The computation is in
[lrc14_tiling_solid_recursion_codex.py](/home/bigo/math/04-computation/lrc14_tiling_solid_recursion_codex.py:1),
with output at
[lrc14_tiling_solid_recursion_codex.out](/home/bigo/math/05-knowledge/results/lrc14_tiling_solid_recursion_codex.out:1).

## Regular Map Split

For a regular `{p,q}` map, the curvature separator is:

```text
delta = 2/p + 2/q - 1.
```

The five Platonic solids are the positive-curvature cases:

```text
tetrahedron   {3,3}  (V,E,F)=(4,6,4)
cube          {4,3}  (V,E,F)=(8,12,6)
octahedron    {3,4}  (V,E,F)=(6,12,8)
dodecahedron  {5,3}  (V,E,F)=(20,30,12)
icosahedron   {3,5}  (V,E,F)=(12,30,20)
```

This matches the existing "five Platonic tournaments" reflection:
tournament object, tournament cube, arc-coordinate space, quotient surface,
and quotient-dual surface.

The Euclidean zero-curvature cases are exactly:

```text
triangular {3,6}
square     {4,4}
hexagonal  {6,3}
```

So the square tiling is self-dual, while triangular and hexagonal tilings are a
dual pair.

## Recursion Numbers

The user numbers split cleanly by type.

Square self-recursion:

```text
1,4,9,16,25,... = m^2
```

Triangle self-recursion:

```text
1,4,9,16,... = integer-scale triangular subdivision
```

Triangle/hexagon recursion:

```text
6 = local dual-incidence count in {3,6} <-> {6,3}
```

This is not a lattice index.  The atlas checks that `6` is absent from both
the Gaussian and Eisenstein norm lists up to `100`.

Hexagon self-recursion has two distinct 7s:

```text
centered flower count: 1 + 6 = 7
Eisenstein index:      7 = N(3 + omega)
```

Repeated hex self-similarity gives:

```text
49 = N((3 + omega)^2) = N(8 + 5*omega).
```

Thus the clean typed reading is:

```text
4,16   square/triangle integer-scale self-recursion
6      triangle-hex local duality
7,49   hex/Eisenstein self-similarity
```

## Norm Atlas Guardrail

The script computes Gaussian and Eisenstein norms through `100`.

Key targets:

```text
4:  Gaussian yes, Eisenstein yes
6:  Gaussian no,  Eisenstein no
7:  Gaussian no,  Eisenstein yes
14: Gaussian no,  Eisenstein no
16: Gaussian yes, Eisenstein yes
22: Gaussian no,  Eisenstein no
27: Gaussian no,  Eisenstein yes
31: Gaussian no,  Eisenstein yes
49: Gaussian yes, Eisenstein yes
```

This matters for the petal material.  `22/7` can remain an angular alias and
`cube_root(31)` can remain a petal scale alias, but neither is automatically an
LRC proof unit.  The atlas says:

```text
22 is not a Gaussian/Eisenstein index;
31 is an Eisenstein index, but still only a typed scale carrier until the
AP/Goddyn-Wong/C27 labels are attached.
```

## Archimedean and Johnson Solids

The script also recomputes the thirteen Archimedean f-vectors from their
vertex configurations by the exact per-vertex curvature:

```text
curv = 1 - d/2 + sum(1/p_i).
```

These are uniform finite defects between Platonic symmetry and arbitrary
regular-faced patches.  The 92 Johnson solids are the non-uniform finite
defects.  Their proof role is therefore not "universal symmetry"; it is
stress-testing whether an AP/Goddyn-Wong/C27 label leaks under finite
curvature perturbation.

## LRC14 Proof Interface

At `n=14`:

```text
apex half-size = 7
C = 2n - 1 = 27 = 3^3.
```

That makes the tri/hex/Eisenstein carrier more relevant to the current LRC14
branch than generic polyhedral symmetry:

```text
C27 shell labels
-> q=3 / F3^3 unital chart
-> Eisenstein triangular/hexagonal lattice carrier.
```

This does not override HYP-2942.  It sharpens it.  HYP-2942 says the q=3 unital
is branch-local because lambda=1 blocks cannot globally carry both the GW and
K33 `H12` branch packets.  HYP-2943 says the surrounding Euclidean geometry
should be typed the same way: use the tri/hex carrier only after C27 and
AP/Goddyn-Wong labels have already fixed the branch.

## Candidate Lemma

```text
After AP/Goddyn-Wong and C27 labels are attached, any low-gap LRC14 row that
uses the triangular/hexagonal carrier must preserve the branch-local C27 marked
transfer chart; otherwise it falls into the HYP-2942 unital pair-repeat
obstruction or the Farey/K33 child branch.
```

This is a proof-interface lemma, not a proof.  It proposes the correct order
of reductions:

```text
exact M/Farey node
> AP/GW C27 marked shell transfer
> q=3 unital branch-local block chart
> Eisenstein tri/hex norm carrier
> triangular-hexagonal local duality
> square/Gaussian self-dual carrier
> Platonic tournament-level carrier
> Archimedean finite defect
> Johnson finite defect
> raw numerical analogy.
```

## Tournament Analysis

Define a binary relation on the ten carriers:

```text
x -> y iff x retains more LRC proof data than y.
```

The resulting tournament is transitive:

```text
vertices = 10
edges = 45
c3 = 0
Hamiltonian paths = 1
```

The point is not the transitivity itself.  The point is the role ordering.  It
prevents a finite-solid analogy from being used before the exact LRC labels
that make the analogy meaningful.

## Guardrail

The most useful negative result is:

```text
6 is not the same kind of number as 7.
```

`6` is local tiling duality.  `7` is both the first centered hex flower and a
hex/Eisenstein self-similarity index.  LRC14 has both a `7` apex and a
`C27=3^3` shell chart, so the live carrier is the typed
`C27 -> q=3/F3^3 -> Eisenstein tri/hex` passage, not undifferentiated
polyhedral numerology.
