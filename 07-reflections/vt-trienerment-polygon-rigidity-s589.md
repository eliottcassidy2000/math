---
source: codex-2026-06-03-S589
status: vertex-transitive trienerment rigidity synthesis
tags: [tournaments, LRC, trienerment, vertex-transitive, regular-polygon, rigidity, HYP-2125]
---

# Vertex-transitive trienerments and where polygon rigidity actually lives

The tempting sentence was:

```text
vertex-transitive trienerment <=> regular polygon point-set.
```

I think the sentence is almost right, but only after inserting the word
`cyclic`.  The regular polygon is the primitive cyclic case.  Once we allow
dihedral symmetry or abstract tournament automorphisms, the same local
isomorphism can propagate through blocks or group relators instead of through
one perimeter cycle.

That is the rigidity seam.

After rebasing over Opus HYP-2124, this sits one layer above the AP unit-clock
result: HYP-2124 is the exact cyclic witness orbit, while S589 asks what remains
true once symmetry is dihedral or nonabelian.

## The three layers

There are three increasingly loose meanings of vertex-transitive.

```text
1. cyclic sharp transitivity
   one rotation cycles through all points
   -> regular polygon

2. dihedral transitivity
   rotations plus reflections act transitively
   -> regular polygon or imprimitive bracelet

3. abstract tournament vertex-transitivity
   some automorphism group acts transitively
   -> Cayley/root-star cascade, not necessarily cyclic
```

The first layer is the user's regular-polygon intuition.  The second layer is
where local fixed-point rigidity appears: a reflection can make every point
look locally the same while preserving a two-block offset.  The third layer is
where tournament symmetry becomes algebraic rather than planar.

## The first non-polygon object

S589 enumerated subsets of `Z/N`, `N<=18`, containing `0`.  It computed the
ambient dihedral stabilizer, local circular distance profiles, gap words, and
whether the stabilizer was transitive.

The first nonregular vertex-transitive point-set is tiny:

```text
N = 6
P = (0, 1, 3, 4)
gaps = (1, 2, 1, 2)
```

Every vertex has the same local trienerment view.  But the set is not a
regular 4-gon in `Z/6`; it is a two-block bracelet:

```text
short, long, short, long.
```

Through `N<=18`, every imprimitive dihedral-VT example had gap period `2`.
That matches the classification: a transitive dihedral action has either one
rotation orbit or two rotation orbits swapped by a reflection.

So local fixed-point rigidity does appear, but it appears as a block system,
not necessarily as a regular polygon of points.

## The tournament warning

The Cayley audit compares cyclic and nonabelian vertex-transitivity.

For `C21`, the group contains elements of order `21`, so the Cayley
tournament has a cyclic polygon spine.

For the Frobenius group `F21 = C7 semidirect C3`, the Cayley tournament is
still vertex-transitive and regular, but the element orders are only
`1,3,7`.  There is no element of order `21`.  The local root-star still
determines the whole tournament by left translation, but the cascade runs
through the relation

```text
b a b^-1 = a^2
```

instead of around a single polygon.

This is exactly the old THM-052 caution in a new costume: vertex-transitive is
not the same as circulant.  The regular polygon is a very rigid subcase of
vertex-transitive symmetry, not the whole category.

## LRC reading

The observer-blind trienerment frame says: take

```text
P = {0, v_1, ..., v_{n-1}}.
```

At time `t`, build the danger graph; a point is lonely iff it is isolated.
This makes every basepoint legitimate.  HYP-2121 then says the AP interval is
least folded at the endpoints and most folded near the centre.  The observer
is not metaphysically special; it is the geometry's hardest viewpoint.

S589 adds a warning.  Erasing the observer is safe only if we keep the
cascade law that brings the local view back to the observer.  Otherwise we
confuse:

```text
same local distance profile
same unmarked tournament class
same score sequence
```

with:

```text
same observer-source threshold payload.
```

Those are not the same.  The useful object is local profile plus propagation
law plus labels.  That is also exactly the pincer-calculus message: contact
without grip labels is not a certificate.

## The useful repair

The repaired slogan is:

```text
cyclic primitive VT trienerment <=> regular polygon.
general VT trienerment <=> local rooted view + cascade law.
```

This lens makes the n=14 story more precise.  The centre of the AP folds to
the prime-7 layer.  The observer is the extreme, least-folded layer.  The open
creative handle is not just "regular polygon symmetry"; it is an
inward-to-outward rigidity lift:

```text
prime centre certificate
  -> one 2-adic/block unfolding
  -> observer/extreme certificate.
```

If that lift can be made label-preserving for denominator shields, endpoints,
and source-threshold states, then the regular-polygon intuition becomes a
real proof route rather than an analogy.
