---
id: HYP-2125
status: OPEN synthesis with proved classification lemma in the cyclic/dihedral point-set model; S589 supplies finite audits and VT tournament counter-pressure
source: codex-2026-06-03-S589
related: [HYP-2134, HYP-2133, HYP-2132, HYP-2131, HYP-2130, HYP-2129, HYP-2128, HYP-2127, HYP-2126, HYP-2124, HYP-2123, HYP-2122, HYP-2121, HYP-2120, THM-401, THM-400, THM-389, THM-381, THM-052, HYP-2113, HYP-2109, HYP-2091, HYP-1977]
---

# HYP-2125: vertex-transitive trienerment rigidity is polygonal only in the cyclic primitive case

## Claim

The slogan

```text
vertex-transitive trienerment <=> regular polygon point-set
```

is correct only after adding the word **cyclic** or **primitive**.

There are three rigidity layers:

```text
cyclic sharp-transitive point-set
  -> regular polygon

dihedral vertex-transitive point-set
  -> regular polygon or imprimitive bracelet

abstract vertex-transitive tournament
  -> local root-star plus a group-law cascade
     (cyclic polygon if the group has an n-cycle,
      nonabelian mesh otherwise)
```

So the real invariant is not "regular polygon" by itself.  It is:

```text
local rooted view + the cascade law that propagates it globally.
```

For a regular polygon the cascade law is one cyclic rotation.  For an
imprimitive bracelet it is a block system plus reflections.  For a nonabelian
Cayley tournament it is a group presentation.  For LRC, it is the
observer/source threshold payload plus endpoint and denominator labels.

Relative to HYP-2124, this is the broader taxonomy layer: HYP-2124 proves the
exact AP/unit-clock cyclic witness orbit, while HYP-2125 asks what remains true
when the symmetry is dihedral, imprimitive, or nonabelian.

## Rigidity Lemma

Let `P` be a finite point-set in a cyclic ambient polygon `Z/N`, and let
`D_N` act by rotations and reflections.

If the **rotation subgroup** of the stabilizer of `P` is transitive on `P`,
then `P` is a regular polygon: it is one orbit of a rotation, hence all
circular gaps are equal.

If the full **dihedral** stabilizer is transitive but the rotation subgroup is
not, then the rotation subgroup has exactly two equal orbits on `P`, and a
reflection swaps them.  Thus `P` is an imprimitive bracelet: two regular
sub-polygons interlaced with an alternating two-gap word.

Reason: the rotation subgroup has normal index at most two in any dihedral
subgroup.  Therefore a transitive dihedral action has either one rotation
orbit or two rotation orbits.  The one-orbit case is cyclic.  In the two-orbit
case, reflections exchange the two orbits and preserve the rotation step,
leaving a block system rather than a primitive polygon.

This is where local rigidity around a fixed point appears: the stabilizer of a
point can force identical local distance/trienerment profiles, but if the
rotation subgroup has two orbits, that local rigidity does not collapse the
global block offset.

## Evidence From S589

S589 enumerates all subsets `P` of `Z/N` with `0 in P` for `3 <= N <= 18`.
It computes the ambient dihedral stabilizer, its orbits on `P`, circular gap
words, distance profiles from every point, and a Cayley-tournament audit.

Point-set audit:

```text
dihedral_vertex_transitive = 83
regular_polygons           = 31
imprimitive_bracelets      = 52
local_distance_profile_equal = 83
local_equal_but_not_dihedral_VT = 0
```

The first nonregular vertex-transitive trienerment is:

```text
N=6, P=(0,1,3,4), gaps=(1,2,1,2).
```

Every vertex has the same local distance/trienerment profile, but the set is
not a regular quadrilateral in the ambient circle.  It is a regular polygon
of two-point blocks, with gap period `2`.

All imprimitive examples through `N<=18` have gap period `2`, matching the
dihedral classification.  This supports the reading that the first loss of
rigidity is not chaos; it is a block system.

Abstract tournament audit:

```text
C21 interval circulant:
  vertex-transitive, has element of order 21, cyclic polygon spine

F21 alternating Cayley:
  vertex-transitive, regular score histogram, c3=385,
  but element orders are {1,3,7}; no element has order 21
```

Thus a vertex-transitive tournament can be globally rigid without being a
regular polygon in the cyclic sense.  The local star at the identity still
determines all arcs, but it propagates through the nonabelian relator mesh
`b a b^-1 = a^2`, not around one perimeter cycle.

This agrees with the repo's MISTAKE-013/MISTAKE-014 correction: non-circulant
vertex-transitive tournaments exist at order `21`, and cyclic/circulant
symmetry is a stricter condition than vertex-transitivity.

## LRC Reading

The symmetric trienerment frame treats

```text
P = {0, v_1, ..., v_{n-1}}
```

as an observer-blind geometry.  A point `p` is lonely at time `t` iff `p` is
isolated in the danger graph.  HYP-2121's basepoint-gradient work shows that
the AP interval is not vertex-transitive as a basepoint geometry: the extreme
observer is least folded and hardest, while the centre is folded and loose.

S589 says the creative route should not be "erase the observer and hope local
symmetry decides."  It should be:

```text
1. identify the local rooted/trienerment profile,
2. identify the cascade law that propagates it,
3. keep the block, relator, threshold, or endpoint labels that survive the
   propagation.
```

This connects directly to:

- HYP-2120: source-perspective recursion is the exact rooted slice;
- HYP-2121: rooted classes still need incident threshold payload;
- HYP-2122: denominator shields are rooted observer-coupled gates;
- HYP-2123: pincer certificates need grip labels and escape ledgers.

Regular polygon AP is the cyclic primitive boundary case.  `V*`, shifted AP,
and nonabelian VT tournaments warn that the next cases are not random; they
are still rigid, but the rigidity is block/relator/labelled rather than purely
polygonal.

## Proof Route

1. Formalize the dihedral point-set classification:
   transitive rotation subgroup gives regular polygon; transitive dihedral
   but nontransitive rotation subgroup gives two-orbit bracelets.
2. Extend the trienerment statement from static distance profiles to
   threshold movies: a dihedral action preserving `P` preserves the full
   danger/safe three-state movie, so vertex-transitive point-sets have
   isomorphic local trienerments at every time.
3. Prove a "local profile is not enough" warning or theorem.  S589 finds no
   local-distance-profile-equal but non-dihedral-VT examples through `N<=18`,
   but LRC needs the stronger threshold-labelled/rooted profile anyway.
4. For tournaments, split vertex-transitivity into:
   cyclic regular action (circulant/polygon),
   dihedral/converse augmentation,
   nonabelian Cayley relator mesh,
   and nonregular automorphism actions.
5. For LRC, route each symmetry type to the right proof object:
   cyclic polygon -> delta-clock/straddle pairs;
   bracelet/block -> basepoint folding and block lift;
   nonabelian mesh -> labelled quotient/automaton;
   observer-rooted source -> endpoint/denominator pincer.

The hoped-for n=14 payoff is an inward-outward rigidity theorem: centre
folding sees the prime-7 layer, while the extreme observer unfolds one
2-adic/block layer.  Proving that this unfolding preserves a labelled pincer
certificate would turn the "regular polygon" intuition into a usable LRC
descent.

## Tournament Analysis

S589 uses rigidity lenses as tournament vertices, not ordinary vertices.

The pair observable is:

```text
(global cascade, local rigidity, LRC payload,
 fixed-point handling, primitive polygonality)
```

The resulting transitive ranking is:

```text
cyclic_regular_polygon_spine
> nonabelian_Cayley_VT_mesh
> dihedral_imprimitive_bracelet
> source_root_threshold_payload
> local_distance_profile_only
> regular_score_sequence_only.
```

Score histogram is `{0:1,1:1,2:1,3:1,4:1,5:1}`, directed 3-cycles `0`, and
Hamiltonian path count `1`.

The hidden no-shortcut fact is that local equality or regular scores cannot
jump directly to LRC safety.  They must pass through a cascade law with the
root/source threshold labels intact.

## Assumption Challenge

Do not assume the vertices are runners, and do not assume
vertex-transitivity means cyclic polygonality.

Candidate vertex sets considered here include circle points, basepoints,
dihedral blocks, group elements, root stars, connection-set clauses,
threshold labels, endpoint owners, denominator shields, and proof obligations.

The chosen quotient preserves the predicate "does a local rooted view have a
specified global propagation law?"  It destroys fine endpoint geometry and
exact LRC time intervals; those re-enter through HYP-2122/HYP-2123 pincer
labels.

The challenged assumption is the strong equivalence
`vertex-transitive trienerment <=> regular polygon`.  The repaired statement
is:

```text
cyclic primitive VT trienerment <=> regular polygon;
general VT trienerment <=> local profile + cascade law.
```

## See

`04-computation/vt_trienerment_polygon_rigidity_s589.py`,
`05-knowledge/results/vt_trienerment_polygon_rigidity_s589.out`,
`07-reflections/vt-trienerment-polygon-rigidity-s589.md`,
HYP-2124, HYP-2123, HYP-2122, HYP-2121, HYP-2120, THM-400, THM-389, THM-381,
THM-052, HYP-2113, HYP-2109, HYP-2091, HYP-1977.
