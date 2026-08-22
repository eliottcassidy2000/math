---
id: THM-3173
title: "Six-state free-factor actions and pointed frame cube"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/codex-thm3088-push-2026-08-02
audit: >
  Two independent hostile audits rebuilt the action trichotomy, order-18
  fibre product and free-kernel rank, all cube automorphisms and translations,
  both S4 complements, orientation cocycle, and the semantic/physical boundary.
  The final companion exhausts all 8! frame permutations and fresh normal and
  optimized runs LF-byte-match the stored transcript and declared hashes.
depends_on:
  - THM-3157-pointed-resolvent-c3-lift-edge-hexagon-and-relative-tournament
  - THM-3145-bass-serre-two-three-tree-and-tetrahedral-congruence-quotient
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
  - THM-2641-modular-abelianization-theta-blindness-and-637-residue-no-go
  - THM-2768-modular-c2-c3-quotients-to-a4-s4-and-bass-serre-cycle-ranks
related:
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
  - THM-2951-fifth-compound-reconstruction-and-v-four-phase-scalarization-boundary
script: 04-computation/six_state_free_factor_pointed_frame_cube_thm3173.py
output: 05-knowledge/results/six_state_free_factor_pointed_frame_cube_thm3173.out
hash_basis: LF-normalized bytes
script_sha256: cad2d746e2a297cd54eee1e0e8907f0d2678fb5875969e06941f71d2041f05bc
output_sha256: 9d0a686466aea53e512657298a75f1fdb5c630981e954068d8d87db0160d2f93
semantic_sha256: 960ba61084e4d5cce401b42f44b6232a958460807da62975dc27310d803f1d0d
---

# THM-3173 -- six-state free-factor actions and pointed frame cube

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## Trigger

Use this card whenever two constructions both produce six states or an
unoriented hexagon and someone proposes identifying them.  Cardinality and
abstract graph type are insufficient.  First record the acting group,
generator images, edge colours, stabilizer, and lost kernel.

## The three objects

| object | state carrier | action | exact survivor | exact loss |
|---|---|---|---|---|
| `C2*C3=PSL2(Z)` | infinite Bass--Serre/Farey ancestry, not six states | noncommuting free factors | reduced-word order and tree position | nothing until a quotient is chosen |
| regular `S3` | six Farey total orders / matching permutations | left-regular `S3`; Cayley `C6` edges are right adjacent transpositions | conjugating frame, reflection parity, alternating `L/R` colours | the index-six free kernel |
| regular `C6` | six exponent classes | commuting half-turn and third-turn; one-step rotation exists | exponents mod `2` and `3` | every commutator and all Farey ancestry |

The regular action cycle inventories distinguish the latter two without any
labels:

```text
regular S3:  1^6 once, 2^3 three times, 3^2 twice, no 6-cycle;
regular C6:  1^6 once, 2^3 once, 3^2 twice, 6 twice.       (1)
```

Thus an abstract `C6` support graph does not identify its symmetry action.
THM-2597/2632 already supply the two sides of this distinction; THM-2641
proves their gluing is noncanonical from the current modular residue data.

## Sharp joint quotient

Fix the standard quotient conventions

```text
s -> (transposition,3),        r -> (three-cycle,2)
             in S3 x C6.                                  (2)
```

The joint image is not the `36`-element product.  The two coordinates share
their binary character:

```text
im(C2*C3) = S3 x_(C2) C6
          = {(g,a):sgn(g)=(-1)^a},             order 18. (3)
```

Both projections are onto and each has a `C3` kernel.  The kernel of the
joint map is torsion-free of index `18`; Euler characteristic
`chi(C2*C3)=-1/6` makes it free of rank `4`.  Hence the sharp common state is
an eighteen-state fibre product, not a preferred bijection between the two
six-sets.

Two minimal hostiles show that neither six-state quotient determines the
other:

```text
[s,r]:   invisible in C6, nontrivial in S3;
(sr)^2:  identity in S3, residue 4 in C6.                 (4)
```

This is the reusable test: find one commutator witness and one quotient-face
relation witness before attempting any six-state identification.

## Oriented pointed quartic frames

THM-2606 gives four unoriented owner cycles `C_o`.  THM-3157 adds a chosen
matching generator and obtains the two orientations at each owner.  The eight
frames form the homogeneous space

```text
F={(o,epsilon):o in X, epsilon=+-1}=S4/C3,               (5)
tau(o,epsilon)=(tau o,sgn(tau)epsilon).                  (6)
```

The stabilizer of one frame is exactly the owner-fixing `C3`.  Give `F` the
adjacency

```text
(o,epsilon) ~ (o',-epsilon)       iff       o!=o'.       (7)
```

Then `F` is `K4,4` minus the owner matching, hence `Q3`.  After matching one
base `C3` frame, this is an `S4`-equivariant realization of the eight-vertex
`S4/C3` side of THM-2768's suppressed cube.  If that construction retains
only the `C3` subgroup and not a generator/base frame, there are two such
identifications, differing by the antipode below.  Intrinsically,
`(o,epsilon)` is the unique
three-cycle `R_o^epsilon`: the eight frames are exactly the conjugacy class
of all three-cycles in `S4`.  Two frames are cube-adjacent exactly when the
product of their three-cycles is a nonidentity `V4` double transposition.
The antipodal frames at one owner are the inverse pair `R_o,R_o^-1`.

The orientation flip

```text
a(o,epsilon)=(o,-epsilon)                              (8)
```

is the cube antipode.  It centralizes the `S4` action but is not induced by a
sheet permutation; adjoining it gives an order-`48` cube action.  Quotienting
by `a` gives the four owners and the `K4` owner quotient.

The cube has a sharper canonical translation structure.  The normal quartic
`V4` already acts regularly on the four owners, and the antipode is the deck
group of orientation-forgetting.  Thus

```text
T=V4 x <a>                                                  (9)
```

acts regularly on the eight frames.  This action is intrinsic; choosing a
base **frame** identifies `F` with `T`, while choosing only a base owner leaves
the two antipodal lifts and hence an `a`-gauge.  For each nonidentity
`v in V4`, put

```text
t_v=(v,a).                                                 (10)
```

The three `t_v` are exactly the three cube-edge translations.  They generate
`T`, two distinct axes multiply to the remaining pure `V4` translation, and
the product of all three axes is `a`.  Hence

```text
Q3=Cay(V4 x C2,{(v,a):v!=1}),                              (11)
Aut(Q3)=T semidirect S3 isomorphic to S4 x C2.             (12)
```

The set of three axes is canonical and `S4`-stable; choosing coordinates only
labels them.  The displayed `S4` factor, however, is not recoverable from an
abstract cube alone.  There are two natural `S4` complements to `<a>`: the
physical sheet action

```text
(o,epsilon) -> (tau o,sgn(tau)epsilon)
```

and the plain owner action `(o,epsilon)->(tau o,epsilon)`.  They differ on
odd permutations by `a` and intersect exactly in `A4`.  The quartic sheet
semantics and descent cocycle select the first; the unlabeled cube does not.

The THM-3157 maps

```text
(o,epsilon) -> directed h_(o,epsilon) hexagon,
(o,epsilon) -> relative tournament                        (13)
```

are `S4`-equivariant bijections onto eight distinct structures.  Forgetting
hexagon orientation is exactly the antipodal quotient and leaves the four
inherited `C_o`.

## Exact descent cocycle

Relative to the set-theoretic section `o -> (o,+)`, the failure of
equivariance is

```text
kappa(tau,o)=sgn(tau) in C2,                              (14)
```

independent of `o`.  It obeys the cocycle law.  Exhaustion of all `2^4`
ownerwise orientation choices gives

```text
S4-equivariant sections: 0;
A4-equivariant sections: 2 (the two global orientations). (15)
```

Equivalently, an odd transposition fixing an owner reverses its two directed
hexagons.  This is the sharp descent statement: restricting the abstract
`S4` action to `A4` splits the frame cover into two four-point orbits, while
the odd deck interchanges them.  An unmarked restriction or retained
discriminant square/class does not choose one orbit.  A chosen discriminant
square root/orientation **does** choose one orbit and the corresponding
matching generator, exactly as in THM-3157; it still supplies no sheet owner.

## Counterindications and physical boundary

- Do not call `C2*C3` a six-state action; only a specified quotient is.
- Do not identify the regular `S3` and `C6` actions from their common
  uncoloured hexagon.
- Do not infer a modular sheet action from `h`: its binary factor is the
  external edge complement, not a quartic `S4` permutation.
- Do not infer a physical owner, current, Keller operation, or LRC move from
  the semantic `S4/C3` cube.  THM-2681 remains the global resolvent hostile.
- A tournament or hexagon selector is physical only after a common carrier
  supplies the owner, matching orientation, and nonzero coefficient/current.

## Exact control

The companion

```text
04-computation/six_state_free_factor_pointed_frame_cube_thm3173.py
```

checks `(1)--(15)`, all `18` joint quotient elements, both projection kernels,
the exact `18-9-6+1=4` Bass--Serre rank census, all `24` sheet actions, all
`16` possible ownerwise sections, the full cube distance table, the central
antipode extension, the regular eight-element translation group, all `8!`
cube permutations and the exact `48` automorphisms, both `S4` complements,
every hexagon cycle type, every tournament pair, and all hexagon/tournament
equivariance cases.
Its output is

```text
05-knowledge/results/six_state_free_factor_pointed_frame_cube_thm3173.out
```

Final LF-normalized hashes are

```text
script cad2d746e2a297cd54eee1e0e8907f0d2678fb5875969e06941f71d2041f05bc
output 9d0a686466aea53e512657298a75f1fdb5c630981e954068d8d87db0160d2f93
```
