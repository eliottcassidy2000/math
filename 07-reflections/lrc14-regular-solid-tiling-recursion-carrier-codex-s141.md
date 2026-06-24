# LRC14 Regular Solids and Tiling Recursion Carrier

The solids prompt could have dissolved into pretty analogies.  The useful move
was to force every geometric object to say what it preserves.

The regular map boundary is the first guardrail:

```text
(p-2)(q-2)<4  spherical finite maps
(p-2)(q-2)=4  Euclidean tiling wall
(p-2)(q-2)>4  hyperbolic maps
```

So Platonic solids are positive-curvature finite skeletons, while the regular
plane tilings `{3,6}`, `{4,4}`, and `{6,3}` are the zero-curvature recursion
wall.  That boundary matters for LRC because finite symmetry and repeatable
recursion should not be fused into one scalar "regularity" word.

The plane tiling recursions split cleanly:

```text
square:      self-dual Gaussian axis scalings, indices 4,9,16,25,...
triangle:    Eisenstein self-scalings, with dyadic spine 4,16,64,...
triangle/hex local bridge: 6 triangles per hexagon / per vertex star.
hexagon:     Eisenstein norm-7 self-recursion, 7,49,343,...
```

The hexagonal guardrail is important: centered rings are a different sequence,
`7,19,37,61,...`.  The prompt's `7,49,...` is better read as a norm-index
chain via `N(3+omega)=7`, not as ordinary centered patch growth.

The main new object is not Platonic, though the Platonic table is a useful
calibration.  It is the prism/antiprism annulus:

```text
n-gonal prism:      (4,4,n)
n-gonal antiprism:  (3,3,3,n)
```

Both have:

```text
kappa = 1/n
V = 2n.
```

So `n=14` gives:

```text
V = 28
kappa = 1/14.
```

That is too aligned with the current LRC14 machinery to ignore.  HYP-2942 gave
a 28-point q=3 unital that preserves pair incidence but not cyclic C27 carry.
The 14-gonal prism/antiprism gives a 28-vertex object that preserves cyclic
annular order and has exactly the LRC threshold curvature `1/14`.

The two objects are complementary:

```text
q=3 unital        -> pair incidence, lambda=1, repeated-pair obstruction
14-annulus        -> two 14-cycles, cyclic order, half-step/twist geometry
```

This suggests the concrete next experiment: build a 28-vertex annular label
model, put `AP,GW,H1..H13,D1..D13` on it, and ask what the HYP-2942 H12
conflict becomes.  Does `GW H12->D3` versus `K33 H12->D9` become a twist?  A
diameter?  A forced two-chart obstruction?  That is much more actionable than
saying "solids are beautiful."

Archimedean solids preserve one vertex-figure word, so they are local-quotient
analogues.  Johnson solids preserve regular faces but lose vertex uniformity,
so they are finite residual-atlas analogues.  That is close to the current LRC
state: exact Farey/C27/unital labels first, uniform channels when they exist,
then bounded AP/GW/petal/K33 residual catalogues.

The final guardrail is unchanged: the hexagonal `7` and `49` are not
automatically the forbidden tournament value story.  They only become relevant
if the LRC residual constructs the correct THM-572 state-lift packet.  Geometry
gives the carrier; the proof still has to preserve the unit.
