# LRC14 Regular Solids and Tiling Recursion Carrier

The solids prompt initially looked like it might become a cloud of nice
analogies.  The computation made it much sharper.

The three regular Euclidean tilings are exactly the zero-curvature vertex
figures:

```text
3^6, 4^4, 6^3.
```

That is a clean way to read the user's recursion list.  Squares carry the
self-dual grid law `n^2`.  Triangles carry dyadic self-subdivision `4^k`.
Triangles and hexagons exchange through six triangles per hexagon.  Hexagons
carry the center-plus-six-petals law: centered patches `1+3r(r+1)` or
hexaflake packets `7^k`.

The main new object is not Platonic, though the Platonic table is useful.  It
is the prism/antiprism annulus:

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

That is too well-aligned with the current LRC14 machinery to ignore.  HYP-2942
gave a 28-point q=3 unital that preserves pair incidence but not cyclic C27
carry.  The 14-gonal prism/antiprism gives a 28-vertex object that preserves
cyclic annular order and has exactly the LRC threshold curvature `1/14`.

The two objects are complementary:

```text
q=3 unital        -> pair incidence, lambda=1, repeated-pair obstruction
14-annulus        -> two 14-cycles, cyclic order, half-step/twist geometry
```

This suggests a concrete next experiment: build a 28-vertex annular label model,
put `AP,GW,H1..H13,D1..D13` on it, and ask what the HYP-2942 H12 conflict
becomes.  Does `GW H12->D3` versus `K33 H12->D9` become a twist?  A diameter?
A forced two-chart obstruction?  That is much more actionable than saying
"solids are beautiful."

Johnson solids also got a useful role.  They are not a global recursion law.
They are finite nonuniform defect packets: exactly the right metaphor for a
finite residual atlas after the uniform and annular carriers have done their
work.

The guardrail is the same as before: the hexagonal `7` and `49` are not
automatically the tournament forbidden value story.  They only become relevant
if the LRC residual constructs the correct THM-572 state-lift packet.  Geometry
gives the carrier; the proof still has to preserve the unit.
