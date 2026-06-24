# LRC14 Regular-Solid And Tiling Recursion Carrier

S141 tested whether Platonic, Archimedean, and Johnson solids, together with
the Euclidean triangular/square/hexagonal tilings, contribute a useful labelled
recursion grammar for LRC14.

The answer is yes, but only after the analogy is made local and labelled.

The useful grammar is:

```text
Euclidean tilings       = flat/tight recursion wall
Archimedean near-flat   = labelled branch-local semiregular chart
Johnson local defects   = finite residual surgery patches
```

The square, triangular, and hexagonal tilings have zero local deficit.  This is
why they are proof baselines rather than endpoint obstructions.  Square
self-recursion (`4,9,16,25,...`) is valuable mostly as a warning: a self-dual
grid can be too scalar and too global.  Triangle/hex duality is more live
because it keeps the C27/F3^3 chart distinct from the six-sector cover.  The
hexagonal norm recursion `N(3+omega)=7`, hence `7,49,...`, is the clean
geometric avatar of the mod-7 petal/coimage seam; centered hex counts
`7,19,37,...` are a different carrier and should not be conflated with it.

The Platonic solids are positive-curvature caps.  Their global symmetry is
beautiful, but for LRC14 it is usually too rigid unless exact `M`/Farey/C27
labels have already been attached.  Archimedean solids are more useful: the
near-flat mixed words such as `4.6.10`, `3.4.5.4`, and `5.6.6` behave like
uniform branch-local charts.  Johnson solids are the closest match to the
actual residual frontier because they are local nonuniform defects rather than
one global atlas.

Tournament Analysis appeared at two levels.  At the global proof scale, the
transitive carrier order keeps exact `M`/Farey and C27/unital-block labels
above all geometric analogies.  At the local-curvature scale, the transitive
order was:

```text
johnson_local_defect
> archimedean_near_flat
> tri_hex_dual
> hex_heptadic
> triangle_self
> platonic_positive
> square_self
> raw_runner_vertices
```

The pairwise observable was channel preservation: exact `M`/Farey label, C27
carry, mod-7 petal seam, pair-incidence unit, branch locality, and resistance
to scalarization.

Assumption challenge: before settling on this quotient, I considered runners,
gaps, fixed circle sections, section boundaries, wall-crossing events,
residues, cover arcs, Fourier modes, matroid circuits, solid vertex figures,
local defect patches, and proof obligations as possible tournament vertices.
The chosen quotient preserves branch labels and proof obligations but destroys
individual runner ownership and detailed time order inside each cell.

Low-gap translation:

```text
AP                    -> flat Euclidean baseline
GW                    -> triangle/hex plus hex-heptadic labelled collision
12->36                -> Johnson-style nonunit K33 defect
10->20, 13->26        -> hex-heptadic unit-visible petals
two-swap S138 rows    -> Johnson-style finite surgery splices
```

The main guardrail is inherited from HYP-2942: a geometric carrier may be
branch-local and still useful.  It should not be treated as a global atlas
unless its preserved unit is declared.  The next serious proof test is to show
that every non-tight residual in the AP/GW/C27 frontier has a unique local
defect address after exact `M`, Farey, and C27 labels are retained.
