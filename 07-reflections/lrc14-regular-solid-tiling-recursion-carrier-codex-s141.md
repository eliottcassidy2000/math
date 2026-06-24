# LRC14 regular-solid and tiling recursion carriers

The prompt's geometry is useful once the exact carrier is named.

Platonic solids give the positive-curvature finite regular maps.  The regular
plane tilings are the zero-curvature wall: `{3,6}`, `{4,4}`, and `{6,3}`.
That boundary is the part that matters for LRC-style recursion, because it
separates finite symmetric skeletons from repeatable Euclidean packets.

The plane tiling recursions split cleanly:

```text
square:      self-dual Gaussian axis scalings, indices 4,9,16,25,...
triangle:    Eisenstein self-scalings, with dyadic spine 4,16,64,...
triangle/hex local bridge: 6 triangles per hexagon / per vertex star.
hexagon:     Eisenstein norm-7 self-recursion, 7,49,343,...
```

The guardrail is that hexagonal centered rings are a different sequence:
`7,19,37,61,...`.  The prompt's `7,49,...` is better interpreted as a norm
index, not as the ordinary centered hex patch count.

Archimedean solids preserve one vertex-figure word, so they are local-quotient
analogues.  Johnson solids preserve regular faces but lose vertex uniformity,
so they are finite residual-atlas analogues.  That is close to the current LRC
state: exact Farey/C27/unital labels first, uniform channels when they exist,
then bounded AP/GW/petal/K33 residual catalogues.

The POKE use is to classify future residual packets by preserved recursion
label:

```text
square self-dual,
triangular dyadic,
triangle-hex support-six,
hex norm-7,
Johnson-like finite residual.
```

The classification is only useful downstream of exact `M/q/Farey` and C27
branch labels.  Raw solid names should stay behind those labels.
