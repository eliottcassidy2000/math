# Unit Spines as Traceable Sections (S628)

S626 left the unit spine as a strong pattern: the Moser `n=14`, `n=21`, and
`n=22` witnesses all had Hamiltonian paths made entirely of unit edges.  S628
makes the reusable part rigorous.

The key move is to stop treating the point set as a flat cloud.  In the Moser
carrier, the good witnesses are layered objects.  Each full layer is an
8-vertex square-prism slab with a Gray-code path through it, and consecutive
layers are glued by the unit bridge `(-1,1,0,0)`.  The caps at `a=0` truncate
the last slab without breaking the path.

This gives THM-408:

```text
P_m^+ has |P_m^+| = 8m+6 and a unit spine.
P_m^- has |P_m^-| = 8m+5 and a unit spine.
```

The exact S626 witnesses sit inside this theorem:

```text
n=14 = P_1^+  with 33 unit edges
n=21 = P_2^-  with 57 unit edges
n=22 = P_2^+  with 60 unit edges
```

So the Moser spine was not a lucky Hamiltonian path found after the fact.  It
is the visible section of a carrier fibration.

Abstractly, a unit spine is a **traceable section**:

- the base is the layer order `m,m-1,...,0`;
- each fiber is a traceable local slab;
- the gluing map is a unit edge between the endpoint of one fiber and the
  start of the next;
- the resulting global section is a Hamiltonian path in the unit graph.

This is the same kind of object the LRC side keeps asking for under other
names: certificate sheaf, endpoint-compatible ear, owner-preserving section,
and quotient with retained side-channel.  The theorem is small, but it clarifies
the architecture.  A tournament built from an arbitrary point order can forget
the section; a tournament whose tie path is chosen from the section records it.

The remaining open question is therefore sharper:

```text
Do all proof-relevant dense Moser extremizers admit a traceable-section
decomposition, or can a true extremizer force three incompatible ears?
```

That is a graph/certificate-gluing problem, not merely a distance-to-tournament
mapping problem.  The next computational proof route should classify `21`-core
extensions by whether the new vertex is an endpoint-compatible ear for some
unit spine, then compare that with the HYP-2176/HYP-2188 degree-`4/5` extension
ledger.

Incoming S627/HYP-2204 sharpens the same warning from the H-gap side: `H=7`
and `H=21` are forbidden scalar-collapse tests, not literal unit-distance
tournament counts.  THM-408 supplies one retained carrier packet for that
grammar, namely the spine section.  The tile/bulk packet and the spine packet
should be kept separate until the proof obligation says which quotient is safe
to forget.
