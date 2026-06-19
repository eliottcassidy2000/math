# Unit-Distance Spine/Ear Recursion (Codex S16)

The useful new object is not another point-to-tournament convention.  It is an
ear certificate.

S625 and S626 already answered the first layer: the "flop" at `n=7` is a
coordinate-order flop.  The unit graph still has Hamiltonian paths.  S16 asks
why those paths keep surviving, and the answer in the tested carriers is very
strong: every vertex is endpoint-compatible.  Delete it, find a Hamiltonian
path in the smaller unit graph, and the deleted vertex touches at least one
possible endpoint.

This reframes the mandatory Hamiltonian path as a section:

```text
smaller unit spine + endpoint-compatible ear -> larger unit spine
```

That is exactly the kind of retained side channel the repo keeps finding in
LRC, OCF, A000568, and unit-distance work.  Raw scalar counts are usually too
flat.  The proof object is the certificate that survives quotienting.

The Moser rows are especially clean.  The exact or frontier beam rows
`n=14,21,22` each have one exact-edge deletion ear:

```text
33 -> 30
57 -> 54
60 -> 57
```

So the visible recursion is not just "there is some unit Hamiltonian path."
It is "there is a degree-compatible boundary cap that preserves the smaller
edge target and the unit spine."  This matches THM-408's slab/cap families
`P_m^-` and `P_m^+`: the path is the traceable section through the carrier
fibration, and the cap is an endpoint ear.

The complement signal is also now cleaner.  The nonunit complement can thread
from `n=6`, fails at `n=7` because the compact hexagon center is isolated in
the complement, and then returns from `n=8`.  That means "unit versus nonunit"
is not a stable dichotomy.  In large rows both graphs may have Hamiltonian
paths.  The invariant is not which color has some path; it is which path is
retained by the construction/proof certificate.

The next real danger is not a lexicographic flip.  It is an incompatible ear
packet.  A true geometric flop should look like:

```text
all candidate deletions either break traceability or attach only to internal
vertices of every smaller spine
```

or like a separator with three path-obligatory branches.  This is the graph
version of an LRC wall packet: all local exits exist in scalar shadows, but no
single retained address can satisfy them simultaneously.

The creative tournament move is therefore to let vertices be proof obligations:
ears, endpoint masks, direction-shell channels, core-extension types, and
S/U/L impurity words.  The point tournament is downstream.  It is a diagnostic
of a chosen gauge, not the source of the proof.

Concrete handoff: run this exact endpoint-mask audit on the stored exact
`n=21` graph6 cores, then classify `n=22` extensions by whether the new vertex
is endpoint-compatible.  A `61`-edge construction, if it exists, becomes much
more interesting if its extra unit edge breaks the unique exact ear pattern.
