# Tournament value-origin ledger - S221

This pass was renumbered after incoming HYP-3054/T1136/LTI-201/LTT-099 landed
as the observer-extension cut payload calculus.  The ledger below is the
count-origin sidecar for that calculus: once an observer/cut payload is named,
we still need to know whether its visible small number came from a class
count, rooted count, fixed branch, deletion fiber, edge sector, path fiber, or
cycle residue.

The session corrected the arithmetic first: `48 + 12 = 60`, so the first
shifted A000568/rooted-perspective failure is actually `48 + 8 = 56`.

That does not make the user's `12` signal disappear.  It makes it more
interesting.  The same value is sitting on three different floors:

```text
R(4)  = 12  rooted/node perspectives on 4 vertices
U(5)  = 12  unlabelled 5-tournament classes
SC(6) = 12  self-converse 6-tournament classes
```

So `12` is not the missing additive mass at the first failure.  It is a
diagonal alignment across quotient types.  The proof-relevant question is
therefore not "why is 12 everywhere?" but "which quotient produced this 12?"

The ledger now gives a useful hierarchy around the first defect:

```text
U(5) parent classes                         = 12
raw incident words                          = 384
parent-automorphism word orbits / R(6)      = 296
unrooted child sinks / U(6)                 = 56
rooted 5-perspective plus incident word     = 1408
ordered-pair perspectives on U(6)           = 1408
directed-edge perspectives                  = 704
```

This reframes the old count coincidence as a typed recursive transport law:

```text
parent class + incident word orbit -> rooted child -> unrooted sink.
```

The cross-sector sidecar from S213 is still the cleanest compact repair at
the first failure.  Sector sizes and internal sector tournaments separate
`55/56`; cross-sector orientation separates `56/56`, with the only
size/internal collision being the converse pair `344/345`.

S217's rectangle/hourglass flow gives the fixed-path analogue.  The redundant
inter-layer lines in a fixed Hamiltonian-path half-tiling decompose into local
rectangle residues plus hourglass residues linking adjacent bridges:

```text
fixed_red(n)=2*C(n-2,3)+C(n-3,2).
```

This feels like the same lesson in a different coordinate system: the count is
not the carrier; the residue explaining duplication is the carrier.

For LRC14, the immediate rule is to tag numerical shadows by origin before
using them.  A quotient should declare whether its small value is a class
count, rooted count, self-converse branch, incident-word orbit, deletion fiber,
edge-sector sidecar, fixed-path fiber, or rectangle/hourglass residue.  If the
origin tag is omitted, a scalar can accidentally merge incompatible proof
obligations.

Assumption challenge: I considered runners, vertices, arcs, ordered pairs,
edge sectors, deletion fibers, self-converse fixed branches, fixed Hamiltonian
paths, rectangle/hourglass cycles, LRC endpoint owners, and proof obligations.
The useful vertices here are value-origin carriers, not tournament vertices.
The preserved predicate is controlled quotient admissibility for LRC-style
proof carriers; the destroyed data is the specific lost coordinate behind a
small numerical coincidence unless the origin type and sidecar are retained.
