# LRC14 q-uniform topology and the q-specific arithmetic break

HYP-3312 named the C2 coordinate as the Borsuk-Ulam antipodal.  The useful
new move is to take that name seriously as a limitation.

The topological degree argument is q-uniform.  It does not know whether the
canonical Goddyn-Wong magnitude break is ON at q=4 and q=7 but OFF at q=5 and
q=6.  The executable HYP-3423 table makes this boring in the best possible
way: the C2/BU residue charge is present on every q-row from 3 through 22,
while the HYP-3413 magnitude switch fires only on q == 1 mod 3.

So the proof split should be policed this way:

```text
C2/BU topology closes residue charge.
C3/real-cubic trace organizes equioscillation.
q mod 3 arithmetic gates the GW magnitude census.
owner-current labels discharge local mixed fibers.
S259/HYP-3418 supplies the two-adic covering-floor warning.
Rprime/decorrelation floor remains the critical inequality.
```

That also clarifies the pulled HYP-3417 result.  The frontier owner-current
cut `{2:g2,11:g1,13:g1}` is one even-cover label plus two binding labels.
That is a beautiful local echo of the global residue/magnitude split, but it
is still a finite packet.  It cannot replace the q-mod-3 switch, and after
HYP-3415 it cannot replace the covering-floor theorem unless it feeds that
floor.  After S259/HYP-3418, the label `2:g2` looks even more important: it is
the live two-adic/even-cover coordinate inside the local owner-current cut.

The new labelled-packet theorem target is a guardrail for quotients:

```text
Any quotient that forgets q-specific magnitude data must either restrict its
claim to the residue/equioscillation half, or carry an arithmetic/floor/owner
sidecar that resurrects the lost coordinate before the magnitude conclusion.
```

This is a small result, but it closes a recurring false route.  We should stop
asking a conserved q-uniform charge to do q-specific work.  The charge names
the seed.  The arithmetic and the floor do the break.
