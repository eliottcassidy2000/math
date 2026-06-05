# Unit Distance Moser Fixed Quantum S648

The LRC/Pillai packet S646 had a clean fixed equation:

```text
A(C) = C
```

for the antipodal gcd-shell mass of `C=2n-1`.  The unit-distance problem does
not hand us the same divisor sum.  Its fixed carrier is dynamical: add one
Moser slab and ask what edge-channel vector repeats.

For the THM-408 Moser spine ladder, the answer is rigid after the cap
transient:

```text
Delta_total = (0,1,8,4,0,4,4,2,4)
Delta_spine = (0,0,0,1,0,1,4,0,2)
Delta_bulk  = (0,1,8,3,0,3,0,2,2)
```

So the stable edge quantum is

```text
27 = 8 + 19.
```

This is the unit-distance version of a fixed point.  The add-one-slab operator
returns the same direction ledger every time:

```text
new slab -> 8 section edges + 19 bulk edges.
```

That matters because `27` is also the LRC `n=14` clock from THM-401 and the
Pillai-fixed clock from HYP-2222.  On the LRC side, `27` is a gcd-shell modulus.
On the unit-distance side, `27` is an edge quantum over the Moser unit shell.
Same number, but the retained side channel is different.

The scalar `27` is also not regular.  There are `9` direction pairs, so
`27 = 3*9`, but the direction vector differs from the uniform `3` vector by

```text
(-3,-2,5,1,-3,1,1,-1,1).
```

That defect is the useful part.  A proof route that remembers only "27 more
edges" throws away the exact direction imbalance that the construction uses.

The `m=2` frontier rows have an extra little balance:

```text
P_2^-: n=21, E=57, spine=20, pure_bulk=20
P_2^+: n=22, E=60, spine=21, pure_bulk=21
```

Pure-bulk direction mass equals the unit-spine length.  That is not true for
all `m`; it singles out the current `21/22` frontier as a section/bulk balanced
pocket.  I do not want to overstate it, but it feels like the right kind of
side-channel invariant: it is invisible to raw edge count and too structured
to be random in the named carrier.

The n=22 proof implication is concrete.  The `P_2^+` row has one degree-3 cap
vertex:

```text
degree_hist = {3:1,4:3,5:8,6:5,7:5}.
```

S614 says a 61-edge 22-point graph cannot have a degree-3 vertex, because
deleting it would leave 58 edges on 21 vertices, impossible if `u(21)=57`.
Therefore a 61-edge improvement near this carrier must repair the cap endpoint
channel.  It cannot simply be "one more anonymous bulk edge."  It must be an
endpoint-compatible ear, or else it must abandon the fixed Moser quantum.

That gives the unit-distance no-leak target:

```text
fixed 27 quantum + unrepaired cap endpoint -> 60
61 -> repaired cap endpoint or different carrier
```

This mirrors S646 exactly enough to be useful:

```text
LRC: mass-changing row -> loose; mass-fixed row -> AP/Vstar or side-label exit
UD:  cap-unrepaired fixed quantum -> 60; 61 -> endpoint repair or new carrier
```

The tournament lesson is also the same.  Vertices should not be points.  For
this proof route, vertices should be carrier lenses: fixed quantum,
spine/bulk split, pure-bulk jackknife, cap endpoint repair, deletion core, and
traceable section word.  Raw edge count belongs at the bottom of that
tournament because it is exactly the quotient that forgot the side channel.
