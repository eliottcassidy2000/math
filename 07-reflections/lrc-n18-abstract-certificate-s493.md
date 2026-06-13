# n=18 LRC abstract certificate session (S493)

The user asked to stay with `n=18` and get more abstract each time the direct
route hit a wall. That turned out to be exactly the right rhythm. Each wall
looked like a failure of the obvious proof object, but each one also exposed a
smaller finite object.

## Wall 1: the two bridges are locally identical

The first instinct was to branch on the two `n=18` bridge choices, `6` and
`12`. Locally, this gives nothing: after the forced fan

```text
(1,5,7,9,11,13,17)
```

the two bridges cover exactly the same six residual owner-18 endpoints:

```text
53/324, 55/324, 161/324, 163/324, 269/324, 271/324.
```

So the bridge variable is not visible at the gate row. The abstraction is to
move from the owner-18 row to global owner-row shadows. There the bridge has a
signed profile across owner rows `7,8,9,10,11,12,16`; rows `8,9,10,12` favor
bridge `6`, while rows `7,11,16` favor bridge `12`.

This is the first finite certificate target: choose row weights so both signs
are expensive.

## Wall 2: ladders only trade gap for debt

The row-parent, gate, and double-gate ladders do not make progress by
shrinking the real gap. They preserve the product:

```text
gap/th * endpoint_debt = 1.
```

At first this feels circular. The abstraction is to quotient by the dyadic
translation. After normalizing the scale `9*2^r`, all three ladders expose
the same depth packet:

```text
96*(0,2) + 16*(0,3) + 64*(4,2).
```

This is the square-core packet. It is a little like a Cayley-Dickson-style
doubling shadow: the obstruction is copied along the 2-adic axis while the
`3^2` core stays fixed. The proof should not chase each copy separately.

## Wall 3: smaller gap can be misleading

A one-slot scan again showed that nonmultiples can shrink the archimedean gap
more than multiples. But those candidates pay extra endpoint debt. The search
objective has to be two-coordinate:

```text
(real gap, endpoint/product debt).
```

This is the same moral as the user's "gap is archimedean size, endpoint debt
is 2-adic size" heuristic. For `n=18`, the conserved product is not a finished
proof, but it is a warning that scalar optimization is looking at only one
projection.

## Wall 4: pressure still peels

The last check reused the S492 pressure tournaments on `n=18` lpd, gate, and
single-gate repair rows. Nearest pressure, two-neighbor pressure, and deficit
pressure all had largest SCC `1` and no directed pressure triangle in the
sampled exact rows.

That is bad for disproof hunting but good for proof building. Pressure leaves
should be assets: a counterexample-like row has to survive endpoint-private
peeling and mobile pressure-leaf peeling simultaneously.

## Current proof shape

The proof shape after the walls is:

```text
1. Local fan lemma:
   owner-18 endpoints force the unit/half fan and one bridge variable.

2. Global bridge-charge lemma:
   owner rows 7,8,9,10,11,12,16 distinguish the two bridge signs.

3. Product-packet lemma:
   row-parent/gate/double-gate ladders are dyadic translates of one
   square-core packet.

4. Peeling lemma:
   pressure leaves remove mobile blockers unless a labelled pressure SCC
   appears.
```

The next concrete move is linear: build a row-weight certificate on the global
bridge rows plus the square-core packet, with pressure SCC as the only escape
hatch.
