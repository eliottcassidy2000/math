# The heptagon node is a basepoint, not the sheet state

*codex-2026-07-14-S6, after THM-769/771, the HYP-6825 atlas, the incoming
tooth-winding guardrail, and the folded unit-grid work.*

The original mapping question now has an exact answer, and the answer exposes
the missing object.

At the prime sheet `c=7`, a seven-owner exact cover can be sent into the
tournament-tiling explorer in two canonical ways.  A marked cut makes the
owner order transitive and lands at `n7-a000`.  Circular three-step precedence
makes the rotational heptagon tournament and lands at `n7-a267`.  Choosing the
least Hamiltonian path then maps every one of the `7!=5040` owner-labelled
sheet assignments onto exactly the 25 staircase masks stored in `a267`'s
inverse fibre.

So a map exists.  What fails is the hope that its target node is the whole
state.

## The exact forgetting ladder

The prime-sheet calculation gives a rare opportunity to count every layer:

| layer | number | retained structure |
|---|---:|---|
| owner-to-sheet bijections | 5040 | every owner label and marked sheet |
| rotation orbits | 720 | labelled circular tournament / cyclic differences |
| dihedral orbits | 360 | circular order without orientation |
| least-path staircase masks | 25 | fixed-Hamiltonian-path presentations of `R_7` |
| merged isomorphism node | 1 | the heptagon class `n7-a267` |

The 25 is not numerology.  `R_7` has 175 Hamiltonian paths and automorphism
group `C_7`, so the path presentations have size `175/7=25`.  HYP-6825's
objective order places the node at rank 267, local depth six, local path word
`SSRRSS`, and blue/black root word `BBBBBB`.  The cut and circular gauges
differ on the same six long chords.

This is the requested objective position of the cyclic tiling inside the
merged metagraph.  It also quantifies every fibre above that position.

## What the prime-sheet object actually is

For owner `a`, away from endpoints,

```text
k_a(x)=-w_a^{-1} round(w_a x) in F_7.
```

Thus the fibre over a core-safe base point `x` is not a bare tournament.  It
is a labelled token configuration in `F_7^r`.  Its exact polynomial is

```text
D_x(X)=product_a (X-k_a(x)).
```

The sheet deck is covered precisely when

```text
X^7-X divides D_x(X).
```

At seven owners this becomes equality, and the six nontrivial power moments
are a complete instantaneous certificate.  At more than seven owners the
quotient

```text
Q_x(X)=D_x(X)/(X^7-X)
```

on a covered chamber is an exact redundancy polynomial.  For eight owners it
is linear: its root is the duplicated sheet.  This looks like a better
`r>=8` sidecar than a scalar overlap count because it retains *where* the
redundancy lives.

An endpoint event is not an edge between two tilings.  The event owner has no
strict bad sheet at the wall, so its factor disappears.  The wall is covered
only if the remaining polynomial is still divisible by `X^7-X`.  In
particular an eight-owner covered wall must have one event owner and seven
remaining owners forming a heptagon state.

This is exactly what happens in HYP-6835's survivor:

```text
W=(108,169,143,213,206,197,30,162),  x=19/216.
```

Owner `108` is absent.  The other tokens are `(6,5,3,1,4,2,0)`, and their
least-path presentation is mask `32153` at `n7-a267`.  The eight-owner wall
contains a codimension-one seven-owner atlas stalk.

## Transport is additional mathematical data

Crossing an event moves a token by

```text
k_a -> k_a-w_a^{-1}.
```

The first moment therefore gives a return holonomy:

```text
sum crossed w_a^{-1}=0 mod 7
```

between exact chambers.  Higher moments give the complete destination test.
But even the complete instantaneous token assignment is not a complete future
state.  THM-773 exhibits the same assignment and the same atlas mask `27833`
in two speed rows.  One next frees sheet 6 at an event of owner 10; the other
next frees sheet 5 at an event of owner 4.  Their inverse steps and endpoint
schedules differ.

There is also a carry around the base circle:

```text
k_a(x+1)=k_a(x)-1.
```

Resetting `x=1` to `x=0` without translating the sheet labels caused the
first implementation failure in this session.  That bug was useful: it is a
literal demonstration that monodromy is not prose decoration.

The continuation-complete local state therefore has the form

```text
(metric core-safe base,
 owner-to-sheet incidence,
 inverse-step vector,
 cyclic endpoint word,
 global sheet carry,
 threshold/strictness flag).
```

The node and moment vector are checksums of this packet.

## How this refines the four-coordinate LRC object

HYP-6815's full carrier was the incidence suspension in `(u,t,c,lambda)`.
THM-773 is the slice `c=7`, `lambda=1/14`.  Even on that slice, the object is
not a two-dimensional picture or a tournament.  It is a ramified finite local
system over the core-safe subset of the `x` circle:

```text
core-safe base x
  -> finite owner-token fibre
  -> polynomial cover predicate
  -> endpoint deletion/update maps
  -> metagraph node and path-gauge shadows.
```

The concurrent tooth-winding atlas supplies the complementary warning.  A
topological winding can be correct while every metric component escapes.  So
the finite-field polynomial solves the *sheet* predicate but cannot replace
the metric base.  The folded unit-grid work likewise supplies obligations on
which base points must be covered; it does not replace persistent ownership
between them.  The three strands fit without being conflated:

```text
unit-grid obligations = required base samples,
metric tooth word     = geometry between samples,
prime token polynomial= exact finite-fibre coverage.
```

## Fold before quotient

THM-774 arrived after the first prime-sheet computation and sharpens the
methodological distinction.  Its passage from two odd owners `(x,y)` to the
half-sum and half-difference modes `(a,b)` is not compression: the signed
eligibility/colour predicate is exactly conjugated to one folded-diamond
inequality.  Compression happens only later, when the whole pointwise diamond
is replaced by its measure.  The former preserves truth; the latter retains
only a necessary scalar shadow.

The seven-token analogue should therefore be sought before another quotient.
Possible coordinates include the discrete Fourier modes of the labelled
occupancy function, signed chord differences between tokens, or the elementary
symmetric coefficients of the redundancy polynomial.  A useful transform
must carry both the endpoint deletion and the translations
`k_a -> k_a-w_a^(-1)` exactly.  Only after finding such a conjugacy should we
ask which coordinates can be discarded under continuation equivalence.  This
reframes “choose better tournament vertices” as the second question; the first
is “which nonlinear change of coordinates makes transport local?”

## The next frontier

The most focused next problem is the eight-owner transport problem:

> Follow the absent owner over the heptagon stalk while `x` moves through one
> core-safe component.  Either a simultaneous event occurs, polynomial
> divisibility fails, or a rigid periodic word of `a267` fibre masks survives.

At eight owners a simultaneous two-owner event cannot stay covered, because
only six tokens remain.  A hypothetical survivor must therefore thread only
simple events, with the seven non-event owners returning to exact heptagon
states each time.  The inverse-step holonomy and the 25-mask word are now
finite exact invariants of that thread.

For more owners, study the redundancy polynomial `Q_x`, not only the number
of overlaps.  Its degree is `r-7`; wall events delete factors and chamber
motion transports roots.  The natural metagraph vertices may then be
redundancy factors, absent-owner obligations, or wall events—not runners.

The broad lesson is precise: an isomorphism class is an excellent base-chart
address.  Essential structure is the labelled stalk plus its transport.
