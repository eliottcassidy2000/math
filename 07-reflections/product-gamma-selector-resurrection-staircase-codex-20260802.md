# Product-Gamma selector resurrection staircase

**Status:** research reflection.  The depth-four/degree-12 cell is proved in
THM-3155.  The depth-five/degree-13 cell is proved and independently audited
in THM-3158.

## The first closed two-axis cell

For the support `(1,3)`, bank-`I2` selector problem, write

```text
C_D^(<=d)={laws on physical prefixes of depth <=d
           that are Hasse-positive in every degree 5,...,D}.
```

The exact proved pattern is

```text
C_11^(<=3)=empty,
C_11^(<=4) nonempty,
C_12^(<=4)=empty.
```

The positive corner has an explicit four-state law.  The two empty corners
have exact positive Farkas certificates.  This is stronger than a list of LP
outcomes: it is a closed birth/death cell in a monotone two-parameter family.

The monotonicities point in opposite directions:

- adjoining prefix states can only enlarge `C_D^(<=d)`;
- adjoining degree constraints can only shrink it.

Depth four crosses the first wall; degree 12 crosses the second.

## Why the cell is structurally small

The depth-three obstruction uses four upset rows.  The depth-four death adds
only the degree-12 singleton-complement row.  Every row factors through the
THM-3147 profile by partition length and singleton count.  Hence the active
dual carrier grows as

```text
four compressed rows  -->  the same four plus one rank tail.
```

The huge integer coefficients are coordinates of a ray in this small carrier,
not evidence that the obstruction intrinsically has huge dimension.  The
correct next algebraic object is the cone-valued pole-subtraction recurrence
on the endpoint jet.

## The second staircase cell

THM-3158 supplies an exact depth-five law on seven states with
denominator `10^6` and numerators

```text
(1)             97856
(1,1,1)         56951
(1,1,1,1)      140643
(1,1,2,2,2)    194498
(1,2,2,3,3)      7398
(2,2,2,3,3)    118572
(5,5,6,7,8)    384082.
```

They sum to `10^6` and have gcd one.  Exact streaming enumeration proves
strict positivity on all `403,539` nontrivial upsets through degree 12.  A
separate primitive positive combination of nine upset facets is strictly
negative on all 682 physical depth-at-most-five states in degrees through 13.
Thus the next proved cell is

```text
C_12^(<=4)=empty,
C_12^(<=5) nonempty,
C_13^(<=5)=empty.
```

The depth-four five-row separator is crossed by exactly two depth-five states,
`(1,2,2,2,3)` and `(5,5,6,7,8)`.  As at depth three, the old obstruction dies
on a tiny new-state locus.  This is the finite-convex analogue of a bistellar
flip: the carrier remains almost valid, two cells change side, and a new
section appears.

The nine-row depth-five separator is crossed by 24 of the 507 new depth-six
states: nine low states over `{1,2,3,4}` and fifteen high-tail states clustered
around `{5,5,6,7,8}`.  Nevertheless a full separation oracle still reports
`C_13^(<=6)=empty`.  Until its new exact certificate is extracted, this last
statement is discovery evidence only.  The 24-state crossing locus shows why
the old pattern of merely adjoining one singleton-complement facet has ended:
the dual wall now mutates along a structured local patch.

## Possible global theorem and its obstruction

The physical pole multiset has only 16 letters.  Its complete submultiset
bank is finite:

```text
(4+1)(3+1)(2+1)^3(1+1)^3-1=4319 states.
```

Therefore this is not an infinite-dimensional escape problem.  If every
finite degree horizon were feasible on the full 4,319-state simplex, compactness
would produce one law satisfying every degree constraint.  Conversely, any
global failure has a finite Farkas certificate.

This makes two ambitious but decidable programmes visible:

1. continue the exact staircase until it stabilizes or dies before depth 16;
2. search directly for a structured law on the full state bank, using the
   finite staircase laws as moment constraints or warm starts.

Neither programme yet connects the averaged virtual-prefix current back to the
original product-Gamma response.  That realization map remains load-bearing.

## Product-law and stopping-law probes

The sparse laws found so far are arbitrary points of a simplex.  THM-3163
shows that bare sequential realization is automatic: sample a terminal set,
uniformly order it, and use posterior transitions.  That construction is
global and response-blind.  Two genuinely structured families should still
be tested before treating the coefficients as accidental:

- **independent deletion:** each pole label is retained/deleted by a Bernoulli
  parameter, possibly depending only on pole value;
- **sequential hazard:** remove one remaining labelled pole at each step with
  a value-dependent transition weight and stop with a depth-dependent hazard.

For either structured family, the source is a small parameter space, the target is the
state-law simplex, and the preserved predicate is exact multiplicity validity.
The destroyed information is arbitrary correlation among deleted poles.  The
cheapest test is whether the depth-four and depth-five rational laws lie in the
closure of either family, or whether their low-order inclusion moments violate
the model identities.

If a structured family survives, it supplies the missing sequential sidecar.
If both fail, their first polynomial relation gives a useful no-go theorem.

## Two Hasse networks and the missing chain map

THM-3163 separates two flows that earlier language blurred together.

- The **prefix network** is the Boolean/submultiset Hasse graph of physical
  pole deletions.  Its augmented unit flows have the terminal law as their
  sink distribution.  Every finite law admits such a flow.
- For each degree `N`, the **partition network** is the refinement Hasse graph
  of `P_N`.  Its nonnegative edge flows have the selector current
  `G_N(lambda)` as boundary.  Existence here is the substantive positivity
  condition.

At vertices the obvious observable is

```text
delta_sigma  |-->  G_N^sigma.                              (A)
```

A genuine sequential proof would require more than `(A)`: it needs an edge
map from physical prefix moves to partition-Hasse transports that commutes
with boundary and preserves the nonnegative cones.  In schematic form,

```text
prefix edge flow --L_N--> partition-Hasse edge flow
       | boundary                    | boundary
       v                             v
 terminal law   --(A)-->       selector current.           (B)
```

THM-3160 identifies information destroyed by the lower arrow when one tries
to evolve it: same-degree selector rows omit cross-degree endpoint Pluecker
minors.  The full exterior tensor restores a signed flat connection, but not
the positivity of `L_N`.  THM-3128's labelled-deletion obstruction is the
canonical hostile against the cheapest positive edge lift.

Thus the next bridge has explicit typing:

- **source:** a state-dependent prefix-flow edge;
- **target:** a partition-refinement transport, simultaneously in the needed
  degree range;
- **preserved predicate:** boundary and nonnegative flow;
- **lost datum:** cross-degree endpoint phase/minors;
- **needed sidecar:** the complete Pluecker pole tensor or a smaller
  fixed-bank quotient;
- **cheapest decisive test:** on the live support-`(1,3)`, bank-`I2` fibre,
  test kernel inclusion for the projected parent-current matrix versus every
  one-pole child matrix, then test cone rather than span if linear transport
  survives.

This is the selector analogue of a chain-level holotopy problem: path
independence exists in the signed enriched object, while positivity is lost
under projection.  It is more precise than seeking another isolated law or
another static common-section certificate.

## Updated concept board

1. **Two-axis selector barcode:** exact birth/death cells in `(D,d)`.
2. **Compressed dual jet:** every current active row uses only `(ell,m_1)`.
3. **Structured crossing loci:** one depth-three crosser, two depth-five
   crossers of the previous wall, then a 24-state low/high depth-six patch.
4. **Finite full bank:** only 4,319 nonempty physical states exist.
5. **Structured stopping laws:** abstract Markov realization is automatic;
   independent-deletion and value/depth-hazard models test real locality.
6. **Original-response map:** still absent; no selector theorem alone proves
   NC2 or GMC.
7. **Positive chain map:** prefix stopping flow and partition Hasse flow are
   different objects; their boundary-compatible intertwiner is the real
   sequential sidecar.

The creative lesson is to study the *motion of the dual face* and the *small
crossing locus* together.  Each alone looks unstable.  Their paired evolution
is rigid enough to suggest a recurrence, and finite enough to settle exactly.
