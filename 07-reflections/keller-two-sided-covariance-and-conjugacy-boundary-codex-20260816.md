# Fixed Keller covariance: the commuting square and the inserted-map boundary

**Historical research reflection.  Current truth is THM-3532, not this file.**

## Inheritance pass

The closest proved mechanism was THM-3528--3531: a fixed raw norm-packet
tower, a finite-sheet unit, prime image edges, an exact Jelonek component
tower, and its intrinsic discriminant square class.  The incoming THM-3533
strengthened parity at the newest prime to normalization-discriminant
multiplicity one and exposed the local index square.

The canonical hostile was MISTAKE-413: rational scalar units cannot be erased
from a field square class.  The corrected near miss was MISTAKE-415: exposed
norm faces do not by themselves make a closed recurrence.  The least-used
sidecar was the coordinate atlas itself: the five weights are filtrations on
a chosen chart, not intrinsic functions of a polynomial self-map.

The live board was:

1. naturality of the cubic field norm;
2. the finite inverse divisor and its unit valuation;
3. prime-image ancestry and generic degree;
4. trace discriminants, integral closure, and the reduced newest different;
5. self-iteration versus left/right equivalence;
6. standard weights and raw scalar normalization.

## The organizing object

For `G=tau o F o sigma`, the useful object is not a comparison of expanded
formulas.  It is the extension isomorphism

```text
K_t  --F^*-->  K_s
 |              |
tau_#         (P -> P o sigma)
 |              |
 v              v
K_t  --G^*-->  K_s.
```

Norms are determinants of multiplication, trace discriminants are
determinants of trace Gram matrices, and integral closure is characterized by
integrality.  All three constructions commute with this square.  Prime
divisors and generic degrees commute because both vertical maps are
isomorphisms.  Properness commutes because polynomial automorphisms and their
inverses are proper.  This packages every one-step covariance statement in
one proof rather than five unrelated substitutions.

| source | target | map | preserved predicate | lost data | required sidecar |
|---|---|---|---|---|---|
| old source equation `P` | new source equation | `P -> P o sigma` | norm input, finite-sheet order | standard source weights | transported source chart |
| old target equation `Q` | new target equation | `Q -> Q o tau^-1` | cleared norm, image prime, newest different | primitive content, standard target weights | transported target chart and scalar |
| `S_F` | `S_G` | `tau` | nonproperness set and component geometry | literal component equations | chosen equations `tau_#P_j` |
| old field extension | new field extension | commuting square above | norm, trace, integral closure | literal primitive coordinate | transported generator/order |

The cheapest decisive tests were correspondingly small: exact inverse checks
for `T_1,T_2`, one pulled cubic discriminant, five weight signatures, and one
second-iterate witness for each W-map.  The finite bank audits the
implementation; it is not the proof of naturality.

## Why one step does not automatically become a tower

The obstruction is visible before any elimination:

```text
(tau F sigma)^r=tau F(sigma tau)F ... (sigma tau)F sigma.
```

The inserted `sigma tau` changes the next input.  In old coordinates the fed
packet operator is

```text
P -> (sigma tau)_#(L^e N_F(P)),
```

not the original raw operator.  Honest conjugacy, `sigma=tau^-1`, kills the
insertion identically and therefore transports every rung and every iterate.
A different two-sided presentation could still work if it has an independently
proved symmetry or intertwiner for the insertion.  Thus conjugacy is the
universal assumption-free boundary, not a classification of every exceptional
equivalence.

## W1/W2 as controls

The maps retained in THM-2465 are

```text
W_1=T_1 F,    T_1(x,y,z)=(x,y,z+x^2),
W_2=T_2 F,    T_2(x,y,z)=(x+yz,y+z^2,z).
```

They are target postcompositions.  Their transformed boundaries, norm
outputs, finite divisor, coordinate-core discriminants, and Jelonek sets are
exact target pullbacks for one application.  They are not conjugates as
self-maps.  Exact hostile witnesses already occur at depth two:

```text
W_1^2(0,0,-1)=(1,-3,0) != (0,0,-2)=T_1F^2(0,0,-1),
W_2^2(-1,0,0)=(62,4,0) != (-2,0,0)=T_2F^2(-1,0,0).
```

Their honest companions `T_i F T_i^-1` transport the complete tower.  This
pair of negative and positive controls is sharper than calling all four maps
“tame variants.”

The standard seed signature also records the chart loss:

```text
L:        (1,-1,-5,2,-8),
L o T_1^-1: (6,-1,-5,2,-8),
L o T_2^-1: (1,-8,-16,8,-40).
```

Transporting the coordinates restores `A(1,0)` exactly.  Total degrees,
multidegrees, term counts, primitive integral contents, and named literal
coordinates are not rescued by that statement.

## Incoming newest-prime result

THM-3533 fits the same square without adding a new computation.  The integral
closure of the old target ring maps to the integral closure of the transformed
target ring, so the reduced newest different and its height-one valuation one
transport.  For a transported integral primitive element, the local index
length is unchanged, hence

```text
v(Disc(theta))=1+2i
```

is covariant too.  A literal coordinate can cease to be the chosen primitive
or locally maximal element; that is destroyed normalization data, not a
failure of the intrinsic theorem.

The later incoming THM-3535 sharpens the observation sidecar.  Every nonzero
constant linear form is primitive at every depth for `F`.  Under conjugacy it
transports as `ell o phi^-1`.  Affine conjugacy keeps these, up to irrelevant
constants, as all standard linear directions; nonlinear conjugacy generally
turns them into polynomial observations.  Thus primitivity itself transports,
while the label “constant-linear coordinate” does not.  The transported local
index is unchanged, but THM-3535 does not determine it.

## Validity-gate repair

The first companion incorrectly used Python `assert` for all truth-bearing
checks.  `python -O` therefore erased the certificate while preserving its
printed transcript.  The repair replaced every such gate by explicit
`require`/raise calls, excluded executable asserts by source inspection, and
replayed normal and optimized modes against the stored output.  This became
part of MISTAKE-419 and of the promotion gate, rather than being hidden as a
mechanical edit.

## Remaining frontier

The theorem classifies transport on one fixed polynomial-conjugacy orbit.  It
does not classify arbitrary Keller maps, distinct tame-equivalence classes,
or maps sharing only a degree or determinant.  The strongest nearby open
questions are now better typed:

- find a genuine non-conjugate intertwiner for the inserted-map recurrence, or
  prove its absence in a specified family;
- compute old-prime positive normalization-discriminant multiplicities;
- decide all-level local maximality for a fixed literal coordinate;
- construct a family parameter space on which packet atlases and scalar
  normalizations are part of the data.

These are classification questions.  W1/W2 alone provide one-step controls,
not evidence that arbitrary tame Keller maps share the fixed tower.
