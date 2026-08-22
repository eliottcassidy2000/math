# The missing level-three coordinate was a separator, not another norm

**Status: independent audit reflection for proved THM-3519.**  The theorem
and exact companion are the proof sources; this note records the mechanism,
the hostile, and the boundary of the result.

## Verdict

The incoming level-three candidate is sound.  The target `(77,62,4)` is a
lawful good-reduction witness over `F_101`; its first two inverse levels split
as `3` and `3+3+3`, and the last level is the product of nine separable cubic
quotients.  In that rank-27 algebra every actual source coordinate generates:

```text
coordinate       power rank       power determinant mod 101
x                     27                       25
y                     27                       64
z                     27                        1
```

The three degree-27 characteristic-polynomial products are independently
squarefree and reproduce the candidate hashes exactly.

## Why the inverse formula is now transparent

The clean-room audit does not treat the quotient formula as an oracle.  Put
`r=1+xy`; eliminating `z` from the map reduces the first two target equations
to two small relations in `r`.  With

```text
K=2+(9ac-b)x,
Q=3c(1+bx)-4x,
D=(12a-b^2)x^2+bx+2,
g=Lx^3+(4-3bc)x-2c,
```

the central identity is

```text
DQ-xK^2=-3g.
```

On the cubic quotient this turns the two relations into the displayed
rational section for `y`, while the third map coordinate gives `z`.  The
symbolic residuals are multiples of `g`; the finite quotient code then
inverts `D` and `x^3` by unrelated matrix solves and replays the map.  The
formula has both a generic derivation and a hostile implementation check.

## The parent target explains why the new target matters

At `(93,28,83)`, every local cubic coordinate polynomial is squarefree, but
the global `y` product is not.  Exactly two blocks share `T+76`, hence the
same geometric `y` value `25`, and the `y` power rank is exactly `26`.

That is a useful boundary case.  Local primitivity in all nine cubic blocks
does not imply primitivity in their product: the coordinate must also
separate different blocks.  The global derivative gcd detects this, while
the full power matrix measures the lost dimension directly.

## Why one modular fibre is enough

The generic power determinant is a rational function on the open target
chart where the monic cubic tower and inverse sections exist.  Every relevant
leading coefficient and denominator is a unit at the new mod-101 point.
Its nonzero reduction therefore proves that the cleared characteristic-zero
numerator is not the zero polynomial.  This is an open-condition argument;
it neither lifts the finite point to a preferred rational target nor claims
that every specialization is primitive.

## The common class is trace geometry, not divisor numerology

Once `x,y,z` are all primitive, their power bases are three bases of the same
degree-27 trace space.  Their Gram determinants differ by squares of basis
change determinants.  THM-3495 already gives the exact `x` class `[-2J]`, so
the other two inherit precisely that class.

The constant `[-2]` is load-bearing.  Dropping it because `J` is the only odd
divisor would repeat MISTAKE-413.  The proof transports the complete class by
literal squares; it never infers the unit from divisor parity.

This also clarifies the three level-one cubic discriminants.  Their common
`-4` factor is a square times the sign unit, and primitivity places the three
coordinate views in one trace-form class.  The fixed tower now reads

```text
level 1, degree  3:  [-L]       (THM-2546)
level 2, degree  9:  [ H]       (THM-3508)
level 3, degree 27:  [-2J]      (THM-3519)
```

This is a dynamical atlas for one fixed Keller map, not an exact
classification of all dimension-three counterexamples.  Equality of square
classes retains only the sign character of the root permutation action; it
forgets the roots, block system, full wreath monodromy, multiplicities, and
the rest of the nonproperness tower.

## What to test next

The next fixed-map question is level-four `y,z` primitivity.  THM-3504 gives
the new image prime and degree-81 `x` information, but the cheapest clean
test is again a simultaneous good-reduction separator, not a discriminant
expansion.  The parent hostile suggests the correct search statistic:

```text
local cubic generation
    + pairwise cross-block coordinate separation
    + full 81-column power rank.
```

For the outside infinite family, none of these fixed-map targets transports
automatically.  A family-level theorem would need a common finite-flat model,
one explicit separator minor as a function of the family parameter, and a
separate audit of its constant discriminant unit.

## Exact artifacts

```text
04-computation/keller_level_three_three_coordinate_primitive_independent_audit_20260816.py
05-knowledge/results/keller_level_three_three_coordinate_primitive_independent_audit_20260816.out
01-canon/theorems/THM-3519-level-three-sporadic-keller-three-coordinate-primitivity-and-common-discriminant-class.md
```

The LF-normalized script/output hashes are

```text
10520f0314a521ca479c0971916b52794084cc287f446bd05bc4eb0ab55e75e2
859104c1aee350220b54142d8be55270ddbb3630e047ce492bf479d3df2e1e91
```

and the semantic digest is

```text
51b5945a5f26c3335b6837b169a916d9562acd8448d2dde4f1fe1449c1e157a3.
```

The boundary remains fixed-map and level-three: no exact discriminant
multiplicities, level-four or all-level primitivity, arbitrary Keller result,
outside-family classification, `JC(2)`, or LRC consequence.
