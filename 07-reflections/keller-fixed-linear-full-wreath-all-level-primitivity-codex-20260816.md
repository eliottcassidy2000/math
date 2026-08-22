# Fixed Keller linear primitivity: the missing coordinate was a block theorem

**Status: proof reflection for PROVED + VERIFIED-EXACT THM-3535.**  The
theorem file is the proof source.  This note records the inheritance pass,
the failed shortcut, and the mechanism that closed the frontier.

## Inheritance and concept board

The anchor was THM-3531's intrinsic all-level discriminant class.  Its
closest proved mechanism was THM-3530's sequence of new absolute image
primes with generic image degree one.  The canonical hostile was THM-3508's
intermediate coordinate: a separable cubic view inside a degree-nine tower
can have a cubed characteristic polynomial and zero degree-nine
discriminant.  The corrected near miss was the fixed-level generic-form
argument in THM-3531: it cannot be iterated by a bare countable-avoidance
claim over `Q`.  The least-used sidecar was the **support** of inertia at the
newest prime, finer than its discriminant parity.

The live board was:

1. newest image prime and unique degree-one ancestry;
2. one supported transposition versus an unlabelled sign character;
3. natural blocks of a full ternary wreath action;
4. one-step descent of a constant linear observation;
5. rational countable avoidance;
6. a finite-field all-direction hostile.

The niche (wreath blocks) overtook the initial coordinate computation.  The
wildcard was the pair of rational collision fibres already present in
THM-2473.

## The closing mechanism

At the generic point of the newest prime `P_(n-1)`, the preceding iterate is
finite etale because all of its nonproperness primes are older.  The
degree-one composite `V(L)->V(P_(n-1))` picks exactly one of its
`3^(n-1)` sections.  On that section the reciprocal cubic

```text
L+Tu^2-2cu^3
```

has one unit root and one tame quadratic pair; every other section is
finite etale.  Thus inertia is not merely odd: it is one transposition
supported in one bottom block.

The group amplification is elementary.  Leaf transitivity makes the
setwise stabilizer of a bottom block transitive on its three leaves.  A
supported transposition therefore grows to the full local `S3`; conjugation
gives every bottom factor.  Induction gives the full iterated wreath group.
Its blocks are exactly ancestor cylinders, so the inverse tower has no
intermediate fields beyond its ancestor chain.

The remaining observation test collapses to one step.  Differences from the
`m=1,2` rational fibres give

```text
(2,-3,0), (4,-3/2,0), (1,-3/2,27/4),
det=243/4.
```

No nonzero constant linear form is constant on both fibres, hence none lies
in the immediate parent field.  With no non-ancestor intermediate field,
every such form is primitive at every depth.

## Connection typing

```text
source:   newest reduced image prime P_(n-1)
map:      strict-henselian inertia -> one bottom transposition
target:   full ternary wreath group -> ancestor block chain
preserved: labelled block support, transitivity, parent field
lost:     affine index, exact discriminant multiplicity, special-fibre data
sidecar:  degree-one V(L) ancestry plus reciprocal Newton polygon
test:     two-fibre difference determinant 243/4
```

This explains why THM-3531's square class alone was insufficient: a sign
character remembers an odd number of transpositions but forgets whether one
transposition is isolated in one block.  THM-3530 restores exactly that
support coordinate.

## Quantifier hostile and exact probe

Enumerating `P^2(Q)={r_1,r_2,...}` and deleting the first `s` points gives a
nested sequence of nonempty rational Zariski opens with empty rational
intersection.  Thus “one open set per depth” never proved a simultaneous
rational form.  THM-3535 bypasses the trap by proving every nonzero constant
direction works.

The companion separately enumerates `F_41^3`.  There are 48 split
degree-nine fibres of `F^2`, and across them every one of the 1,723
directions in `P^2(F_41)` separates a full fibre.  This is a useful hostile
probe, not the all-level proof.

## Boundary

The result is fixed-map and constant-linear.  It proves full-degree
separable minimal polynomials and, with THM-3531, their all-level square
class.  It does not compute their affine order indices or exact old-prime
multiplicities; THM-3533 separately fixes the newest normalization
coefficient at one.  It does not transport to tame conjugates automatically,
classify other Keller maps, or imply properness, invertibility, `JC(2)`,
`DC(2)`, or LRC(14).

The research move is a candidate meta-pattern: **an intrinsic parity class
becomes a full permutation carrier only after a prime-ancestry sidecar
localizes its support to one block**.  It should be promoted to
`META-PATTERNS.md` only after a second, distinct thread exhibits the same
upgrade.
