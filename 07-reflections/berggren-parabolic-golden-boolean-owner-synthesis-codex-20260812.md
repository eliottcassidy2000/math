# Parabolic cusp, golden ray, Boolean planes, and the missing owner

**Synthesis date:** 2026-08-12.  This reflection routes proved objects and
open transfers; theorem files, not this narrative, are the truth sources.

## Current board

| lane | exact carrier | status | theorem boundary |
|---|---|---|---|
| parabolic branch | consecutive parameters `(t,t+1)` and the fixed Farey cusp `(1,1)` | **PROVED** | [THM-3334](../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md) |
| golden branch | Fibonacci parameters, three ancestry rays, affine cocycle/current no-go | **PROVED** | [THM-3339](../01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md) |
| equal-plane collisions | Gaussian factor-choice torsors `F_2^r/<1>` | **PROVED** | THM-3334 |
| internal four-box | `W(m,n)=(n-m,m,n,n+m)` with a signed matching current | **PROVED** | THM-3339 |
| square/triangular selector | intended Pell/Markov compiler | **RESERVED / OPEN** | [THM-3335](../01-canon/theorems/THM-3335-square-triangular-pell-markov-pythagorean-selector.md) is an empty stub |
| Gaussian content curvature | intended multiplication/content cocycle | **RESERVED / OPEN** | [THM-3336](../01-canon/theorems/THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md) is an empty stub |
| LRC spectral closure | owner, phase, clock, endpoint word, global exit | **OPEN** | neither THM-3334 nor THM-3339 removes a row |

## One interval, two genuinely different rays

Both constructions live on reduced fractions `m/n` in `(0,1)`, but their
dynamics have different Jordan type.

The THM-3334 spine uses

```text
g=[2 -1;1 0],          g(t,t+1)=(t+1,t+2),
(g-I)^2=0.
```

It is parabolic.  Its fractions `t/(t+1)` approach the rational fixed cusp
`1`, and consecutive vertices form a fan of Farey triangles about `(1,1)`.
Gaussian squaring turns linear parameter drift into quadratic triple growth.

The THM-3339 ray uses

```text
G=[0 1;1 1],           G(F_k,F_(k+1))=(F_(k+1),F_(k+2)).
```

It is hyperbolic up to orientation, with eigenvalues `phi` and `-1/phi`.
Its fractions approach the quadratic irrational `1/phi`.  Three Fibonacci
steps equal two Berggren parameter steps,

```text
G^3=AB,
```

except that the odd/odd residue lane must first pass through the primitive
normalizer `T`, where `T G^3 T^(-1)=CB`.  This produces the three ancestry
rays

```text
(BA)^r,             A(BA)^r,             C(BC)^r.
```

The parabolic and golden loci meet only at `(1,2)` and `(2,3)`, giving
`(3,4,5)` and `(5,12,13)`.  The next parameters `(3,4)` and `(3,5)` already
separate.  The initial coincidence is not a common infinite branch.

## The four-versus-six puzzle is a typed diagram

There are two unrelated `K4` vertex sets.

1. The **internal** `K4` has the four entries of one window.  Its six edge
   products encode one Pythagorean triple only through asymmetric operations:
   one edge, twice an opposite edge, a sum, and a signed difference.
2. The **external** `K4` at `c=1105` has four distinct equal-hypotenuse
   parents.  Its six edges are coloured by the three Gaussian prime-XOR
   directions.

Both have three perfect matchings because every `K4` does.  That shared
incidence does not identify their vertices, operations, or ancestry.

On the Fibonacci window three finite shadows coexist:

- the four vertex values give the same transitive `T4` for every strict
  window and therefore forget Cassini sign, content, and ancestry ray;
- the six edge products give a transitive `T6` whose middle arc alternates;
- oriented Farey reduction gives a `C6` of all six transitive `T3` orders on
  the three matching channels.

The alternating `T6` arc is the isolated edge swap `(03 12)`.  It is odd,
whereas every edge permutation induced by `S4` is even.  So it is a ranking
sidecar, not a relabeling of the underlying `K4`.  At the quotient level each
adjacent matching reflection has four lifts.  Exactly four lift pairs fix an
owner, one for each vertex; the other twelve generate transitive `S4`.  The
six-state loop sees none of this `V4` coordinate.

## The two appearances of `-4` point in different directions

For the parabolic spine,

```text
c_t=t^2+(t+1)^2=2t^2+2t+1,        disc_t(c_t)=-4.
```

This fixed quadratic discriminant says that a prime divides some `c_t`
exactly when `-1` is a square modulo that prime.  CRT then builds arbitrarily
large Gaussian factor-choice ranks and hence unbounded equal-plane fibres.

For the known three-dimensional Keller counterexample, THM-1310 has

```text
disc_x(N)=-4 Q^2 L.
```

There the variable square class `-L`, not the visible constant `-4`, carries
the cubic `S3` monodromy.  The square factor `-4Q^2` cannot classify the
counterexample family or identify it with the Pythagorean fibres.  The golden
Pell form has discriminant `5`, adding a third, distinct arithmetic carrier.

## What the LRC lane actually inherits

THM-3334 contributes a canonical fixed-cusp Farey fan, exact Lorentz lifts,
and large Boolean side fibres.  THM-3339 contributes a branch word, a
six-state channel order, and a complete finite owner-lift obstruction.  These
are useful coordinates, but the LRC target needs a physical signed owner and
phase on one typed ancestry base.  Reduction modulo two retains the order of
three channels and destroys the integral Gram, height, endpoint history, and
owner.  The exact theta-slaved `r=5` contraction is positive evidence on a
typed row, not a bridge that restores those lost coordinates.

The branch-word question now has an exact answer.  After a lawful channel
frame and one displacement calibration, the branch matrices define a
`V4`-valued cocycle with affine image `D4`.  Its owner sequence is
`0,p,p,r,q,q`, and a moving gauge flattens it.  But the full signed current
has trivial translation stabilizer: every nonzero correction changes either
`(a,b)` to `(b/2,2a)`, flips the antisymmetric `c`, or both.  The current's
four signatures decode the missing translation rather than erase it.

The next lawful question is therefore not “is the puzzle a tournament?” or
“can the owner be flattened?” It is:

```text
Can the 24-state residue-order x signed-current-orbit bundle be transported
to one physical LRC ancestry base without forgetting owner, phase, or sign?
```

A static gauge is insufficient, and THM-3339 proves that even the derived
branch correction cannot preserve the fixed current while flattening the
owner.  A physical use must retain the current orbit as part of its state.

## Research queue opened by the synthesis

1. Classify square, triangular, and square-triangular values on `c_t` and
   `Q_t=2c_t+1`; separate genuine Pell selectors from finite coincidences.
2. Derive the Gaussian multiplication content cocycle before quotienting by
   primitive normalization; determine whether it acts on, merely grades, or
   destroys the fixed-hypotenuse Boolean torsor.
3. Build the exact `6*4` residue/order-current bundle, including the odd/odd
   leg-order normalizer, and test whether its `D4` holonomy has a physical
   LRC owner/phase interpretation.  The flatten-and-preserve route is closed.
4. Test whether depth dispersion inside the rank-`r` Gaussian fibre is
   unbounded and whether any prime-XOR toggle has a finite Berggren word
   transducer.  The scalar plane label alone cannot see this.
5. For LRC, retain the fixed-cusp integral ancestor and ask for one explicit
   owner/phase/current transport on the canonical typed row.  The Boolean
   rank or residue hexagon by itself is not a spectral-closure certificate.

The useful meta-pattern from incoming THM-3340 is methodological only: a
coarse quotient defect may be repaired by a delayed donor when an explicit
injection proves capacity.  In the present lane the odd/odd normalization
ray is a tempting donor, but THM-3339's stabilizer no-go says no donor can
flatten the owner while leaving the same signed current fixed.  A different
target bundle would need its own capacity theorem.  This analogy is not a
theorem dependency.
