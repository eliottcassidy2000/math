# Planar Jacobian counterexample design after three false survivors

**Status:** research synthesis.  The obstruction theorems cited below are
proved in their stated scopes; the proposed construction lanes are OPEN.
Nothing in this note is a counterexample to the characteristic-zero planar
Jacobian conjecture.

The useful outcome of this session was not a candidate pair.  It was a more
adversarial specification of what a candidate must retain.  Three objects
passed the usual cheap plausibility tests:

- a many-term, positive-area, high-genus first coordinate;
- a quotient of a genuine three-dimensional Keller collision; and
- polynomial graphs interpolating that collision.

All three fail for different structural reasons.  Together they say that
degree, support size, fibre genus, collision data, and long formal response
jets are weak evidence unless the construction also retains an exact second
Hamiltonian channel.

## 1. Inheritance pass

### Closest proved mechanisms

- THM-2045, THM-2063, THM-2071 and THM-2118 close broad low-fibre and
  factorized first-coordinate families.
- THM-2102, THM-2127, THM-2132 and THM-2134 force proper-power top faces and
  isolate cross-factor compatibility as the remaining Newton problem.
- THM-2230 quotients mates by target shears; it does not manufacture a mate.
- THM-3543 proves that the full torus-invariant quotient of the THM-1300
  collision pays a ramification square.
- THM-3544 forces every nonzero target-pencil member of a hypothetical
  nonautomorphic planar Keller pair to have total degree at least six.
- THM-3545 turns the separated self-intersection thickening into a
  nonterminating Catalan series.
- THM-3546 gives exact coordinate-graph descent gates from dimension three.
- THM-3548 separates conductance magnitudes from the missing phase and Segre
  sidecars; THM-3549 closes one-sided, total-degree-at-most-three, and every
  transverse-affine correction of the ramified quotient.
- THM-3550 excludes prime target-pencil degrees and raises the pencil-height
  floor to eight.
- THM-3551 closes three all-degree one-invariant response families.
- THM-3552 closes a two-channel cyclic family by holomorphic nonexactness.
- THM-3553 closes the most literal polynomial graph chart of THM-1300 in all
  degrees.
- THM-3554 identifies the natural curved collision surface as
  `G_m x A^1`, with an exact degree-two etale Kummer restriction whose
  affine-plane filling must ramify.  THM-3555 identifies the Catalan branch
  with the universal cubic marked-root cover and proves that a correction
  fixing its ramification line pointwise cannot be Keller.

### Canonical hostile and corrected near misses

The canonical positive collision source is still the explicit THM-1300 map:
it has constant three-dimensional Jacobian and a three-point collision.  Its
torus quotient is a hostile, not a descent: it remembers a collision but
forgets the transverse character, and the planar Jacobian acquires the
square of the missing scale.

The initial quotient probe also overclaimed its factor-stripping search.  It
enumerated only powers of the ramification coordinate `gamma`.  Factoring

```text
P=gamma A,             Q=gamma^2 u B
```

and stripping `A` and `gamma^2 B` instead recovers the affine pair
`(gamma,u)` with Jacobian `-1`--but precisely by destroying the collision.
Thus the correct statement is a full-*invariant-quotient* obstruction, not
an obstruction to every polynomial recombination.

The apparent hyperelliptic scaffold

```text
x+(xy)^a+lambda x^(2a-1)y^(2a)
```

was another coordinate illusion.  With `T=x^(a-1)y^a` it is just
`x(1+T+lambda T^2)`, so THM-3551 kills it in one weight sector.  Its high
genus does not create an independent response direction.

### Least-used sidecar

THM-3064's sheetwise cofactor decoder remains underused.  It can certify a
supplied pointed pair and detect tame valuation, but it needs ownership and
does not reconstruct the missing transverse coordinate.  Combined with
THM-3546, it suggests searching for a coordinate hypersurface that preserves
both collision ownership and the cofactor direction, rather than quotienting
either away.

## 2. Assumptions that failed

| tempting assumption | hostile verdict | first failed implication |
|---|---|---|
| `(P_x,P_y)=(1)` is nearly a mate | false | submersion is only local solvability of the Hamiltonian vector field |
| more monomials or positive Newton area means more channels | false | powers of one primitive invariant occupy one exact weight/residue sector |
| high generic-fibre genus helps produce noninjectivity | usually backwards | the forced response may become a nonzero holomorphic, hence nonexact, differential |
| a collision surviving a quotient is enough | false | the quotient can lose scale and acquire ramification in its Jacobian |
| collision interpolation makes graph descent plausible | false in the fixed chart | the top tangential minor survives every interpolation constraint |
| a long formal solution should terminate after a clever perturbation | false as evidence | one-ray responses and the Catalan branch have arbitrarily long nonterminating jets |
| filling the natural punctured collision surface should be harmless | false | the missing unit is the Kummer square root, and every affine-plane filling of that quadratic field extension ramifies |
| nonconstant Jacobian factors can cancel under polynomial composition | false | `J(B o A)=(JB o A)JA`; a constant unit product in a polynomial ring makes both factors units |
| Laurent reciprocal cancellation can be cleared later | false in basic models | clearing the pole reinstates a nonconstant Jacobian divisor |
| internal cyclic monodromy is map monodromy | category error | the cyclic cover can be the fibre normalization of a coordinate with no mate at all |

The practical correction is to label every signal by what it actually
measures: gradient unimodularity, formal solvability, fibre topology,
collision retention, response exactness, or global polynomiality.  None may
stand in for another.

## 3. Exact hostile objects from this session

### 3.1 One invariant, arbitrarily many terms: no mate

THM-3551 proves all-degree no-go statements for

```text
x phi(x^p y^q),       x+h(xy),       x+h(x^p y^q).
```

The first and third have a unique constant-producing sector with a terminal
nonzero leading coefficient.  The diagonal family is stronger: in generic
coordinates `z=P,v=xy`, a mate would require

```text
partial_v Q=kappa/(z-h(v)),
```

whose nonzero simple residues forbid even a rational primitive.

### 3.2 A four-vertex genus-ten fibre: still no rational mate

THM-3552 takes `T=x^2y^5` and

```text
P=x(1+T+T^2)+T+T^2/2+T^3/3.
```

It has unit gradient ideal, Newton area `20`, and generic genus `10`.  Yet on
the normalized fibre `P=z`, every mate would force

```text
dQ=y dT/[5(1+T+T^2)T].
```

That differential is nonzero and holomorphic, with divisor orders

```text
0; 1,1; 2,2,2; 10.
```

So it cannot be exact.  This is a particularly useful hostile because it
survives the gradient, support, degree, and topology filters simultaneously.

### 3.3 The obvious graph descent closes in every degree

For the THM-1300 map and graphs `z=h(x,y)`, THM-3553 proves

```text
top Jac(F1|_h,F2|_h)=-3S partial_x S,
S=x^3y^2 h_d.
```

It is nonzero for every positive-degree leading form `h_d`; constants retain
the coefficient `-89`.  Affine, quadratic, and cubic collision interpolation
were therefore not unlucky searches.  The whole fixed chart is closed.

### 3.4 Small ansatzes that collapse before geometry

Exact hostiles also closed the following common constructions:

- independent crossed shears give `J=1-f'(y)g'(x)` and are forced one-sided;
- compositions cannot cancel nonunit polynomial Jacobian factors;
- a shared inner polynomial evolves affinely and conjugates to a triangular
  shear;
- symmetric nilpotent-Hessian attempts through quartic order collapse to an
  isotropic one-axis shear;
- radial multiplicative pairs `(xA(xy),yB(xy))` force `AB` constant, while
  the Laurent escape pays a pole.

These are stopping reasons, not evidence that a degree threshold alone is
the missing ingredient.

## 4. Live concept board

| lane | object | preserved datum | likely loss | cheapest decisive test |
|---|---|---|---|---|
| Anchor | coordinate-hypersurface descent of THM-1300 | genuine Keller determinant and collision | transverse cofactor or a puncture under chart choice | coordinate test + tangential minor + Kummer divisor + pencil floor |
| Niche | mixed cubic-root-cover surgery | explicit three-sheet fibre and formal constant-J branch | polynomial termination or ramification | move the ramification line, then run defect-by-defect Groebner/Hensel recurrence |
| Niche | two nonparallel toric invariants | two response sectors and Newton faces | gradient unitness or exactness at infinity | gradient resultant, weight-face gate, fibre differential divisor |
| Wildcard | cross-factor compatibility scheduling | proper-power top face with several irreducible owners | eventual common approximate root | staged compatibility ranks and factor-invoice ledger |
| Wildcard | Laurent cancellation plus boundary modification | exact symplectic identity off a divisor | polynomiality at the divisor | valuation/residue audit on a compactification |

The board deliberately separates representation from invariant.  A graph,
quotient, Newton polygon, and fibre curve are not competing pictures of the
same information: each forgets different coordinates.

## 5. Proposed construction lane A: change the graph chart, not the degree

### Contract

- **Source:** the THM-1300 map `F:A^3->A^3` and one of its collision pairs.
- **Target:** a planar restriction to a polynomial coordinate hypersurface
  `rho_s=0`, mapped into a coordinate target hypersurface `rho_t=0`.
- **Map:** `rho_s | rho_t(F)`, followed by restriction and two target
  coordinates.
- **Preserved predicate:** the ambient Keller determinant and the chosen
  collision, provided both points lie on `rho_s=0`.
- **Destroyed information:** potentially the normal cofactor direction.
- **Needed sidecar:** coordinate verification for `rho_s,rho_t` and the
  sheetwise cofactor/tangential minor.

THM-3553 says increasing the degree of `h` in the fixed equation `z=h(x,y)`
cannot work.  The next move is therefore a polynomial source coordinate
change, a different target coordinate pair, or a nongraph coordinate
hypersurface--not degree four in the same chart.

The natural nongraph surface is no longer a candidate either.  THM-3554
puts the curved component `C=0` into the exact form

```text
G_m x A^1 -> A^1 x G_m,          (s,b) -> (b,4s^2).     (1)
```

It is etale and retains the collision precisely because `s` is a
nonconstant unit.  Filling `s=0` restores the quadratic branch divisor, and
no etale `A^2` filling can induce the same function-field extension.  Thus a
live hypersurface must change the extension, not merely compactify this one.

### Cheapest search

Enumerate tame coordinate changes of small length around the three collision
points.  For each, require:

1. both source points satisfy `rho_s=0`;
2. the common target satisfies `rho_t=0`;
3. `rho_s` divides `rho_t(F)`;
4. both `rho_s` and `rho_t` pass a coordinate test;
5. the induced tangential Jacobian is constant;
6. every target-pencil member meets the degree-six floor.

A leading-form obstruction should be extracted before a coefficient search.
The fixed-chart identity in THM-3553 is the model: if every chart in a tame
orbit has a nonzero top cofactor, close the orbit rather than extending the
degree cap.

## 6. Proposed construction lane B: mixed cubic-root-cover surgery

THM-3545 solves the separated collision thickening

```text
P=v^2+A(w),             Q=v^3-v+vB(w)
```

exactly.  Constant Jacobian forces an infinite Catalan series.  The live
repair is not another choice of `A`; it is to add mixed corrections

```text
P=v^2+A(w)+C(v,w),
Q=v^3-v+vB(w)+D(v,w),                                  (1)
```

with `C,D` vanishing at the chosen collision points.

THM-3555 now explains the Catalan square root geometrically.  After adjoining
`r=sqrt(1-3 kappa w)`, the map is affinely equivalent to the universal
marked-root cover

```text
(t,p) -> (p,-t^3-pt).                                  (2)
```

The formal constant Jacobian cancels the cubic cover's ramification factor
with `dr/dw`.  A polynomial correction that fixes the ramification line
pointwise still has zero Jacobian at its cusp.  Therefore the mixed terms in
`(1)` must move that line already at order zero while retaining a chosen
simple three-point fibre.  They may redirect the Catalan defect into a finite
nilpotent block, or merely hide the same algebraic tail.  The decisive
observable is the transition operator on successive defects, not the number
of solvable jets.

### Cheapest search

- impose the collision and constant-J equations weight by weight;
- require the correction to move the cubic ramification line, not merely its
  higher jets;
- quotient target shears immediately using THM-2230;
- compute the affine space of corrections at each defect;
- search for a periodic or nilpotent transition on the obstruction module;
- reject a lane if every branch has a nonzero scalar Catalan component.

A positive signal is a finite invariant subspace on which the defect operator
is nilpotent.  A solution through degree `N` without that structure is only a
formal near miss.

## 7. Proposed construction lane C: two nonparallel invariants

Take primitive monomials

```text
T=x^a y^b,             S=x^c y^d,       ad-bc!=0,        (2)
```

preferably with `|ad-bc|=1`, and seek

```text
P=x Phi(T,S)+Psi(T,S).                                  (3)
```

This is the smallest honest escape from THM-3551 and THM-3552: `S` cannot be
a power or rational function of `T`.  The unimodular determinant minimizes
ramification introduced by changing toric coordinates.

### A genuinely nonparallel hostile

The first exact implementation still fails, but for a new reason.  For
`m>=2`, put

```text
T=x^(m-1)y^(2m-1),             S=x^m y^(2m+1),
P_m=x[1-4mT-(8m^2/(2m-1))T^2+S(1+T)].                 (3)
```

The primitive exponent determinant is `-1`, and an exact support test rules
out every representation `x phi(W)+Psi(W)` with monomial `W`.  For `m=4`,

```text
P=x-16x^4y^7-(128/7)x^7y^14+x^5y^9+x^8y^16.          (4)
```

It has Newton area `3` and the polynomial Bezout certificate

```text
[1-Cx(9+16T)/9]P_x+[Cy(5+8T)/9]P_y=1,
C=(P_x-1)/x.                                           (5)
```

Thus it is a genuine two-direction polynomial submersion.  In the
unimodular Laurent coordinates `(T,S)`, however, its generic fibre is

```text
(1+T)S^8+(1-16T-(128/7)T^2)S^7-zT^9=0.               (6)
```

The Newton polygon has four interior points and generic genus `4` (excluding
`z=0` and the face-degeneracy value `z=-16^8/7`).  A mate
would force

```text
dQ=T^3S^3 dT/F_S,                                      (7)
```

the nonzero holomorphic Newton adjoint indexed by the interior point `(4,4)`.
It is not exact, so even a rational mate is impossible.  Exact coefficient
systems remain inconsistent through mate degree `24`, where the kernel has
already grown from constants to `<1,P>`.

This is `SCRATCH-PROVED + FINITE-EXACT + INDEPENDENTLY DERIVATION-AUDITED`,
not yet a canonical theorem: the companion still needs explicit four-face
and connectedness gates before promotion.  It is genuinely two-directional
but only affine-linear in `S`, with two parallel `T` rows, so it is not a
general mixed two-invariant theorem.  It sharpens rather than invalidates the
lane.  Nonparallel invariants must make the response cohomology classes
cancel; a unimodular exponent matrix alone is not enough.

### Hostile-first protocol

1. Compute `Res(P_x,P_y)` or a Groebner basis before solving for a mate.
2. Check every positive-weight top face against the proper-power walls.
3. Decompose the mate equation into the two-dimensional charge lattice; do
   not project it to one residue class.
4. On the generic fibre, derive the actual response differential and audit
   all poles, residues, and holomorphic components.
5. Only then solve bounded coefficient systems for `Q`.

The failure modes are crisp.  If the response is holomorphic and nonzero,
THM-3552's exactness wall returns.  If it has a nonzero residue, THM-3551's
rational wall returns.  A genuine positive signal is a differential of the
second kind with zero residues whose cohomology class can be cancelled by
the `S`-channel.

## 8. Proposed construction lane D: alternate cross-factor invoices

The remaining multi-face gap in THM-2134 suggests starting with several
smooth factors rather than one invariant.  An exact probe used weights
`(2,3)` and

```text
f=y^2-x^3,       g=y^2-2x^3,       h=f^2g^3,
P0=h^2,          Q0=h^3.                                  (4)
```

At the first correction, compatibility forced only `g|S`, not `h|S`; the
next defect had rank `15` and augmented rank `16`, killing the simplest seed.
This is a useful negative signal: one stage invoices one irreducible factor,
and the next invoices the debt it left.

The live variant uses three factors or a different multiplicity vector so
successive compatibility conditions can alternate owners without ever
forcing a single global approximate root.  Record at each stage:

- which factor is invoiced;
- the solution-space dimension;
- divisibility of the residual by each factor;
- whether the owner sequence becomes periodic;
- whether a target shear removes the apparent freedom.

This lane is FINITE-EXACT exploratory.  A few successful stages do not imply
a polynomial Keller pair.

## 9. Decision rules for future candidates

A proposed planar counterexample should be rejected or promoted using the
following order.

1. **Refactor first.**  Identify primitive invariants and remove fake
   multi-face structure.
2. **Gradient gate.**  Prove `(P_x,P_y)=(1)` exactly; do not infer it from
   numerical sampling.
3. **Target-pencil gate.**  Enforce composite degree at least six in every
   nonzero direction and pencil height at least eight.
4. **Response-sector gate.**  Identify every sector capable of producing the
   scalar Jacobian term.
5. **Compactified exactness gate.**  Compute the divisor and cohomology class
   of the forced fibre differential.
6. **Collision gate.**  Verify an actual pair of distinct source points, not
   only a singular or contracted quotient.
7. **Polynomiality gate.**  Audit poles introduced by every division,
   factor stripping, or Laurent conjugacy.
8. **Only then increase a degree cap.**

The most promising object is no longer “a high-degree polynomial with many
terms.”  It is a construction that can name two independent response
channels and explain, before coefficient search, how their pole divisors or
cohomology classes cancel without losing the collision.

## 10. Bottom line

No planar characteristic-zero counterexample was found.  The session did
produce three exact all-degree walls and a narrower target:

```text
one invariant          -> one terminal response sector;
controlled second term -> holomorphic but nonexact response;
fixed graph descent    -> nonzero top tangential minor.
```

The honest frontier is the intersection of their complements: a construction
with genuinely nonparallel invariants, a chart that retains transverse
cofactor data, and mixed corrections whose response class cancels for a
structural reason rather than for many initial coefficients.
