---
id: THM-4146
title: "Rational three-cycle order-six lift, horizontal-section divisor collision, and fibre firewall"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED; EXTENDS THM-4139; JC(2)
  OPEN. Beyond THM-4139's rational graph, AP rigidity, specific lift,
  zero-fibre/horizontal-section separation and finite-ring controls, every
  exact quadratic three-cycle has a canonical trace-one SL_2 lift with cube
  -I and sixth power I. For x^2-29/16 the lift is the complete six-point
  integral hexagon on X^2+2XY+13Y^2=48, while the signed Pythagorean template
  forces the 3:4:5 similarity class and 29. THM-4134's horizontal section has
  the same translated root divisor at a=-48 but splits it into two q-fibres;
  this is compatible with THM-4139's stronger target-point separation. The
  quadratic sixth dynatomic is irreducible of degree 54 with 18 real roots.
  The order-six action belongs only to the comparison divisor. THM-4138 has
  already closed the degree-16/15 Delta_V wall by independent mechanisms; no
  new monodromy exclusion or planar-Jacobian case is claimed here.
source: codex-planar-jc-dynamics-20260825
depends_on:
  - THM-4139-rational-three-cycle-order-six-lift-and-horizontal-carrier
  - THM-4134-delta-v-collision-wall-strict-transform-and-high-branch-exclusion
related:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-3465-nonreal-cyclic-character-keller-rigidity-and-hfc-separation
  - THM-3742-square-triangular-pell-mod13-central-sign-projective-cycle
  - THM-4057-stern-brocot-depth-pullback-and-rational-edge-tournament-gauge
  - THM-4138-delta-v-horizontal-carrier-monodromy-exclusion
script: 04-computation/rational_three_cycle_order_six_horizontal_carrier_thm4146.py
output: 05-knowledge/results/rational_three_cycle_order_six_horizontal_carrier_thm4146.out
independent_audit_script: 04-computation/rational_three_cycle_order_six_horizontal_carrier_thm4146_independent_audit.py
independent_audit_output: 05-knowledge/results/rational_three_cycle_order_six_horizontal_carrier_thm4146_independent_audit.out
script_sha256: 599bd10756498f27f9c528e71829673a14406c60a4c3bb70ff4ac559eecc4c19
output_sha256: 0f4c26ef092827ce54afcc9cf215e0c350af51554c7c9ac6b6affe44128311e2
independent_audit_script_sha256: 81969e3d145f5fd2a83c794e2871d5d901c5c5b8641f9a6c9a3f11c81e9ef1d1
independent_audit_output_sha256: 0ac6ff15ba0bd55fb4474fa3f298b59ad24391d814e96b4a6f0f8360d9e9c1c0
semantic_sha256: 014f11897d8ca216faa897f117db118a9760b791d393d664289b621d99456d3c
hash_basis: raw LF bytes
audit: >
  PASS. The primary certificate checks the complete bounded rational graph,
  universal trace parameter, two matrix models, invariant norm, complete
  integral level set, AP and Pythagorean forcing, THM-4134 section identity,
  uniform fibre split, a Rabin irreducibility certificate modulo 11, exact
  Sturm count, derivative-degree/Mersenne identities and the full scalar
  mod-63 functional graph. A separate script imports no primary code and
  rebuilds the results using Fraction matrix arithmetic and a different
  exact polynomial-division path. Normal and optimized streams byte-match.
---

# THM-4146 -- universal order-six lift and horizontal fibre firewall

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED; JC(2) OPEN.**

THM-4139 already proves the complete rational graph, the unique AP-supported
cycle, the specific projective lift, and the separation between its zero-fibre
normalization incidence and THM-4134's horizontal section. This extension
proves the universal three-cycle lift, the complete integral hexagon,
Pythagorean forcing, the stronger two-fibre/root-divisor firewall, and the
quadratic map's algebraic period-six census. The scalar map has real algebraic
six-cycles but no affine rational one. Its canonical length-six object over
`Q` is instead the linear lift of its projective three-cycle. Projectivization
retains the model divisor and forgets the central sign of its defining section.

## 1. Exact rational dynamics

Put

```text
f(x)=x^2-29/16.                                          (1)
```

Then

```text
-7/4 -> 5/4 -> -1/4 -> -7/4.                            (2)
```

More strongly, on the affine rational line,

```text
PrePer(f,Q)={m/4:m in {+-1,+-3,+-5,+-7}},               (3)
```

and `(2)` is its unique affine rational cycle. If `f` is instead extended to
`P^1(Q)`, infinity is of course one additional fixed point.

### Proof of the denominator and height gates

Let `x in Q` be preperiodic. For every odd prime `p`, if `v_p(x)<0`, then
the two terms in `(1)` have unequal valuations and

```text
v_p(f(x))=2v_p(x).
```

That valuation falls strictly and can never enter a finite orbit. Thus the
denominator of `x` is a power of two.

At `p=2`, the constant in `(1)` has valuation `-4`. If `v_2(x)<-2`, then
the valuation again doubles downward. If `v_2(x)>-2`, the next valuation is
`-4`, after which it falls. If `v_2(x)=-2`, write `x=u/4` with `u` a
2-adic unit. Since `u^2-29=4 mod 8`, one has
`v_2(f(x))=v_2(u^2-29)-4=-2` exactly. Hence every affine rational
preperiodic point has exact denominator four.

Write `x=y/4`, with `y` odd. The induced integer recurrence is

```text
g(y)=(y^2-29)/4.                                        (4)
```

Odd squares are `1 mod 8`, so `(4)` preserves odd integers. For `|y|>=9`,
`g(y)>|y|`, and the positive orbit then grows strictly. The remaining graph is

```text
y     -7  -5  -3  -1   1   3   5   7
g(y)   5  -1  -5  -7  -7  -5  -1   5.                  (5)
```

This proves `(3)`, the uniqueness of `(2)`, and the absence of a rational
exact six-cycle. QED.

## 2. Every quadratic three-cycle has an order-six determinant-one lift

Let `F_c(X)=X^2+c` over a characteristic-zero field and let

```text
alpha -> beta -> gamma -> alpha                           (6)
```

be a cycle of three distinct points. Put

```text
sigma=alpha+beta+gamma,
p=alpha beta+beta gamma+gamma alpha,
r=alpha beta gamma.                                      (7)
```

Subtracting consecutive recurrence equations and writing

```text
u=alpha+beta,       v=beta+gamma,       w=gamma+alpha
```

gives

```text
u(w-v)=u-w,       v(u-w)=v-u,       w(v-u)=w-v.          (8)
```

The solution `u=v=w` would make the points equal. Solving the other branch
of the first two equations gives

```text
v=-(u+1)/u,             w=-1/(u+1),             uvw=1. (9)
```

Substitution into `(7)` and one recurrence equation yields

```text
c=-(sigma^2+sigma+2),
p=-(sigma^2+2sigma+3),
r=-(sigma^3+2sigma^2+3sigma+1).                       (10)
```

Thus the cycle polynomial is

```text
P_sigma(X)=X^3-sigma X^2-(sigma^2+2sigma+3)X
             +(sigma^3+2sigma^2+3sigma+1),             (11)

disc(P_sigma)=(4sigma^2+6sigma+9)^2.                   (12)
```

Define

```text
M_sigma(X)=((sigma+1)X-(sigma^2+sigma+1))/(X-sigma),   (13)

L_sigma=[[sigma+1, -(sigma^2+sigma+1)],
         [1,                    -sigma]].               (14)
```

Direct division gives

```text
F_c(X)-M_sigma(X)=P_sigma(X)/(X-sigma).                (15)
```

Consequently `M_sigma` agrees with the quadratic map exactly on the cycle
scheme. Moreover

```text
det L_sigma=1,       tr L_sigma=1,
L_sigma^2-L_sigma+I=0,
L_sigma^3=-I,        L_sigma^6=I.                       (16)
```

The induced `PGL_2` element has order three, while the trace-one `SL_2` lift
has order six. The central sign is not decorative: its negative is the other
lift, of order three.

## 3. The specific integral hexagon

For `(2)`, `sigma=-3/4`. In `y=4x` coordinates, `(13)` becomes

```text
M(y)=(y-13)/(y+3),
A=[[1,-13],[1,3]],                 B=A/4 in SL_2(Q).     (17)
```

The exact identities are

```text
det A=16,       tr A=4,       A^3=-64I,
B^2-B+I=0,      B^3=-I,       B^6=I,                    (18)

g(y)-M(y)=((y+7)(y-5)(y+1))/(4(y+3)).                  (19)
```

The six-vector orbit is

```text
(-7,1) -> (-5,-1) -> (2,-2) -> (7,-1)
       -> (5,1) -> (-2,2) -> (-7,1).                   (20)
```

Projective slopes collapse `(20)` to `-7 -> 5 -> -1`. The form

```text
q_B(X,Y)=X^2+2XY+13Y^2=(X+Y)^2+3(2Y)^2                (21)
```

is `B`-invariant, and `(20)` is the complete set of integral solutions of

```text
q_B(X,Y)=48.                                           (22)
```

Indeed `(21)` bounds `|Y|<=2` and `|X+Y|<7`, after which the six-case check is
exact. Thus the lift is a six-vertex Eisenstein hexagon. The requested
arithmetic is structural:

```text
13=29-16,                    disc(q_B)=-48=-3*16.       (23)
```

Homogenize the scalar cycle polynomial:

```text
L(X,Y)=(X+7Y)(X-5Y)(X+Y).                              (24)
```

Then

```text
L(B(X,Y))=-L(X,Y).                                     (25)
```

The three-line divisor is invariant, but its defining cubic section has the
nontrivial central character. More precisely, `B` cyclically permutes the
three factors with multipliers `2,-1,1/2`, whose product is `-1`. This is a
global defining-cubic eigencharacter for the introduced `C_6` action; it is
not yet a local normal-bundle or meridian character on THM-4134. Equation
`(25)` is the exact sidecar that a divisor-only or projective passport loses
inside this comparison model.

## 4. Why the arithmetic progression and `3:4:5` are forced

### 4.1 Unique arithmetic-progression three-cycle

Suppose a rational exact three-cycle of `X^2+c` is an arithmetic progression.
Choose its cyclic order as

```text
m-r -> m+r -> m -> m-r,                    r!=0.        (26)
```

Subtracting recurrence equations gives

```text
r(-4m-1)=0,                 r(2m+r-1)=0.               (27)
```

Therefore

```text
m=-1/4,             r=3/2,             c=-29/16.       (28)
```

So `(2)` is the unique rational quadratic three-cycle whose point set is an
arithmetic progression.

### 4.2 Scale-covariant Pythagorean characterization

Let positive `a,b,h` satisfy `a^2+b^2=h^2`, and ask for

```text
G_D(y)=(y^2-D)/b
```

to have the signed cycle

```text
-(a+b) -> h -> -(b-a) -> -(a+b).                        (29)
```

The second and third arrows give

```text
D=a^2+2b^2-ab.                                          (30)
```

Compatibility with the first arrow is `h=3a-b`. Hence

```text
0=(3a-b)^2-a^2-b^2=2a(4a-3b).                          (31)
```

It follows that

```text
(a,b,h)=(3k,4k,5k),                 D=29k^2.             (32)
```

Conversely `(32)` satisfies `(29)`. At `k=1`, the points are

```text
-(3+4)=-7,          5,          -(4-3)=-1.              (33)
```

This characterizes the `3:4:5` similarity class; it is not a numerical
analogy. It also types the remaining integers:

```text
7=3+4,       16=4^2,       29=5^2+2^2=2^2+3^2+4^2,
29-16=13=2^2+3^2.                                       (33a)
```

## 5. The THM-4134 horizontal-section collision

THM-4134's polynomial horizontal section, subsequently used in THM-4138's
proved exclusion of the old degree-`16/15` residual, is

```text
U=a/2+16rho^2/(9a^2),
V=-rho-64rho^3/(27a^3),
q=a^3/2+rho^2,                                          (34)
```

satisfying

```text
V^2-U^3+(3/4)a^2U+a^3/4=q.                             (35)
```

Set `t=4rho/(3a)`. Then

```text
U=a/2+t^2,           V=-t^3-(3a/4)t.                   (36)
```

At `a=-48`,

```text
V=-t(t-6)(t+6).                                         (37)
```

Under `t=y+1`, the cubic is exactly the negative of `(24)` in affine
coordinates:

```text
t(t^2-36)=(y+7)(y-5)(y+1).                             (38)
```

The value `a=-48` is also the characteristic discriminant of

```text
A_t=[[2,-12],[1,2]],       tr(A_t)^2-4det(A_t)=-48.    (39)
```

Indeed, with `S=[[1,1],[0,1]]`, one has `A_t=SAS^-1`, and its projective
action

```text
t |-> 2(t-6)/(t+2)
```

cycles `-6 -> 6 -> 0 -> -6`. Thus the root-divisor and discriminant collision
is exact in these fixed coordinates. This Möbius action is introduced from the
rational-cycle model on the zero divisor. It is **not** an automorphism of the
horizontal elliptic surface, the BC cover, or a Keller map, and it supplies no
cover monodromy.

### The fibre firewall

The bridge does **not** put the roots on one elliptic fibre. If `a!=0` and
`r^2=-3a/4`, the roots of the cubic in `(36)` have

```text
t=0:       (U,V,q)=(a/2,0,a^3/2),
t=+-r:     (U,V,q)=(-a/4,0,5a^3/64).                  (40)
```

Their fibre difference is `27a^3/64`, which never vanishes. At `a=-48`,

```text
t=+-6:     (U,V,q)=(12,0,-8640),
t=0:       (U,V,q)=(-24,0,-55296).                    (41)
```

The outer pair coincides in the full `(U,V,q)` triple; only its `t`, equivalently
`rho`, sheet sign distinguishes the two roots. Thus `(38)` preserves the
section polynomial, root divisor, AP spacing and discriminant, but destroys
the common-fibre predicate, cyclic dynamics, injectivity and `rho`-sheet sign.

### Compatibility with THM-4139's stronger separation

There are three distinct maps on the same parameter roots, and conflating
them was the original near miss. THM-4139 feeds `t=-6,6,0` into the nodal
zero-fibre normalization

```text
nu_-48(t)=(t^2-48,t(t^2-72)),
```

obtaining `(-12,216),(-12,-216),(-48,0)`, all at `q=0`. The actual horizontal
section instead gives the two-fibre values `(41)`. The Möbius map in `(39)`
cycles only the source root divisor. Thus THM-4139's refutation of target-point
equality and this theorem's parameter-divisor collision are simultaneously
true:

```text
same three source roots; different target maps; no common target fibre;
no descended automorphism and no BC-cover monodromy.                 (41a)
```

Any future reuse of this comparison must retain `a`, `t`/`rho`, `q`, the
chosen normalization and target map, and the actual labelled meridians of the
new carrier under study.

## 6. Exact algebraic and real six-cycles

For `(4)`, form

```text
Phi_6(y)=((g^6(y)-y)(g(y)-y))/
         ((g^3(y)-y)(g^2(y)-y)).                       (42)
```

Its primitive integral numerator `H_54` has degree `54`. Reduction modulo
`11` passes the Rabin certificate

```text
y^(11^54)=y mod H_54,
gcd(H_54,y^(11^27)-y)=gcd(H_54,y^(11^18)-y)=1.         (43)
```

Hence `H_54` is irreducible over `Q`. It is coprime to the period factors for
`1,2,3`, and exact Sturm counting gives eighteen real roots. Consequently
`(1)` has

```text
three real exact six-cycles,
six nonreal exact six-cycles,
zero rational exact six-cycles.                         (44)
```

No numerical approximation is load-bearing.

## 7. Why `63` belongs to the lift, not to a scalar six-cycle

There are four exact appearances of `64=2^6`:

```text
(-2)*4*(-8)=64             (the y-cycle pair sums),
L(-3,1)=64,
A^3=-64I,
deg((g^6)')=2^6-1=63.                                  (45)
```

The last follows from

```text
(g^6)'=2^(-6) product_(j=0)^5 g^j,                     (46)
```

whose factor degrees sum to `1+2+4+8+16+32=63`. Also

```text
63=2^6-1=3^2*7,       ord_3(2)=2,       ord_7(2)=3.    (47)
```

Thus exponent six introduces no primitive prime divisor. The old prime `3`
lifted to `3^2` nevertheless gives `2` exact order six modulo `63`, and `A`
also has exact order six modulo `63`.

As a hostile control, exhaustive iteration of `(4)` modulo `63` has only

```text
(5,62,56),       (14,26,20),       (35,47,41),          (48)
```

all of length three. Thus `63` is the sixth-iterate derivative degree and a
modulus supporting both the multiplicative order-six clock of `2` and the
central-sign lift `A^3=-I mod 63`; it is not the affine three-cycle multiplier
(`35/8`). Nor is it the exact-period-six point count, which is `54`; it hides
no scalar period-six orbit modulo `63`.

## 8. Planar-Jacobian compatibility and stopping boundary

The matrices `B` and `C=[[0,-1],[1,1]]` are determinant-one linear
polynomial automorphisms, so the lift lies on the trivially invertible side
of `JC(2)`. After real linear conjugacy it is a sixty-degree rotation;
THM-3465 already closes the real dagger-paired pure-character Keller sector.
Arbitrary complex cyclic-equivariant pairs are a different carrier and can
contain nonlinear triangular automorphisms.

THM-3742 already proves, on a different Pell/conic carrier, that ordinary
projectivization kills a central sign and that a stereographic sidecar can
restore the full clock. Relative to THM-4139, the new content here is the
universal quadratic-cycle interpolant, complete integral hexagon,
Pythagorean forcing, exact algebraic period-six census, and the
horizontal-section root-divisor/two-fibre decomposition `(40)`--`(41a)`.

The intrinsic horizontal-section information is `(40)`: its three zeroes form
a set-theoretic singleton `q`-fibre plus a two-reduced-point outer `q`-fibre
for every `a!=0`; the latter maps to the node of the horizontal image `S`.
Equation `(25)` is instead a hostile comparison model after the chosen
root-divisor identification `(38)`. It does not place a `C_6` action on the
horizontal surface or BC cover. THM-4138's actual exclusion uses its
independently proved Mordell--Weil, vanishing-loop and orbit-merger mechanisms,
not this model.

For an argument that uses this particular comparison route, the retained data
before quotienting must include at least

```text
root divisor + model defining-section character + a + q-fibre + t/rho sheet,
actual BC normalization and target map + labelled meridians when constructed. (49)
```

The model divisor has projective period three while its defining cubic sees
the order-six central sign; the actual horizontal response then refutes the
common-fibre inference. This is not asserted to be a universal minimal carrier
for other Jacobian cells. In particular, no local normal character or
labelled-meridian system for a new open cell has been constructed here.

This theorem does not furnish an entry theorem placing a hypothetical Keller
counterexample in this symmetry class, construct a Keller map on the
horizontal section, reprove or strengthen THM-4138, address the now-closed
exact-`M=8` trichotomy, enter an `M>=9` cell, or prove `JC(2)` or `DC(2)`.

## 9. Replay

Run

```text
python3 -B 04-computation/rational_three_cycle_order_six_horizontal_carrier_thm4146.py
python3 -B -O 04-computation/rational_three_cycle_order_six_horizontal_carrier_thm4146.py
python3 -B 04-computation/rational_three_cycle_order_six_horizontal_carrier_thm4146_independent_audit.py
python3 -B -O 04-computation/rational_three_cycle_order_six_horizontal_carrier_thm4146_independent_audit.py
```

The normal and optimized streams byte-match. Both exact paths have semantic
digest

```text
014f11897d8ca216faa897f117db118a9760b791d393d664289b621d99456d3c.
```

**QED.**
