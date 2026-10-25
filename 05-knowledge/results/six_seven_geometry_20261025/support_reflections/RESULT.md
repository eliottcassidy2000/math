# A geometric repair of the Boolean seven, and a Cremona symmetry of the actual hourglass K3 quotient

2026-10-25. **PROVED**, with [independent audit: PASS](../audits/SUPPORT_REFLECTION_AUDIT.md). The
elementary coordinate and combinatorial claims below have FINITE-EXACT
controls. Classical double-cover geometry is identified explicitly; no
historical-priority claim. All surface statements are over **C**, except
where finite-field or pointed-monoid conventions are explicitly named.

## 0. Inheritance, portfolio, and the arriving connection

Closest inherited mechanisms are the modified-complement seven of
[reflection_bridge/seven.md](../../reflection_bridge_20261025/seven.md),
the honest projective carrier in
[geometric_carrier.md](../../reflection_bridge_20261025/geometric_carrier.md),
and Session59's [root-cover theorem T1](../../session59_f1_geometry_and_information_20261025/toric-shared-roots/TORIFICATION_OBSTRUCTION.md).
The canonical hostile is the lost order under identifying empty and full;
the corrected near miss is treating a seven-point bijection as a scheme
automorphism. The least-used sidecar is the **exceptional direction** above
each coordinate point. We use the META-PATTERNS cards “Refine and saturate
before transporting a factor or shadow” and “Controlled forgetting and
unlabeled quotients require a sidecar.”

The anchor is a geometric realization of that seven-point reflection; the
niche is an actual Cremona symmetry of the central hourglass six-line
double plane; the wildcard is its six transitive-tournament flags and seven
anticanonical sections. The board is support strata / exceptional directions /
line-pair products / root characters / fixed loci / coordinate flags.

During this work, live `origin/main` acquired `f7041582a8`. Its audited
[Session60 hourglass-symmetry packet](../../session60_six_seven_geometry_20261025/hourglass-symmetry/RESULT.md)
already proves that the central six-line double plane is Kummer, identifies
the genus-two curve with `v²=u^6-1` over C, and obtains Picard rank twenty.
Those are **recovered results**, not discoveries of this lane. The present
addition is an explicit quadratic Cremona involution, its two lifts, and the
valuation obstruction to lifting it back to all individual roots.

Conventions: “degree-six del Pezzo” means `K²=6`, the blowup of **three**
points of P². This avoids the degree/blowup-index switch in MISTAKE-643,
[MISTAKES.md](../../../../01-canon/MISTAKES.md).

## 1. The seven-point reflection becomes geometric after resolving directions

Let `P=P²`, with coordinates `[x0:x1:x2]`. Consider the standard quadratic
Cremona map

```
c(x0:x1:x2)=(x1*x2:x0*x2:x0*x1).
```

Its graph closure has an especially small presentation:

```
Gamma = {(x,y) in P² x P² : x0*y0=x1*y1=x2*y2}.       (1)
```

Projection to the first factor is the blowup of the three coordinate
points. On the dense torus the second factor is coordinatewise reciprocal.
Above a support-two point the image is its missing coordinate point; above
the coordinate point `e_i`, the second factor is the entire opposite line
`y_i=0`. In local coordinates at `e_0`, the equations reduce to the usual
blowup incidence `x1*y1=x2*y2` with `y0=x1*y1`; the other charts are the
same. Thus Gamma is smooth. Swapping the two factors is an honest
involution `c_tilde`, resolving the original rational map.

There are six coordinate torus-fixed points `(e_i,e_j)` with `i!=j`.
Delete those six points and call the resulting invariant open surface `V`.
Use exactly the pointed-monoid `A1={0,1}` convention of the inherited
reflection packet. A monoid point of (1) has nonempty supports `(S,T)`
such that either `S=T={0,1,2}`, or `S and T` are disjoint. There are thirteen
such pairs: one dense-torus stratum, six one-dimensional boundary strata,
and six zero-dimensional strata. On V only seven remain:

```
(S,S^c),  empty != S != {0,1,2};    (full,full).        (2)
```

Projection to `S` is a bijection `V(A1) -> P²(A1)`. Under that bijection,
the geometric factor swap is **exactly** the old modified complement:
singletons and opposite pairs exchange, while full support is fixed.
This is a positive geometric realization, not a rejection of the old seven.

The price is measurable. All six proper-support labels in V now have torus
rank one. In P² the three singleton labels had rank zero. The blowdown is
continuous but its inverse seven-point bijection is not: for example
`{0} subset {0,1}` in P², whereas `( {0},{1,2} )` and
`( {0,1},{2} )` are incomparable in V's product-support order. The old
failure to preserve order has become a specific exceptional-direction cost.

Over F_q, blowing up each rational point adds q points, so

```
#Gamma(F_q)=q²+4q+1,    #V(F_q)=q²+4q-5,
#P²(F_q)=q²+q+1.                                      (3)
```

The blowdown fibres above the three coordinate points in V are G_m,
not points. At q=2 they each have one point: the seven-point bijection also
holds on F_2-points. At q=3, V has sixteen points and P² has thirteen.
There is no scheme isomorphism hiding in the `7=7` count. Its standard
torification has six G_m strata plus G_m², so the formal q=1 value is zero,
while its literal A1 point count is seven. The original carrier
`U=(P²xP¹xP¹) minus two sections` still has ten rank-zero and sixty-one
literal A1 points; replacing its P² factor by Gamma changes that carrier.
The original U reflections are coordinate permutations, not this new
Cremona inversion.

## 2. The hexagon joins six transitive tournaments and seven sections

The six boundary components of Gamma are the three exceptional divisors
and the strict transforms of the opposite coordinate lines. Their incidence
graph is a hexagon. They are labelled by the six proper nonempty subsets
of `{0,1,2}`; opposite components receive complementary labels.

The **six zero-dimensional torus orbits**, a different six-element set,
are indexed by full coordinate flags

```
{i} subset {i,j},        (i,j,k) a permutation of (0,1,2).
```

Give the three labels the transitive tournament order `i -> j -> k`.
Two such torus-fixed points lie on one invariant boundary curve iff the
orders differ by an adjacent transposition, equivalently their transitive
tournaments differ in exactly one arc. This recovers the `n=3` case of
[THM-207, permutohedron-transitive](../../../../01-canon/theorems/THM-207-permutohedron-transitive.md)
as the actual toric one-skeleton. Inversion sends the flag to
`{k} subset {j,k}`, hence sends the tournament to its converse. The two
cyclic tournaments are absent, exactly as in
[THM-2632, Farey V4 theta channels, §1](../../../../01-canon/theorems/THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar.md).
No arbitrary tournament predicate or path count is transported.

There is also an exact lattice connection to the **actual hourglass**.
Its six peripheral normalized forms are `E+a*x+b*z` with
`(a,b) in R={±(1,0),±(0,1),±(1,1)}`. In cyclic order

```
(1,0),(1,1),(0,1),(-1,0),(-1,-1),(0,-1),
```

these are precisely the primitive rays of Gamma's fan: start with the P²
rays `(1,0),(0,1),(-1,-1)` and insert the sum of each adjacent pair.
Every adjacent determinant is one. Thus normalized peripheral branch
labels identify with boundary-ray labels in a specified integral lattice,
and the central label E identifies with the zero cone/dense-torus label.
The D12 action found in the incoming Session60 packet acts on the same fan;
its central `-I` is opposite-ray exchange, hence the same modified
complement on the seven labels of V. This connects the two sevens through
an actual action and lattice, not only their size. It does **not** turn
branch-line incidence into boundary-curve incidence: opposite hourglass
lines meet the central branch line at a triple point, whereas the two
opposite boundary components are disjoint. The affine offset E and the
branch incidence are the retained sidecar.

Adjunction in the two `(1,1)` equations of (1) gives

```
-K_Gamma=O_Gamma(1,1),    K_Gamma²=6.                  (4)
```

The square is the coefficient of `h1²*h2²` in `(h1+h2)^4`, namely six.
Equivalently `-K=3H-E0-E1-E2`. Its seven sections are exactly the plane
cubics vanishing at the three coordinate points:

```
x0*x1*x2;  x0²*x1, x0²*x2, x0*x1²,
x0*x2², x1²*x2, x1*x2².                             (5)
```

After division by `x0*x1*x2`, these are the center and six vertices of an
integral hexagon. The dimension seven is thus a lattice-point computation,
and inversion exchanges three pairs of outer weights while fixing the
center. In this representation the 3+1+3 action has become linear.

Seven sections do not identify this surface with the hourglass: the
anticanonical Gamma has degree six and nine independent quadratic
relations (`28-19`); the seven-root hourglass has degree sixteen and four
quadratic equations. The same ambient P^6 and dimension of sections have
different multiplication laws. The sidecar is the graded ring, not another
cardinality label. Also boundary divisors and torus-fixed points both
number six but have different S3 actions; their incidence hexagon, not an
unmarked bijection, is the connection.

## 3. An explicit Cremona involution of the actual central hourglass quotient

Delete the central root E from the hourglass. Write its six branch forms in
the fixed order

```
(L0,L1,L2,L3,L4,L5)=(E+x,E-x,E+z,E-z,E+x+z,E-x-z).
```

The corresponding square-root surface `S0` is the smooth six-root K3.
Let `D` be the normal double plane

```
w²=product_i L_i=(E²-x²)(E²-z²)(E²-(x+z)²),           (6)
```

and `D_tilde` its minimal resolution. Its fifteen branch crossings are
double, with no triple; D has fifteen A1 singularities and D_tilde is K3.
The morphism from S0 to D is `[r_i] -> [E:x:z:product_i r_i]`, of degree
sixteen. It is the quotient by the even-sign subgroup of order sixteen.
The induced map `S0 -->> D_tilde` is rational, not claimed to be a finite
morphism between minimal K3 surfaces. Session60 already establishes this
quotient and the stronger Kummer statements.

Choose the three paired intersections

```
(L0,L2), (L1,L4), (L3,L5).
```

They are noncollinear. Explicitly set

```
[E:x:z] = B[X:Y:Z],
B = [ 1 -1 -1 ]
    [-1 -1  0 ],       det B=5.
    [-1  2 -1 ]
```

Then the forms, still in order, are

```
-2Y-Z, 2X-Z, Y-2Z, 2X-3Y, -X-2Z, 3X-2Y.             (7)
```

Let F be their product. The twisted Cremona involution

```
c_B(X:Y:Z)=(YZ:XZ:-XY)                              (8)
```

swaps each selected pair. Pair products acquire factors `-X²,-Y²,+Z²`,
respectively, giving the exact identity

```
F(c_B(X,Y,Z))=(XYZ)² F(X,Y,Z).                        (9)
```

This is the operative bridge: the three lost directional factors cancel
to a square in the **total product**.

Blow up the three chosen intersection points. They are the coordinate
vertices in (7), so the new surface is the same degree-six del Pezzo Gamma.
The strict transform of the sextic has divisor class `-2K_Gamma`. It has
twelve remaining transverse nodes. The normalized double cover of Gamma
branched in that divisor has twelve A1 singularities. Its minimal resolution
is D_tilde: locally at each blown-up branch crossing, writing `w²=xy` and
then `y=xt`, normalization uses `w=xv`, giving `v²=t`, the standard
resolution of the A1 point. At the remaining crossings one still resolves
`w²=xy` in the usual way.

On the base torus put `U=Y/X,V=Z/X` and `eta=w/(XYZ)`. Equation (9) gives
two involutions of the cover:

```
sigma_+(U,V,eta)=(1/U,-1/V, eta),
sigma_-(U,V,eta)=(1/U,-1/V,-eta).                     (10)
```

They extend over the normal double cover: (9) is a linearization of its
branch-square-root line bundle. They then lift uniquely over its A1
minimal resolution. Equivalently their weighted-plane formulas send w to
`-XYZ*w` and `+XYZ*w`, respectively; the sign change arises because the
product of the three coordinates in (8) is `-(XYZ)²`.

The base fixed locus is precisely the four torus points
`U=±1,V=±i`. No boundary component has a fixed point, since opposite
hexagon components are exchanged and disjoint. None of these four torus
points lies on the branch sextic: the four values of `F/(XYZ)²` are
`25,25,-175-600i,-175+600i`. The twelve remaining branch nodes form six
exchanged pairs, so their exceptional curves create no new fixed points.
Consequently **sigma_+ has eight fixed points and sigma_- is free** on
D_tilde. A global two-form is, on this torus,

```
omega=(dU/U wedge dV/V)/eta.
```

Thus sigma_+ preserves omega and sigma_- negates it. The first is a
symplectic (Nikulin) involution. The second yields an Enriques quotient:
the smooth quotient has an unramified K3 double cover, equivalently it has
`q=p_g=0` and nontrivial canonical bundle of order two. This is the
classical del Pezzo/Cremona/K3 mechanism of Peters--Sterk, here with a fully
specified nodal branch arrangement from the actual hourglass.

## 4. The symmetry does not lift to all six individual roots

The positive quotient bridge has a precise failure boundary. Let
`K=C(P²)` and let M be its six-root field, generated by square roots of
the ratios `L_i/L_0`. All square classes defining M have even valuation
at the new coordinate axis `X=0`: none of the six lines (7) is that axis.
But direct substitution gives

```
c_B^*(L0/L1) = X*(Y-2Z)/(Y*(X+2Z)),                 (11)
```

which has odd valuation at `X=0`. In a multiquadratic extension, the base
elements becoming squares are exactly the span of its defining square
classes. One proof is that each quadratic subfield corresponds to an
index-two character of the elementary abelian Galois group, and is
generated by a product of defining square roots. Equation (11) therefore
cannot become a square in M.

Any automorphism of M covering c_B would carry the existing square root
of `L0/L1` to such a square root, a contradiction. Hence there is **no
birational lift covering c_B to the full six-root surface**. This does not
deny other automorphisms of that K3, including the projective D12 symmetries
already recovered by Session60. The failed implication is “a symmetry of
the product square root preserves all individual square roots.” Its
strongest survivor is the two exact quotient involutions (10); the missing
coordinate is the labelled square-class extension. The same new-axis test
also obstructs this c_B on the seven-root field after adjoining sqrt(E)
in the usual ratio presentation, since `E=X-Y-Z` is not `X=0` either.

## 5. A smooth auxiliary family explains the seven-section K3 polarization

This auxiliary family is separate from the actual line arrangement above.
It gives a nonsingular nearby model of the geometric mechanism. Consider

```
F_lambda = sum_{i!=j} xi^4*xj² + lambda*x0²*x1²*x2².  (12)
```

Its strict transform on Gamma belongs to `|-2K_Gamma|`. At each coordinate
vertex its quadratic tangent cone is `y²+z²`, so the plane branch has a
node there for every lambda. On the torus, divide by `(x0*x1*x2)²` and put
`a=(x1/x0)²,b=(x2/x0)²`. The branch equation becomes

```
g(a,b)+lambda=0,   g=a+a^-1+b+b^-1+a/b+b/a.
```

Its logarithmic critical equations sum to
`(a+b)(ab-1)/(ab)=0`. If `b=1/a`, the first equation is
`(a-1)(a+1)(a²+a+1)/a²=0`; if `b=-a`, it is `a-a^-1=0`.
The only critical values of g are `6,-2,-3`. On each open coordinate line
the strict-transform equation has the simple roots of `1+t²`, and it avoids
the six torus-fixed corners. Consequently its strict transform is smooth
**iff lambda is not -6,2,3**. At lambda=2 it factors into six lines:

```
F_2=(x0²+x1²)(x0²+x2²)(x1²+x2²).
```

For a permitted lambda the double cover `pi:T->Gamma` is smooth. Its
canonical bundle is `pi^*(K_Gamma-K_Gamma)=O`, and
`pi_*O_T=O_Gamma direct-sum K_Gamma` gives `h¹(O_T)=0`, using rationality
and Serre duality on Gamma. It is K3. The pullback polarization
`A=pi^*(-K_Gamma)` has `A²=12`, so genus `1+A²/2=7`, and

```
h⁰(T,A)=h⁰(Gamma,-K_Gamma)+h⁰(Gamma,O)=7+1=8.        (13)
```

The six-to-seven numerical resonance here has a concrete mechanism:
degree-six anticanonical geometry supplies seven sections, and a double
cover supplies a genus-seven K3 polarization. This is different from the
canonical seven sections of the seven-root general-type hourglass.
Peters--Sterk already give lambda=1 in their Lemma 3.1; the explicit
three-parameter-value exclusion above is our direct computation, with no
priority claim.

## 6. Connection contracts and reproducibility

| Source -> target | Map and preserved data | Lost information / sidecar | Decisive test |
|---|---|---|---|
| Boolean modified seven -> Gamma open V | (2), blowdown support bijection, exact involution | three rank-zero labels become rank-one; exceptional direction | F3 counts 16 versus 13 |
| Six transitive tournaments -> Gamma fixed points | total order to complete coordinate flag; one-arc adjacency and converse | two cyclic tournaments omitted; labels retained | all 15 pairs of six orders |
| Actual six-root K3 -> double-plane K3 | total root product; parity character and degree16 | individual roots, arithmetic lifting, marked extension | odd axis valuation (11) |
| Double-plane K3 -> Enriques quotient | free lift sigma_- | anti-invariant form and sheet label | four base fixed points avoid branch; twelve nodes pair |
| Degree-six Gamma -> auxiliary K3 | double cover of a -2K branch divisor | branch coefficients indispensable | lambda=-6,2,3 hostiles |

Run from this directory:

```
python3 verify.py > normal.json
python3 -O verify.py > optimized.json
```

The checker uses explicit exceptions and imports no inherited producer.
It exhausts the 49 pairs of nonempty three-bit supports; all projective
point pairs over the prime fields q=2,3,5,7; six coordinate flags and their
15 unordered pairs; all 15 intersections of the actual six branch lines;
all twelve remaining nodes and four base fixed points. Symbolic identities
verify the branch pairing, involution square, odd-axis witness, seven and
nineteen anticanonical monomial dimensions, and the auxiliary critical
split. It prints 126 successful checks. These exact controls supplement
the geometric proofs; they do not infer a complex theorem from a finite
sample or certify a lift of rational points to individual integer squares.

Primary source checked 2026-10-25: Chris Peters and Hans Sterk,
[*On K3 double planes covering Enriques surfaces*](https://link.springer.com/article/10.1007/s00208-018-1718-4),
§3.1 Lemma 3.1 and §3.2 Theorem 3.4. The standard smooth
Cremona-invariant sextic construction and two lift types are recovered;
our actual nodal-branch instance is justified locally above rather than
silently substituted for their smooth-branch hypotheses. No LRC, arbitrary
tournament classification, or integer-square existence theorem follows.
