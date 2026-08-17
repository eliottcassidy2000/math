# The newest different is intrinsic; a literal coordinate still needs every resultant

**Status: proved synthesis for THM-3538.**  The all-level result is an exact
local index criterion.  Index zero is verified only at newest levels two,
three, and four.  Equality at every depth remains open.

## Inheritance pass

The closest proved mechanism is THM-3533: at `p=(P_(n-1))` the preceding
cover is finite etale, exactly one last-step block meets `L=0`, and the
normalization discriminant exponent is one.  THM-3535 supplies full-degree
primitivity of each literal coordinate.  THM-2546 supplies constant-leading
cubics for `y,z` and the exact square factors in their discriminants.  The
canonical hostile is now THM-3537: a primitive literal `x_2` observation on
the old-`L` transversal has index length two because its quartic and
quadratic packets collide in residue.

The least-used sidecar was the complete product-discriminant formula.  The
normalization sees ramification inside blocks; a prescribed coordinate also
sees coincidences between different blocks.

## The carrier

After strict-henselian splitting, write the `3^(n-1)` predecessor sections as
`q_i`, with `q_0` the unique section satisfying `v(L(q_0))=1`.  For
`theta=y,z`, or `theta=u=1/x` when every preceding `c(q_i)` is a unit, let
`f_(theta,i)` be the monic last-step coordinate cubic.  Its exact
discriminant is

```text
Disc(f_(theta,i))=-L(q_i)h_theta(q_i)^2.
```

The complete square/index carrier is

```text
C_(theta,n)
 = product_i h_theta(q_i)
   product_(i<j) Res(f_(theta,i),f_(theta,j)).
```

Because THM-3535 makes the product of the blocks the minimal polynomial,

```text
v Disc(m_theta)=1+2v(C_(theta,n)),
length(B_n/R[theta_n])=v(C_(theta,n)).
```

This turns local maximality into an iff test: every unramified block must be
squarefree, the ramified block may have only its forced double shadow, and
all distinct block reductions must be coprime.  Equivalently, the derivative
gcd of the total reduced polynomial is exactly the one forced linear shadow.
Primitivity guarantees degree but does not imply this test.

## Exact finite closures

The companion closes exactly three newest levels.

```text
n=2 over Q:    degree 9,  q0=(2/27,1,1),
               gcds Y-1/3, Z, U;
n=3 mod 53:    degree 27, q0=(2,23,4),
               gcds Y-7, Z, U;
n=4 mod 41:    degree 81, q0=(7,26,3),
               gcds Y-18, Z, U.
```

The level-two calculation factors the other two predecessor sheets as one
exact quadratic algebra over `Q`.  The level-three row uses nine completely
split predecessor blocks and independently agrees with quotient norms.  The
level-four row uses nine split parents followed by exact cubic quotient
norms, covering all 27 algebraic predecessor blocks.  Each row checks the
reciprocal chart, unique `L` block, full degree, internal squarefreeness, and
cross-block noncollision.  Smoothness of `q_0`, `det J_F=-2`, and the checked
unit denominators Hensel-lift each row to characteristic zero, so these are
nonidentity certificates for the generic divisor, not an all-level pattern
extrapolation.

## Reciprocal versus literal x

At the newest prime, under the reciprocal chart gate,

```text
v Norm(x_n)=-1.
```

Thus `x_n` is not integral while `u_n=1/x_n` is.  If `N=3^n`, `m_u` is the
monic reciprocal polynomial, and `I` is its local index, then

```text
E_x(X)=X^N m_u(X^(-1))=m_u(0)m_x(X),
v Disc(m_u)=v Disc(E_x)=1+2I,
v Disc(m_x)=3-2N+2I.
```

The exponent-one object is therefore either the integral reciprocal order or
the primitive nonmonic raw-cleared reversal.  Calling the monic literal-x
polynomial an integral power order would erase its pole and give the wrong
normalization.

## Why the old-L hostile does not conflict

THM-3537 has normalization exponent four, inertia `(4)(2)(1)^3`, literal
`x_2` order exponent eight, and index length two.  Its quartic/quadratic
packets share the residue value zero.  The newest-prime criterion does not
apply there: before the last inverse step, the old-prime cover is already
ramified and nonproper, so there is no finite-etale split with a unique
last-step `L=0` block.  The hostile validates the need for resultants rather
than refuting the criterion.

## Connection contract and boundary

| field | exact answer |
|---|---|
| source | newest-prime finite-etale predecessor algebra |
| target | split last-step coordinate cubics |
| map | block product followed by discriminant/resultant factorization |
| preserved | primitive degree, internal discriminants, residue collisions, local index |
| lost by primitivity alone | every internal square factor and cross-block resultant valuation |
| reciprocal sidecar | all preceding `c`-coordinates are units |
| exact positive rows | newest levels `2,3,4` have index zero for `y,z,u` |
| hostile | old-`L` canonical `x_2` has index two |
| open | unit carrier for every `n>=5` |

THM-3532 adds a separate covariance warning: nonlinear conjugacy transports
the observation together with its coordinate chart; it does not preserve the
standard literal-coordinate family automatically.  Nothing here classifies
arbitrary Keller maps, proves an all-level old-prime law, or connects the
carrier to LRC, tournaments, or a Jacobian-conjecture classification.
