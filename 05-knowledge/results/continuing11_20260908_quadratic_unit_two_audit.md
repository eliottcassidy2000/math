# Independent audit: a quadratic unit of order two generates both torsion arms

**Verdict: PASS.** This is an independent analytic and exact reconstruction
of `continuing11_20260908_quadratic_unit_two.md`. Work over C on the actual
surface W2. The conclusion is about the ORIGINAL response module
C[x,t]/D_F C[x,t], with D_F=F_x partial_t-F_t partial_x. It is not a
quotient of O(W2), and it is not a polynomial Jacobian-one pair.

## Actual function and complete submersion check

Use u=x-1, w=1+u^2*t and L=(u+1)w-4. The proposed first function is

    F=u*w*L=x(x-1)^5*t^2+2(x-1)^3(x-2)*t+(x-1)(x-4).

Its degree in t is exactly two. Set y=u*w. On u!=0 the change to (u,y)
has Jacobian u^3 and F=(1+1/u)y^2-4y. A simultaneous zero of the two
derivatives would have y=0 from F_u=-y^2/u^2, whereas F_y=-4 there.
At u=0, direct source differentiation gives F_x=-3. This proves that
the entire affine source is smooth, with no parameter or height bound.

In the actual other chart x=1/r,t=-r^2-r^4*b put
zeta=r[1+b(1-r)^2]. Direct substitution yields

    F_inf=(1-r)(zeta^2-4),   y_inf=(1-r)(2-zeta).

Both are polynomials. On the whole boundary r=0, F=-4 and F_r=4.
Thus F is globally regular and has no critical point anywhere on W2.
No proposed local inverse is used to infer global submersion.

## Rational primitive and every component

The exact primitive is

    G=((1-2u)w+2)/(3u^2 w^2)
      =F/(3y^3)+2/y^2-1/y.

The source Jacobian J(F,G)=1 is checked by direct differentiation.
Independently,

    u=y^2/(F-y^2+4y),   t=(y/u-1)/u^2,
    D_F y=-y^4/(F-y^2+4y).

These identities pay C(x,t)=C(F)(y) and the exact differential constant
field C(F). Differentiating the displayed expression for G while holding
F fixed also gives D_F G=1.

The fibre F=0 has precisely the factors u,w,L. They are pairwise comaximal:
w mod u=1, L mod u=-3 and L mod w=-4. Each occurs once. The latter two
are primitive linear polynomials in t with coprime coefficients, hence
irreducible. No component label can be discarded. For c!=0, the
discriminant of u(u+1)w^2-4uw-c is

    4u[(c+4)u+c].

Its zero at u=0 is simple for every c!=0, including c=-4, so it is not
a square in C(u). The quadratic is irreducible after inverting u;
F(0,t)-c=-c rules out an extra factor lost in that localization. Thus
all other ORIGINAL-SOURCE fibres are irreducible. The fibre at -4 on
the whole surface also contains the boundary and is a different object.

## Exact order, full scalar parts and the canonical connection

G has no source poles outside u=0 and w=0. The complete scalar principal
parts in g=F, ordered by the components (u=0,w=0,L=0), are

    (9/g^2+4/g, (32/3)/g^2+4/g, 0).

The checker subtracts the first expression in (u,t), and the second in
(u,w), and cancels denominators to verify regular remainders. The third
component has neither denominator zero generically. In particular,

    F^2 G=L^2[(1-2u)w+2]/3

is polynomial, with derivative F^2. A polynomial primitive of q(F), with
ord_0 q<2, would be q(F)G+H(F), because the constant field is C(F).
The required pole correction on w=0 creates an uncancelled pole on L=0.
Hence the complete C[F]-annihilator of theta=[1] is (F^2). This rules
out polynomial mates of every degree and proves the exact unit order.

The proved component-jet theorem in
`planar_jc48_sep06_torsion.md` applies: the affine gradient is a unit
ideal, its rational constants are C(F), and the fibres are fully accounted
for. It identifies all torsion with two copies of
Pr_0=C[g,g^-1]/C[g], and the canonical connection nabla with d/dg on
the scalar parts. This imported statement is used in its original ring.

Modulo the common component diagonal, write A=(9,32/3), B=(4,4).
Their determinant is -20/3, and theta=A/g^2+B/g. Then

    g theta=A/g,             (g nabla+2)theta=B/g

in the principal-part module. Derivatives of 1/g give every negative
power in characteristic zero. Therefore theta generates BOTH full
torsion arms under multiplication and the canonical connection. Its
j-th derivative has exact scalar pole order j+2 with nonzero leading
vector (-1)^j (j+1)! A.

There is an exact cyclic presentation. Let A1 be the Weyl algebra with
nabla*g-g*nabla=1. The left annihilator of theta is exactly A1*g^2,
and Tors(C[x,t]/D_F C[x,t]) is A1/(A1*g^2), sending 1 to theta.
For an elementary proof, put every operator in normal order with g on
the right. Modulo the left ideal A1*g^2 it is a(nabla)+b(nabla)g.
Its B component on theta is a(nabla)/g, which vanishes only when a=0:
the derivatives of 1/g are nonzero multiples of distinct negative powers.
Then its A component is b(nabla)/g and forces b=0. Surjectivity was
proved by the two displayed extraction identities. This proves all
operator degrees, not only the tested derivative orders.

## Scope, prior gap, and reproduction

The earlier `continuing10_20260908_dg_linear_all_m.md` excludes unit order
two for globally smooth rationally integrable source-linear functions.
The actual quadratic above therefore proves that degree two is the
smallest possible source t-degree for this phenomenon on W2. The
response unit is being classified; other classes and connection
derivatives are not covered by the earlier gap.

The standalone audit checks 76 always-active identities and predicates,
including all-point algebraic critical-point tests, whole generic-fibre
identities, exact full scalar parts, and derivative controls j=0..12.
The finite c controls supplement the uniform coefficient 4c proof;
they are not a census of fibres. No producer is imported or executed.
Run the source normally and with -O. The frozen raw LF stdout agrees,
and its certificate regenerates unchanged. The all-degree conclusions
are proved above and in the inherited component-jet theorem.

General JC(2), polynomial pair entry, and the separate quartic closure
are not consequences or dependencies of this example. The incoming
[quartic theorem](planar_jc48_sep08_quartic_closure.md) was accepted before
filing this audit; its polynomial-pair scope is distinct.
