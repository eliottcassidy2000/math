# The fixed fifth-order quadratic pencil has an exact unit-cyclicity wall

**Status: PROVED + FINITE-EXACT; independently audited.**
On the fixed surface W2, this exhausts the globally submersive members
of one precisely specified discriminant pencil and classifies their
distinguished unit's Weyl module. Their scalar repair order is uniformly
two, while the number of unit-generated torsion arms changes on a genuine
parameter wall. No broader quadratic classification or JC(2) conclusion
is claimed.

## 1. Object, inheritance and statement

Keep the actual source and second chart

    W2=(P1_x x P1_Z) minus {Z=x^2},
    U0=A2_(x,t), Uinf=A2_(r,b),
    x=1/r, t=-r^2-r^4b, omega=dx wedge dt=r^2 dr wedge db.

The complete quadratic section space is inherited from
[planar_jc48_sep08_dg_quadratic.md](planar_jc48_sep08_dg_quadratic.md).
The fixed exact space `u^5 span{1,u}` is one of the whole-pencil spaces in
[planar_jc48_sep08_exact_pencils.md](planar_jc48_sep08_exact_pencils.md).
The closest realized point is
[quadratic_unit_two](continuing11_20260908_quadratic_unit_two.md), with its
independent [Weyl audit](continuing11_20260908_quadratic_weyl_audit.md).
The complete affine torsion and canonical connection are inherited from
[planar_jc48_sep06_torsion.md](planar_jc48_sep06_torsion.md), retaining
**THM-3412**,
`01-canon/theorems/THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms.md`,
and **THM-3770**,
`01-canon/theorems/THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate.md`.
No external-priority claim is made for exactness or Weyl module theory.

For a fixed h in C put u=x-h, w=1+u^2t. Consider genuine source-quadratic
global functions `F=N(x)t^2+P(x)t+Q(x)` whose **whole discriminant pencil**

    span{N,P^2-4NQ}=u^5 span{1,u}                       (1)

is exactly this fixed two-dimensional space. An affine target shift is
allowed. No unrelated exact pencil or arbitrary global quadratic is part
of the exhaustion claim.

**Theorem.** The functions satisfying (1) and having no critical point
anywhere on W2 are precisely

    F=s+(a u^2+d u)w^2+k u w,
    k=-4ah,
    a*d*h*(d-4ah)!=0.                                 (2)

Every such F has the rational mate

    G=[2(d-2au)w-k]/[6d^2u^2w^2], J_(x,t)(F,G)=1.      (3)

Set g=F-s and use only the fixed affine response quotient

    C_F=C[x,t]/D_F(C[x,t]),
    D_F=F_x partial_t-F_t partial_x, theta=[1].

The complete torsion is two principal-part arms supported at g=0, and
the distinguished unit has exact C[g] annihilator `(g^2)` throughout
(2). Its j-th canonical derivative has exact order j+2.

Let `D_alg=C<g,nabla>/(nabla g-g nabla-1)`. Then:

* If d!=6ah, theta generates the entire two-arm torsion and its exact
  left annihilator is `D_alg g^2`.
* If d=6ah, theta generates one full arm inside that unchanged two-arm
  torsion. Its exact left annihilator is

      D_alg g^2 + D_alg [g nabla+2-(6a/k^2)g].        (4)

The wall d=6ah is included in the submersion domain: there d-4ah=2ah
is nonzero. Thus this is not an apparent wall created by singular or
inadmissible parameters.

The five live concepts are full global coefficient matching, boundary
normal derivatives, labelled scalar principal parts, scalar unit order,
and the canonical Weyl action. The canonical hostile is the one-arm
order-two model from the earlier Weyl audit; here it is realized by an
actual globally submersive family member. The corrected near miss is to
infer unit generation of all available arms from the repair order alone.
The least-used sidecar is the lower scalar coefficient together with the
regular third component's zero principal part.

## 2. Complete coefficient matching within the declared pencil

Write `N=u^5(a u+d)`, initially with a,d not both zero. The global
quadratic numerator criterion gives

    deg Q<=4,
    deg(P-2Qx^2)<=4,
    deg(N+Qx^4-Px^2)<=4.                              (5)

Equivalently `N=Qx^4+B x^2+C0`, `P=2Qx^2+B`, where each of Q,B,C0
has degree at most four. Since deg N<=6, the first equality forces
deg Q<=2: a term of degree more than six in Qx^4 cannot cancel against
Bx^2 or C0. Consequently deg P<=4.

The pencil condition makes `u^5` divide `P^2-4NQ`, hence `u^3` divides
P. Write `P=u^3(pu+q)`. The two highest equations in the last bound
of (5) then force

    Q=(p-a)u^2+[q-d+(4a-2p)h]u+s.                    (6)

The coefficient of u^8 in `P^2-4NQ` is `(p-2a)^2`. Its vanishing
over C forces p=2a. Put k=q-2d; equations (6) become exactly

    N=u^5(a u+d),
    P=u^3(2a u+2d+k),
    Q=a u^2+(d+k)u+s,

which is the family in (2), before the submersion conditions are imposed.
Conversely these coefficients satisfy all three degree boxes and

    P^2-4N(Q-c)
      =u^5[((k^2+4a(c-s))u+4d(c-s))].                 (7)

The two coefficient vectors of N and `P^2-4NQ` have determinant
`-d k^2`. Thus the pencil is the full two-dimensional space (1)
exactly when dk!=0. This pays both necessity and sufficiency, with no
search over bounded mate degrees. A target shift changes s and leaves
the same pencil space.

## 3. Global submersion: the sharp iff and the actual rational field

All displayed family functions are global for arbitrary a,d,k,h,s;
the relation k=-4ah is a submersion condition, not a condition for
global regularity. Put

    M=2h-r[h^2+b(1-hr)^2], Y=(1-hr)M.

Literal substitution gives the full polynomial second-chart expression

    F_inf=s+aY^2+d r(1-hr)M^2+kY.

On the entire added boundary D,

    F=s+4ah^2+2kh,
    F_b=0,
    F_r=4dh^2-(4ah+k)(3h^2+b).                       (8)

If 4ah+k is nonzero, the last expression vanishes at one actual value
of b. If it is zero, it is nonzero everywhere precisely when dh!=0.
This accounts for all boundary points.

On u!=0 use the actual rational source coordinate y=uw; its source
Jacobian with u is u^3, and

    F=s+(a+d/u)y^2+k y,
    F_u|y=-d y^2/u^2,
    F_y|u=2(a+d/u)y+k.                               (9)

Under dk!=0 these two derivatives cannot vanish together. On the
omitted source line u=0 the actual derivative is F_x=d+k. Thus the
full affine source is submersive exactly when `dk(d+k)!=0` in this
full-pencil family. Combining this with (8) gives precisely (2).
For example d=4ah is excluded because the original source line becomes
critical; d=6ah is not excluded.

The rational field and nonzero derivation are

    u=d y^2/(F-s-a y^2-k y),
    C(x,t)=C(F)(y),
    D_F y=-d u y^2
         =-d^2 y^4/(F-s-a y^2-k y).                  (10)

Their rational constants are exactly C(F). Holding F fixed, the function

    G=(F-s)/(3d^2 y^3)-a/(d^2 y)-k/(2d^2 y^2)

has D_F derivative one, and substitution of y=uw gives (3). Thus the
mate is one rational source function, not independently chosen fibrewise
primitives.

## 4. Every component, both scalar coefficients, and exact unit order

Set `L=(a u+d)w+k`. Then g=uwL. The zero fibre has exactly the three
reduced irreducible affine components E_u, E_w, E_L. Their separation
values are w|u=0=1, L|u=0=d+k and L|w=0=k, all nonzero in (2).
For irreducibility, w is primitive linear in t, and

    L=(a u+d)u^2t+a u+d+k

is primitive linear as well: its constant term is nonzero at u=0 and
at u=-d/a. This retains the whole residual component, even though G
is regular there.

For every c!=0, g-c has no original-source u factor. In the valid
(u,w) chart its quadratic discriminant is

    u[(k^2+4ac)u+4dc].                               (11)

It has an odd-order zero at u=0 because dc!=0, including the value
at which the bracket becomes a nonzero constant. Hence it is nonsquare
in C(u), and the polynomial is irreducible. Passing back to the source
introduces no component supported on u=0 because g-c=-c there. Equivalently,
gcd(N,P) is supported only at u=0; its multiplicity can jump to four
when k=-2d, and no constant multiplicity claim is needed. Therefore
zero is the only reducible fibre on the original affine source.

The full scalar parts of G in the same target parameter g are

    E_u: alpha_u/g^2+beta/g,
    E_w: alpha_w/g^2+beta/g,
    E_L: 0,

where

    alpha_u=(2d-k)(d+k)^2/(6d^2),
    alpha_w=-k^3/(6d^2),
    beta=-ak/d^2.                                    (12)

At E_u, impose w=1+u^2t before expanding. Then
`g=(d+k)u+a u^2+O(u^3)` and
`G=(2d-k)/(6d^2u^2)-2a/(3d^2u)+O(1)`.
At E_w, u is a unit and
`g=kuw+u(au+d)w^2`. These expansions give (12); the exact source
also checks the complete remainders after subtracting both scalar terms.
There are no affine poles besides E_u and E_w.

In particular alpha_w is always nonzero, so the unit has order two.
Explicitly

    P2=g^2G=L^2[2(d-2au)w-k]/(6d^2)

is polynomial with D_F P2=g^2. Conversely gG has a genuine simple
pole on E_w and is regular on E_L. Any other rational primitive of
g differs by H(F), using (10); cancelling the E_w pole forces a pole
on E_L. Hence no polynomial primitive of g exists. This proves the
exact C[g] annihilator `(g^2)` without a degree bound.

At the admitted parameter d=-2ah, equivalently 2d-k=0, alpha_u vanishes.
G then has only a simple pole on E_u, since beta remains nonzero.
The order-two conclusion survives because alpha_w is nonzero. It would
be false to assert that G always has order two on both polar components.

The boundary value in (8) becomes F=s-4ah^2, different from the
special value s. Its fibre on W2 contains D; no claim that s is the
only reducible global fibre is made. The Hamiltonian response quotient
continues to use only C[x,t]. Since omega vanishes at D but F_r=4dh^2
does not, the source derivation is not a regular derivation of O(W2).

## 5. The full torsion and the exact Weyl wall

The affine gradient and rational constant-field hypotheses of the
component-jet theorem are now paid. With three components at the only
reducible affine fibre, the complete torsion is

    W tensor J,
    W=C^3/C(1,1,1), J=g^-1 C[g^-1].

Multiplication by g removes nonnegative powers; nabla acts as d/dg.
Choose representatives with last coordinate zero and write

    A=(alpha_u,alpha_w,0), B=(beta,beta,0),
    theta=A/g^2+B/g.

The regular E_L is not disposable: without its zero principal part,
the B vector would incorrectly become diagonal and disappear. In the
actual two-dimensional component quotient,

    det(A,B)=(alpha_u-alpha_w)beta
            =-ak(2d+3k)/(6d^2).                      (13)

All factors except 2d+3k are nonzero in (2). Since k=-4ah, this last
factor vanishes exactly at d=6ah.

There is a general elementary mechanism within this inherited model.
For any finite packet `v=sum_(j=1)^N A_j g^-j` in `W tensor J`, one has

    D_alg v=span{A_1,...,A_N} tensor J.                (14)

Indeed E=g nabla has eigenvalue -j on the j-th summand. The interpolation
operator `product_(ell!=j)(E+ell)/(ell-j)` isolates A_j g^-j. Multiplication
by g^(j-1) brings this to A_j/g, and derivatives generate its whole arm.
Conversely Weyl operators act only on the scalar factor, so cannot leave
the displayed coefficient span. This proves (14) for arbitrary finite
packets, including dependent or zero coefficients. It is an explanation
of the existing principal-part model, not an external-priority claim.

Thus a distinguished unit generates the whole torsion exactly when its
scalar coefficient vectors span W; its exact pole order sees only its
last nonzero coefficient. If A_1,...,A_N are all independent, the exact
left annihilator is D_alg g^N: in a PBW remainder
`sum_(j=0)^(N-1) p_j(nabla)g^j`, the A_1 coefficient forces p_0=0,
then the A_2 coefficient forces p_1=0, and so on. Each step uses
independence of the derivatives of g^-1. The two cases below are
the independent N=2 case and an actual dependent-coefficient wall.

Off the wall, the actual operations

    g theta=A/g, (g nabla+2)theta=B/g

recover both independent first-level vectors. Differentiation recovers
every higher level, proving theta generates the complete two-arm torsion.
By PBW normal ordering, every Weyl operator has unique remainder
`a0(nabla)+b0(nabla)g` modulo the left ideal D_alg g^2. Applied to
theta, its B coefficient is `a0(nabla)(g^-1)`. The derivatives of
g^-1 are nonzero scalar multiples of distinct negative powers, so a0=0;
the A coefficient then forces b0=0. Thus the exact off-wall left
annihilator is D_alg g^2 and

    D_alg theta ~= D_alg/(D_alg g^2) ~= J direct_sum J.

On the wall d=-3k/2, the leading coefficients agree:

    alpha_u=alpha_w=alpha=-k^3/(6d^2)!=0,
    lambda=beta/alpha=6a/k^2.

Both A and B lie on the single component direction `(1,1,0)`. Since
g theta is a nonzero multiple of that direction times g^-1, the unit
generates its entire arm J and none of a complementary arm. Its scalar
order nevertheless stays two, and the whole module still has two arms.

For the exact wall annihilator, the same PBW remainder acts as

    alpha [b0(nabla)+(lambda-nabla)a0(nabla)] g^-1

on that direction, because g^-2=-nabla(g^-1). It vanishes if and only
if `b0=(nabla-lambda)a0`. Therefore its remainder is a left multiple of

    1+(nabla-lambda)g = g nabla+2-lambda g.

Together with g^2 this proves both containments in (4). This is an
exact left-ideal statement for arbitrary finite Weyl operators, not a
comparison of finitely many jets. The module generated by theta is
isomorphic to J, though theta itself is not the lowest-level generator.

For every parameter, nabla^j theta has exact scalar order j+2 because
the leading A vector is nonzero. Thus this integer statistic is constant
across the wall and cannot detect the change in unit-generated arms.

## 6. Literal controls, mechanism and scope

The point a=d=h=1 gives the earlier explicit quadratic example. The
wall point a=h=1,d=6,k=-4,s=0 is globally submersive, with

    A=(8/27,8/27,0), B=(1/9,1/9,0), lambda=3/8,
    F_x|u=0=2, F_r|D=24.

Its additional unit annihilator is `g nabla+2-(3/8)g`.
The full-pencil point a=d=h=k=1 is global and source-submersive but
has its boundary critical point at b=-11/5. It tests the necessity
of k=-4ah without dropping the independent-pencil hypothesis.
The admitted point a=h=1,d=-2 loses only the leading Eu pole and
remains off-wall. The admitted point a=h=1,d=2 has gcd(N,P)=u^4;
it remains off-wall with determinant -4/3. These are hostile controls
against treating a pole-order drop or a source gcd jump as the cyclicity
wall.

The source checks the complete coefficient identities, full boundary
chart, all-parameter derivative formulas, actual rational mate and
constant-field inverse, complete scalar remainders and the determinant.
Seven exact parameter controls include two wall points, pole drops,
the gcd jump and ordinary off-wall points. Their entire PBW remainder
spaces through derivative degree six have rank deg+2 on the wall and
2(deg+1) off it. These finite ranks verify the formulas; the complete
Weyl proofs are the PBW arguments above.
Euler interpolation additionally recovers every coefficient of all
declared two-component packets of orders one through six, checking the
general coefficient-span mechanism separately from the geometric family.

The successful connection preserves the fixed global pencil, actual
source ring, target shift, three component labels and both scalar
coefficients. Replacing the tuple by its highest pole order loses the
directional alignment measured by (13). Restoring that information
turns a formal one-arm hostile into a genuine smooth member of the
same geometric family.

Reproduce normally and with `-O`:

    python -B 04-computation/continuing11_20260908_quadratic_family.py
    python -B -O 04-computation/continuing11_20260908_quadratic_family.py

The standalone source imports no producer; before installation it writes
its certificate beside itself. Normal and optimized runs pass 173
always-active exact gates with identical raw LF output and regenerated
certificate bytes. The independent analytic audit below is accepted.

This family contains rational mates and globally submersive first
functions but no polynomial mates, since theta is nonzero throughout.
It does not classify other exact pencils, arbitrary DG quadratics,
or any Keller entry mechanism. JC(2) remains OPEN.

## Accepted independent audit

The [independent family referee](continuing11_20260908_quadratic_family_audit.md) accepts the complete
coefficient exhaustion, global iff, full principal parts, both annihilator
cases and the general coefficient-span lemma. No mathematical repair was
required; source, output and certificate remain frozen.
