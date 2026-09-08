# The last two independent quadratic pencils: a smooth elliptic survivor

**Status: PROVED ANALYTICALLY + FINITE-EXACT; INDEPENDENTLY AUDITED.**
This closes the remaining two entries of the six-entry independent exact radical
pencil table on the fixed surface W2. The elliptic entry has a globally smooth
survivor with a rational mate and source unit response of exact order two. Its
unit generates a pointed Weyl module isomorphic to that of an earlier rational
generic-fibre example. Its generic completed fibre has genus one. No dependent
one-dimensional pencil, global regular mate, full response-module isomorphism,
or case of the general Jacobian conjecture is claimed.

## 1. Inheritance, conventions, and exact scope

The inherited whole-table supplier is
[All exact radical pencils through degree eight](planar_jc48_sep08_exact_pencils.md),
and the exact global coefficient criterion is
[Global DG quadratics](planar_jc48_sep08_dg_quadratic.md).
The response supplier is the component-jet theorem in
[Vertical component jets](planar_jc48_sep06_torsion.md), with
[the coefficient-span lemma](continuing11_20260908_quadratic_family.md).
These are proved local dependencies; this report makes no literature-priority
claim. The first four pencil entries were paid in
[the fifth-order family](continuing11_20260908_quadratic_family.md) and
[the other three translated pencils](continuing12_20260908_quadratic_pencils.md).

The closest mechanism is coefficient matching followed by both-chart gradient
tests. The canonical hostile is a source-submersive function with a critical
point on the added divisor. The corrected near miss is treating a rational
mate as a rational generic-fibre coordinate: the elliptic survivor below has
a rational mate but a cubic function-field extension. The least-used sidecar
is the discarded exact elliptic radical pencil. The live concepts are whole
pencil matching, boundary jets, geometric integrality, labelled principal
parts, and the pointed Weyl module. The connection maps the last two objects
to the unit's operator response while destroying generic-fibre genus; the
missing sidecar is the function field of the generic fibre.

Use the two charts

    x=1/R,  t=-R^2-R^4 B,
    omega=dx wedge dt=R^2 dR wedge dB.

The added divisor is D={R=0}. For F in C[x,t], use
D_F G=F_x G_t-F_t G_x. Every response statement refers to the original ring
C[x,t]/D_F C[x,t], not a substituted affine chart or the ring O(W2).

For a source quadratic F=N(x)t^2+P(x)t+Q(x), globality is equivalent to

    deg N<=8,
    P=2N_8 x^6+2N_7 x^5+sum_(j=0)^4 P_j x^j,
    Q=N_8 x^4+N_7 x^3+(P_4-N_6)x^2+(P_3-N_5)x+Q_0.       (1)

Here N_j means the coefficient of x^j. We classify the cases in which the
*whole* two-dimensional space V=span{N,P^2-4NQ} equals either span{1,x} or
u^5 span{1,u^3}, with u=x-h. In particular N and P^2-4NQ are independent.
The inherited table and rational-mate criterion imply that these are the
only remaining independent cases, not that dependent pencils are closed.

## 2. The affine-linear pencil has no globally smooth member

Let V=span{1,x}. Then N=ax+d. Formula (1) initially allows P of degree at
most four and Q=P_4 x^2+P_3 x+s. Successively the coefficients of x^8,
x^6, x^4, and x^2 in P^2-4NQ force P_4=P_3=P_2=P_1=0. Thus exactly

    F=(ax+d)t^2+p t+s,  a*p!=0.                         (2)

The last condition is the exact pencil rank condition, since the coefficient
determinant of N and P^2-4NQ is -a p^2. In the source chart F_x=a t^2,
so F_x=0 forces t=0 and then F_t=p!=0. However the complete boundary
expression is

    F_inf=s-pR^2+aR^3+(d-pB)R^4+2aBR^5
            +2dBR^6+aB^2R^7+dB^2R^8.                 (3)

Both first derivatives vanish everywhere on D. Consequently every member
of this whole independent pencil fails global submersion, despite source
submersion.

## 3. Complete matching in the elliptic pencil

Let V=u^5 span{1,u^3}. Necessarily N=u^5(a u^3+d). Since P^2-4NQ is
divisible by u^5, P is divisible by u^3. Its two highest coefficients are
fixed by (1). Write

    P=u^3(2a u^3-4ah u^2+p u+2d+k),
    Q=s+a u^4-4ah u^3+(4ah^2+p)u^2+(d+k-2hp)u.         (4)

The constant s is Q(h). Conversely (4) satisfies every global coefficient
condition (1). The *complete* discriminant calculation is

    (P^2-4NQ)/u^5
      =-4ds+(k^2+8dhp)u+2(pk-8adh^2)u^2
         +(p^2-8ahk-4as)u^3.                          (5)

Therefore exact whole-pencil membership, including rank, is precisely

    k^2+8dhp=0,  pk-8adh^2=0,
    d*(p^2-8ahk)!=0.                                  (6)

For actual gradient tests put w=1+u^2t and Y=u(uw-2h). The full source
expression is

    F=s+aY^2+pY+d u w^2+k u w.                        (7)

To retain every added-chart point, define

    M=2h-R[h^2+B(1-hR)^2],
    X_inf=(1-hR)M,
    Y_inf=-(1-hR)[h^2(3-hR)+B(1-hR)^3].

Since w=RM and u=(1-hR)/R, the complete polynomial extension is

    F_inf=s+aY_inf^2+pY_inf+d R(1-hR)M^2+kX_inf.       (8)

At D set C=3h^2+B. Direct differentiation gives

    F_inf|D=s+a C^2-p C+2kh,
    F_B|D=2aC-p,
    F_R|D=(-2aC+p)4h(C-2h^2)+4dh^2-kC.               (9)

If h!=0, equations (6) force

    p=-k^2/(8dh),  a=-k^3/(64d^2h^3),
    p^2-8ahk=9k^4/(64d^2h^2)!=0.

In particular k!=0. At the actual point C=4dh^2/k, both derivatives in
(9) vanish. Thus every translated h!=0 member is globally critical.

If h=0, equations (6) give k=0 and dp!=0. Now (8) simplifies completely to

    F_inf=s+a B^2-pB+dR^3B^2.                         (10)

For a!=0 the point (R,B)=(0,p/(2a)) is critical. For a=0 the boundary
tangential derivative is the nonzero constant -p. The source test below
shows that this is also sufficient. Hence the exact global-submersion
locus of the entire elliptic pencil is

    h=a=k=0,  d*p!=0,  s arbitrary,
    F=s+x(1+x^2t)[d(1+x^2t)+p x].                    (11)

This includes all degree drops consistent with whole-pencil rank. The
discarded conditions d=0 or p=0 are not part of this independent pencil.

## 4. The survivor is globally smooth and has a rational mate

Put u=x, w=1+u^2t, L=dw+pu, and g=F-s=u w L. At u=0, the actual
source derivative g_x is d!=0; w cannot be held independent of u there.
On u!=0, (u,w) are valid coordinates, and

    g_u=w(dw+2pu),  g_w=u(2dw+pu).

If g_u=0 and w=0, then g_w=pu^2!=0. Otherwise dw=-2pu and
g_w=-3pu^2!=0. Equation (10), with a=0, pays all points of D. This
proves global submersion on W2 for every parameter in (11).

Let v=w/u. Although v is rational, its source Jacobian with u is u:
J_(x,t)(u,v)=u. Since g=u^3v(dv+p), we obtain

    D_g v=3g,
    G=v/(3g)=1/[3u^2(dw+pu)],  D_g G=1,
    g^2 G=w^2(dw+pu)/3 in C[x,t].                    (12)

These are identities in the original function field. G is a rational
mate, not a polynomial or globally regular mate. Indeed the unit order
computed below rules out any polynomial correction.

## 5. Geometric integrality, constants, and genus

The full field is

    C(x,t)=C(g)(v,u),  u^3=g/[v(dv+p)].               (13)

It is not C(g)(v). Over the algebraic closure of C(g), the right side
has a pole of order one at v=0 and therefore is not a cube. The cubic
is irreducible, as the field contains third roots of unity. This proves
geometric integrality of the generic curve; equivalently C(g) is
algebraically closed in C(x,t).

Any D_g-constant must be algebraic over C(g). Otherwise it would be
transcendental over C(g) in a field of relative transcendence degree one;
the full field would then be algebraic over a field killed by D_g.
Characteristic zero would force D_g to vanish, contradicting D_g v=3g.
Geometric integrality consequently proves the required exact constants:

    ker(D_g:C(x,t)->C(x,t))=C(g).                    (14)

For each c!=0, the complete source quadratic g-c has discriminant

    Delta_c=u^5(p^2u^3+4dc).                         (15)

Its valuation at u=0 is five, so it is not a square in C(u). The
polynomial g-c is primitive in t, since its leading and middle
coefficients have only u as a possible common root, whereas its constant
coefficient at u=0 is -c. Gauss's lemma proves that every nonzero source
fibre is irreducible.

Removing the square u^4 from (15), its smooth projective completion is
the double cover

    eta^2=u(p^2u^3+4dc).                              (16)

The quartic has four distinct roots: its discriminant is
-27 p^4(4dc)^4!=0. The double-cover genus formula gives genus one for
every c!=0. Adding the W2 boundary point does not change the function
field or its smooth projective genus. Geometric integrality, not
rationality, was the correct hypothesis for (14).

## 6. Complete labelled principal parts and the exact Weyl module

The source fibre g=0 has exactly three reduced, irreducible, pairwise
disjoint components

    E_u: u=0;  E_w: w=0;  E_L: dw+pu=0.               (17)

The last two equations are primitive linear polynomials in t; their
constant values at u=0 are respectively 1 and d, proving irreducibility
and separating them from E_u. On E_w, u is invertible and L=pu!=0,
separating E_w from E_L. All other source fibres are irreducible by
(15). Thus the inherited component-jet theorem, using global source
smoothness and (14), gives exactly two full torsion arms.

For the rational primitive in (12), the complete scalar principal parts
in the common parameter g, with component order (E_u,E_w,E_L), are

    pp(G)=(d/(3g^2)+p/(3d g),  0,  -p/(3d g)).        (18)

At E_u, the actual source relation w=1+u^2t gives
g=du+pu^2+O(u^3) and
G=1/(3du^2)-p/(3d^2u)+O(1). These determine the first entry of (18);
the remaining difference is regular. At E_w, G is regular. At E_L,
gG=w/(3u) evaluates to -p/(3d), giving the complete simple pole.

The regular E_w component is essential and cannot be discarded before
quotienting by the diagonal. With E_w as the zero reference, write

    theta=[1] <-> A_2/g^2+A_1/g,
    A_2=(d/3,0),  A_1=(p/(3d),-p/(3d)).              (19)

Their determinant is -p/9!=0. In particular theta has exact scalar
annihilator (g^2). This also follows directly from (12): gG=w/(3u)
has a genuine E_u pole and is regular at E_w, so no rational correction
in C(g) can make it polynomial. The constants field (14) pays every
possible rational correction.

Let A1=C<g,nabla>/(nabla*g-g*nabla-1) be the first Weyl algebra, with
the canonical connection nabla acting by principal-part differentiation.
Equations (19) give

    g theta=A_2/g,
    (g*nabla+2)theta=A_1/g.

Differentiation now generates every negative power in both directions.
The unit generates the *entire* two-arm torsion module. Moreover

    Ann_left(theta)=A1*g^2,
    (tors C_g,theta) isomorphic to (A1/A1*g^2,[1]).    (20)

For the exact converse, normal-order an arbitrary operator modulo A1*g^2
as a(nabla)+b(nabla)g. Its A_1 component is a(nabla)(1/g); the distinct
derivatives of 1/g are linearly independent, so a=0. Then its A_2
component forces b=0. This proves the full operator statement, not only
the scalar order. The j-th canonical derivative of theta has exact
scalar order j+2 for every j>=0.

## 7. What the completed table and connection do, and do not, establish

The six independent exact pencils now have the following global outcome.

| Exact whole pencil | Globally smooth members | Supplier |
| --- | --- | --- |
| span{1,x} | None: whole added divisor critical | Section 2 |
| u^3 span{1,u} | None: source critical line | Continuing12 quadratic pencils |
| u^4 span{1,u^2} | None: source or boundary critical | Continuing12 quadratic pencils |
| u^5 span{1,u} | Classified fifth-order family | Continuing11 quadratic family |
| u^7 span{1,u} | Classified seventh-order family | Continuing12 quadratic pencils |
| u^5 span{1,u^3} | Exactly (11) | Sections 3-6 |

The inherited rational-mate criterion makes this a complete classification
of globally smooth quadratic first functions with rational mates *whose
discriminant pencil has dimension two*. Dimension-one pencils remain a
separate obligation. The assertions are in the fixed original source
coordinates, with their actual W2 completion.

For a concrete genus-zero comparison, the already audited example
F=(x-1)(1+(x-1)^2t)[x(1+(x-1)^2t)-4] in
[the unit-order-two theorem](continuing11_20260908_quadratic_unit_two.md)
has rational generic fibre and the same pointed module (20). Thus even
all legal Weyl operations on the distinguished unit do not determine the
genus of the generic completed fibre. This is a comparison of pointed
unit-generated torsion modules, not an identification of whole response
modules or a claimed map between the underlying fibrations.

## 8. Exact control protocol

The paired source independently expands every matching coefficient, full
added-chart expression, declared critical point, source gradient,
primitive and repair, all labelled principal parts, and the quartic
discriminant. It checks the Weyl generators and derivatives j=0,...,6,
plus four declared nonzero rational parameter choices. This is a finite
control universe supporting the all-parameter proofs above; neither
parameter exhaustion nor arbitrary Weyl operators are inferred from
finite samples. All checks are always active under Python optimization.

Reproduction after relocation:

    python 04-computation/continuing12_20260908_last_quadratic_pencils.py
    python -O 04-computation/continuing12_20260908_last_quadratic_pencils.py

The runner writes its deterministic certificate next to the paired report
when located in 04-computation, or beside itself in the external packet.
Normal and optimized standard output and certificate bytes are compared
before freeze. No producer is imported, and no global theorem is inferred
from the computation alone.

## Independent acceptance

The [independent audit](continuing12_20260908_last_quadratic_pencils_audit.md) accepts
this scoped result and its analytic proof, with 174 exact gates per
normal and optimized pass. Its producer-independent controls retain the
declared positive and hostile boundaries. The parent replays both modes
against frozen raw LF outputs and certificates.

Subsequently, the [dependent-pencil supplement](continuing12_20260908_dependent_pencil_closure.md)
closed the separate dimension-one obligation. The [current synthesis](continuing12_20260908_synthesis.md)
routes the resulting complete submersive quadratic rational-mate classification.
