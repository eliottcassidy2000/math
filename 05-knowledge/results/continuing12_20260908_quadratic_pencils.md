# Three further exact quadratic pencils: two obstructions and a cubic unit

**Status: PROVED ANALYTICALLY + FINITE-EXACT; INDEPENDENTLY AUDITED.**
This classifies three specified whole discriminant pencils of global
quadratics on the fixed surface W2. Two admit no globally submersive
member. The third has an explicit complete submersion family with unit
order three, exactly two full source torsion arms, and two complex
coefficient-alignment branches. No classification of the remaining exact
pencils or general Jacobian conjecture is claimed.

## 1. Inheritance and exact scope

Use the actual charts x=1/r,t=-r^2-r^4b on
W2=(P1_x x P1_Z) minus {Z=x^2}, with omega=dx wedge dt=r^2 dr wedge db.
Put u=x-h, h arbitrary until submersion is imposed. For
F=N(x)t^2+P(x)t+Q(x), N nonzero, the pencil is the entire space
span{N,P^2-4NQ}; it is required to equal the specified two-dimensional
space, rather than merely contain one exact member.

The closest mechanisms are the complete six-pencil classification in
`planar_jc48_sep08_exact_pencils.md`, the full global coefficient criterion
in `planar_jc48_sep08_dg_quadratic.md`, and the original-source component
torsion theorem in `planar_jc48_sep06_torsion.md`. The recent
`continuing11_20260908_quadratic_family.md` supplies the coefficient-span
interpretation of the canonical Weyl action. Its fifth-order pencil is
already closed and is not reexamined here.

The retained hostile is a source-submersive rational-mate function with a
critical point on the added chart. The corrected near miss is to discard
the u=0 line when changing fibre coordinates. The least-used sidecars are
the two-moving-root fourth-order pencil, the translated seventh-order
root, and every scalar coefficient on every labelled fibre component.
Concept board: complete coefficient boxes; affine critical lines;
boundary normal derivatives; complete fibre factors; scalar pole order;
coefficient rank. The decisive tests below are whole-pencil matching and
actual gradients before any primitive search.

## 2. The whole third-order pencil is source-critical

Suppose span{N,P^2-4NQ}=u^3 span{1,u}. Then N=u^3(au+d), deg N<=4.
Globality gives deg P<=4 and Q=P_4 x^2+P_3 x+constant. Since u^3 divides
P^2, write P=u^2(pu^2+qu+v). The u^8 coefficient of the discriminant is
p^2, so p=0. Its subsequent u^6 coefficient is q^2, so q=0. Thus exactly

    F=s+u^3(au+d)t^2+v u^2t,
    d*v!=0

has the required full pencil. The nonzero rank condition is d*v^2!=0;
the converse follows directly from the displayed discriminant.
Both source derivatives vanish on the entire u=0 line. Therefore this
whole pencil contains **no** globally submersive first function.

## 3. The whole fourth-order pencil fails in one of the two charts

Suppose the pencil equals u^4 span{1,u^2}. Write N=u^4(au^2+d).
The global boxes and u-adic divisibility give P=u^2(pu^2+qu+v) and

    Q=(p-a)u^2+[q+(4a-2p)h]u+s.

The u^8 discriminant coefficient forces p=2a. The remaining pencil
condition and independence are precisely

    q(v-2d)=0,
    a(v-2d)^2-dq^2!=0.

These split into two disjoint exhaustive branches. Put w=1+u^2t.

If q=0, write k=v-2d; then a*k!=0 and

    F=s-d-k+(au^2+d)w^2+kw.

The entire original u=0 line is critical. If v=2d, then d*q!=0 and

    F=s-d+(au^2+d)w^2+quw.

This second branch is source-submersive. Off u=0, its derivatives in
(u,w) are w(2auw+q) and 2(au^2+d)w+qu. If the first vanishes, the second
cannot: for w=0 it is qu; otherwise it is -dq/(au) when a!=0, while
a=0 leaves the first derivative qw. On u=0 the actual source derivative
is q.

Nevertheless, in the full added chart, set
M=2h-r[h^2+b(1-hr)^2]. Direct substitution gives

    F_infinity=s-d+a(1-hr)^2M^2+d r^2M^2+q(1-hr)M,
    F_b|D=0,
    F_r|D=-(4ah+q)(3h^2+b).

The actual point r=0,b=-3h^2 is always critical. If 4ah+q=0, the whole
boundary is critical. Thus this second whole pencil also contains
**no** globally submersive first function. Source submersion alone
would have missed the obstruction.

## 4. The complete seventh-order global submersion family

Suppose span{N,P^2-4NQ}=u^7 span{1,u}. Then N=u^7(au+d).
Globality and u^7 dividing P^2 force

    P=u^4[2au^2+(2d-4ah)u+p].

The complete four lower global equations determine Q up to its value s
at u=0. Put k=p+4dh and z=u-2h+u^3t. The resulting entire family is

    F=s+u z[(au+d)z+k],
    N=u^7(au+d),
    P=u^4[2au^2+(2d-4ah)u+k-4dh],
    Q=s+u(u-2h)[a u(u-2h)+d(u+2h)+k-4dh].

Conversely these satisfy every global equation and

    P^2-4N(Q-c)=u^7[(k^2+4a(c-s))u+4d(c-s)].

The whole pencil has rank two exactly when dk!=0. There is no additional
matching equation. In particular a=0 is retained.

The exact global-submersion locus is

    d*k*h*(k-2dh)!=0,    a and s arbitrary.

Indeed off u=0, y=uz has J(u,y)=u^4 and
F=s+(a+d/u)y^2+ky. Its u derivative at fixed y is -dy^2/u^2; vanishing
forces y=0, where the other derivative is k!=0. On u=0 the actual
source derivative is -2h(k-2dh).

For the complete second chart set R=h^2(3-hr)+b(1-hr)^3. Then z=-rR and

    F_infinity=s+(1-hr)R[(a(1-hr)+dr)R-k].

On D the tangential derivative is 2aR-k. If a=0 it is -k!=0. If a!=0,
on its zero locus the normal derivative is dR^2, nonzero since
R=k/(2a), d,k!=0. Thus all boundary points are paid; no exceptional
boundary matching is required in this pencil.

## 5. Rational mate, all source fibres and the exact unit order

Put g=F-s and L=(au+d)z+k. The rational field and derivation are

    u=dy^2/(g-ay^2-ky),   C(x,t)=C(g)(y),
    D_F y=-du^2y^2=-d^3y^6/(g-ay^2-ky)^2.

Their rational constants are exactly C(g). Integrating with g fixed gives
one actual rational mate

    G=g^2/(5d^3y^5)-gk/(2d^3y^4)
      +(k^2-2ag)/(3d^3y^3)+ak/(d^3y^2)+a^2/(d^3y),
    D_F G=1.

Equivalently, with

    V=2(8a^2u^2-4adu+3d^2)z^2+k(7au-3d)z+k^2,

one has G=V/(30d^3u^3z^3). Thus g^3G=L^3V/(30d^3) is a polynomial
in the original source ring.

The special source fibre g=0 has exactly three reduced, disjoint,
irreducible components u=0,z=0,L=0. For z and L, primitivity as linear
polynomials in t follows from the nonzero constants -2h, k-2dh at
u=0 and k at u=-d/a when a!=0. Their pairwise separation values are
-2h,k-2dh,k. On a nonzero source fibre g=c, the discriminant has odd
u-valuation seven, and the only possible common root of N,P is u=0,
where the constant coefficient is -c. Hence every other source fibre is
irreducible. No corresponding statement about every global fibre is
needed or inferred.

At z=0 the highest scalar coefficient of G is k^5/(30d^3), nonzero.
At L=0 it is regular. Therefore g^2G has a genuine pole on the first
component and is regular on the second; no common rational constant
correction H(g) can make it polynomial. Since g^3G is polynomial, the
distinguished unit in C[x,t]/D_F C[x,t] has exact annihilator **(g^3)**.
The source component theorem gives exactly two full torsion arms. Thus
this whole pencil cannot supply more than two arms, despite its larger
unit pole order.

## 6. Complete principal parts and a complex alignment locus

Represent the three-component quotient by setting the L component to
zero. Write theta=C3/g^3+C2/g^2+C1/g, with vectors in C^2 corresponding
to (u=0,z=0). Their coordinates are

    C3_u=(k-2dh)^3(k^2+6dhk+24d^2h^2)/(30d^3),
    C3_z=k^5/(30d^3),
    C2_u=(k-2dh)(4ad^2h^2+2adhk+ak^2-6d^3h)/(3d^3),
    C2_z=ak^3/(3d^3),
    C1_u=C1_z=a^2k/d^3.

At u=0 the actual jet is z=-2h+u+u^3t; it must be imposed before
expansion. The source verifies the complete polynomial remainders to
order three in each actual local parameter. There are no other affine
poles. The zero principal part along L remains essential.

The unit-generated module is span{C1,C2,C3} tensor J, where
J=g^-1 C[g^-1], by the inherited Euler-interpolation lemma. If a=0,
C2_u=-2h(k-2dh)!=0, C2_z=0, and C3_z!=0, so the rank is two.
If a!=0, C1 is nonzero diagonal, and rank one holds exactly when

    4ah^2-6dh+3k=0,
    12d^2h^2-15dhk+5k^2=0.

Indeed C2_u-C2_z=-2h(4ah^2-6dh+3k)/3, while
C3_u-C3_z=-8h^3(12d^2h^2-15dhk+5k^2)/15. These are two admissible
complex branches: k/(dh)=(15 plus-or-minus sqrt(-15))/10 and
a=3(2dh-k)/(4h^2). They lie in the exact smooth locus. There is no
real alignment locus, since the second quadratic has discriminant
-15d^2h^2<0 for nonzero real d,h.

Thus the unit generates both full arms except on these two complex
branches, where it generates one. Its scalar order stays three
everywhere; the j-th canonical derivative has order j+3.

## 7. Exact order-three Weyl relation compiler

The order-three unit does not have left annihilator Dg^3 when its
coefficient rank is two or one. Retain D=C<g,nabla>/([nabla,g]-1).
For each constant relation lambda1*C1+lambda2*C2+lambda3*C3=0, set

    R_lambda=lambda1+(lambda2+lambda1*nabla)g
       +(lambda3+lambda2*nabla+lambda1*nabla^2/2)g^2.

The exact left annihilator is generated by g^3 and R_lambda for a
basis of the constant relation space. Thus one extra generator suffices
off the alignment locus; two suffice on it. Minimality of these chosen
generator lists is not claimed.

To prove equality, write a PBW remainder p0(nabla)+p1(nabla)g+
p2(nabla)g^2 modulo Dg^3. Its coefficient of 1/g is

    C1*p0+C2*(p1-nabla*p0)
       +C3*(p2-nabla*p1+nabla^2*p0/2).

The kernel of this constant coefficient matrix over C[nabla] is its
constant relation space tensored with C[nabla]. Inverting the displayed
triangular substitution gives exactly the R_lambda above. This proves
both containments for arbitrary finite Weyl operators. The relation
compiler retains information that scalar order alone discards.

## 8. Reproduction and stopping boundary

Run `04-computation/continuing12_20260908_quadratic_pencils.py` normally
and with `-O`. Outside installation its certificate is written beside
the source; when filed under `04-computation`, it uses
`05-knowledge/results`. The verifier imports no producer and passes
82 always-active exact gates, with matching LF output and certificate
bytes in both modes.

The finite universe is explicit: complete symbolic matching identities,
two full chart formulas, the rational mate and polynomial witness,
complete local cubic remainders, the exact quotient by
5q^2-15q+12 for both complex alignment branches, four stated real
parameter controls, and the symbolic PBW triangular compiler. These are
controls for the analytic proofs, not a finite sample asserted to prove
an all-parameter classification.

The two eliminated pencils and the complete seventh-order family are
paid stopping points. More than two full torsion arms must come from
another exact pencil or a higher source degree. Every response-module
statement here retains the original ring C[x,t]; no global regular
mate, unrestricted quadratic classification, or general JC result follows.

## Independent acceptance

The [independent audit](continuing12_20260908_quadratic_pencils_audit.md) accepts
this scoped result and its analytic proof, with 316 exact gates per
normal and optimized pass. Its producer-independent controls retain the
declared positive and hostile boundaries. The parent replays both modes
against frozen raw LF outputs and certificates.
