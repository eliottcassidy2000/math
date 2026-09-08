# Independent analytic referee: the fixed fifth-order quadratic pencil

**Status: PASS / PROVED analytical statements within the declared pencil and
original response ring.** The complete frozen primary report, including its
final Euler-interpolation addition, has been compared with the independent
derivation below. No mathematical repair is requested. This referee supplies
an analytic proof audit, separate from the producer's finite exact controls.

Audited primary: `continuing11_20260908_quadratic_family.md`, raw SHA256

    33bdd575435ef450e879ff008c97a80ceab3a97290255749e60caab61c088591.

The frozen companion source has SHA256
`360c01c346d92a28b5f4a53a05ece1eeb833ec7bcc603a71c2a133ac4cd949ef`;
its producer reports 173 always-active gates with identical normal/optimized
LF output. This referee does not substitute that finite gate count for the
all-parameter arguments below.

## 1. Exhaustion of the declared pencil

Write u=x-h and F=N(x)t^2+P(x)t+Q(x), with N nonzero. Assume precisely that
span{N,P^2-4NQ}=u^5 span{1,u}; no other exact pencil is included in this
classification. Thus N=u^5(au+d). The complete global quadratic criterion
in `planar_jc48_sep08_dg_quadratic.md` forces deg P<=4. Since u^5 divides
both the discriminant and N, it divides P^2, so P=u^3(pu+q).

The same global criterion, translated carefully from x to u, gives

    Q=(p-a)u^2+[q-d+(4a-2p)h]u+s.

The u^8 coefficient of P^2-4NQ is (p-2a)^2. It must vanish because the
declared pencil has degree at most six. Hence p=2a, and putting k=q-2d
gives the entire family

    w=1+u^2t,
    F=s+(au^2+du)w^2+kuw.

Conversely this expression is globally regular for every parameter value
and its discriminant is u^5[(k^2-4as)u-4ds]. The determinant of the N and
discriminant coefficient vectors is dk^2. Thus the exact two-dimensional
pencil assumption is equivalent to d*k!=0. This is a full coefficient
exhaustion within that specified pencil, not a search over a sample or a
classification of the other five exact pencils.

## 2. Every affine and boundary critical point

On u!=0, put y=uw. The actual source Jacobian J(u,y)=u^3 is nonzero and

    F=s+(a+d/u)y^2+ky,
    F_u=-dy^2/u^2,
    F_y=2(a+d/u)y+k.

Within the declared pencil d,k are nonzero, so F_u=0 forces y=0 and then
F_y=k!=0. On the omitted affine line u=0, the source gradient has
F_x=d+k and F_t=0. Thus source submersion is equivalent to d+k!=0.

For the complete second chart x=1/r,t=-r^2-r^4b set

    z=2h-r[h^2+b(1-hr)^2].

Direct substitution gives the polynomial identity

    F_infinity=s+(1-hr)[(a(1-hr)+dr)z^2+kz].

It includes all b-values and all points of this chart. Along D={r=0},

    F=s+4ah^2+2kh,
    F_b=0,
    F_r=4dh^2-12ah^3-3kh^2-(4ah+k)b.

If 4ah+k!=0 there is exactly one boundary critical point. If k=-4ah,
the normal derivative is the constant 4dh^2. Combining both charts and
the pencil condition proves the exact global-submersion locus

    k=-4ah,    a*d*h*(d-4ah)!=0,    s arbitrary.

No localization has removed an affine critical line or a boundary point.

## 3. Rational primitive, all source fibres and unit order

Put g=F-s and L=(au+d)w+k. Then g=uwL and

    G=[2(d-2au)w-k]/(6d^2u^2w^2)
     =g/(3d^2y^3)-k/(2d^2y^2)-a/(d^2y).

In the rational function field,

    u=dy^2/(g-ay^2-ky),
    D_F y=-duy^2=-d^2y^4/(g-ay^2-ky).

Consequently D_F G=1 and the rational constants are exactly C(g): the
field is C(g)(y), and a nonzero multiple of differentiation in y has
precisely this constant field in characteristic zero.

On the global-submersion locus, g=0 has the three reduced, disjoint,
irreducible source components u=0,w=0,L=0. The last factor is primitive
linear in t because its coefficients u^2(au+d) and au+d+k have no common
root: their potentially common values are excluded by d+k!=0 and k!=0.
The other disjointness checks are w|u0=1, L|u0=d+k and L|w0=k.

For a nonzero fibre value v of g, the discriminant is
u^5[(k^2+4av)u+4dv], so its u-valuation is odd and it is not a square.
The only possible common root of N and P is u=0, where Q-s-v=-v!=0.
Hence all other **original source** fibres are irreducible. The exact
gcd power can rise from u^3 to u^4 at k=-2d; the support argument pays this
allowed parameter case without an unjustified fixed gcd exponent.

The boundary lies in the different global fibre g=-4ah^2, which is
nonzero. Therefore the assertion about other source fibres must not be
transferred to all fibres on W2.

The complete source-polynomial witness is

    g^2G=L^2[2(d-2au)w-k]/(6d^2).

After multiplication by g, G still has a nonzero simple pole along w=0
but is regular along L=0. A correction from C(g) that cancels the former
would create a pole on the latter. Thus, in the original response module
C[x,t]/D_F C[x,t], the unit has exact scalar annihilator (g^2).

## 4. Principal parts, the cyclicity wall and exact left annihilators

Retain w=1+u^2t when expanding at u=0; use the actual (u,w) chart at w=0.
The complete scalar principal parts along (u=0,w=0,L=0) are

    theta=A/g^2+B/g,
    A=((2d-k)(d+k)^2/(6d^2), -k^3/(6d^2), 0),
    B=(-ak/d^2, -ak/d^2, 0).

The zero entry along L=0 is essential before taking the quotient by the
diagonal. In the quotient identified using this third entry, the two
coefficient vectors have determinant

    -ak(2d+3k)/(6d^2).

The inherited torsion theorem applies because the source is a submersion,
the rational constants have been established and exactly one source
fibre has three components. Thus the full torsion has two complete arms.

Off 2d+3k=0, g theta=A/g and (g nabla+2)theta=B/g recover independent
arms. The unit generates the entire torsion under the Weyl algebra
D=C<g,nabla>/([nabla,g]-1). Modulo the left ideal Dg^2, every operator has
form a(nabla)+b(nabla)g. Its B coefficient acting on theta forces a=0,
because distinct derivatives of 1/g are independent; its A coefficient
then forces b=0. The exact left annihilator is Dg^2.

On the global-submersion locus the determinant wall is d=6ah. It is
allowed, since it does not equal the excluded d=4ah. Here A1=A2!=0 and
B=lambda A with lambda=6a/k^2. The unit still has scalar order two and
generates one full torsion arm, since g theta=A/g, but it does not generate
the second arm. Its exact Weyl left annihilator is

    Dg^2 + D(g nabla+2-lambda g).

For completeness, the second generator equals 1+(nabla-lambda)g after
normal ordering. Modulo these two left-ideal generators, an arbitrary
a(nabla)+b(nabla)g reduces to [b-a(nabla)(nabla-lambda)]g. This acts as the
same polynomial in nabla on A/g, so it vanishes only when that polynomial
is zero. This proves equality of the annihilator, not just two exhibited
relations.

The wall is therefore a concrete globally submersive hostile to inferring
two-arm cyclicity from scalar unit order two alone. All module statements
refer to the original response ring. G is rational with the stated poles;
no global regular mate or general Jacobian-conjecture conclusion follows.

A different allowed locus, d=k/2=-2ah, kills only A's first entry. The
w=0 leading entry and the shared simple coefficient remain nonzero, and
the determinant remains nonzero. Thus dropping one component's leading
pole does not drop the unit order or two-arm cyclicity. This is distinct
from both the d=6ah cyclicity wall and the d=2ah gcd-power jump.

## 5. General coefficient-span explanation

The final addition in the frozen primary is also valid in the inherited
principal-part model. For v=sum_(j=1)^N A_j g^-j, the Euler operator
E=g nabla acts diagonally with eigenvalue -j on level j. Consequently

    product_(ell!=j) (E+ell)/(ell-j)

isolates exactly A_j g^-j, including when that coefficient is zero.
Multiplication by g^(j-1) and successive derivatives then generate the
entire A_j direction tensored with J. Conversely the Weyl action never
changes the constant coefficient span. Thus Dv=span{A_j} tensor J for
every finite packet, with no independence assumption on its coefficients.

If all N coefficients are independent, g^N kills v and its exact Weyl
left annihilator is Dg^N. In a PBW remainder sum_(j=0)^(N-1)
p_j(nabla)g^j, the A_1 coefficient first forces p_0=0. After that the
A_2 coefficient forces p_1=0, and so on; each step is a polynomial in
nabla acting on 1/g, whose distinct derivatives are independent. This
is an all-order proof, not an extrapolation of the finite interpolation
controls. It is correctly scoped as an elementary structural explanation
inside the inherited model.

For the geometric family the leading vector A is always nonzero, so
nabla^j theta has leading coefficient (-1)^j(j+1)! A at level j+2.
This verifies the claimed exact scalar order j+2 on and off the wall,
and explains why that statistic cannot detect the change in arm count.
