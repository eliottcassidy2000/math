# Independent audit: the complete 5+1+1+1 quartic boundary exclusion

**Audit status: PASS — complete analytic/source review, independent
section-space and inverse-coefficient reconstruction, and normal,
optimized and frozen replay agreement.** September 8, 2026. Root owns
primary status promotion.

Primary: [fivefold plus three simple roots](planar_jc48_sep08_five_one_one_one.md).
Source: [exact controls](../../04-computation/planar_jc48_sep08_five_one_one_one.py).
Output: [frozen replay](planar_jc48_sep08_five_one_one_one.out).

No producer file was modified during this audit. The final audit includes
the authorized nine-gate extension giving rational mates inside the same
boundary partition. All earlier proof arguments remain valid.

## 1. Accepted statement and precise boundaries

In the fixed DG surface

    W=(P1_x x P1_z) minus S,  S={z=x^2},
    t=1/(z-x^2),  omega=dx wedge dt,
    x=1/r,  t=-r^2-r^4 b,  omega=r^2 dr wedge db,

let H in L2, L in L1, deg_t H=2 and F=H^2+L. Suppose H's nonzero
complete binary-octic boundary section has four distinct zeros of
multiplicities 5,1,1,1. The accepted conclusion is that no G in C[x,t],
of any degree, has nonzero constant Jacobian with F. G need not extend
to W.

If all four points are finite, the proof first shows that any **rational**
mate would force L constant, then uses polynomiality for the final
contradiction. The new explicit rational family proves that this
distinction is essential even inside the exact global coefficient class.
If one boundary point is infinity, the separate leading differential
already excludes rational mates. There is no arbitrary projective
normalization or unproved surface automorphism.

The three location cases are exactly the following, because the missing
degree of the finite polynomial N is the multiplicity at infinity:

| Boundary placement | Finite polynomial degree and partition |
|---|---|
| All four finite | degree 8, 5+1+1+1 |
| Fivefold point at infinity | degree 3, 1+1+1 |
| One simple point at infinity | degree 7, 5+1+1 |

These exhaust all distinct-point placements. JC(2) is not claimed solved.

## 2. Complete global section spaces

I independently recovered the full coefficient spaces from their
literal bidegree section boxes. For H use every numerator monomial
x^i z^j with 0<=i<=4 and 0<=j<=2, divide by s^2, and substitute
z=x^2+1/t. The fifteen columns give precisely

    H=N t^2+P t+Q,
    deg N<=8,
    P=2N8*x^6+2N7*x^5+sum_{i=0}^4 p_i*x^i,
    Q=N8*x^4+N7*x^3+(p4-N6)*x^2+(p3-N5)*x+q0.

Independent Fraction row reduction gives total rank fifteen and
boundary-restriction rank nine. The nine N coefficients, five p_i and
q0 are freely variable. Every induced coefficient relation was checked
on each literal section basis vector, establishing completeness without
a selected-generator or fixed-prefix assumption.

For L, the six numerator monomials x^i z^j with 0<=i<=2 and 0<=j<=1,
divided by s, give total rank six and restriction rank five. The complete
form is

    L=M t+R,
    M=sum_{i=0}^4 m_i*x^i,
    R=m4*x^2+m3*x+r0.

The second-chart calculation agrees: its only possible negative powers
are r^-2 and r^-1, cancelled by R2=m4 and R1=m3. In particular
**deg M<=4 and M=0 forces the entire L to be constant**. This last
implication uses globality; it is false for an arbitrary polynomial
that merely has t-degree zero. The primary keeps that distinction.

These calculations agree with the current proved
[DG filtration](planar_jc48_sep08_dg_quadratic.md). They pay every
coefficient used later without placing a bound on the unknown mate.

## 3. The second inverse coefficient and rational trace descent

Normalize the nonzero Jacobian to one by scaling G. Choose a^2=N and
K=C(x)(a). The currently PROVED
[leading-exactness theorem](planar_jc48_sep08_leading_exactness.md)
pays the formal inverse, field injection and coefficientwise derivation.
I checked the additional coefficient directly.

Set U=vT and solve

    (N U^2+P U v+Q v^2)^2+M U v^3+R v^4=1,
    U(0)=1/a.

The derivative in U at v=0 is 4a, which is a nonzero element of K.
The unique formal coefficients are

    T_-1=1/a,
    T_0=-P/(2a^2),
    T_1=(P^2-4a^2Q)/(8a^3),
    T_2=-M/(4a^2)=-M/(4N).

All P,Q terms are retained before their cancellation from T2 is proved.
R first enters the scaled equation at order four. I independently
reconstructed the expansion through v^3 using sparse Laurent-polynomial
arithmetic over fractions.Fraction, treating a,P,Q,M as independent
symbols and allowing negative a powers. No producer or SymPy definitions
were imported in that reconstruction.

For any rational G, substitution t=T is legitimate: a nonzero polynomial
denominator has a unique lowest Laurent term, from its highest t-degree,
and cannot vanish identically. The derivation d/dx extends uniquely to
the characteristic-zero algebraic field K and then coefficientwise to
K((v)). The rational chain rule gives

    partial_x G(x,T)=(1/4)v^5 partial_v T.

At v^6 only T2 contributes, with multiplier 2/4. Thus

    g6'=T2/2=-M/(8N),
    d(2g6)=-M dx/(4N),  2g6 in K.

This initially provides an algebraic primitive, not automatically a
rational function of x. The primary correctly pays descent:

    B=Tr_{K/C(x)}(2g6)/[K:C(x)],
    dB=-M dx/(4N),  B in C(x).

Every embedding into a normal closure commutes with the uniquely
extended derivation, so the derivation commutes with trace. The
normalized trace therefore retains the displayed differential. If N
is square, the extension degree is one; no disconnected presentation
is substituted for the field.

There is also an independent local check of the same simple-root
consequence. At a simple root q of N the rational residue is

    -M(q)/(4N'(q)).

On the radical normalization x-q=z^2 it doubles to -M(q)/(2N'(q)).
It cannot disappear under ramification. Either trace descent or this
doubled residue proves **M(q)=0** at each simple root. Each labelled
residue is used individually, not through its sum.

## 4. All local branches at the finite fivefold point

Assume all four boundary points are finite, with fivefold point p.
The full graph-coordinate equation is

    Ncal=N+sP+s^2Q,  Mcal=M+sR,
    E=Ncal^2+s^3 Mcal-cs^4=0,
    eta=s^2 dx/E_s,  dF wedge eta=omega.

The sign and actual multiplier agree: t=1/s gives
omega=ds wedge dx/s^2, and on E=0 one has F_s=E_s/s^4.
At a finite point the remaining coordinate changes multiply eta by
a unit. The infinity zero of the volume is not substituted here.

I re-read the relevant proved local suppliers,
[common-root necessity](planar_jc48_sep08_quartic_common_root.md)
and [shared-root first jets](planar_jc48_sep08_shared_roots.md).
They apply with the following complete scope.

If M(p)!=0, set u=x-p and j=ord_u P, allowing j=infinity. Since
ord_u N=5, the balanced equation 5=3j is impossible.

* For j=0, local Weierstrass degree in s is two. Both cancellation
  determinations have s-order five and eta a unit times u^(5/2)du.
  Their actual quadratic normalization is regular. The apparent third
  Newton root has s-order zero and does not pass through this point.
* For j=1, the Weierstrass degree is three. One simple branch has
  s-order two and regular form. The other two determinations have
  s-order four and eta a unit times u du. These exhaust the degree.
* For all j>=2, including infinity, the initial balance is
  N(u)^2+M(p)s^3=0. Its three nonzero Puiseux determinations have
  s-order 10/3 and simple leading roots. The term 3M(p)s^2 dominates
  E_s, so eta is a unit times du and remains regular on its cubic
  normalization. This is an unbounded analytic case, not a finite-j
  extrapolation.

If M(p)=0 but P(p)!=0, the implicit centre Ncal=0 has s=psi(u) of
order five. For generic c,

    ell=ord_u(Mcal(psi(u),u)-c psi(u))

lies between one and five. A lower order remains if present; otherwise
the nonzero leading term of psi supplies a nonzero generic c-slope at
order five. The two determinations have eta a unit times
u^((5-ell)/2)du, regular after either unramified or ramified
normalization. Weierstrass degree two pays exhaustion. No extra bound
on the higher order of M is assumed.

Every other boundary point is simple and regular for arbitrary Mcal.
At a shared simple point Ncal itself is an invertible tangential
coordinate; the generic leading equation is Ncal^2+unit*s^4=0 and
s^2 ds/E_Ncal is regular on its two smooth branches. Nonshared simple
points are covered by the M-unit supplier. The simple-root residue
gate already makes them shared for a hypothetical rational mate, but
the stronger stated local regularity is also correct.

There is no boundary zero at infinity in this all-finite case, since
N has degree eight with nonzero leading coefficient. Therefore the
compact generic fibre has no branch meeting S there.

## 5. Generic components, first jets and the five counted zeros

If either regular alternative at p holds, eta is holomorphic on every
compact normalized generic component. In W, choose a smooth generic
fibre and avoid the finitely many vertical exceptional values of a
proposed rational denominator. Horizontal pole divisors may meet that
fibre at points; they do not make the rational restriction undefined
identically on a generic component.

The relative form is regular in W, including D where the actual volume
vanishes to order two. Each generic component meets W away from D:
a component contained in D would be the fixed D at a fixed fibre value,
and S is not a component of the cleared generic equation because N|S
is nonzero. Thus eta is nonzero on every component. A rational primitive
cannot have a pole whose derivative is holomorphic; it would be a
global holomorphic function on the compact normalization, hence constant,
contradicting eta!=0. Neither geometric irreducibility nor nonconstant
F|D is required.

Consequently a rational mate would require

    M(p)=0, P(p)=0.

The tangential derivative of N already vanishes to order at least four,
so this is a singular shared point. The finite first-jet theorem requires
P'(p)=M'(p)=0. This use retains the exact finite-point geometry and
imposes no unproved higher-jet condition. Directly, the tangent equation
after s=uZ is

    E=u^4[(bZ+eZ^2)^2+dZ^3+(f-c)Z^4+...],
    b=P'(p), d=M'(p).

If b!=0, after factoring Z^2 a quadratic with nonzero constant term
b^2 has two simple nonzero roots for generic c. If b=0,d!=0, after
factoring Z^3 a linear factor has one simple nonzero root. Each gives
an actual branch with nonzero residue Z0^2/P_c'(Z0). One such branch
suffices; no exhaustive assertion about the Z=0 branches is needed.
Genericity excludes only finitely many exceptional c values.

Let C(x)=prod_i(x-q_i) for the three simple roots. The formal residues
give C|M, and the finite main-point conditions give (x-p)^2|M. Their
coprime product has degree five, while deg M<=4. Therefore M=0 and
the entire L is constant.

Equivalently, write M=C(Ax+B). The two p-jet equations in A,B have
matrix

    [ C(p)*p,          C(p)  ]
    [ C'(p)*p+C(p),    C'(p) ],

with determinant -C(p)^2!=0. The independent confluent Vandermonde
calculation has determinant

    d1^2*d2^2*d3^2*(d2-d1)*(d3-d1)*(d3-d2),

where d_i=q_i-p. Distinctness pays every factor. The affine shift used
to count these polynomial zeros is not asserted to extend to W.

For a polynomial G, the final identity is

    J(F,G)=2H J(H,G).

H has t-degree two and is a nonconstant nonunit of C[x,t]. A product
with this factor cannot be a nonzero scalar. This excludes polynomial
mates of every degree, even if G is not globally regular on W.

For rational G the necessary conclusion stops at L constant. A
rational derivative may contain a denominator H, so the factor argument
does not apply. The final sharp family below proves that this is an
essential distinction within the same partition.

## 6. Both infinity placements

The current leading-exactness and
[boundary-exactness classification](planar_jc48_sep08_boundary_exactness.md)
are PROVED and independently audited. For the actual leading coefficient
N^2, the degree-four formal theorem requires dx/sqrt(N) exact in
C(x)(sqrt(N)).

If the fivefold point is infinity, N has degree three and three
distinct simple finite roots. Its connected normalized radical curve
is elliptic and dx/sqrt(N) is nonzero holomorphic everywhere, including
its sole infinity point. It cannot be a rational derivative.

If a simple point is infinity, the finite partition is 5+1+1 in
degree seven. A possible primitive has only one pole, of order three,
above the fivefold point. The differential has order four at infinity,
so the primitive would have local degree five there, larger than its
total pole degree at most three. This is impossible.

Thus both infinity cases exclude rational mates, and together with the
all-finite proof they establish the claimed all-location polynomial
exclusion. These calculations do not use an unweighted projective
coordinate change or the delicate elliptic coefficient condition in
degree eight. The all-finite proof itself applies to every root position
in its partition, whether or not its first leading gate passes.

## 7. Literal controls and the same-partition rational family

The producer's first global control is

    N=x^5(x^3+1), P=2x^6, Q=x^4-x,
    H=N t^2+P t+Q.

Its exact second chart is b^2+2br+b^2r^3. All induced source rows are
present. The two global corrections

    L1=(x^3+1)t+x,
    L2=x(x^3+1)t+x^2

both have leading coefficients vanishing at the three simple roots.
The first still has M(0)=1; the second has M(0)=0 and M'(0)=1.
Thus neither three simple-root evaluations nor those evaluations plus
one main-point zero suffice. For L2 the literal tangent polynomial is
Z^3-cZ^4; the actual nonzero branch Z=1/c has residue exactly -1.
This explicitly tests the remaining first-jet obligation.

I independently verified the final, stronger rational control in the
same partition. For delta!=0 and arbitrary beta,gamma,c0, put

    h=x^2+x^4t,
    H=(x^4+delta*x)(1+x^2t)^2+beta*h+gamma,
    F=H^2+c0.

The actual second chart is

    H=(1+delta*r^3)b^2-beta*b+gamma.

The leading coefficient N=x^5(x^3+delta) has a fivefold zero at zero
and three distinct nonzero simple roots, since the cubic discriminant
is -27delta^2. Thus this is in the exact global boundary class.

Set u=x^-3 on the function field and keep h. Directly,

    J(u,h)=-3,
    H=(1+delta*u)h^2+beta*h+gamma.

Hence

    J(H,1/(3delta*h))=1,
    J(H^2+c0,1/(6delta*H*h))=1.

The second identity also follows by dividing the first primitive by
2H; the derivative of the factor 1/H contributes no bracket with H.
The displayed reciprocal denominators are nonconstant polynomials, so
the mates are rational and not polynomial. No cancellation converts
them into a polynomial.

The field map is genuinely degree three. The original field is
C(x,h), since t=(h-x^2)/x^4, while the subfield is C(u,h). The equation
x^3=u^-1 is irreducible there: u^-1 has valuation -1 at u=0 and cannot
be a cube, and over a field containing the cube roots of unity a
reducible cubic X^3-u^-1 would have a root. Thus the extension has degree
three. The volume is -du wedge dh/3. Neither a birational coordinate
change nor a surface automorphism is claimed.

This proves actual rational-mate existence inside the same partition,
while the theorem excludes every polynomial mate. It is stronger than
the earlier outside-class example F=t^4, G=-x/(4t^3), which remains
a valid simple test of the final factor inference.

## 8. Source, replay and final pins

I read the complete source, including all nine new sharp-family gates.
It uses exact SymPy arithmetic, no other producer imports, and explicit
exception gates that remain active under Python optimization. The
symbolic universe is the full fifteen-dimensional H space, the full
six-dimensional L space, all distinct complex positions and all three
location cases. No bound on a prospective mate or sampled coefficient
census substitutes for a proof.

The source checks the full inverse through the required order, complete
second-chart cancellations, all labelled simple-root residues, both
five-zero determinants, normalization order controls, the actual
linear-jet tangent, the infinity cases and the rational sharp family.
The bounded local controls illustrate the proved unbounded branch
arguments. In particular the j>=2 case and the degree-three subfield
claim are analytic, not inferred from a finite list or from a polynomial
degree check alone.

An independent standard-library computation rebuilt every section-box
column and both restriction ranks, checked all induced coefficient
identities, and reconstructed the formal inverse by Fraction Laurent
arithmetic. The determinant, trace, local residue and sharp-family
calculations were independently checked analytically above.

Reproduction from the worktree root:

    python3 04-computation/planar_jc48_sep08_five_one_one_one.py
    python3 -O 04-computation/planar_jc48_sep08_five_one_one_one.py

Both final independent replays completed successfully. Their **62-gate**,
**394-byte** outputs are byte-identical to the final frozen output.
The semantic gate hash is
de0cbef63edcb5a98cf214274d92dd91417c5d4bef999d43548d5dd5cd36fecc.

| Accepted artifact | Bytes | SHA256 |
|---|---:|---|
| Final source | 8,859 | 029193e4fc766ebbfc7bf49a8af1bbc928da6cd88ea04f77dd1067d784e9a6fd |
| Final frozen output and each final independent replay | 394 | c6ef3ec938cffd949ad4a755ef8732c52c0a5bafce736415eafb3992ae90f776 |
| Final primary before status promotion | 18,508 | 37bc9e846cca9b7855dd34794a54fb2bb87ac1055480b117b25a88502c93b0ad |

These pins supersede the initial 53-gate versions for this audit.
No mathematical correction remains. This audit is frozen for root's
promotion/checkpoint with the polynomial/rational distinction and all
location hypotheses preserved.

