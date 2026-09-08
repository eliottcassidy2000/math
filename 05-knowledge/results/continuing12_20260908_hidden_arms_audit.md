# Independent referee: fixed unit futures with unbounded ambient torsion

**Status: PASS / analytical proof and independent exact controls.** The
complete frozen primary report has been compared with the independent
derivation below. No mathematical repair is requested. The independent
verifier imports and executes no producer. Its finite controls are
separate from the all-parameter proof below.

Audited primary: `continuing12_20260908_fixed_unit_hidden_arms.md`, raw SHA256

    1f272fdc64da49cf965c5ce4584257eb9f09c6699b92aac7ddd9461b37cbb2c4.

Its source has SHA256
`d347e5a50e2247fcf2934d52252a0eaeae466c1c7b2d74951910fa0b207776bd`.
This audit targets that frozen candidate before any status-only promotion.

## 1. Actual global function and every critical point

Fix W2 with x=1/r_D,t=-r_D^2-r_D^4 b_D. For an integer r>=3 and nonzero
h,b,d, put

    u=x-h, z=u-2h+u^3t, y=uz, W=uz^2,
    K=b+dW^r, g=F=yK.

Here r is the family index and r_D is the added-chart coordinate. They
must not be identified. On the entire added chart set
R=h^2(3-h r_D)+b_D(1-h r_D)^3. Literal substitution gives

    z=-r_D R,
    y=-(1-h r_D)R,
    W=r_D(1-h r_D)R^2.

Thus F is globally regular. On D its tangential derivative is -b, so
every boundary point is noncritical. On u=0 its source derivative is
-2bh, also nonzero. On u!=0 the actual coordinate pair (u,y) has
Jacobian u^4 and

    F=by+d y^(2r+1)/u^r,
    F_u|y=-r d y^(2r+1)/u^(r+1).

If the last derivative vanishes, y=0 and F_y=b. This pays all source
points. No discarded critical locus or global derivation hypothesis is
being smuggled into this coordinate change.

## 2. Complete fibres, rational constants and the primitive

The special fibre is

    g=d u z product_(rho^r=-b/d) (W-rho),

where the product is the complete factorization of K. There are r+2
reduced disjoint source
components: E_u, E_z and the r curves L_rho:W=rho. Their disjointness
follows from z|u0=-2h, W=0 on E_u and E_z, distinct nonzero roots rho,
and K|W0=b.

The curve E_z is primitive linear in t. Every L_rho is primitive
quadratic with discriminant 4rho*u^7, nonsquare in C(u); its constant
coefficient at u=0 is -rho, so no vertical factor has been omitted.
Equivalently L_rho is isomorphic to G_m via
u=rho/z^2 and t=(z-u+2h)/u^3. Hence all the claimed factors are whole
irreducible components, not formal branches.

For every nonzero g=c, u,z,K are units. The inverse formulas are

    u=c^2/(W K^2),  z=W K/c,
    t=(z-u+2h)/u^3.

They give the complete fibre coordinate ring
C[W,W^-1,K^-1]. Therefore every other source fibre is irreducible,
and C(x,t)=C(g)(W). Computing in the original source gives

    J_(u,z)(y,W)=W,   J_(x,t)(u,z)=u^3,
    D_g W=K u^3W=g^6/(W^2K^5).

This verifies both the source orientation and the rational constant
field: a nonzero multiple of differentiation in W has constants exactly
C(g) in characteristic zero.

Define

    H(W)=sum_(m=0)^5 binom(5,m)d^m b^(5-m) W^(rm+3)/(rm+3).

Every denominator is nonzero and H'=W^2K^5. Consequently the single
rational function G=H(W)/g^6 satisfies D_gG=1. Also H=W^3S, so

    G=S(W)/(u^3K^6),    g^6G=H(W) in C[x,t].

This is a literal original-source polynomial response witness and a
global first function, not a pair of global regular functions.

## 3. Every scalar principal part and the exact unit order

At L_rho, W-rho is a local parameter and K has a simple zero.
Since H'=W^2K^5, H-H(rho) is divisible by (W-rho)^6. The local
parameter g differs from W-rho by a unit, so the complete principal
part is just H(rho)/g^6. All lower scalar coefficients vanish.
Furthermore

    H(rho)=b^5 rho^3 c_r,
    c_r=sum_(m=0)^5 binom(5,m)(-1)^m/(rm+3)
       =120r^5 / product_(j=0)^5 (3+jr) !=0.

The all-r identity follows by integrating t^2(1-t^r)^5 on [0,1]
and substituting t^r. Thus none of these leading coefficients vanishes.
The root-character rho^3 can repeat values when r is divisible by three;
distinct coefficients are neither assumed nor needed.

At E_z, u is a unit and K=b, so G is regular. This zero principal
part must be retained in the full component quotient.

At E_u, W=O(u), and S-K^6/(3b)=O(W^r). Hence
G=1/(3bu^3)+O(u^(r-3)), with a regular error precisely in the claimed
r>=3 range. Write g0=bu(u-2h). The actual source jet gives
g-g0=O(u^4), including the u^3t term in z. The exact identity is

    1/(3bu^3)
      -[-8b^2h^3/(3g0^3)-2bh/g0^2]
      =1/[3b(u-2h)^3].

The right side is regular along E_u. Replacing g0 by g changes these
principal parts only by regular terms, because g-g0=O(u^4). Thus the
complete E_u part is -8b^2h^3/(3g^3)-2bh/g^2. This uses the actual
source jet, not an invalid independent-coordinate expansion at u=0.

There are no other affine poles: the denominator u^3K^6 lists all of
them. Since g^6G is polynomial but g^5G has a genuine pole on L_rho
and is regular on E_z, no correction from the rational constants C(g)
can make g^5G polynomial. Therefore the unit has exact C[g]
annihilator (g^6).

## 4. Full ambient module versus every legal unit future

The original affine gradient is unimodular, the rational constants have
been proved, and exactly one source fibre is reducible with r+2
components. The inherited component torsion theorem therefore gives
r+1 complete torsion arms.

Set the E_z coordinate to zero in the quotient by the diagonal.
Let A have coordinates H(rho) on the r L components and zero at E_u;
let B be supported on E_u with value -8b^2h^3/3. Then A,B are
independent, and

    theta=A/g^6+B/g^3+lambda B/g^2,
    lambda=3/(4bh^2).

The coefficient-span theorem gives D theta=(C A+C B) tensor J,
where J=g^-1 C[g^-1] and D is the Weyl algebra [nabla,g]=1.
It has exactly two complete arms for every r>=3, even while the full
ambient torsion has r+1 arms. The j-th derivative retains scalar order
6+j, but adds no coefficient direction.

For fixed b,h, sending A and B to the corresponding vectors of any
other member identifies the **pointed unit-generated modules**. This
identifies the formal target parameter g and the distinguished unit.
It is not an isomorphism of the full ambient modules, their labelled
component spaces or the source surfaces with their first functions.
The parameter d and the integer r disappear from this pointed module,
although they change the surrounding geometric component data.

## 5. Exact annihilator of the pointed unit

Over C[nabla], use the free basis e_A=A/g,e_B=B/g. Let

    A_j=(-1)^(5-j)nabla^(5-j)/(5-j)!,   0<=j<=5,
    B_0=nabla^2/2-lambda*nabla,
    B_1=-nabla+lambda, B_2=1, B_3=B_4=B_5=0.

Then g^j theta=A_j e_A+B_j e_B. The pair g^2theta,g^5theta is a
C[nabla] basis: its matrix has determinant -1. In particular it is
unimodular, not merely generically independent over C(nabla).

For j=0,1,3,4 the exact relations are

    R_j=g^j-B_j g^2-(A_j+B_j*nabla^3/6)g^5.

The full Weyl left annihilator is generated by g^6 and these four
R_j. All kill theta. Conversely normal ordering modulo Dg^6 gives
sum_(j=0)^5 p_j(nabla)g^j. Reducing by the four displayed relations
leaves only coefficients of g^2 and g^5. The unimodular basis makes
both coefficients zero whenever the operator kills theta. This proves
the reverse containment for arbitrary finite operators.

The ideal depends only on lambda. It therefore certifies equality of
all legal unit Weyl futures across the family, not just equality of
scalar order or coefficient rank. No claim about every ambient response
follows. The unbounded invisible complement is exactly the missing
coordinate in that stronger proposed observation principle.

## 6. Independent finite controls

The separate verifier tests r=3,...,8 with full symbolic parameters;
whole original/added chart maps; complete K-root quotient algebras;
all five vanished lower L derivatives; the exact Beta coefficient;
the Eu correction threshold and rational remainder; complete nonzero
fibre inverses; and the entire polynomial PBW relation matrix.
It imports or executes no producer. These 149 always-active gates
support the mechanisms above; the all-r theorem is analytic.

Reproduce with the adjacent `continuing12_20260908_hidden_arms_audit.py`
normally and under `-O`. The certificate is written beside an external
source or into `05-knowledge/results` when filed in `04-computation`.
Both normal and optimized runs pass **149 always-active gates**, with
byte-identical LF output and regenerated certificate bytes. Frozen pins:

* source: `94b1e350a9beb12dd5fc23812b3174db2f762e34c14f25cdd78920fdf965258c`;
* output: `eac123f9c0824aef0daf6e575ec1f917df0389b7ff77a1e524c4455be392a70b`;
* certificate: `f7811e5b4b863c9912347103b2f33933de64f2c0d6f25546bb93ef68a5666ea7`.
