# Independent audit: the all-m alternative chart and exact response order

**Verdict: PASS, with no mathematical repair.** The frozen primary
`continuing10_20260908_dg_unbounded_torsion.md` and its source/certificate
have been compared against the independent argument and controls below.
The unbounded proof is independent of the finite m bank. Work over C, with
m>=1 and a*h!=0. General JC(2) remains OPEN.

The comparison also checks the primary's m=2 constant convention: its full
polynomial part contains6h^2, so B0_old=B0+6ah^2 recovers the previous
source-linear family. The special value and boundary intercept must use
the same convention. No source translation is asserted to be a global
automorphism of W_m.

The chart supplier is the proved and independently audited
`planar_jc48_sep08_dg_genus.md`. The original m=2 source family and its two
special components are in `continuing10_20260907_dg_linear_carrier.md`.
The existing unit-response mechanism is retained from
`planar_jc48_sep06_torsion.md` and `planar_jc48_sep08_unit_torsion.md`.
The direct annihilator proof below does not need the component-jet isomorphism.

## Definitions, global source and the first jet

Put N=2m, e=N-1 and z=x-h. The two charts of W_m are

    U0=A2_(x,t), Uinf=A2_(r,b),
    x=1/r, t=-r^m-r^(2m)b,
    omega=dx wedge dt=r^(2m-2) dr wedge db.

Define Q_m to be the polynomial part of (x-h)^(2m)/x^m, namely

    Q_m(x)=sum_(j=0)^m binom(2m,j)(-h)^j x^(m-j),
    A=a z^(2m), B=a Q_m+B0, F=A t+B,
    c0=B(h), g=F-c0.

The complete second-chart expression is

    F=B0-a(1-hr)^(2m)b
      -a sum_(j=m+1)^(2m) binom(2m,j)(-h)^j r^(j-m).

It is polynomial on the entire chart. In particular F_b at r=0 is -a.
On U0, F_t can vanish only at x=h. There F_x=B'(h), and

    Q_m(h)=(-1)^m binom(2m-1,m) h^m,
    Q_m'(h)=(-1)^(m-1) binom(2m-2,m-1) h^(m-1) !=0.

These identities are valid for every m>=1. For completeness,
sum_(j=0)^k (-1)^j binom(n,j)=(-1)^k binom(n-1,k). Apply this first with
n=2m,k=m. For the weighted derivative sum, use

    (m-j)binom(2m,j)
      =m[binom(2m-1,j)-binom(2m-1,j-1)]

and the same partial-sum identity twice. The difference of the two adjacent
coefficients is binom(2m-2,m-1)/m. The m=1 boundary has derivative1, so it
requires no separate limiting argument. Thus F is a submersion on all W_m.

## The actual alternative affine plane

Let E1={x=h}, a closed copy of A1 in U0, and put u=1/(x-h). On Uinf,
u=r/(1-hr), so u is regular on W_m minus E1, including the boundary r=0.
The morphism (F,u) identifies this open surface with A2_(c,u).

Here is the full inverse, with no torus-only omission. On u!=0 use

    x=h+1/u,
    t=u^(2m)[(c-B0)/a-Q_m(h+1/u)].

On 1+h u!=0 use

    r=u/(1+h u),
    b=(B0-c)(1+h u)^(2m)/a
      -sum_(j=m+1)^(2m) binom(2m,j)(-h)^j
         u^(j-m)(1+h u)^(3m-j).

Both formulas are regular on their declared principal opens, including u=0
in the second chart. The opens cover A2 because (1+h u)-h u=1. Substitution
recovers F=c, and the two expressions satisfy x=1/r and
t=-r^m-r^(2m)b on their overlap. Conversely u=1/(x-h) or r/(1-hr) recovers
the original chart point. This is an isomorphism of open surfaces, not merely
a field identification. The old boundary D is precisely u=0 in this plane.

Every fibre F=c with c!=c0 is now an A1, since it misses E1 and is a coordinate
line in the alternative plane. The special fibre has exactly two components:

    F-c0=a z K,
    K=z^(2m-1)t+[Q_m(h+z)-Q_m(h)]/z.

E1 and E2={K=0} are disjoint because K at z=0 is Q_m'(h)!=0. Both are reduced.
E1 is A1_(t); E2 is the full c=c0 coordinate line in the alternative plane,
also A1. In the original U0, E2 and each nonspecial fibre are G_m; the u=0
boundary point fills each to A1. Omitting that point would change the result.

## Exact annihilator in the original polynomial response module

Let D_F=F_x partial_t-F_t partial_x and

    C_F=C[x,t]/D_F(C[x,t]), theta=[1].

The action of C[F] is well-defined because D_F(F)=0. The claim is the exact
ideal equality

    Ann_(C[F])(theta)=((F-c0)^(2m-1)).

Indeed C(x,t)=C(F,x), because t=(F-B(x))/(a z^(2m)). On this rational field
D_F=-a z^(2m) partial_x with F fixed, so its constant field is exactly C(F).
An explicit rational primitive is

    G0=1/[a e z^e], D_F G0=1.

Therefore every rational solution of D_F H=P(F), for P polynomial, has the
form H=P(F)G0+R(F), with R rational. Suppose H is polynomial in x,t.

First R can have no finite pole. At c!=c0 the divisor F=c is irreducible and
has x-h invertible at its generic point, since A and B-c are relatively prime.
A pole of R at c cannot cancel against P(F)G0 there. At c=c0 use E2: it is a
distinct reduced component, z is invertible there, and again G0 is regular.
Thus no pole at c0 can cancel either. Hence R is polynomial in F.

Now use E1. Its valuation of g is exactly1 because K(0)!=0. Consequently

    v_E1(P(F)G0)=ord_(c0)(P)-e.

Since R(F) is regular there, polynomiality of H forces ord_(c0)(P)>=e.
Conversely

    g^e G0=(a K)^e/(a e)

is polynomial, and its D_F derivative is g^e. This proves both inclusions
in the annihilator ideal, for every m, with no degree bound on a putative
primitive and no finite search. In particular theta!=0 and no polynomial
unit-Jacobian mate exists on the original source U0. A global mate on W_m
would restrict to such a polynomial, so none exists there either.

This quotient is explicitly the original C[x,t] response module. For m>=2,
omega vanishes on D, and D_F need not preserve O(W_m); the statement does not
silently replace its source module by the global coordinate ring.

The primary's entire canonical-derivative assertion is also accepted. Its
supplier `planar_jc48_sep06_torsion.md` requires a unit polynomial gradient
and rational constants C(F); both have been proved above. The unit primitive
G0 has affine poles only on E1 at c0. Its top scalar coefficient in g is
alpha=(a Q_m'(h))^e/(a e), which is nonzero; its principal part on E2 is zero.
Thus its component tuple is nonzero modulo the common diagonal. The inherited
canonical connection differentiates these scalar principal parts. At derivative
order j the highest coefficient is

    alpha*(-1)^j*e*(e+1)*...*(e+j-1),

nonzero in characteristic zero, and its pole order is e+j. Lower powers
cannot cancel it; no other fibre support is introduced. Consequently the full
annihilator of nabla^j theta is exactly((F-c0)^(e+j)) for every j>=0. The
singular h=0,m>=2 control below is not fed into this smooth-source theorem.

## The power map and the equality of orders

On the alternative plane, the same rational primitive is the polynomial
G0=u^e/(a e). A direct source Jacobian calculation gives

    omega=(u^(e-1)/a) dF wedge du,
    dF wedge dG0=omega.

Thus (F,G0) on W_m minus E1 is literally the finite flat power map of degree e.
The source ring is free over C[F,G0] with basis1,u,...,u^(e-1), using the
monic relation u^e=a e G0. If m>=2 it has ramification index e along D={u=0};
the ramification divisor and omega each have vanishing order e-1. These two
integers must not be confused. The exact response order is e, the same as the
degree and ramification index in this family.

For m=1, e=1: the map is an affine coordinate rescaling and is unramified,
although the original response class is nonzero and has order1. Therefore
nonzero response torsion does not require positive ramification. In every case
G0 has its order-e pole on the omitted E1, so the power map is not a morphism
on all W_m or a polynomial map from the original U0.

The unbounded orders vary the actual surface parameter m. They do not show
unbounded order inside a fixed W2 or classify every linear source function.

## Sharp multiplicity boundary and exact controls

The h!=0 condition has a decisive boundary when m>=2. Set h=0. Then

    F-B0=a x^m(1+x^m t).

The source has a critical line x=0, and this special component has multiplicity
m. The alternative chart and degree-e power map still exist. But the same
valuation proof gives

    Ann_(C[F])([1])=((F-B0)^2), since ceil((2m-1)/m)=2.

The residual component1+x^m t=0 still forces R(F) to be polynomial. Along x=0
the obstruction is now m*ord_(B0)(P)-(2m-1), rather than ord(P)-(2m-1).
This proves the claimed order2 and explains the failure of a degree-only
transfer. For m=1,h=0, the first jet remains a and the order remains1.

The independent source reconstructs Q_m by polynomial division, then checks
the complete symbolic source, both chart expressions, both inverse patches,
the cover and transition, first jet, special factorization, form identity,
rational primitive, polynomial annihilator witness and h=0 boundary for
m=1,...,7. A separate integer reconstruction checks both binomial identities
through m=100. The finite bank validates formulas; the arguments above prove
the unbounded quantifiers and absence of every possible polynomial primitive.

After filing, run:

    python 04-computation/continuing10_20260908_all_m_torsion_audit.py
    python -O 04-computation/continuing10_20260908_all_m_torsion_audit.py

The source imports SymPy but no other mathematical producer and writes its
certificate to results. Both modes pass445 always-active exact gates and emit
identical raw LF stdout; certificate regeneration is unchanged.

    Audit source 5d9dcf4ec3d3db6ceac66de14a436e292b5dd66be5a7441f19c6407a7d3076fc
    Audit output 540dd0a4d8a1cb45f43c7b8b6aebeeda9d1df0f3928ba28ace2255c6cf7c12b3
    Audit certificate 87c8f51d8a39ba3dacba9e7f7ace515f370ca969ee9589ee9f50065dc8ab925c
    Primary source 2c6e575521e44b971877455d1984e4286611c0df2e15bec73ef17941d6a6bd49
    Primary output 35272f408c20715cd0b84df71ecc6a48bf92c841a5a2c103f0077ed0755cfd53
    Primary certificate 3517b35f819399b7aaca886bc58cbac6653506f9b387c2409d708c3772fa6ef9

Root owns final promotion, maintained routing and Git.
