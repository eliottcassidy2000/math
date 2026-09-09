# Squarefree quadratic logarithmic pencils cannot produce a global W2 submersion

**Status: PROVED ANALYTICALLY + FINITE-EXACT; INDEPENDENTLY AUDITED.** No bound is imposed on
the degree of the first function or on the polynomial coefficients of the
witness. The hypothesis is squarefreeness of its generic discriminant in
the original x-coordinate. Repeated discriminants remain outside this theorem.

The [independent referee](continuing15_20260909_squarefree_quadratic_audit.md) accepts this result;
this filed copy changes only audit status and adds this link.

## 1. Statement and inheritance

Let F,H belong to C[x,t], F nonzero, and put

    H=A(x)t^2+B(x)t+C(x), A!=0,
    J(F,H)=F_x H_t-F_t H_x=F.

Write v for an independent variable, E=B^2-4AC, and D=E+4Av. Suppose D
is squarefree as a polynomial in x over C(v). Then F cannot extend to an
everywhere submersive regular function on the fixed W2 with charts

    x=1/r, t=-r^2-r^4 b.

In particular, a quadratic polynomial repair witness for a nonzero order-one
unit of a global W2 submersion must have a repeated generic discriminant.
This is a necessary condition, not a classification of that remaining case.
The first function can have arbitrarily high source t-degree.

The inherited closest mechanism is [the source-linear witness theorem](continuing14_20260908_logarithmic_linear_witness.md).
The hostile is the globally regular but boundary-critical F=t(1+xt), for
which adding lambda*F to H=-xt produces a quadratic witness. Another hostile,
F=(x+2t)(t^2-x^2/4-1), has a genuinely quadratic minimum witness, is a
source submersion, and has a pole on W2. Source submersion alone therefore
does not prove the theorem. The unused sidecar is the divisor of dx/y on
the normalized generic H-level, followed by the actual W2 boundary chart.

## 2. From the logarithmic equation to a complete conic alternative

On H=v, let y=2At+B. The function field has y^2=D and

    dF/F=dx/y.                                           (1)

Indeed H_t=y and differentiation at fixed H converts J(F,H)=F into
y partial_x F=F. Work on a smooth projective completion over the algebraic
closure of C(v), or separately on its components in the degree-zero case.
A logarithmic differential has only simple poles with integer residues:
if F=u^m f(u), f(0)!=0, then dF/F=m du/u+df/f.

Because D is squarefree, dx/y has no finite poles. If deg_x D>=3, it has
no poles at infinity either: direct local substitution at each infinity
gives order at least zero. A rational function with holomorphic logarithmic
derivative on a proper curve has neither zeros nor poles, so is constant;
this contradicts (1). Degree zero or one gives a pole of order two at
infinity and is also impossible. Thus deg_x D=2.

Write D=d2(v)x^2+d1(v)x+d0(v). The two infinity residues are
plus/minus 1/sqrt(d2(v)). They must be nonzero integers, so d2 is a
nonzero constant square 1/n^2, n a nonzero integer. Consequently

    deg A<=1, deg E=2, leading coefficient(E)=1/n^2.      (2)

This exhausts the squarefree case. It did not assume the normalization
had positive genus: the conic case is the necessary surviving alternative.

## 3. Constant A: the complete source-submersive normal form

Suppose A=a is a nonzero constant. Put s^2=leading coefficient(E),
E=s^2 x^2+e1 x+e0, and

    U=2at+B+s x+e1/(2s),
    V=2at+B-s x-e1/(2s),
    h0=(e1^2/(4s^2)-e0)/(4a).

These are polynomial coordinates: U-V recovers x and U+V then recovers t.
One has

    H=h0+UV/(4a), J(U,V)=4as,
    J(U,H)=sU, J(V,H)=-sV.

Comparing monomials in C[U,V], every nonzero eigenpolynomial is
U^n P(H) after possibly interchanging U,V and changing the sign of s,
with n a positive integer and ns=1. If F is a source submersion, its
vertical factor U^n has exponent one, hence n=s=1. It is then exactly

    F=U P(H), P nonzero and squarefree, P(h0)!=0.          (3)

Necessity of the last condition follows at U=V=0, and repeated roots of
P give repeated components. Sufficiency follows directly: on U=0 the
gradient is P(h0)dU; at another zero H=gamma the coordinates U,V are
nonzero and that level is smooth. Away from F=0 the bracket prohibits a
critical point. Scalar factors of F are included in P.

This classification is polynomial, so no unverified cancellation of a
rational coefficient in H is involved.

## 4. Linear A: recover the old theorem in a new coordinate

Suppose A=a(x-alpha), a!=0. Divide B uniquely as B=2AL+b with L polynomial
and b constant, and set X=x-alpha, tau=t+L(x). This shear preserves the
original bracket, although it need not preserve the W2 completion.
Polynomial divisibility and (2) give

    H=aX tau^2+b tau+cX+d,
    c!=0.

At fixed tau, this is linear in X. Equivalently -H is a source-linear
witness for the ordered coordinates (tau,X). The complete source-linear
classification therefore applies. For clarity its short argument here is

    F=Y(tau)P(H), Y'/Y=-1/(a tau^2+c).

The two root residues must be +1 and -1 by source submersion and
polynomiality, as in the inherited theorem. Label the positive root rho.
Then

    c=-a rho^2, rho!=0, a rho=-1/2,
    h_plus=b rho+d, h_minus=-b rho+d,
    F=(tau-rho)/(tau+rho) P(H),                          (4)
    P squarefree, P(h_minus)=0, P(h_plus)!=0.

These conditions force b!=0 and deg P>=1; the expression is polynomial
because H-h_minus is divisible by tau+rho. They are also sufficient for
source submersion by the same two-branch check as the inherited theorem.
The reversed residue choice is included by relabelling rho. In particular
source submersion strengthens (2) to leading coefficient(E)=1.

## 5. Pay the actual W2 boundary

In (3), if deg B>=2, then U,V have the same leading term B at x=infinity,
and U P(H) has a pole of degree (2 deg P+1)deg B. If B is affine, U,V are
affine in x,t with different x coefficients. When U has nonzero x
coefficient, it has order -1 in r. Either H has a pole, or V is a nonzero
multiple of t and H tends to h0; P(h0)!=0 prohibits cancellation in the
latter case. Thus F has a pole.

The remaining case is U=2at+u0. If u0!=0 and P is nonconstant, H has a
pole and so does F. If P is constant, F is regular but its derivatives
in r,b vanish along r=0. If u0=0, then U has order two, H tends to h0,
and P(h0)!=0 gives F order two. Again the entire added divisor is critical.

In (4), if deg L=ell>=1, then tau has order -ell and H has a pole of
order 2ell+1, while the ratio tends to one. Hence F has a pole. For L
constant not equal to plus/minus rho, H has a simple pole and the ratio
has a nonzero finite limit; again F has a pole. For L=rho, the ratio has
order two and P(H) tends to the nonzero P(h_plus). Thus F is regular
with the entire added divisor critical. For L=-rho, the ratio has order
-2, whereas H-h_minus has order one: its leading coefficient contains
2a rho, which is nonzero. The simple root of P at h_minus therefore
leaves F with a pole of order one.

These cases are exhaustive and prove the theorem.

## 6. Failure boundary and exact controls

Repeated D can introduce finite logarithmic poles and is not covered.
For example F=x(xt+1), H=xt+lambda F^2 has source submersion and
D=x^2[1+4lambda(v+1)x^2]. This witness is gauge-equivalent to xt and
F is not globally regular. It disproves deleting repeated factors of D
before computing dx/y. In contrast the genuine quadratic family in
Section 1 disproves assuming every witness can be reduced to source degree
one by adding a polynomial in F. Neither hostile is a global submersion.

The companion exact program checks both coordinate normal forms, their
literal original brackets, full affine critical ideals, all boundary cases,
generic discriminants, integer residues and the two distinct hostiles.
It uses raising gates in normal and optimized Python. Its declared finite
sample supports the proof and is not an all-coefficient enumeration.
Run `python 04-computation/continuing15_20260909_squarefree_quadratic_witness.py`
and the same command with `python -O` after relocation into the repository.

General order-one obstruction, repeated quadratic pencils, higher witness
degree, and planar JC remain **OPEN**. No external priority claim is made.
