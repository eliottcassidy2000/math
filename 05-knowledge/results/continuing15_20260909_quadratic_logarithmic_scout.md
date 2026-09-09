# Quadratic logarithmic witnesses: exact pole corridor and sharp affine hostiles

**Status: PROVED ANALYTICALLY + FINITE-EXACT; INDEPENDENTLY AUDITED.**
There is no bound on the degree of the first function. The squarefree
generic-discriminant case has a complete raw polynomial-eigenfunction
existence criterion. The nonsquarefree case has exact local obstructions
and an explicit elliptic example showing why they are not sufficient.
There are source-submersive order-one examples whose minimum original
witness t-degree is exactly two; all examples here fail global W2
submersion. No global classification for the nonsquarefree case is claimed.

The [independent referee](continuing15_20260909_quadratic_scout_audit.md) accepts this result;
this filed copy changes only audit status and adds this link.

## 1. Inheritance and fixed coordinates

The closest proved mechanism is
[the all-degree source-linear witness obstruction](continuing14_20260908_logarithmic_linear_witness.md).
Its rational logarithmic derivative turns polynomial eigenfunctions into
signed root valuations. The least-used extension is normalization of the
generic quadratic level curve before testing its logarithmic poles.
The recovered Euler viewpoint is elementary Hamiltonian eigenfunction
separation; it does not assume a global semisimple action. The canonical
hostile is an affine source submersion acquiring an added critical divisor
or a pole under the actual fixed W2 transition

    x=1/R, t=-R^2-R^4 B.

The live concepts are witness gauge, fixed versus moving discriminant
roots, local residues, normalized curve genus, principal divisors, and
original-chart boundary poles. The source response ring is always
C[x,t]/D_F C[x,t], with D_F=J(F,-). No global-ring preservation is assumed.

We seek a nonzero polynomial F satisfying

    J(F,H)=F, H=A(x)t^2+B(x)t+C(x), A!=0,
    J(F,H)=F_x H_t-F_t H_x.                              (1)

F is automatically nonconstant. Set v=H, K=C(v),

    E=B^2-4AC, y=2At+B, D(x,v)=E(x)+4A(x)v.

The generic level has y^2=D. Along it, (1) is exactly

    dF/F=dx/y.                                          (2)

All pole tests below take place on the complete normalized geometric
generic curve, not on its singular affine equation. On each geometric
component F is nonzero and nonconstant, by (1). A rational logarithmic
derivative has only simple poles with nonzero integer residues, and a
nonconstant rational function on a complete curve has a zero or pole.
These observations also apply if one initially allows a disconnected
geometric generic level.

## 2. Exact fixed-root and infinity tests

Every repeated generic root of D is a fixed common root of A and B.
Indeed, away from A=0, simultaneous equations D=D_x=0 imply
(E/A)'=0. Unless E/A is constant, their x-coordinate is a constant root
of a nonzero polynomial over C, incompatible with v being transcendental.
If E/A is constant, all roots of D are roots of A and hence fixed anyway.

At a fixed common root alpha, write a=ord_alpha A and b=ord_alpha B.
Because v is transcendental,

    m=ord_alpha D=min(a,2b).                            (3)

Thus the full generic square factor is obtained by taking exponent
floor(m/2) at each such root; every remaining moving root is simple.

For m=2e even, there are two normalized points, with local coordinate
z=x-alpha. The differential dx/y has pole order e when e>0. For
m=2e+1 odd, use x-alpha=u^2: then y is a unit times u^(2e+1), so

    dx/y = unit * u^(-2e) du.

The pole order is 2e, not 2e-1. In particular m=3 already gives a
double pole, which cannot be logarithmic. Therefore m can only be
0,1,2. The cases m=0,1 have no finite pole.

If m=2 and a=2, the residue contains a nonconstant square root of an
affine function of v, so it cannot be an integer. The only possible
finite pole is consequently

    ord_alpha A>=3, ord_alpha B=1,
    residues = +/-1/B'(alpha),
    1/B'(alpha) in Z minus {0}.                         (4)

These are necessary conditions, not a sufficient integrability criterion.

Put d=deg_x D=max(deg A,deg E). At infinity, if d=2k, use x=u^-1;
the differential has order k-2. If d=2k+1, use x=u^-2; its order is
2k-2. Hence d=0 or 1 gives a forbidden double pole. For d=2, infinity
gives two simple poles with residues +/-1/sqrt(leading_x D). For d>=3
there is no infinity pole. Therefore:

* d=0,1 is impossible;
* if d=2, necessarily deg A<=1, deg E=2, and its leading coefficient
  is 1/n^2 for an integer n!=0;
* if d>=3, there must be at least one admissible fixed double root (4).

The last point follows because otherwise dx/y has no poles anywhere,
whereas (2) for nonconstant F must have a pole. It recovers the familiar
holomorphic obstruction when D is squarefree, and also shows that a
proportional pencil E=cA cannot give a nonzero polynomial eigenfunction.

If F is additionally a source submersion, (4) strengthens to

    B'(alpha)=+1 or -1.                                (5)

To prove this, the polynomial vector field delta=J(-,H) vanishes at
p=(alpha,-C'(alpha)/B'(alpha)). Its linearization has eigenvalues
B'(alpha), -B'(alpha): in (x,t) its diagonal entries are these values
and the upper off-diagonal entry is zero. Equation delta F=F gives
F(p)=0. Differentiation and dF(p)!=0 force 1 to be an eigenvalue of the
linearization. This proves (5); it does not exclude such a root.

## 3. Complete raw existence iff when generic D is squarefree

Assume D is squarefree over C(v). Then a nonzero polynomial solution F
of (1) exists **if and only if**

    deg A<=1, E=e2*x^2+e1*x+e0, e2=1/n^2,
    n a nonzero integer.                               (6)

Necessity is Section 2. For sufficiency choose n>0 with e2=1/n^2,
kappa=1/n, write A=a1*x+a0, and put

    U=y+kappa*x+(e1+4a1*H)/(2kappa).

This is a nonzero polynomial in x,t and satisfies

    J(U,H)=kappa*U.

Thus U^n is a nonzero polynomial solution. No smoothness or globality
conclusion is attached to this raw existence theorem.

The rational classification is also exact. Set

    V=y-kappa*x-(e1+4a1*v)/(2kappa).

On the generic level,

    U V = e0+4a0*v-(e1+4a1*v)^2/(4kappa^2) = N(v).

The quadratic polynomial D cannot have identically zero discriminant in
x: its v^2 coefficient would force a1=0, and its v coefficient then
forces a0=0, contradicting A!=0. Thus N(v)!=0. The generic conic is
geometrically integral and K(x,y)=K(U), by V=N(v)/U and the displayed
linear formulas for x,y. Its rational constant field is exactly C(H).
Every rational eigenfunction of eigenvalue one is therefore

    F=K(H) U^n,  K(H) in C(H)^*.                        (7)

Which expressions (7) are polynomial and source-submersive is an
additional question; denominators in K(H) may not simply be ignored.
The separately audited global squarefree-pencil exclusion addresses
that next step. The present iff is about raw polynomial existence.

## 4. Identical local residues hide an unbounded exact eigenvalue lattice

For every integer m>=1 consider the whole explicit family

    H_m=x^(m+2)t^2+xt, D=x^2(1+4v*x^m),
    Y=1+2x^(m+1)t, omega=dx/(xY).

The exact rational AND polynomial eigenvalue spectrum is

    {lambda in C: J(F,H_m)=lambda F for some F!=0}
       = [m/gcd(m,2)] Z.                               (8)

Here rational and polynomial functions are taken in the original source
field and ring, respectively; the zero eigenvalue includes constants.
All members have the same local finite residues +/-1 and no infinity
poles for omega. Their least positive permitted eigenvalue nevertheless
tends to infinity.

We prove necessity even after adjoining algebraic constants
L=overline(C(v)). Set R=(Y-1)/(Y+1). Then

    d log R=m omega,
    x^m=R/[v(1-R)^2].                                  (9)

The extension L(R,x)/L(R) has degree m by Eisenstein at R=0 and cyclic
automorphism sigma(x)=zeta*x, with R fixed. An eigenfunction of
eigenvalue lambda has residues +/-lambda, so lambda is an integer.
By (9), F^m/R^lambda has zero differential and is a nonzero element
of L. Also sigma(F)/F is constant, so is a character of the cyclic
group. Thus F=x^j f(R) for 0<=j<m and f in L(R). Comparing valuations
of F^m=cR^lambda at R=0 and R=1 gives

    lambda congruent j (mod m),  2j congruent 0 (mod m).

Hence m divides 2lambda. Conversely the least positive eigenvalue
ell=m/gcd(m,2) is attained by the following polynomials:

    m odd:  F_+=x^(m+2)t^2,   ell=m;
    m even: F_+=x^((m+2)/2)t, ell=m/2.

Let K0=1+x^(m+1)t. Negative ell is also attained polynomially:

    m odd:  F_-=t^m K0^(m+2)=H_m^(m+2)/F_+;
    m even: F_-=t^(m/2) K0^((m+2)/2)=H_m^((m+2)/2)/F_+.

Powers realize all positive and negative multiples, proving (8).

In particular eigenvalue one is possible only for m=1,2. In those
two cases every rational eigenfunction of eigenvalue one has the form
K(H_m)F_+, with K in C(H_m): the generic normalized curve is
geometrically integral and has rational constant field C(H_m).
Polynomiality forces K to be polynomial. A pole at any nonzero level
cannot be canceled by F_+, which vanishes only over H_m=0. A pole at
zero cannot be canceled on the component K0=0, which F_+ omits.
Every polynomial eigenfunction of eigenvalue one therefore retains
the repeated x factor of F_+, and is source-critical. No source or
global submersion follows from the eigenvalue spectrum alone.

The m=3 case supplies an elementary genus-one explanation of the same
obstruction, useful without the cyclic-cover argument:

Consider

    H=x^5 t^2+xt,
    D=x^2(1+4v*x^3).

Its only finite pole candidate has A order five, B simple and B'(0)=1.
There are no infinity poles. Thus every necessary local condition above,
including the source-linearization test, is satisfied.

Nevertheless there is no rational eigenfunction F of eigenvalue one.
The normalized generic curve is

    Y^2=1+4v*x^3,

which is smooth of genus one: the double cover of P1 has the three
simple finite branch points and infinity. The differential dx/(xY)
has only two simple poles, at P_+=(0,1) and P_-=(0,-1), with residues
+1 and -1. If it were dF/F, the divisor of F would be P_+-P_-.
That would give a degree-one map from the genus-one curve to P1,
which is impossible.

The near miss has a precise torsion mechanism:

    d log((Y-1)/(Y+1)) = 3 dx/(xY).

The ratio has divisor 3(P_+-P_-); it is not a rational first root of
that divisor. In the original source one indeed has

    J(x^5 t^2,H)=3*x^5 t^2.

This example identifies the missing global principal-divisor/Pell
coordinate after all local poles and residues have been paid. Local
integral residues alone do not prove a logarithmic primitive. Formula
(8) makes the size of this missing divisor obstruction explicit and
unbounded within one quadratic-witness family.

## 5. Gauge thickening does not remove every admissible double root

Let lambda!=0 and set

    F=x^2t+x=x(xt+1),
    H=xt+lambda*F^2.

Then J(F,H)=F, and F is a source submersion. Its two zero components
x=0 and xt+1=0 are disjoint and reduced; the unit has exact order one.
Here

    A=lambda*x^4,
    B=x+2lambda*x^3,
    C=lambda*x^2,
    D=x^2[1+4lambda*(1+v)*x^2].

This is an actual source-submersive example with an admissible double
root and B'(0)=1. It refutes eliminating every such root merely from
source submersion or polynomial divisibility. The witness is a gauge
thickening: subtracting lambda F^2 leaves the t-linear witness xt.
The first function has an x pole under the W2 transition and is not
globally regular.

There is also a squarefree gauge hostile:

    F=t(1+xt), H=-xt+lambda F.

This satisfies (1), with A=lambda*x, and extends globally as a first
function; its entire added divisor is critical. Raw witness degree
alone is therefore not invariant under adding polynomials in F.

## 6. Genuinely minimal quadratic witnesses with source unit order one

The gauge issue does not imply all quadratic witnesses reduce to
t-linear witnesses. Put

    H=t^2-x^2/4, U=x+2t, V=x-2t, H=-UV/4,
    F=U P(H),

where P is any nonconstant squarefree polynomial with P(0)!=0.
Then J(U,H)=U, so J(F,H)=F. The degree of F in the ORIGINAL t coordinate
is 2 deg P+1, with no upper bound.

Every such F is a source submersion. On F=0 the component U=0 has
H=0 and nonzero normal derivative P(0). Every other component H=gamma
comes from a simple nonzero root gamma of P, is a smooth hyperbola,
and is disjoint from all other zero components. Away from F=0, equation
(1) itself prevents a critical point.

The exact rational field is C(F)(H), since

    U=F/P(H), V=-4H P(H)/F,
    x=(U+V)/2, t=(U-V)/4.

Its rational constant field for D_F is C(F). The primitive H/F has
complete simple principal parts zero on U=0 and gamma/F on each
component H=gamma. At least one gamma is nonzero. A scalar function
of F cannot cancel these different labelled parts, so there is no
polynomial mate. Thus [1]!=0 and its exact scalar annihilator is (F).

Every other polynomial repair H' satisfies H'-H in C(F) intersect
C[x,t]=C[F]. The intersection equality follows from Bezout for a
reduced rational expression p(F)/q(F): polynomiality forces q(F) to
be a unit. Any nonconstant polynomial in F has t-degree at least three,
so it cannot cancel the quadratic leading term of H. Therefore the
minimum original t-degree among ALL polynomial repair witnesses is
exactly two. This is not merely a high-degree representative obtained
by adding a polynomial in F.

For P(H)=H-1 the first function is a concrete cubic, with unit parts
(0,1/F). The linear change to (U,V) recovers an Euler eigenfunction,
but changes the source coordinate used by the fixed W2 charts. Indeed,
if m=deg P, the original first function has leading boundary pole

    F_infinity ~ leading(P)*(-1/4)^m R^(-(2m+1)).

Thus this entire family is non-global on W2. The example proves that
the quadratic-witness operation is real while leaving the intended
global unit-order-one problem open.

## 7. Verification and stopping scope

The independent controls check the exact quadratic-level identity,
normalized pole exponents, complete conic generator formulas, the
elliptic divisor hostile, both gauge families, original critical ideals,
the nonlinear witness family and its exact boundary pole orders.
Finite degree samples support the unbounded proofs; they do not supply
the quantifiers. No producer engine is imported.

Reproduce with `python continuing15_20260909_quadratic_logarithmic_scout.py`
and the same command with `python -O`. Certificates relocate from
04-computation into 05-knowledge/results; outside sources write beside
themselves. The nonsquarefree global W2 corridor remains unresolved.

The frozen engine passes 520 always-active gates. Normal and optimized
Python runs have byte-identical LF output and certificate. The finite
spectrum controls cover m=1..8 polynomial generators of both signs and
all character congruences for m=1..12, -2m<=lambda<=2m.
