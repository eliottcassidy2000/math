# A polynomial-mate first-jet gate for every parabolic source chart

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**

The result is uniform in every integer m>=1 and every polynomial
degree. It concerns an exact subring of the original polynomial
source. It does not assert a rational-mate exclusion, a formal
completion theorem, or a general planar Jacobian conclusion.

## 1. Exact statement and inheritance

Let K be any field of characteristic zero. In the original source
coordinates (u,t), put w=u^m t for an integer m>=1. Suppose
Fhat in K[u,w], set F(u,t)=Fhat(u,u^m t), and assume a polynomial
G in K[u,t] satisfies J_(u,t)(F,G)=1. Then

    Fhat=f0+beta*u+u^2 A(u,w),    f0 in K, beta in K*.

Moreover, after expressing the actual polynomial mate as a finite
Laurent polynomial

    Ghat=G(u,w/u^m)=sum_l u^l g_l(w),

its least nonzero u-exponent is exactly -m, and

    g_(-m)=w/beta.

In particular G(0,t)=t/beta+gamma for some gamma in K.
A nonzero constant Jacobian other than one is handled by a constant
rescaling of G. These are necessary conditions, not a characterization
of all polynomial mates.

The closest inherited mechanism is weighted initial-bracket
cancellation, as in [THM-2102, power-free-weight-face-and-first-defect-descent](../../01-canon/theorems/THM-2102-power-free-weight-face-and-first-defect-descent.md).
That theorem uses positive weights and does not directly supply the
present assertion: the original weights here are (1,-m), and the
chart has a nonconstant Jacobian factor. The proof below is independent.
The [moving-root residue suppliers](planar_jc48_sep08_five_three.md)
and [their infinity continuation](planar_jc48_sep08_five_three_infinity.md)
retain more information appropriate to rational primitives; their
rational conclusions are not replaced by this polynomial gate.

The corrected near miss is to allow arbitrary Laurent coefficients
for a polynomial mate. The missing sidecar is the exact image of
K[u,t] inside K[u,u^{-1},w]. The canonical rational hostile
F=u*w, G=1/(m*u^m) shows why that sidecar is essential for every m.
The five live concepts are the original source ring, signed weights,
finite Laurent rows, the nonzero leading bracket, and a complete
coefficient image. No degree census or descent of a formal solution
is used.

## 2. Complete source-image lemma and the actual bracket

The substitution into the rational function field is injective.
An original monomial u^i t^k becomes u^(i-mk) w^k. Thus a finite
Laurent polynomial sum_l u^l g_l(w) is the image of a polynomial
in K[u,t] if and only if

    w^max(0,ceil(-l/m)) divides g_l(w) for every l.      (1)

Indeed the inverse exponent of u for a term u^l w^k is l+mk.
It is nonnegative exactly under (1). This is a complete monomial
bijection, so no unproved polynomiality criterion is used.
In particular every nonzero negative row g_l is divisible by w
and has positive degree.

The coordinate map is birational on u!=0, with actual Jacobian

    J_(u,t)(F,G)=u^m J_(u,w)(Fhat,Ghat).              (2)

It is not a polynomial automorphism of the whole source. Every
use of (2) retains the image condition (1). The row-to-row contribution
for Fhat=sum_i u^i f_i(w), i>=0, and Ghat=sum_l u^l g_l(w) is

    u^(i+l+m-1) [i f_i g_l' - l f_i' g_l].          (3)

All these expressions are finite Laurent polynomials. Equality to
one is therefore coefficientwise, with no convergence, formal-field
constant extension, or rational-time convention involved.

If every nonzero row of Ghat has exponent l>=0, then every nonzero
contribution in (3) has exponent at least m. The only possible pair
with i+l=0 is i=l=0, and its bracket coefficient is zero.
Hence a polynomial mate must have a negative least exponent L.
This explicitly handles leading constants and the m=1 corner.

## 3. The constant row of Fhat must be constant

Write f(w)=f_0(w). Suppose f is nonconstant. Let L<0 be the least
nonzero row of Ghat. The first nonzero term of the actual bracket
is

    -L f'(w) g_L(w) u^(L+m-1).                      (4)

Its coefficient is nonzero: characteristic zero gives L!=0,
f'!=0, and g_L!=0. Every higher row of Fhat or Ghat contributes
strictly greater u-order. Thus no higher correction can cancel (4).

For the bracket to be one, its exponent must be zero. If m=1,
that would require L=0, contrary to L<0. If m>=2, it requires
L=1-m, but g_L is divisible by w. Its product with the polynomial
f' cannot be a nonzero constant. This contradiction proves f in K.

The argument covers every degree of f, including linear f, and
does not remove its constant term before checking the source image.

## 4. The first row is an arbitrary-degree polynomial, not only affine

Subtract the constant f0 from F, which preserves the Jacobian, and
write

    Fhat=u g(w)+u^2 A(u,w).

If g=0, the original F-f0 is divisible by u^2 in K[u,t]. Both
original partial derivatives then vanish on u=0. A polynomial
Jacobian mate is impossible there. This is an original polynomial
critical-line argument; it is not applied to rational mates.

Suppose next that g is nonconstant, of degree h>=1. Again L<0.
By (1), j=deg g_L>=1. The first possible bracket term is

    [g g_L' - L g' g_L] u^(L+m).                   (5)

The coefficient in brackets has exact degree h+j-1. Its leading
coefficient is the product of the nonzero leading coefficients
of g and g_L times j-Lh. Since j>=1 and L<0, this positive
integer is nonzero in K. Thus (5) really is the first nonzero
term; all omitted rows have strictly greater u-order.

It could equal a constant bracket only when L=-m. Even then its
coefficient has degree h+j-1>=1, so it cannot be one. Therefore
every nonconstant polynomial g is excluded, not merely a row of
the form alpha*w+beta.

It remains that g=beta in K*. For any negative leading row the
coefficient in (5) is beta*g_L', nonzero because g_L is divisible
by w and is not zero. Equality to one forces

    L=-m,       beta*g_L'=1.

The only polynomial solution with g_L(0)=0 is g_L=w/beta.
This proves the stated full initial-row conclusion.

Finally a source monomial surviving on u=0 is t^k, of Laurent
weight -mk. The minimum weight -m forbids k>=2; its unique
k=1 coefficient is 1/beta. The remaining constant gives
G(0,t)=t/beta+gamma. No higher or negative row can alter it.

## 5. Exact scope, positives, and hostiles

For every m, F=f0+beta*u and G=t/beta+B(u), with B any polynomial,
give genuine polynomial positive controls and the asserted initial row.

The necessary rows are not sufficient: F=u+u^2 satisfies them
but has an original critical line u=-1/2, so it cannot have a
polynomial mate. No stronger sufficiency assertion is made.

For every m>=1,

    F=u*w=u^(m+1)t,       G=1/(m*u^m)

satisfies J(F,G)=1 rationally. Its first row g(w)=w is nonconstant.
The putative mate's leading coefficient is the constant 1/m in
weight -m, which fails (1). Thus the polynomial/rational boundary
is sharp for every m, including m=1.

For m>=2 there is also the first-row-independent hostile

    F=w=u^m t,       G=u^(1-m)/(m-1).

Its Jacobian is one, although f(w)=w is nonconstant. Its constant
leading Laurent coefficient again fails (1). No rational analogue
for m=1 is inferred from this second example.

These are exact original-source computations. Neither a pole at u=0
nor the birational chart is silently treated as a polynomial map.

## 6. Useful square-prefix consequences at m=2

Consider explicit polynomials

    H=N(u)t^2+P(u)t+Q(u), L=M(u)t+R(u), F=H^2+L,

with ord_u N>=4 and ord_u P,ord_u M>=2. Then Fhat belongs to K[u,w]
for w=u^2t. Write N_j=[u^j]N, and similarly for the other rows.

If N_4!=0, the constant row of Fhat has w^4 coefficient N_4^2.
It is nonconstant, so there is no polynomial mate. This is an
active-fourfold source consequence under the displayed value/jet
hypotheses. It does not assert that every geometric fourfold boundary
point automatically satisfies those hypotheses.

When ord N>=5, constancy of the initial row

    (P2*w+Q0)^2+M2*w+R0

forces P2=M2=0. If ord N=5, the next row is

    2Q0*N5*w^2+(2Q0*P3+M3)w+2Q0*Q1+R1.

Since N5!=0, the theorem forces Q0=M3=0 and R1!=0.
If ord N>=6, its complete necessary equations instead are

    P2=M2=0,
    M3=-2P3*Q0,
    beta=2Q0*Q1+R1!=0.                             (6)

These statements are about the declared source polynomials. Any
geometric entry argument that gives their value/jet hypotheses is
separate; no unproved geometric entry is promoted here.

For a concrete all-p DG application, set u=x-p and retain the
explicit global rows

    H=u^6t^2+[(A+2)u^4+b*u^3+c*u^2]t
            +(A+1)u^2+(b-2pA)u+d,
    L=(lambda*u^4+mu*u^3+nu*u^2)t
            +lambda*u^2+(mu-2p*lambda)u+e.          (7)

The numerator-degree conditions for the actual fixed surface
z=x^2 hold for every parameter p,A,b,c,d,lambda,mu,nu,e.
The source verifies all three degree-four H numerator rows and
both degree-two L numerator rows directly with x=u+p. Thus this
is not an assumption that translating u extends to that surface.

The first two exact rows of Fhat are

    f=(c*w+d)^2+nu*w+e,
    g=2(c*w+d)[b*w+b-2pA]+mu*w+mu-2p*lambda.

The polynomial-mate gate therefore gives

    c=nu=0,       mu=-2bd,
    beta=-2p*(lambda+2dA)!=0.                       (8)

In particular p=0 is excluded for polynomial mates of this displayed
family. The theorem does not infer a rational exclusion from (8).
It does not depend on a pending six-plus-two geometric classification.
In the usual notation A=a-2, (8) is the direct coefficient consequence
reported to that independent lane.

## 7. Exact reproduction and finite-control universe

The proof is algebraic for all m and all degrees. The standalone source
imports no inherited mathematical implementation. It checks the
general monomial coefficient/weight identities, the complete
source-image criterion, actual full bracket transformations with
higher rows retained, leading-degree noncancellation, the m=1 case,
both rational hostiles, polynomial positives, and the complete
square-prefix/declared-DG row identities.

Its bounded monomial controls use m=1..6, Laurent exponents -12..5
and w-exponents 0..6. Its degree controls use h,j=1..6 and negative
L=-12..-1. These bounds are explicitly controls of the proved
identities, not a census that establishes the all-m conclusion.
Every gate is an explicit exception check, active under optimization.

    python3 04-computation/planar_jc48_sep08_polynomial_weight_gate.py
    python3 -O 04-computation/planar_jc48_sep08_polynomial_weight_gate.py

Both normal and optimized producer replays are byte-identical to the
frozen output: **1,947 gates, 553 output bytes**.

* Source: 5,982 bytes, SHA256
  885b2c01b2498c4dc991b91796ab1ee9365fc190cbd5699bd2ab76676c012509.
* Output SHA256:
  a0a99f9fb55e5e5bfe1075e60d1f800e1a86ef02bb5cfc86aeb052e7f766a6e7.
* Semantic gate trace:
  7c6566cea7919b18480be705f73791af41e44f876d1d78d66837ab208b277b79.

The source/output are frozen. The proof and applications are PROVED
following the independent audit and parent-owned promotion below.


## Accepted independent audit

The [complete independent audit](planar_jc48_sep08_polynomial_weight_gate_audit.md)
accepts the analytic proof, full coefficient scope and actual source
coordinates. Both normal and optimized replays reproduce all1947
gates and the frozen output. The candidate is now in the proved graph.
