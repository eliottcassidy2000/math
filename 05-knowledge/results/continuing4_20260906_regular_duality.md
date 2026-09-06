# Reflection turns a carried boundary into regular trinomial families

Status: **PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The paired
factor law below is analytic for every height. Positivity is certified only
at heights1 through6; general doubled-return noncancellation remains OPEN.
No global novelty or solution of the general Laurent problem is claimed.

## Inheritance and connection

The closest proved mechanism is the all-height singular-block divisor in
`third_20260906_laurent.md` and its independent audit. The real-root supplier
is THM-4436, `01-canon/theorems/THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall.md`,
with its explicit dependence on a cited finite-preserver theorem. The
canonical hostile is a doubled row with a lost inverse carry: the literal
support(-15,1,9) in `trinomial_doubled_row_literature_empty_core_three_ray_sep06.md`.
The corrected near miss is the non-positive-semidefinite coefficient of an
otherwise positive Bezout matrix in the incoming factor proof. The sidecar
recovered here is the complementary factor at a negative parameter, which
had supplied a valuation bound but was not used as an actual Laurent family.

Source: the original carried rows at x=-r. Target: regular rows with reflected
exponents. Map: specialize coefficients, remove their zero-root factors, and
reflect u to1/u. Preserved: the literal complementary polynomial and response
at every nonzero complementary root. Lost: the singular block and, without
a separate gcd condition, the first support-return mass. Needed sidecars:
the canceled inverse carry, monomial powers, normalization, and primitive
charge congruence. Cheapest tests: exact coefficient identities and a gcd
hostile, both included in the source. This map generates a new operation to
study, rather than inferring a cross-row theorem from individual real roots.

[Independent proof and exact audit](continuing4_20260906_regular_duality_audit.md) passes.

## 1. An actual family with a complete regular doubled row

Let H>=1 and z>=1 be integers with gcd(z,3H)=1. Put g=3H+z and

    f(u)=alpha*u^(-3H)+beta*u^z+gamma*u^(6H+3z),
    alpha*beta*gamma != 0,  tau=alpha^2*gamma/beta^3,
    X=alpha^z*beta^(3H).

Coefficients may be arbitrary complex numbers. The charge equation at mass m
is -3H*m+g*n_beta+3g*n_gamma=0, so g divides m. The complete count triples at
g and2g are respectively

    (z+2j,3H-3j,j), j=0,...,H,
    (2z+2j,6H-3j,j), j=0,...,2H.

Both contain their j=0 term. Define monic polynomials, with rising products,

    Phi(t)=sum_j H! (z+2j+1)^[2H-2j]/[j!(3H-3j)!] *t^j,
    Psi(t)=sum_j (2H)! (2z+2j+1)^[4H-2j]/[j!(6H-3j)!] *t^j.

The sums have the respective ranges above; a zero-length product is1.
Their coefficients are in Q[z]. Literal multinomial coefficients give

    CT(f^g)=X*binom(g,H)*Phi(tau),
    CT(f^(2g))=X^2*binom(2g,2H)*Psi(tau).

THM-4436 applies with A=1,B=3,x=z,r=0,z_base=0, and gives H distinct
strictly negative roots of Phi. In particular it is not necessary to assert
the inapplicable A=2 parameterization with an unbounded base remainder.
There is no inverse carry in this family. Dropping gcd(z,3H)=1 is invalid:
at z=3H the support is(-3H,3H,15H), and the first support return is2,
whereas the proposed g is6H. With all coefficients1 its CT at mass2 is2.

## 2. Universal paired divisor for the response characteristic polynomial

Let M be multiplication by Psi modulo the monic Phi on Q[z][t]/(Phi),
in the basis1,t,...,t^(H-1), and write

    det(wI-M)=w^H+sum_(k=1)^H c_(H,k)(z)*w^(H-k).

For every H>=1 and1<=k<=H,

    D_(H,k)(z)=product_(ell=1)^H
                 [(z+2ell-1)(z+2ell)]^max(0,k-H+ell)

divides c_(H,k) in Q[z]. Its degree is k(k+1), and

    deg(c_(H,k)/D_(H,k)) <= 4Hk-k(k+1).

Proof of divisibility. Fix m in{1,...,2H}, write ell=ceil(m/2), and set
delta=z+m. At delta=0 the first ell coefficients of Phi have a simple
zero and the coefficient of t^ell is nonzero. Thus Phi=t^ell*a(t) with
a(0)!=0 at this parameter, and its constant coefficient has a simple
delta zero. The coefficients of Psi below degree m have a simple delta
zero; those at degree at least m are regular. These statements follow
directly from the consecutive factors in the two rising products.

The two relatively prime factors t^ell and a lift formally over C[[delta]].
After delta=epsilon^ell, the small roots have the form
t_i=epsilon*v_i(epsilon), where v_i(0) are the distinct nonzero roots of
a(0)*v^ell+partial_delta Phi(0,0). Formal lifting is legitimate since the
derivative at each such v is nonzero. At a small root, every term of Psi
below degree m has order at least epsilon^ell, and every higher term has
order at least epsilon^m. Since m>=ell, all ell small response eigenvalues
have order at least delta. Their jth elementary symmetric function has
delta order at least j. These symmetric functions lie in C[[delta]], by
the lifted degree-ell factor and its multiplication operator.

The Chinese remainder decomposition leaves a regular complementary
operator of size H-ell. A term contributing to c_(H,k) must take at least
max(0,k-H+ell) small eigenvalues. Therefore (z+m) to this power divides
c_(H,k). Grouping m=2ell-1 and m=2ell proves the stated divisor. This is
a lower valuation bound; possible higher-order cancellation does not harm it.

Proof of degree. The coefficient of t^j in Phi has z-degree2(H-j), and
that in Psi has degree4H-2j. Monic reduction preserves the weighted bound
deg_z coefficient of t^i <= weight-2i. Hence M_(i,j) has degree at most
4H+2j-2i. In a principal minor the row and column indices cancel, giving
deg c_(H,k)<=4Hk. Subtract the divisor degree. This also bounds every
coefficient used in the independent interpolation audit.

## 3. Exact reflection of the carried boundary

For the inherited original height h and1<=r<=h put H=h-r,z=2r+1.
Writing p_h,x and q_h,x exactly as in `third_20260906_laurent.md`,

    p_h,-r(t)=t^r*Phi_(H,2r+1)(t),
    q_h,-r(t)=t^(2r)*Psi_(H,2r+1)(t)/(4h+2)!.

These are coefficient identities, including H=0 with Phi=Psi=1. In the
first row the falling product vanishes below j=r; above it, set j=r+k
and use (h-r)_(h-j)=H!/k!. In the doubled row it vanishes below e=2r;
above it, (2h-2r)_(2h-e)=(2H)!/(e-2r)!. The displayed identities follow
from the factorial forms of Phi and Psi.

The old support at x=-r is(-6h-3,-(2r+1),3(h-r)). Reflecting the exponents
and reordering gives(-3H,z,6H+3z). For H>=1 and the separate primitive
condition this is the actual family of Section1. The factor identities
themselves require no gcd restriction.

The normalization also survives exactly. If the reflected coefficients are
(alpha_new,beta_new,gamma_new)=(gamma_old,beta_old,alpha_old), the phase tau
is unchanged. The old formal mass is g=3h-r+1=3H+z, and
X_old*tau^r=X_new. Moreover binom(g,2h+1)=binom(g,H),
(2g)_(4h+2)/(4h+2)!=binom(2g,2H), and
gcd(g,6h+3)=gcd(z,3H). Thus the complementary polynomial map is compatible
with the entire literal constant-term normalization and primitive condition.

At a nonzero complementary root the inherited polynomial response agrees
with this specialized q. At the small zero-root block it need not: the
inverse-carry coefficient and p_0 both vanish, so their ratio must first
be canceled in the generic quotient. Naive specialization at t=0 discards
the very term responsible for the old carried divisor. No singular-block
equivalence is asserted here.

## 4. Certified heights and actual first detection

The companion exact JSON retains all residual coefficients, in ascending
powers of z, for H=1,...,6 and k=1,...,H. Every one is strictly positive.

    H       residual degrees                         coefficient count
    1       2                                        3
    2       6,10                                     18
    3       10,18,24                                 55
    4       14,26,36,44                              124
    5       18,34,48,60,70                           235
    6       22,42,60,76,90,102                       398

There are833 coefficients in21 residual polynomials. For z>=1 every paired
factor is positive, so all nonleading characteristic coefficients are
positive. The polynomial det(wI-M) has no real root w>=0. At every root
rho of Phi, Psi(rho) is a real eigenvalue of M, hence is strictly negative.
This proves noncancellation of the first two moments for every actual
parameter and all nonzero complex coefficients in these six families.

Consequently the first nonzero constant term is exactly g or2g. Each is
attained at every admissible z: a generic tau avoids the H first roots;
choosing any one of those H distinct negative roots cancels the first
moment and leaves the doubled one nonzero. Every tau is attainable by
setting alpha=beta=1 and gamma=tau. Both masses are strictly below the
endpoint width3g. H=1,2 have negative endpoints3,6 and serve as inherited
small-endpoint controls. H=3,4,5,6 give the explicit negative-endpoint
families9,12,15,18. A bounded repository search recovered no previous
uniform certificate for this exact family; that is a retrieval statement,
not a claim of literature novelty.

The producer performs monic polynomial division and Faddeev--LeVerrier over
Q[z], verifies the full Cayley--Hamilton remainder, divides the proved
factors exactly, and checks each of the833 rational coefficients. Literal
complete fibres and multinomial normalizations, all boundary identities
through h=6, and the nonprimitive earlier-return hostile are additional
controls. Run `python 04-computation/continuing4_20260906_regular_duality.py`
with an optional output JSON path, or with `python -O`; its gates are active
in both modes. The independent audit must reconstruct literal integer
multinomial rows, rather than trusting the producer's matrix routine.

## 5. Stopping boundary and next test

The all-height paired factor law is proved analytically; the all-height
positivity of the residual coefficients is OPEN. Its failure would not
by itself refute the weaker response-sign claim. The new family arose by
retaining a discarded complementary object, not by extrapolating the old
height table. A next structural test should compare the regular paired
factors with the carried singular block while retaining the response
operator. An orthogonality label or an unmarked quotient loses the required
cross-row sign. Neither LRC(14) nor general trinomial endpoint detection is
closed by this result.
