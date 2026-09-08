# An all-degree obstruction for a source-linear logarithmic witness on W2

**Status: PROVED ANALYTICALLY + FINITE-EXACT; INDEPENDENTLY AUDITED.**
The theorem bounds the original source degree of the repair witness, not the
degree of the first function. It does not exclude unit order one when every
polynomial repair witness has source t-degree at least two.

## 1. Statement and inherited mechanism

Use the fixed surface charts

    x=1/r, t=-r^2-r^4 b, dx wedge dt=r^2 dr wedge db,
    J(F,H)=F_x H_t-F_t H_x.

**Theorem.** There is no polynomial F in C[x,t] which extends to an
everywhere submersive regular function on W2 and satisfies

    J(F,H)=F,       H=A(x)t+B(x) in C[x,t].                 (1)

There is no bound on deg_t F, deg A, or deg B. In fact all source-submersive
pairs (1) are classified below, and all their globally regular first
functions have the entire added divisor r=0 critical.

For the original response module C[x,t]/D_g C[x,t], where g=F-c and
D_g=J(g,-), scalar unit order one requires [1]!=0 and g[1]=0. The second
condition is exactly J(g,H)=g for a polynomial H. Thus the theorem proves
that any such repair witness for a global W2 submersion must have original
t-degree at least two. The converse is not claimed. No change of response
ring to O(W2), nor preservation of that global ring by D_g, is assumed.

The closest proved mechanisms are the original-source component-jet
theorem in [Vertical component jets](planar_jc48_sep06_torsion.md) and the
fixed W2 boundary checks in [the univariate supplier](continuing13_20260908_univariate_fixed_unit.md).
The recovered elementary viewpoint is the Hamiltonian eigenfunction
equation from [THM-1345 / jc2-equivariant-category-poisson-reframing-dc1-shadow](../../01-canon/theorems/THM-1345-jc2-equivariant-category-poisson-reframing-dc1-shadow.md).
Only its elementary eigenfunction viewpoint is used; no semisimple action
or equivariant Jacobian theorem is needed here. Our sign convention is
explicit in (1).

The canonical hostile is F=t(1+xt): its original-source unit has exact
order one and it is globally regular, but its added divisor is critical.
The corrected near miss is mistaking an affine logarithmic eigenfunction
for a global submersion. The least-used operation is taking the rational
logarithmic derivative, then comparing its residues by Newton identities.
The live concepts are witness degree, signed root valuations, eigenfunction
separation, boundary normal order, and source versus global response.

The one-arm construction [THM-3975 / danielewski-one-arm-modification-cubic-control-and-hyperelliptic-no-mate](../../01-canon/theorems/THM-3975-danielewski-one-arm-modification-cubic-control-and-hyperelliptic-no-mate.md)
uses a different completion and its no-mate family starts at a different
exponent. It is not a prohibition on the present n=1 hostile. The
[prime-leading filter](planar_jc48_sep08_prime_leading.md) requires a
polynomial mate and cannot be applied to arbitrary rational mates here.

## 2. Rational separation and the residue obstruction

Suppose F is everywhere submersive on the original source C^2 and (1)
holds. Then F is nonconstant and nonzero. If A=0, the left side is
-B' F_t, whose t-degree is smaller than that of a nonzero F, or is zero.
This is impossible.

For A!=0 put v=H. The rings C(x)[t] and C(x)[v] agree. At fixed v, (1)
becomes A partial_x F=F. Every nonzero coefficient of F as a polynomial
in v solves the same first-order equation, and the ratio of any two
solutions has derivative zero. Consequently

    F=Y(x) P(H),   P in C[v] nonzero,   Y in C(x)^*,
    Y'/Y=1/A.                                           (2)

A rational logarithmic derivative has only simple finite poles, with
integer nonzero residues at its poles, and has no polynomial part.
Therefore A cannot be a nonzero constant, A is squarefree, and at every
root alpha of A

    m_alpha=ord_alpha Y=1/A'(alpha) in Z minus {0}.

There are no other zeros or poles of Y. Its nonzero scalar is absorbed
into P.

First P must be squarefree (a nonzero constant is allowed). Indeed, a
multiple root gamma would give a repeated nonvertical irreducible factor
of H-gamma in F. Such a factor exists because H-gamma is linear with
nonzero t coefficient over C(x); vertical rational factors Y cannot
cancel it. It has actual source points with A!=0. Both derivatives of F
vanish there, a contradiction.

At a root alpha of A, the order of H-B(alpha) at the vertical divisor is
one: its first coefficient is A'(alpha)t+B'(alpha), not the zero
polynomial. Hence ord_alpha P(H) is r=0 or 1, depending on whether
P(B(alpha)) is nonzero or zero. If m_alpha>0, the vertical factor in F
has exponent m_alpha+r; submersion forces m_alpha=1 and r=0. If
m_alpha<0, polynomiality requires r>=-m_alpha, so m_alpha=-1 and r=1.
We have proved

    m_alpha=+1 => P(B(alpha))!=0;
    m_alpha=-1 => P(B(alpha))=0.                        (3)

These signs force deg A<=2. If N=deg A>=2, expansion of (2) at infinity
gives

    sum_alpha m_alpha alpha^j=0,  0<=j<=N-2.             (4)

The j=0 equation makes N=2n, with n positive and n negative roots. If
n>=2, then N-2>=n. The power sums of the two n-element root sets agree
through degree n. Newton identities make their monic degree-n
polynomials identical, contradicting that the root sets are disjoint.
Thus N=2. This is an all-degree argument, not a finite search.

The use of source submersion is necessary. For example

    A=x(x^2-1)/2, H=A t, Y=(x^2-1)/x^2,
    P(v)=v^2, F=(x^2-1)^3 t^2/4

satisfy (1)-(2) with deg A=3 and residues (+1,-2,+1). The repeated factor
t^2 makes F critical. Rational logarithmic exactness alone does not
force deg A<=2.

## 3. Complete source-submersive normal forms

In each row P is squarefree and nonzero; all nonzero scalar factors of F
are incorporated in P. The displayed conditions are necessary and
sufficient for source submersion and equation (1).

**One positive root.** For arbitrary alpha and B,

    A=x-alpha, H=(x-alpha)t+B(x),
    F=(x-alpha)P(H), P(B(alpha))!=0.                    (+)

P may be constant. The vertical zero component has nonzero normal
derivative P(B(alpha)); all other zero components are simple levels of
H away from A=0.

**One negative root.** For arbitrary beta, b0, C and R,

    L=-t+C(x), H=b0+(x-beta)L,
    P(v)=(v-b0)R(v), R(b0)!=0, P squarefree,
    F=L R(H), A=-(x-beta), B=b0+(x-beta)C(x).            (-)

The component L=0 has derivative L_t=-1 and R(b0)!=0. Other zero
components are simple H-levels away from x=beta, and are disjoint from
L=0.

**One root of each sign.** For arbitrary distinct alpha,beta and B,

    A=(x-alpha)(x-beta)/(alpha-beta), H=A t+B(x),
    F=(x-alpha)/(x-beta) P(H),
    P(B(beta))=0, P(B(alpha))!=0.                      (+-)

This expression is polynomial. More explicitly, with

    b0=B(beta), P(v)=(v-b0)R(v),
    C=(B-b0)/(x-beta),
    L=(x-alpha)t/(alpha-beta)+C,

one has F=(x-alpha)L R(H). The vertical component has nonzero normal
derivative P(B(alpha))/(alpha-beta). The level L=0 avoids x=alpha
because C(alpha)!=0, and its t derivative is nonzero there. Remaining
simple H-levels avoid both roots of A and all other zero components.

For all three rows there can be no critical point with F!=0, directly
from (1). The checks on F=0 therefore prove sufficiency everywhere.
Conversely Sections 2 and (3) give exactly these three rows, including
all allowed constant-P and constant-R boundaries.

## 4. Complete W2 boundary obstruction

Take x to infinity with b fixed under the actual chart substitution.

In (+), if P is constant then F has a pole. If B is nonconstant and P
is nonconstant, H has the leading polynomial term B(x), and F has a
pole of positive order 1+deg P deg B. If B is constant beta0, then
H=beta0+O(1/x), and F~x P(beta0) has a pole by the source condition.
Thus no source-submersive (+) member is globally regular.

In (-), put M=deg R>=0. If C has positive degree k, then
F=L R(H) has a nonzero leading term of degree k+M(k+1), so has a pole.
If C=kappa!=0 is constant and M>=1, the pole degree is M. The remaining
members, and only these, are globally regular:

    C=0, any R;   or   C constant, R constant.           (5)

In the first case

    H=b0-(x-beta)t,
    F_infinity=r^2 R(b0)+O(r^3),  R(b0)!=0.

Both H and -t are polynomial in r,b. In the second case F is a nonzero
scalar multiple of kappa-t, and its boundary expression is that scalar
times (kappa+r^2+r^4b). In both cases both first derivatives vanish
along the entire added divisor r=0.

In (+-), P has positive degree and the two values B(alpha),B(beta)
must differ. Hence B is nonconstant. But Y=(x-alpha)/(x-beta) tends
to 1, and A t tends to -1/(alpha-beta). Thus F=Y P(H) has the nonzero
leading polynomial term of P(B(x)), a pole. Equivalently, forcing B
constant would contradict the two conditions in (+-).

This proves the theorem and identifies the strongest survivor: every
globally regular source-submersive pair with source-linear witness is
in the negative-root row (5), and its whole added divisor is critical.

## 5. A genuine order-one hostile and a unit-zero boundary

Set

    F=t(1+xt), H=-xt, G=H/F=-x/(1+xt).

Then J(F,H)=F and J(F,G)=1. The two reduced zero components t=0 and
1+xt=0 are disjoint; F is a source submersion. The function F is global
on W2, but Section 4 shows that its added divisor is critical.

Its original-source scalar unit has exact order one, rather than zero.
Indeed C(x,t)=C(F)(H), with t=F/(1-H) and x=-H(1-H)/F. Since
D_F H=F, the rational constant field is C(F), and every rational mate
is G+R(F). On t=0 the mate G is regular and H restricts to 0. On
1+xt=0 it has simple scalar part +1/F. A rational scalar R(F) cannot
cancel the latter without introducing a pole on the former. Thus
there is no polynomial mate; F[1]=0 and [1]!=0.

By contrast F=-t+C(x), H=-(x-beta)t+(x-beta)C(x)+b0 are in the same
negative-root classification with R=1, but F has the polynomial mate x
and its source unit is zero. This prevents identifying every
logarithmic identity with a nonzero order-one unit.

## 6. Scope, reproduction and next obligation

The proof is over C and covers all coefficients and degrees of the
first function and all source-linear witnesses. Exact symbolic controls
check the bracket separation and normal forms, positive and negative
root conditions, complete Laurent globality in finite declared parameter
grids, the actual added critical divisor, the nonzero-unit hostile,
the coordinate zero-unit boundary, and the cubic residue counterexample.
They support the analytic proof; the finite grids do not supply its
quantifiers.

The frozen engine passes 467 always-active gates, including six separate
Groebner ideals of both source partial derivatives. Normal and optimized
Python runs have byte-identical LF output and identical certificates.

Reproduce with `python continuing14_20260908_logarithmic_linear_witness.py`
and `python -O continuing14_20260908_logarithmic_linear_witness.py`.
The engine locates its certificate output beside the source, or in the
repository results directory after relocation from 04-computation.

The next unresolved operation is a polynomial witness H with genuine
source t-degree at least two. The residue separation used here relies
on H being linear over C(x), and is not an obstruction for that case.
No general cubic or quintic classification, no general unit-order-one
impossibility, and no Jacobian-conjecture conclusion follows.

Independent acceptance: [referee and reproducible controls](continuing14_20260908_logarithmic_linear_audit.md). Only this audited status and routing paragraph were added after the referee-pinned producer report; the proof and computational bytes are unchanged.
