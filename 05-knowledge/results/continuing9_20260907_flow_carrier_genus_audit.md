# Independent audit of invariant genus and rational flow exclusion

**Status: AUDIT PASSED / PROVED with cited curve inputs; FINITE-EXACT controls.**
This audit covers the complete affine-multiplier classification, every
univariate multiplier, primitive monomial invariants, and polynomial
compositions in the same invariant in
[the primary proof](continuing9_20260907_flow_carrier_genus.md).
It does not exclude general multipliers or compositions of different flows.

The auditor independently read the inherited completed-carrier and
nonrationality proofs, the current curve-input sidecar, and the primary
Stacks pages. The source imports no research producer. In particular it
keeps an arbitrary nonzero coefficient of y in its discriminant calculation
and reconstructs the rational inverse of the actual source coordinates.

## 1. Generic curve, discriminant, and boundary audit

For `I=p²(p³-y²)(alpha+beta p+gamma y)`, with gamma nonzero,
the cubic equation after normalizing gamma has three roots of p-valuation
`-2/3` at zero. A rational function of p has integral valuation, so the
cubic has no rational root even over an algebraic closure of K(c).
This is an irreducibility proof, not just an assertion about the displayed
affine equation. The generic curve is geometrically integral.

The independent unnormalized computation of `disc_y(I-c)/p⁴` gives
degree 13, leading coefficient `4gamma⁴`, and value
`-27gamma²c²` at p=0. The independent gradient resultant has degree nine
and leading coefficient `169gamma⁵`. Thus normalization has not hidden a
parameter specialization where the critical-locus argument disappears.

The primary proof correctly separates triple roots from double roots.
At a nontriple double root the first variation of the cubic discriminant
is a nonzero multiple of the variation of its value at that root. A
repeated finite discriminant zero would therefore be a critical point
of I with value c. The exact resultant forces its coordinates to be
algebraic over K, incompatible with transcendental c. The separate triple
root resultant is a nonzero cubic in p. All 13 nonzero finite branch
points are consequently simple, uniformly in alpha and beta.

The audit reconstructs the local base changes directly from I-c:
at zero, `p=r³,y=r^-2Y`; at infinity, `p=r^-2,y=r^-3Y`.
The zero chart gives a single index-three point. The infinity chart gives
one index-two point from the branches near Y=+1,-1 and one unramified
point from Y=0. The deck involution fixes the original p,y, so the three
upstairs branches are not incorrectly counted as three downstairs points.
Total ramification is `13+2+1=16`; the degree-three cover has genus six.

## 2. General univariate multiplier and monomial checks

For arbitrary nonzero `R(p)=A(p)²B(p)`, take A monic and B squarefree,
leaving any nonsquare scalar in B. The exact birational coordinate
`v=pAB y` gives

`v²=B(p)(p⁵R(p)-c)`.

The second factor has no generic repeated root: such a root must also
annihilate the nonzero derivative of p⁵R, and hence be algebraic over K,
forcing c algebraic. The factors are coprime because at a root of B the
second factor is -c. The branch polynomial is therefore squarefree.
With m=deg R and b=deg B, m-b is even, so its degree m+b+5 is odd and
its genus is `(m+b+4)/2 >=2`. This proves geometric integrality as well.
The repeated-factor controls in the audit independently check this
square-removal step; no rational square root of a scalar is assumed.

For primitive `J=p^aD^b`, the Bezout torus map has a rational inverse
over K(c) itself. The audit reconstructs it for 63 coprime exponent pairs
and checks the double-cover branch count independently. The analytic
formula is

`g=(a+3b+(a mod2)+(b mod2)-2)/2 >=2`.

The distinction between the primitive invariant J and the Hamiltonian
f(J) is essential: the latter can have disconnected generic fibres.
The proof retains J and correctly tests the exact carrier entry
`a*ord_0(f-f(0)) >=2`, since b and the order are positive.

## 3. Why the genus really rules out the proposed flow

The actual function fields agree: from `p=t(1+x²t), y=xtp` one recovers
`s=y/p=xt`, `tau=p-s²=t`, and `x=s/tau`. The source verification checks
these identities directly. The fixed-input comparison of the two formal
models is inherited with its injectivity and convergence hypotheses;
it is not an identification of arbitrary completed elements with
rational Laurent expressions.

For a univariate or affine multiplier, the invariant has leading form
`I=tau W(s)+O(tau²)` with nonzero W. If the first nonconstant outer term
is `f_k I^k`, then

`exp(lambda delta)(p)-p = 2lambda k f_k s W(s)^k tau^k + O(tau^(k+1))`.

For a monomial invariant replace k by bk and use its displayed W.
Every positive iterate has a nonzero leading displacement because the
field has characteristic zero and lambda is nonzero. Thus the proposed
rational map could not have finite order.

A rational specialization would give a nonconstant selfmap of the
smooth proper curve over the retained invariant field. The function-field
correspondence has the correct contravariant direction.
[Stacks, Theorem 53.2.6](https://stacks.math.columbia.edu/tag/0BY1)
For genus greater than one, Riemann--Hurwitz forces this selfmap to have
degree one. [Stacks, Section 53.12](https://stacks.math.columbia.edu/tag/0C1B)
Its automorphism group is finite: the finite-type group scheme has no
tangent sections since `deg T_C=2-2g<0`.
[Stacks, Section 109.7](https://stacks.math.columbia.edu/tag/0DST)
This contradicts the explicit infinite-order displacement.

The canonical rational non-LND example `H=x²t`, constant Hamiltonians,
zero time, and composition with inverse time all remain outside the
claimed exclusion or are explicitly exempted. This result does not
resolve the planar Jacobian conjecture.

## 4. Reproduction and explicit universe

The independent engine checks the arbitrary-coefficient cubic identities,
both boundary charts, nine univariate multipliers including repeated and
square factors, 63 primitive exponent pairs with `1<=a<=12,1<=b<=8`,
the actual source inverse, and outer leading orders 1 through 4.
These controls corroborate the analytic statements; they are not a
finite substitute for the all-parameter proofs.

    python3 04-computation/continuing9_20260907_flow_carrier_genus_audit.py
    python3 -O 04-computation/continuing9_20260907_flow_carrier_genus_audit.py

Both modes pass 313 always-active exact gates with identical raw LF output.
The session manifest pins the final source, transcript, and proof bytes.
