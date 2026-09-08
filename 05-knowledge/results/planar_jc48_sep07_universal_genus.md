# A fixed exceptional curve excludes rational time in the entire universal carrier

**Status: PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [independent analytic and source audit](planar_jc48_sep07_universal_genus_audit.md)
accepts the complete all-parameter theorem and both frozen replays.
JC(2) remains OPEN. The new argument concerns the specified universal
Hamiltonian carrier, not arbitrary Keller endomorphisms or compositions
of flows with different invariants. No external priority claim is made.

## The precise stronger result

Write Delta=p3-y2. For every nonzero polynomial I in C[p,y] divisible by
p Delta, **every geometric component of the generic fibre of I has genus
at least two**. Geometric integrality is not assumed and need not hold.
The statement includes arbitrary repeated factors, arbitrary coefficients,
all degrees, and specializations of those coefficients.

Consequently, for every nonconstant universal carrier

    S in C+p2 Delta C[p,y],

its actual completed Hamiltonian flow at every nonzero scalar time fails
to have both p and y images rational in C(p,y)=C(x,t). In particular it
cannot be a polynomial source automorphism. This settles the rational-time
question throughout this carrier, including every polynomial lift of the
finite source response. It leaves later finite compatibility and changed
carriers as separate problems. It also leaves compositions retaining
different invariants outside the argument.

## Inheritance, hostile, and changed representation

The closest proved mechanism is the generic-curve/finite-automorphism
argument of [the incoming nonrational-time theorem](planar_jc_long_20260906_nonrational.md).
The [variable y-linear genus theorem](planar_jc48_sep07_carrier_frontier.md)
and the [actual supplier genus27 theorem](planar_jc48_sep07_supplier_genus.md)
retain exact genus by calculating whole discriminants. The
[completed source-image theorem](planar_jc48_sep07_completed_response.md)
leaves other carrier lifts, including their six terminal directions, free.
The present result handles all these remaining polynomial carriers at once.
The exact local-cusp question was already explicitly proposed, but left
**OPEN / unused**, in the [incoming carrier-genus stopping boundary](continuing9_20260907_flow_carrier_genus.md), Section4: multiplicities and arbitrary extra factors needed a separate
local/global argument. The contribution here pays that stated obligation;
it does not claim to originate the question.

The canonical hostile to using non-local-nilpotence alone is the rational
flow of x2t, outside the carrier. The corrected near miss is to assume
I is primitive: I=(p Delta)2 is a valid universal carrier with a geometrically
reducible generic fibre. We retain its relative constant field below.
The least-used sidecar is just one exceptional curve over the fixed
intersection p=Delta=0. Its nearby covering surface already has genus two;
additional global ramification is unnecessary for the lower bound.

The five live concepts are weighted initial forms, actual embedded
bordered surfaces, connected components of cyclic covers, generic constant
fields, and scalar iteration. The map is the third ordinary blowup
(p,y)=(v2w,v3w). It retains an actual open subset of each nearby source
fibre. Passing to its exceptional leading equation loses the rest of the
fibre; the analytic coordinate change below pays that loss exactly.
The cheapest hostiles remove p or Delta, or take an outer power.

## 1. The weighted initial form and the actual blowup chart

Take the least weight of I for wt(p)=2,wt(y)=3. Since I is divisible by
p Delta, its nonzero weighted initial form factors over C as

    I_min=C p^A y^B product_(i=1)^r (p3-lambda_i y2)^e_i,   (1)

where C!=0, A>=1, B>=0, r>=1, e_i>=1, and the lambda_i are distinct
nonzero constants. The factor lambda=1 occurs because Delta divides I.
Factors with lambda=0 are absorbed into p^A, and a y factor is absorbed
into y^B. The factorization follows directly by solving 2a+3b=n among
the monomials of a weighted homogeneous polynomial: after removing
common p and y powers, its terms are a homogeneous binary polynomial
in p3,y2. Repeated binomial factors are explicitly allowed.

Put E=sum e_i and M=2A+3B+6E. Three ordinary blowups at the origin and
the successive cusp intersection give the chart

    p=v2w, y=v3w,
    I=v^M(Phi(w)+v Psi(v,w)),
    Phi(w)=C w^n0 product_i(w-lambda_i)^e_i,
    n0=A+B+2E,       ninf=A+B+3E.                       (2)

Here Psi is polynomial in this chart: every higher weighted monomial
has v-exponent at least M+1. The new exceptional component is v=0;
w=0 and w=infinity are its intersections with the two older exceptional
components. The strict Delta branch meets at w=1. Formula (2) is obtained
by literal substitution, independent of a completed resolution of I.

Let K be the compact sphere with disjoint small open disks removed around
0,infinity and all lambda_i. It lies in the displayed finite w chart,
and Phi is nonzero on a neighborhood of K. For sufficiently small v,
|v Psi/Phi|<1 uniformly there. The single-valued binomial root gives the
analytic coordinate

    vtilde=v(1+v Psi/Phi)^(1/M).

It has derivative one at v=0 and is fibrewise invertible on one uniform
neighborhood over K. Hence for every sufficiently small c!=0, the actual
piece of I=c above K is isomorphic to

    vtilde^M Phi(w)=c,      w in K.                    (3)

All M roots have uniformly small modulus. No numerical root proposal or
Milnor-fibre approximation is being substituted for this analytic
identification. Because v!=0 and w stays away from0, blowdown is an
isomorphism on this piece: it is a compact bordered surface actually
embedded in the original affine curve I=c.

## 2. Every connected piece has genus at least two

The monodromy increments for (3), on M cyclically labelled sheets, are
-n0 around zero, -e_i around lambda_i and ninf around infinity. Their
sum is zero. The number of connected components is

    d=gcd(M,n0,e_1,...,e_r)=gcd(A,B,e_1,...,e_r).         (4)

Indeed E is divisible by every common divisor of the e_i, and
M-2n0=B+2E then recovers B, while n0-B-2E recovers A. Conversely all
of M,n0,e_i are divisible by the right-hand gcd.

Capping the boundary circles of (3) gives the compact cyclic cover with
these branch exponents. Each of its d components has the same genus g.
For an exponent n the total number of points over that puncture is
gcd(M,n), and its total ramification contribution is M-gcd(M,n).
The characteristic-zero [Riemann--Hurwitz formula](https://stacks.math.columbia.edu/tag/0C1B)
therefore gives exactly

    2d(g-1)=rM-gcd(M,n0)-gcd(M,ninf)-sum_i gcd(M,e_i).    (5)

Dividing A,B and all e_i by d divides M,n0,ninf and the right side of
(5) by d without changing the component genus. We may thus assume d=1.
The normalized A is still at least one and E is still positive.

If B>0, use

    gcd(M,n0)<=B+2E,       since M-2n0=B+2E,
    gcd(M,ninf)<=B,        since M-2ninf=B,
    sum_i gcd(M,e_i)<=E.

Then the right side of (5) is at least

    rM-2B-3E >= 2A+B+3E >= 6.

If B=0, then ninf=M/2, while gcd(M,n0)<=2E. The right side is at least

    rM-2E-M/2-E=(2r-1)A+6(r-1)E >= 1.

It is even by (5), so it is at least two. In both cases g>=2.
This argument uses no generic-coefficient assumption and makes no claim
that adding factors preserves an exact genus. It establishes the lower
bound from the retained exceptional covering piece alone.

The bound is sharp: p Delta and p2 Delta give connected genus-two pieces.
The power (p Delta)2 has two components, each of genus two. Removing p
gives Delta, whose smooth proper generic curve has genus one; removing
Delta gives the rational generic fibre of p. Thus both mandatory factors
in the statement have a visible role.

## 3. From the embedded piece to all generic components

Choose a small nonzero c outside the finite exceptional set for a smooth
proper model of the polynomial pencil, allowing a finite base change
which separates the geometric components. Such a model is obtained by
resolving the rational map from P2 to P1 and normalizing its finite
Stein base; generic smoothness and properness give an open base of smooth
proper fibres. Only finitely many values are removed from a curve base.

Every connected surface in (3) embeds in one affine component of I=c,
and therefore in its compactification. Capping the boundary is used only
to measure the genus of the bordered surface; the capped surface itself
is not asserted to embed. An embedded genus-g bordered surface forces
ambient genus at least g: its g handle pairs retain their nondegenerate
intersection form in the ambient oriented surface. Hence at least one
component of a general I=c has genus at least two.

All geometric generic components have the same genus. Algebraically,
I-c is irreducible over C(c) before extending constants, since it is
irreducible and linear in the variable c in C[p,y,c]. The relative
algebraic closure E0 of C(I) in C(p,y) is a finite extension; the
geometric components are conjugate under the embeddings of E0 into an
algebraic closure. Their smooth proper genera agree. Alternatively the
finite Stein base is connected and the smooth proper component genus
is constant on its nonempty smooth open. The preceding embedded piece
therefore proves genus at least two on every geometric generic component.
This is deliberately a componentwise theorem, not an incorrect geometric
integrality assertion for every I.

## 4. All nonzero scalar times in the universal carrier

Take S=c0+I with 0!=I in p2 Delta C[p,y]. The completed scalar-time
automorphism and its group law are supplied by
[the actual completed-carrier theorem](planar_jc_long_20260906_hamiltonian.md).
Its fixed-input rationality comparison is supplied by the independently
audited [nonrationality note](planar_jc_long_20260906_nonrational.md).
We do not identify the whole completion with a rational Laurent chart.

Factor I=Delta^e T with e>=1 and Delta not dividing T. In logarithmic
coordinates p=s2+tau,y=sp, the actual bracket is

    {F,G}=tau(F_s G_tau-F_tau G_s).

Since Delta=tau(s2+tau)2,

    I=tau^e W(s)+O(tau^(e+1)),
    W(s)=s^(4e) T(s2,s3) != 0.                         (6)

The last nonvanishing follows from ker[C[p,y]->C[s],(p,y)->(s2,s3)]
being exactly (Delta). The derivation raises tau order by at least e.
For the actual source coordinate p=s2+tau,

    exp(lambda delta_S)(p)-p
      =2lambda e s W(s) tau^e+O(tau^(e+1)).             (7)

No later iterate can cancel that coefficient. The same nonzero
coefficient occurs with lambda replaced by n lambda for every positive
integer n in characteristic zero.

Suppose both images of p,y were rational at lambda!=0. The faithful
fixed-input comparison and formal invertibility make them algebraically
independent; substitution defines an injective endomorphism sigma of
L=C(p,y). It fixes C(I). Let E0 be the finite relative constant extension
from Section3. Then sigma(E0) is contained in E0 and is an automorphism
over C(I), because E0 is finite over that field. Some positive power
sigma^n fixes E0 pointwise. On the smooth proper genus-at-least-two curve
with function field L/E0, it is a nonconstant selfmap.

Riemann--Hurwitz forces that selfmap to have degree one, and its geometric
automorphism group is finite. Thus a further positive power of sigma is
the identity. The scalar-time group law would give identity at a positive
integer multiple of lambda, contradicting (7). This also handles every
outer-composite Hamiltonian without pretending its generic fibre is
geometrically integral or requiring an explicit polynomial primitive.

The source identities x=yp/Delta and t=Delta/p2 show that simultaneous
rationality of p,y and simultaneous rationality of x,t are equivalent.
The primary [curve/function-field correspondence](https://stacks.math.columbia.edu/tag/0BY1)
and [finite automorphism group](https://stacks.math.columbia.edu/tag/0DST)
are the same already-read inputs as in the earlier genus theorems.

The same two statements hold over every characteristic-zero coefficient
field k. For the genus assertion, descend the finitely many coefficients
of I to a field k0 finitely generated over Q and embed k0 into C.
Geometric components and the genera of their smooth proper models are
preserved by algebraically closed field extension, so the complex result
gives the result over k. For the rational-time assertion, also include
lambda and the coefficients of both proposed rational images in k0.
After clearing their denominators, the completed-flow equalities are
coefficientwise identities over k0: their coefficients are obtained from
the polynomial derivation by field operations and division by factorials.
The embedding into C preserves these identities, the nonzero scalar time,
and nonzero denominators, contradicting the complex result. No embedding
of the entire possibly larger field k into C is required.

The proof therefore applies after extending a coefficient field by an external
boundary parameter sigma and specializing it: j(sigma) is a scalar time,
not a function of the source's logarithmic variable s. Time zero and
constant Hamiltonians are exempt. Nothing here prohibits inverse-time
cancellation or a rational composition retaining no common invariant.

## 5. Exact controls and completed audit

The [standalone verifier](../../04-computation/planar_jc48_sep07_universal_genus.py)
checks13104 complete multiplicity tuples: A=1..12,B=0..12, one to three
distinct binomial factors, each exponent1..4. Every tuple checks the true
component gcd, integral genus, primitive normalization and both inequalities.
A separate180-case literal permutation-cover computation reconstructs
connected sheet orbits and each puncture's ramification without using the
gcd genus formula for those counts. Actual polynomial substitutions include
repeated initial factors not repeated in the full polynomial. Sharp equality,
power and missing-factor controls retain their distinct meanings.

    python3 04-computation/planar_jc48_sep07_universal_genus.py
    python3 -O 04-computation/planar_jc48_sep07_universal_genus.py

The finite controls do not prove the analytic embedding or the generic
constant-field consumer. The independent audit checks both maps and the
complete source, and its normal and optimized replays match the frozen
395-byte output with105404 always-active gates in each mode. The source
SHA256 is `ffed4b59e5a5832e43aa74fd325e8589c6973cc008066276be702f8f0b00f872`;
the output SHA256 is `910002d880a2849c8509c00cddb577a87155e719c48a88c66a5b79747abe31d6`.
The exact-genus formulas of the earlier
restricted families remain valid and now become quantitative refinements
of this universal lower-bound and nonrational-time theorem.
