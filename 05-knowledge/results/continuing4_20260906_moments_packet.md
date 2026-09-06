# The eighth moment recovers the two-interlacer geometry

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**

The separate bounded-support moment packet through degree seven still admits
an exact positive full response at an original phase. Adding the eighth
moment repairs the lost geometry completely: for the prescribed coefficient
pencil, positivity of the two ordinary 5-by-5 Hankel matrices is equivalent
to the full weak two-interlacer model. This is an exact model decoder and a
sharp obstruction to the lower moment relaxation. It does not prove the
remaining finite-phase sign theorem.

[Independent proof and exact audit](continuing4_20260906_moments_packet_audit.md) passes.

## 1. Inheritance and scope

The closest proved mechanism is the two-measure residue floor and original
tail2500 in `05-knowledge/results/continuing3_20260906_residue_floor.md`.
The canonical hostile is its inherited coefficient surrogate at
(x,y,s)=(104,50,15/2): lower Hankel positivity lost the support endpoint.
The corrected near miss is identifying necessary residue inequalities with
the actual beta-root geometry. The least-used sidecar here is the exact
degree-five denominator recurrence, including factors canceled from a ratio.

The concept board is: separate residue measures; bounded-support localizers;
the finite denominator quotient; common canceled factors; and the original
phase response with every carry. The new connection keeps the native
coefficient pencil and recurrence, rather than treating the moments as freely
specified sequences.

Let x,y,z be nonnegative real numbers and put

    B(v)=v^5-13v^4+55v^3-xv^2+yv-z,
    C(v)=v^4-12v^3+45v^2-(2/3)xv+(3/7)y,
    D(v)=v^4-11v^3+36v^2-(5/12)xv+(1/7)y.

The model consists of B having five nonnegative real roots, counted with
multiplicity, and both C,D weakly interlacing B. Its beta roots then have
sum13 and square sum59. Zero roots, repeated roots, and common factors are
allowed. Define the formal moments at infinity by

    C(v)/B(v)=sum_(j>=0) mu_j v^(-j-1),
    D(v)/B(v)=sum_(j>=0) nu_j v^(-j-1).

For either sequence m the initial entries are

| j | mu_j | nu_j |
|---|---|---|
|0|1|1|
|1|1|2|
|2|3|7|
|3|x/3-16|7x/12-19|
|4|16x/3-373-4y/7|115x/12-632-6y/7|

Every later entry, including m_5,...,m_8, is defined exactly by

    m_j=13m_(j-1)-55m_(j-2)+xm_(j-3)-ym_(j-4)+zm_(j-5), j>=5. (1)

No independence between these later moments is introduced.

## 2. Exact degree-eight criterion

**Theorem.** On x,y,z>=0 the full model is equivalent to

    H_C=(mu_(i+j))_(i,j=0)^4 >=0,
    H_D=(nu_(i+j))_(i,j=0)^4 >=0.                         (2)

These are positive-semidefinite conditions, not positive-definite ones.
They use only moments0 through8. No separate endpoint or shifted condition
is needed in this equivalence.

Here is a finite-recurrence lemma with its proof. Suppose A,B are real
polynomials, B monic of degree d and deg A=d-1. Write m_j for the formal
moments of A/B. If H=(m_(i+j))_(i,j=0)^(d-1) is positive semidefinite, then
A/B is a finite nonnegative-residue Cauchy transform on real nodes, after
cancellation. For m_0=1 its weights sum to one.

Work in the d-dimensional algebra E=R[v]/(B). Define the functional L by
L(v^j)=m_j for 0<=j<d. The denominator recurrence makes
L(v^j mod B)=m_j for every j. Thus H represents the form

    <p,q>=L(pq).

Multiplication T by v is selfadjoint for this form because
L((vp)q)=L(p(vq)). The radical of a positive-semidefinite form is its
nullspace, and selfadjointness makes that radical T-invariant. On the
quotient by the radical the form is positive definite and T is a real
selfadjoint operator. The spectral theorem gives

    L((w-T)^(-1)1)=sum_i w_i/(w-a_i), w_i>=0, a_i real.

Its expansion at infinity is exactly the series of A/B, hence the rational
functions are equal. The argument applies to singular H and removes any
nilpotent or nonreal denominator part invisible to the functional. Values
beyond degree2d-2 follow from the recurrence; they are not extra positivity
hypotheses.

For the prescribed C,D pencil, those invisible nonreal parts cannot coincide.
Indeed, at any common zero r!=0, solving C(r)=D(r)=0 gives

    x=(24/7)r^3-36r^2+108r,
    y=3r^4-28r^3+63r^2.                                  (3)

If r=a+ib with b!=0, reality of x implies

    b^2=3a^2-21a+63/2.

Reality of y, after this substitution, implies

    0=-6(2a-7)(2a^2-14a+21).

At a=7/2 the required b^2 is -21/4. At either root of the quadratic it is
zero. Thus **C,D have no common nonreal zero for any real x,y**.

Apply the lemma to both matrices in (2). Every nonreal B zero would have to
cancel from both ratios, contradicting (3). Thus all B roots are real.
For u>0,

    B(-u)=-(u^5+13u^4+55u^3+xu^2+yu+z)<0,

so none is negative. Each positive-residue reduced ratio has one numerator
zero between consecutive poles. Restoring the common real factors to its
numerator and denominator preserves weak interlacing: their root-counting
functions acquire the same increments. Hence both C,D weakly interlace B.
Conversely, weak interlacing gives nonnegative residues after cancellation
and therefore (2), exactly as in the inherited residue proof.

This proves the equivalence, including common and repeated roots. An exact
nonvacuous singular control is

    (x,y,z)=(648/7,54,27/7),
    B=(v-3)(7v^4-70v^3+175v^2-123v+9)/7.

Both Hankels have rank4 and are positive semidefinite; all B roots are
positive and gcd(B,C,D)=v-3. Replacing (2) by strict positivity would wrongly
discard this valid model point.

### Corollary: every fixed-(x,y) admissible fibre is a compact interval

The degree-eight criterion also gives a general version of the incoming
four-anchor fibre geometry. For fixed x,y>=0, each Hankel in(2) is affine
in z. Indeed the moments0 through4 are independent of z, and the z slopes
of moments5 through8 are respectively

    C: (1,14,130,904+4x/3),
    D: (1,15,147,1067+19x/12).

These follow directly from recurrence(1) and are independently checked
symbolically by the referee. Thus the admissible set of z>=0 is an
intersection of two closed convex positive-semidefinite matrix pencils
and a half-line. By the theorem it is exactly the full model fibre.
The decoded B roots are nonnegative with sum13, so AM-GM bounds their
product z by(13/5)^5. Each fibre is therefore empty, a singleton, or a
compact nondegenerate interval. This supplies an exact convex feasibility
problem at every fixed(x,y); it neither gives the endpoints explicitly nor
proves the carried response sign. The incoming special fibre x=84,y=35
has its stronger explicit domain[0,14sqrt2-16] in
`long_frontier_sep06_fibre_domain.md`.

## 3. A positive response survives the entire degree-seven packet

The following is an exact coefficient surrogate:

    x=78071/1000, y=601/50, s=57/2,
    z=127473806477/203203019250.                           (4)

Retain the original source and its variable t, including all Laurent shifts:

    beta=t^-1(1+13t+55t^2+xt^3+yt^4+zt^5),
    C_raw=t^-1(1+12t+45t^2+(2x/3)t^3+(3y/7)t^4),
    D_raw=t^-1(1+11t+36t^2+(5x/12)t^3+(y/7)t^4),
    O=sum_j binom(14,2j+1)t^j, E=sum_j binom(14,2j)t^j,
    P=O star beta,
    Q=(O^2+t^-1 E^2) star(beta^2+2t C_raw D_raw).

Star denotes coefficientwise multiplication. Formula (4) uses exactly

    z=12y/(7s)-x/s^2+10/s^3-1/(11s^4),

so P(-s)=0. The complete original response is strictly positive:

    Q(-57/2)=1473043757735348617617612017/28089600000 >0.   (5)

The source independently forms the literal convolutions, including the
factor2t. Every Q exponent from -1 through8 is retained and q_(-1)=28.
The original first polynomial has four simple positive phases. The beta
polynomial has exactly three real roots, all positive, and a nonreal pair.
These root counts are exact rational Sturm certificates, not floating signs.

Nevertheless put M_0=707/100. For BOTH m=mu and m=nu, the matrices

    H=(m_(i+j))_(i,j=0)^3,
    K=(m_(i+j+1))_(i,j=0)^3,
    U=(M_0 m_(i+j)-m_(i+j+1))_(i,j=0)^3                  (6)

are strictly positive definite. Every leading principal minor is recorded
as an exact positive fraction in the JSON and transcript. The degree-six
matrices (M_0 m_(i+j+1)-m_(i+j+2))_(i,j=0)^2 are also positive definite.

This support cap is stronger than the exact anchor endpoint:

    M_0 < (13+6 sqrt(14))/5 <71/10.

The first inequality follows by evaluating 5M_0^2-26M_0-67<0 with
M_0>13/5. The surrogate additionally satisfies all four strict Newton
inequalities, the proved residue floor and slope, and 5z<M_0 y. Thus this
is not a replay of the old support-endpoint violation.

In fact (6) supplies genuine probability measures matching ALL the moments
through degree seven, not only a list of determinant checks. On polynomials
of degree<=3 give the inner product Gram matrix H and let T=H^(-1)K.
It is selfadjoint and 0<T<M_0 I. For the monomial basis e_0,...,e_3,

    T e_j=e_(j+1), j=0,1,2.

Consequently e_0 is cyclic. The spectral measure of e_0 has four distinct
nodes in (0,M_0), strictly positive weights and total mass1. It reproduces
m_j for j<=6 by splitting j into two exponents<=3. For j=7 use

    <e_0,T^7 e_0>=<e_3,T e_3>=K_(3,3)=m_7.

The source constructs each rational quartic denominator and cubic numerator
for this four-atom transform, verifies its four interior roots, and checks
formal equality of moments0,...,7. Both transforms disagree with the native
moment8. Their canonical quadrature supports are disjoint. This last fact
is a property of these two constructions, not a claim that every pair of
representing measures must have disjoint supports.

Therefore every necessary separate-measure inequality involving only
moments0 through7 and nonnegative polynomials on the anchor interval is
satisfied by (4). Such a packet cannot prove nonpositivity of the response.

## 4. The exact failed gate and the repaired question

At (4) the first omitted determinants are already negative:

    det H_C =
     -4589218342226103628274440890013740853428732686393332483
      /95479173982126390784148568451886227718750000000000000,

    det H_D =
     -78037916605582141428435773185728629446352083104557342123
      /1251464629218527029285992116412563163955200000000000000.

The source-to-target map is now precise. The original model maps to two
degree-seven positive moment functionals with a common support interval.
That map preserves positivity and the initial recurrence identities but
loses positivity of the degree-four square, the last norm needed to close
the degree-five quotient algebra under multiplication. The hostile proves
that the loss affects the actual carried original-phase response.

Adding the two degree-eight Hankels restores the model exactly, by Section2.
Thus degree8 is sharp in this hierarchy of separate initial bounded-support
moment packets. This is not a minimality claim among arbitrary inequalities
or cross-measure sidecars. The full finite-phase sign question is equivalently
the three-coefficient semialgebraic problem (2), x,y,z>=0, P(-s)=0, with
the original response and s in (1/80,2500). Solving that problem remains OPEN.

## 5. Reproduction and controls

The standalone [producer](../../04-computation/continuing4_20260906_moments_packet.py),
[transcript](continuing4_20260906_moments_packet.out), and
[exact certificate](continuing4_20260906_moments_packet_certificate.json)
use only Python's standard library. Run from the repository root (or substitute the full source path):

    python 04-computation/continuing4_20260906_moments_packet.py
    python -O 04-computation/continuing4_20260906_moments_packet.py

The explicit universe is one exact hostile, two full-model controls (one
strict, one with a common canceled node), the inherited C-only repeated
boundary, all stated matrix minors, and all named exact polynomial identities.
The two four-atom transforms are reconstructed from the moment
linear systems and checked by formal division. The strict model control has
literal beta roots (1,3,9,22,30)/5 and direct positive-residue checks. The
common-node control has a separate rational-interval proof of both reduced
interlacings. The C-only boundary (0,0,3,5,5) has a
positive C transform but D's ordinary3 determinant is -37/16, so it remains
rejected. Numerical optimization was only a locator for (4); no optimization
output is used as evidence. All verification gates remain active under -O.

Normal and optimized runs both pass **128 always-active exact gates** and
produce byte-identical LF transcripts. Frozen identities:

- Source SHA256: `1b7d02c20c631f95b8a539af8948bf92a9b6dc1504ddc05d5a00b6b13892a9bc`.
- Output SHA256: `54e693133895b249e31aeeed325d786ee4ed6c46d8343fa082e802adcbf6932e`.
- JSON SHA256: `e396e221c4cd8e4b668f0f07a572ed5314b25cb7c66a8ac68eb439d84979170f`.

The independent proof audit passes, including the compact-interval corollary.
The producer phase made no maintained-file or Git changes.
