# Independent audit: constant leading coefficients and the conductor-six boundary

20 October 2026. **PASS for the precise theorem below.**

All reviewer writes are confined to `session13-constant-leading-review`.
No repository or remote write was performed. The parent owns publication.

## Verdict

The proposed globalization mechanism is valid. Its hypotheses must include the
actual nonzero constant Jacobian, not merely immersion or saturation. All
ordinary degree reductions can be made over C, so the reduced pair remains a
Keller pair globally. The exceptional reduced (4,6) embedding has enough
coefficient rigidity to recover the canonical cubic over C[x].

A stronger independent verification is available: once reduced to (3,4) or
(4,6), direct Jacobian coefficient equations give contradictions. This avoids
using the moving-period theorem for the final exclusion and avoids needing
saturation for the globalization statement itself.

**Precise theorem.** Put R=C[x], K=C(x). There do not exist P,Q in R[y]
with one coordinate monic of positive y-degree, J(P,Q) in C*, and
Gamma(K[P,Q])=<3,4>.

In particular, under the already proved actual-source monic/primitive
hypotheses of the low-delta packet, generic conductor degree six is excluded:
the quadratic alternative was excluded there, and the remaining <3,4>
alternative is excluded here. The resulting bound is at least eight for a
noninvertible source in that stated scope. This is not an exclusion of all
source semigroups, not a bound on the function-field degree of the planar map,
and not a resolution of JC(2).

## Inheritance, sidecars, and controls

Closest proved mechanisms: the exact reduced-pair classification in
`fixed_cubic_quartic_embedding_classes_20261020`, the moving monic entire
obstruction in `moving_monic_cubic_quartic_obstruction_20261020`, and the
saturation/filtration statements in the low-delta and globalization-boundary
packets. All four requested proof packets and the fixed-packet audit were read.

Canonical hostile: R[y^4, x y^3+y^2+y], whose Jacobian is -4y^6.
Corrected near miss: global monicity does not follow from generic degree data
and saturation alone. The new sidecar is the *top Jacobian coefficient after
every reduction*, not just the initial leading coefficient.

Concept board: constant-leading gate; constant target reductions; one-place
characteristic bounds; exceptional embedding coefficients; polynomial units;
direct terminal Jacobian coefficients.

A second hostile starts with two monic generators,
(y^4, y^8+x y^3+y^2+y), but cancellation exposes the nonconstant coefficient x.
Its Jacobian is still -4y^6. Thus initial constant leading coefficients alone
do not prove their persistence. The Keller equation is essential at every step.

## 1. Constant-leading gate and reduction over C

Write p=y^n+lower terms, q=b(x)y^m+lower terms, where n,m>0. Then

    [y^(n+m-1)] J(p,q) = -n b'(x).

Indeed p_x has y-degree at most n-1, while the only degree n+m-1
contribution comes from -p_y q_x. Since n+m-1>0 and J is constant,
b'=0. Thus b is a nonzero complex constant and q can be made monic
by constant scaling.

More generally, after both coordinates have constant nonzero leading
coefficients, whenever one y-degree divides the other, subtract the unique
constant multiple of the corresponding power of the smaller coordinate.
This is an ordinary triangular target automorphism over C; it preserves
R[p,q] and the nonzero constant Jacobian. One unchanged coordinate remains
monic after constant scaling, so the same gate makes the new positive-degree
coordinate's leading coefficient constant. The sum of the positive degrees
strictly decreases.

A constant terminal coordinate is impossible in the <3,4> sector. Over K,
that would leave K[r] as the entire algebra. Since its fraction field is K(y),
r would have degree one, contradicting its semigroup. Equivalently, a direct
Jacobian calculation with a y-constant coordinate forces the other degree
to be one. A zero coordinate would have zero Jacobian. Hence termination gives
two positive degrees, neither dividing the other.

There is no x-dependent target transformation in this argument. Later
x-dependent affine expressions are used only to compare generators or compute
Jacobians; they are never asserted to preserve the Keller determinant.

## 2. Why the only terminal degrees are (3,4) and (4,6)

The degree semigroup has precisely three gaps, 1,2,5. Choose monic elements
h,k of degrees three and four. Coprime degrees give K(h,k)=K(y), and the
K[h]-basis 1,k,k^2 has distinct residues modulo three. Therefore K[h,k] has
exactly those three gaps and equals the original algebra by equal codimension.
After a translation t=y+s over K and changes of these chosen generators,

    h=t^3+a t,   k=t^4+b t^2+c t,   a,b,c,s in K.

A plane presentation of this algebra has normalization K[y] and exactly one
place at infinity. Base extension to an algebraic closure preserves the degree
semigroup: the finite-degree coefficient kernels and degree filtration commute
with extension of the scalar field. The nonprincipal characteristic-sequence
bound therefore applies to the terminal pair.

The two initial pole orders really are its coordinate degrees. One can see
this directly by taking the formal parameter tau=p(y)^(-1/n): q has leading
term of order -m, so the first characteristic pole orders are n,m. The generic
coordinate algebra is geometrically integral and its infinity place is
residually rational, as required. No abstract classification of unrelated
semigroups is substituted for the given plane presentation.

The cited bound gives initial orders at most conductor+1=7. Thus each is in
{3,4,6,7}. Nondividing unordered pairs are (4,3), (6,4), (7,3), (7,4),
and (7,6). A pair involving seven already has gcd one, so the characteristic
sequence cannot subsequently add the missing primitive generator three or
four. The pair (6,4) has gcd two and its next, final order must be three.
Thus the terminal degree pairs are exactly (3,4) and (4,6).

Primary-source check: Sathaye--Stenerson, *On Plane Polynomial Curves*, author
text dated 10 September 2014, was freshly opened at
https://www.ms.uky.edu/~sohum/ma561/affgeom/papnew.pdf. Section 2.2 gives the
field and one-place conventions; Section 3.1, printed pages 8--9, gives the
characteristic gcd axioms, nonprincipal bound and primitive-generator rule.
The frozen packet's Abhyankar--Sathaye source also identifies the coordinate
orders and elementary reductions. Fresh DOI and AMS-PDF retrievals failed;
no fresh downloaded copy of that paper is claimed. The indexed original and
the frozen independent source audit agree with the use here. No general
embedding-uniqueness theorem is used.

## 3. Independent audit of the exceptional (4,6) classification

For A=K[h,k], its low-degree filtration is

    A_<=4 = K + K h + K k,
    A_<=6 = K + K h + K k + K h^2.

A monic degree-(4,6) generating pair therefore has the form

    p=k+beta h+d,
    q=h^2+eta k+gamma0 h+e.

Subtracting d and then eta times the first coordinate, and removing a
constant, gives a pair k+beta h, h^2+gamma h over K. This comparison preserves
the K-algebra, but it is not used as an x-dependent Keller transformation.

The following is a direct free-module calculation, valid over any
characteristic-zero field. If a monic quartic p and sextic q generate the
field K(y), their algebra is free over K[p] with basis 1,q,q^2,q^3. Replacing
the last two entries by r,qr, where r=q^2-A(p)q-B(p), is unimodular. When
its four basis degrees d_i have distinct residues modulo four, the exact
codimension is sum_i(d_i-(d_i mod 4))/4.

- If beta is nonzero, r=q^2-p^3 has leading term -3 beta t^11.
  The basis degrees 0,6,11,17 give codimension seven, not three.
- Hence beta=0. After subtracting (4a-3b)pq from q^2-p^3, the degree-nine
  coefficient is 2 gamma-3c. If nonzero, degrees 0,6,9,15 give codimension six.
- Hence gamma=3c/2. If c=0, both coordinates are even, so they cannot generate
  A. If c is nonzero but 2a differs from 3b, subtracting -2a(a-b)p^2 gives a
  degree-seven remainder with leading coefficient -c(2a-3b)/2. The degrees
  0,6,7,13 give codimension five.

Consequently beta=0, gamma=3c/2, b=2a/3, c nonzero. The original terminal
pair is exactly

    p=k+d,
    q=q0+eta k+e,   q0=h^2+(3c/2)h.                 (1)

No intrinsic normalization automorphism, ambiguous embedding label, or
uniform relative lift has entered this conclusion. It is an exact identity
for the given parametrization.

## 4. Global coefficient recovery, with two valid endings

Since p is a global monic quartic and p=k+d, its y^3 coefficient is 4s.
Thus s belongs to R. The coordinate t=y+s is therefore a polynomial source
shear of determinant one. In t-coordinates, the coefficients of t^2,t,1 in
p are b,c,d; hence they belong to R. Since a=3b/2, a also belongs to R.
Thus the canonical cubic h=t^3+a t is itself a global polynomial.

The coefficient of t^4 in q is 2a+eta, so eta belongs to R; its constant
coefficient is e, so e belongs to R. Formula (1) therefore holds globally.
At t=0, its Jacobian is

    J(p,q)|_(t=0) = c ((3a/2+eta)d' - e').          (2)

The nonzero constant Jacobian makes c a unit of C[x], namely c in C*.
This is also the fiberwise common-derivative argument: if c vanished, both
p_y and q_y would vanish at t=0 on that fiber.

### 4A. Direct algebraic globalization, without saturation

The universal identity

    q0^2-k^3-2a kq0+(2a^2/3)k^2
      +(2a^3/27+3c^2/4)q0-(a^4/9+a c^2)k=(c^3/8)h

recovers h globally because c is a unit and k=p-d,
q0=q-eta k-e. Therefore R[p,q]=R[h,k]. This proves the proposed monic
cubic/quartic globalization, without any saturation hypothesis.

### 4B. The originally proposed saturation argument is also sound

In the actual-source setting, the inherited all-fiber primitivity statement
makes S saturated in R[y], so h in S_K intersect R[y] implies h in S.
The algebra T=R[h,p] has monic degrees three and four on every fiber. Those
fibers are primitive because their field index divides both three and four;
the inherited monic saturation lemma makes T saturated. Generic equality
T_K=S_K then implies S=T by clearing a base denominator and saturating.
One must use saturation of T for this last step, not just saturation of S.
The original sketch did not state that distinction, but all required
hypotheses hold. Thus this is a repair of exposition, not a mathematical gap.

## 5. Independent endpoint: direct exclusion of terminal (4,6)

Work in the global t-coordinate of Section 4, with c in C*. Put A=a'.
For the exact pair (1), the coefficient of t^7 gives eta'=-A. The coefficient
of t^5 then gives d'=2aA/9. The coefficient of t^3 gives

    e'=a^2 A/3+2a eta A/9.

Substituting these three identities into the *entire* Jacobian gives exactly

    J(p,q) = -(a' c^2/2)t.                         (3)

This cannot be a nonzero constant. This independently excludes the exceptional
terminal pair without the global recovery identity, saturation, conductor
monicity, a fiber singularity classification, or moving collision periods.

## 6. Independent endpoint: direct exclusion of terminal (3,4)

Constant rescaling and swapping put the terminal pair in monic degrees
three and four. A polynomial source shear centers the cubic. Write

    p=t^3+a t+d,
    q=t^4+eta t^3+b0 t^2+c0 t+e0,

with all coefficients in R. The t^5 coefficient of J is -3 eta', so eta
is constant. A constant target shear q -> q-eta p gives

    p=t^3+a t+d,
    q=t^4+b t^2+c t+e.                            (4)

For (4), the positive-degree coefficient equations are

    4a'-3b'=0,
    -3c'+4d'=0,
    2b a'-a b'-3e'=0,
    c a'-a c'+2b d'=0.

Writing b=4a/3+B with B in C, they imply

    d'=3c'/4,
    e'=(4a/9+2B/3)a',
    ((2a+3B)c)'=0.                                (5)

The constant coefficient gives

    j=(3c^2/8-4a^3/27-Ba^2/3)'.                   (6)

Let (2a+3B)c=D in C. If D is nonzero, both factors are polynomial units,
so a,c are constant, contradicting j nonzero. If D=0 and c is not identically
zero, then a=-3B/2 is constant, and (6) requires c^2 to have degree one.
This is impossible for a polynomial c. Finally, if c=0 identically, (6)
requires 4a^3/27+Ba^2/3 to have degree one. If a is constant its derivative
is zero; if a is nonconstant that polynomial has degree 3 deg(a), again
impossible. Thus no terminal (3,4) Keller pair exists.

The argument does not assume that immersion freezes a general moving
cubic/quartic family. It uses the actual terminal pair's exact coefficient
equations, including the scalar terms d,e and the initially present cubic
term eta. Omitting those terms would be an invalid normalization.

## 7. Reproduction and status boundaries

Run from this review directory:

    python3 independent_checks.py > results.json
    python3 -O independent_checks.py > results-optimized.json
    cmp results.json results-optimized.json

All 77 checks pass, with identical outputs. They import no producer code and
use explicit exceptions, not optimization-sensitive assertions. The tests
include 36 bounded leading-gate coefficient probes, a genuine Keller
positive control reducing to the affine-line sector, both globalization
hostiles, every universal exceptional-class coefficient and recovery
identity, polynomial normalization readouts, both full terminal-Jacobian
eliminations, and the finite characteristic-pair enumeration.

The controls corroborate but do not replace the universal arguments above.
The only classical theorem needed for the independent exclusion is the
one-place characteristic-semigroup theorem with its nonprincipal initial
order bound. The direct (4,6) classification is completely replayed in
Section 3. The moving monic cubic/quartic theorem is a separate, consistent
endpoint for the valid globalization proof, but is not a dependency of
Sections 5--6.

No counterexample to the stated mechanism survived. No issue remains in the
normalizations, constant-leading invariance, field descent, global coefficient
recovery, or generic-to-global equality once the target saturation distinction
in Section 4B is made explicit. The theorem remains sharply scoped to generic
semigroup <3,4> and a monic original coordinate.

The producer checker was also replayed independently: its 130 checks pass both
normally and with `-O`, with byte-identical outputs saved as
`producer-replay.json` and `producer-replay-optimized.json`.
`DEPENDENCY_HASHES.json` pins the four frozen proof dependencies and the
reviewed producer snapshot; it is a review-time snapshot, not a restriction
on subsequent explicitly reviewed edits. The only requested producer
clarification was to name subtraction of the quartic constant d before the
normalized generic basis calculation.
