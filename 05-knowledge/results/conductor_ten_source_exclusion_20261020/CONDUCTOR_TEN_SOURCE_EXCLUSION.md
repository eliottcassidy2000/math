# Genus-five Keller sources and the conductor-twelve bound

20 October 2026. **PROVED; independently audited PASS in the exact primitive monic generic-source scope.**

Prepared publication packet; not yet published.

## Result and scope

The quartic/sextic coefficient obstruction was proposed and independently
checked. Every integration and both first integrals were rederived from
the whole Jacobian, retaining scalar terms. The complete genus-five
classification and dependency synthesis were also independently audited.

The genus-five classification is as follows. The symmetric
numerical semigroups are exactly

    <2,11>, <4,6,7>, <5,6,7,8>.

The last is not planar. The only nonprincipal characteristic sequences of
the first two, including both orientations, are

    (2,11), (11,2), (4,6,7), (6,4,7).

Both plane cases have explicit immersed polynomial realizations. Thus every
actual primitive monic Keller source of generic codimension five would have
a terminal pair of degrees (2,11) or (4,6). The elementary actual-quadratic
lemma excludes the former; the audited coefficient argument excludes the
latter. Therefore **generic codimension five is impossible** under these
hypotheses.

Combining precisely this new exclusion with the independently audited
bound-ten synthesis gives generic codimension at least six and generic
normalization-conductor degree at least twelve for a noninvertible primitive
monic Keller source. Section 9 spells out that dependency; the
[independent review](INDEPENDENT_REVIEW.md) records the audit scope. This does not
resolve the planar Jacobian conjecture, does not handle arbitrary nonmonic
source algebras, and does not assert a bound on every special fiber or on the
function-field degree of the planar map.

## Inheritance and decisive distinctions

The inherited mechanism is the full degree-separated module basis, not merely
one low-degree polynomial in the algebra. The main hostile is the moving
family

    P=t^4+x t,
    Q=t^6+(3/2)x t^3+(3/8)x^2,
    J(P,Q)=-(3/4)x^2.

It satisfies every positive-degree Jacobian equation and both constant first
integrals but is not stationary. Its generic codimension is three, not five: it is a
hostile for the whole delta<=5 coefficient mechanism, not a genus-five
realization. The final requirement that the scalar
Jacobian be a nonzero constant cannot be discarded.

The corrected near miss is to omit x-dependent constant terms from P,Q while
computing J. They can be subtracted for a comparison of K-algebras, but not
for the actual Jacobian. The least-used sidecar is polynomial factor degree:
its lower bound closes the H nonzero branch after the two first integrals
have been combined into a difference of squares.

Concept board: symmetric gaps; characteristic gcd drops; polynomial
realizability; full free-module remainders; differential first integrals;
polynomial units and factor degrees. The anchor is genus-five exclusion;
the niche is classification of all t-independent Jacobians in the forced
normal form; the wildcard is its nonstationary non-Keller boundary.

## 1. Exhaustive symmetric genus-five classification

Let Gamma be symmetric with five gaps. Its Frobenius number is nine, so
exactly one of n and 9-n belongs to Gamma for 0<=n<=9. Since 1 is a gap,
8 belongs. The number 3 cannot belong, since then 3+3+3=9 would belong;
hence 6 belongs.

If 2 belongs, then 4 belongs, and symmetry excludes 5 and 7. The gaps are
{1,3,5,7,9}, and Gamma=<2,11>.

If 2 is absent, then 7 belongs. Exactly one of 4,5 belongs. Choosing 4 gives
{1,2,3,5,9}, the gaps of <4,6,7>; choosing 5 gives {1,2,3,4,9}, the gaps
of <5,6,7,8>. Closure imposes no further choice. Each displayed monoid
contains every integer at least ten, as can also be read from its Apéry set.
This is an exhaustive proof, not an extrapolation from a finite census.

For a primitive two-generated polynomial algebra A=k[p,q] in k[t] in
characteristic zero, the inherited trace/residue proof identifies its
codimension with the number of degree gaps and proves this symmetry. The
same symmetry also follows from the classical characteristic-sequence
conductor formula. No immersion hypothesis is needed for this step.

## 2. Plane characteristic sequences and terminal coordinate degrees

We use the source-checked characteristic input recorded in [PRIMARY_SOURCES.md](PRIMARY_SOURCES.md).
For a nonprincipal characteristic sequence delta_0,...,delta_h with
conductor ten, its initial entries are at most eleven and later entries
at most nine. Writing d_i=gcd(delta_0,...,delta_(i-1)), the final gcd is
one, each later entry strictly decreases the gcd, n_i delta_i belongs to
the previous monoid for n_i=d_i/d_(i+1), and

    delta_i < n_(i-1) delta_(i-1)   for i>=2.

Every primitive generator occurs among the entries. The first two entries
are the actual coordinate pole orders of a nonprincipal plane presentation.

For <2,11>, primitive 11 cannot occur later than the initial pair because
later entries are at most nine. Its partner is in {2,4,6,8,10}; 11 is prime,
so every such pair immediately has gcd one. No later entry can occur, and
the primitive 2 must already be the partner. Thus only (2,11) and (11,2)
remain. In particular the composite-nine alternate (6,9,2) from genus four
has no analogue here.

For <4,6,7>, the positive allowed initial entries are {4,6,7,8,10,11}.
All three of 4,6,7 are primitive, so an initial coprime pair is impossible.
The noncoprime nondividing pairs, ignoring orientation, are

    (4,6), (4,10), (6,8), (6,10), (8,10).

Each has gcd two, so exactly one later odd entry is possible before the
gcd reaches one. To include all three primitives, the initial pair must
therefore be 4,6 and the later entry 7. Both orientations satisfy the
membership and strict inequality axioms: 2*7=2*4+6 and 7<12.

For <5,6,7,8>, there are four primitive generators. Its possible initial
entries are {5,6,7,8,10,11}. The only noncoprime nondividing pairs are
(6,8), (6,10), (8,10), all with gcd two, permitting only one later entry.
A coprime pair permits none. No sequence can include all four primitives.
A principal initial pair can first be reduced by elementary target changes,
without changing the algebra or semigroup. Hence this symmetric semigroup
is not the degree semigroup of any such plane polynomial algebra. This
also follows from the source's primitive-prime criterion, but that extra
criterion is not needed here.

All finite bounds are cited structural bounds, not arbitrary search cutoffs.
The independent checker enumerates the entire bounded initial integer
universe, then separately verifies a target-filtered recursion.

### Base-field and normalization hypotheses

For k=C(x), extend constants to an algebraic closure only to apply the
one-place theorem. Filtered dimensions and degrees are unchanged: flat base
change preserves the kernel of A -> k[t]/k[t]_<=n. The identity recovering
t rationally from p,q also survives. The integral closure is k[t], with a
single rational infinity point. Thus the classical characteristic input
has the correct hypotheses. The terminal operations themselves take place
over k and, for an actual Keller pair below, over C.

## 3. Polynomial realizations are explicit

The quadratic-undecic case is

    p=t^2, q=t^11+t.

The power basis 1,q over k[p] has degrees 0,11 in distinct residues modulo
two, so its degree semigroup is exactly <2,11> and its codimension is five.
The derivatives have no common zero.

A nontrivial quartic/sextic realization is

    p=t^4+t,
    q=t^6+(3/2)t^3+t,
    r=q^2-p^3=2t^7-(3/4)t^6+3t^4-t^3+t^2.

The degrees 4 and 7 show that k(p,q)=k(t): the intermediate field index
divides both. The monogenic basis 1,q,q^2,q^3 over k[p] transforms
unimodularly to 1,q,r,qr, with degrees 0,6,7,13. Their residues modulo
four are distinct. Hence the exact degree semigroup is

    (4N) union (6+4N) union (7+4N) union (13+4N)=<4,6,7>,

and its codimension is (0+6+7+13-0-1-2-3)/4=5. Direct differentiation gives
gcd(p',q')=1. Swapping p,q realizes the other orientation. Thus no abstract
planar-to-polynomial realizability implication is being assumed.

## 4. Why delta<=5 forces the quartic/sextic normal form

Work first over any characteristic-zero field k. Center a monic quartic,
and for this algebra comparison only subtract scalar terms. Write

    p0=t^4+a t^2+b t,
    q0=t^6+A t^5+c t^4+B t^3+e t^2+h t.

Assume k(p0,q0)=k(t). Since [k(t):k(p0)]=4, the minimal polynomial of q0
has degree four over k(p0). Its coefficients are in k[p0] by integrality
and integral closedness. Thus 1,q0,q0^2,q0^3 is a free k[p0]-basis.
Replacing q0^2 by r=q0^2-U(p0)q0-V(p0) gives the unimodular basis
1,q0,r,q0r.

If deg r=m is odd, its basis degrees 0,6,m,m+6 have distinct residues
modulo four. No leading cancellation is possible in a nonzero linear
combination over k[p0]. Consequently the codimension is exactly (m+3)/2.

The degree-eleven coefficient of q0^2-p0^3 is 2A. If A is nonzero the
codimension is seven, so delta<=5 forces A=0. Subtracting
(2c-3a)p0q0 from q0^2-p0^3 then gives degree-nine coefficient 2B-3b.
If nonzero the codimension is six. Thus B=3b/2. No condition on h is
imposed at this stage. The next odd coefficient is

    3ab/2-2bc+2h,

and a nonzero value gives the genus-five boundary.

Restoring the actual scalar terms d,f gives exactly

    P=t^4+a t^2+b t+d,
    Q=t^6+c t^4+(3b/2)t^3+e t^2+h t+f.            (N)

For an actual monic Keller pair over C[x], the centering t=y+s(x) uses
the quartic cubic coefficient divided by four, so s and every coefficient
in (N) belong to C[x]. This polynomial source shear has determinant one.
Subtracting d,f was used only in the K-algebra argument, not as a claimed
Jacobian-preserving transformation.

## 5. Every coefficient equation, including the scalar equation

Use J=P_x Q_t-P_t Q_x, and let primes mean x derivatives. The complete
coefficient list for (N) is

    J7 = 6a'-4c',
    J6 = 0,
    J5 = -2ac'+4ca'+6d'-4e',
    J4 = -3ab'+(9/2)ba'-bc'+4cb'-4h',
    J3 = -2ae'+3bb'+4cd'+2ea'-4f',
    J2 = -2ah'+(9/2)bd'-be'+2eb'+ha',
    J1 = -2af'-bh'+2ed'+hb',
    J0 = hd'-bf'.                                  (J)

The vanishing of J7,J5,J4 integrates respectively to

    c=3a/2+C,
    e=3d/2+3a^2/8+Ca+E,
    h=b(3a/4+C)+H,                                  (I)

where C,E,H are constants. The J3 equation gives

    f'=ea'/2+3bb'/4+cd'-ae'/2.

Its full scalar primitive is

    f=C d+3ad/4-a^3/16+Ea/2+3b^2/8+F,               (S)

with F constant. No terms depending on d or f have been discarded.

Set

    v=d-a^2/4+2E/3,
    I1=Ha+3bv,
    I2=Hb+(3/2)v^2-(3/4)ab^2.

After (I),(S), the WHOLE Jacobian is the universal identity

    J = I1' t^2 + I2' t
        + H d'-(3/4)(b^2 b'+b v a').               (U)

In particular a t-independent Jacobian forces I1=L and I2=M for constants
L,M. These are the two remaining first integrals. Formula (U) was verified
both by differentiation of symbolic
functions and independently in a polynomial jet ring.

## 6. The H nonzero factor-degree obstruction

Suppose H is nonzero. The first integral yields a=(L-3bv)/H. Substituting
in the second and completing a square with

    w=v+3b^3/(4H)

gives

    (w-3b^3/(4H))(w+3b^3/(4H))
        =(L/(2H))b^2-(2H/3)b+2M/3.                (F)

If b is nonconstant, let r=deg_x b>0. The polynomial on the right is
nonzero: substitution z -> b(x) is injective on C[z], and its linear
coefficient -2H/3 is nonzero. Its degree is at most 2r.

The factors on the left are therefore both nonzero. Their difference is
3b^3/(2H), of degree 3r, so at least one factor has degree at least 3r.
The other has nonnegative degree. Since degrees of nonzero polynomial
products add, the product has degree at least 3r, contradicting 2r<3r.
This handles cancellations in w correctly; it never assumes both factors
have degree 3r and it separately rules out a zero factor.

Hence b is constant. Equation (F) now implies that w^2, and thus w and v,
are constant. The first integral makes a constant. Then (I),(S) and the
definition of v make d,c,e,h,f constant. Therefore J=0.

The argument depends on polynomiality; a factor with a pole need not have
nonnegative degree. No rational-coefficient generalization is claimed.

## 7. The H=0 unit argument, with its exact failure boundary

When H=0, the original scalar coefficient in (J) is

    j=b((3a/4+C)d'-f').

If j is a nonzero constant and all coefficients are polynomials, b must be
a nonzero polynomial unit. Hence b is a nonzero constant. Then I1=3bv=L
makes v constant; I2=(3/2)v^2-(3/4)ab^2=M makes a constant. It follows
again that all coefficients are constant, contradicting j nonzero.
The branch b=0 is already impossible because the displayed scalar j is zero.

There is a stronger complete description of the moving boundary, useful as
a hostile control. Under only J being independent of t:

- H nonzero forces every coefficient constant, as above.
- H=0 and b a nonzero constant also force all coefficients constant.
- H=0 and b=0 allow a to move, with v constant, but J=0.
- H=0 and b nonconstant force v=0, since bv is constant; then ab^2 is
  constant, forcing a=0. Formula (U) becomes

      J=-(3/4)b^2 b'=-(1/4)(b^3)'.

  This is nonconstant whenever b is nonconstant, since its degree is
  3 deg b-1>=2. Such families genuinely exist, including the x example
  in the introduction. These nonconstant-b families have generic codimension
  three, not five. Indeed, after algebra-only scalar subtraction and removal
  of the C multiple of P, write p0=t^4+b t and q0=t^6+(3b/2)t^3. Then

      q0^2-p0^3+(3b^2/4)q0=(b^3/8)t^3.

  The full module degrees 0,6,3,9 give codimension three and primitive
  fraction field. This is a hostile for first-integral stationarity, not a
  genus-five realization.

Thus the constant-nonzero contradiction is exact, and its strongest
survivor is a classification of these t-independent Jacobians, not a claim
that every first-integral family is stationary.

## 8. From terminal degrees to the complete genus-five exclusion

Let R=C[x], K=C(x), and let P,Q belong to R[y], with one coordinate monic
of positive y-degree, J(P,Q) in C*, and K(P,Q)=K(y). Suppose delta=5.

If one coordinate has degree n and leading coefficient alpha in C*, and
the other has degree m>0 and leading coefficient beta(x), then

    [y^(n+m-1)]J=-n alpha beta'.

Thus beta is a nonzero constant. Whenever one y-degree divides the other,
a constant multiple of a power of the smaller-degree coordinate can be
subtracted from the larger. This is a target automorphism over C, strictly
reduces the degree sum, and preserves a nonzero constant Jacobian and the
generic algebra. The unchanged coordinate keeps constant leading
coefficient; the same top-coefficient argument applies again.

A y-constant terminal coordinate r(x) would give J=r' Q_y in C*, making
Q linear in y. Then K[P,Q]=K[y] and delta=0, a contradiction. Hence the
process terminates in a nonprincipal pair with positive degrees. Constant
rescaling makes it monic. Sections 1–2 give only terminal (2,11) or (4,6).

For terminal (2,11), complete the square by a polynomial source shear:
P=t^2+a(x). At t=0, a' Q_t(x,0)=j in C*, so a'=alpha in C*. Thus
P=alpha x+beta+t^2, and (P,t) are polynomial coordinates. In these coordinates
J=alpha times the t derivative of Q with P held fixed. Integrating gives

    Q=(j/alpha)t+G(P),    G in C[P].

Consequently K[P,Q]=K[t] and delta=0, contradicting delta=5. This is an
actual-quadratic-coordinate proof, not globalization of a generic quadratic
semigroup generator. The all-m quadratic theorem is unnecessary at genus
five because 11 is prime and there is no alternate terminal embedding.

For terminal (4,6), Sections 4–7 apply because delta=5<=5. They give
J=0, the other contradiction. Thus no generic codimension-five source
satisfies the stated actual primitive monic Keller hypotheses.

## 9. Exact dependency synthesis: conductor at least twelve

The independently audited predecessor is the
[primitive monic source-conductor bound ten](../conductor_eight_source_exclusion_20261020/CONDUCTOR_TEN_SYNTHESIS.md).
Its SHA-256 is recorded in [DEPENDENCY_HASHES.json](DEPENDENCY_HASHES.json),
and its [exact-dependency PASS audit](../conductor_eight_source_exclusion_20261020/review/SYNTHESIS_PUBLICATION_AUDIT.md) was checked. It says that under precisely the hypotheses of
Section 8, a nonautomorphism has delta>=5 and deg_y c=2 delta, where cB is
the conductor of A=K[P,Q] in B=K[y]. Its proof covers delta=0 by one-line
injectivity and every positive delta<=4 by individually audited exclusions.
It does not assume or use the present genus-five result.

Section 8 excludes delta=5, so delta>=6. The inherited exact conductor
identity then gives

    deg_y c >= 12.

This implication has no new geometric or specialization hypothesis and no
circular dependency. This audit checks the predecessor theorem's statement,
its independent audit, and its exact file hash. It does not claim to repeat
the entire all-m quadratic proof or every older exclusion. This packet is
prepared for publication, not a claim that publication has
already occurred.

## 10. Reproducibility and boundaries

Run

    python3 independent_checks.py > results.json
    python3 -O independent_checks.py > results-optimized.json
    cmp results.json results-optimized.json

Both executions pass **88** explicit checks and have byte-identical output.
No bare Python assertion enforces a load-bearing check. The script imports
no producer code and uses both symbolic functions and an independent jet
ring. It tests the full finite characteristic universe, a second recursion,
all scalar equations, the full primitive for f, the completing-square
identity, zero-factor and moving-family hostiles, exact module remainders,
both polynomial realizations, source-centering chain rule, and actual
quadratic Keller positive controls.

The written arguments prove the universal exclusions. Computation checks
identities and explicitly bounded numerical universes; it is not a finite
search for Keller maps. No statement excludes arbitrary terminal (4,6)
pairs with delta six or seven: the wrong sextic t^3 and t^5 coefficients
are excluded only by the delta<=5 hypothesis. No semigroup generator is
silently globalized, and no scalar subtraction is silently treated as a
Jacobian-preserving target operation.
