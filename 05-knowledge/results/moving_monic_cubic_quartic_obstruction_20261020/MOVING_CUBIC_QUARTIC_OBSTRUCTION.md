# Entire monic cubic/quartic source algebras have no generating Keller pair

20 October 2026. **PROVED; independently reviewed PASS within the stated entire global-monic hypotheses.**

See [the independent audit](review/AUDIT.md) for the proof, source, and
exact-control review.

## Statement and exact boundary

Let R=C[x], K=C(x), and

    h=y³+a(x)y,   k=y⁴+b(x)y²+c(x)y,   a,b,c∈R,
    E=R[h,k] ⊂ R[y].

There do not exist P,Q∈R[y] such that

    R[P,Q]=E  and  P_x Q_y−P_y Q_x=j∈C*.

Here “entire” means polynomial over the whole affine x-line, and equality
means equality of the entire R-algebras. Generic generation over K alone
is not substituted for it. The coefficients of the chosen cubic and
quartic are globally polynomial and their leading coefficients are units.
Any global monic degree-three/degree-four presentation can be put in this
form: a polynomial triangular source translation removes the cubic
quadratic term, and subtracting a base multiple of h and base constants
normalizes the quartic and both constant terms. The source translation
has Jacobian one, so the conclusion covers that unnormalized form too.

This is a result about globally presented source algebras. It neither
proves that every generic conductor-six source has such generators nor
excludes the full conductor-six sector or settles the planar Jacobian
conjecture. The final section gives the precise conditional application to
an actual source.

## Inheritance and concept board

Closest proved mechanism: monic conductor division plus singular-fiber
cotangent descent in the moving quadratic-cylinder obstruction. Closest
hostile: generic cubic/quartic presentation need not globalize, as proved
in the separate cubic/quartic globalization boundary packet. Corrected
near miss: finite embedding classes do not identify their periods; the
second class must be retained separately. Least-used sidecar: an exact
polynomial primitive reduced in the gap basis 1,2,5.

Anchor: the entire moving obstruction. Niche: generic ambient reduction
over a non-algebraically closed coefficient field. Wildcard: replacing
an elliptic-curve/Mason input by coprime differentiation. Concept board:
immersion resultant; degeneration of collision pairs; monic conductor;
ambient embedding label; global cotangent units; period gap remainder.

## 1. Uniform fibers and the exact conductor

Over R[h], the elements 1,k,k² are a free basis of E. Indeed y is integral
over R[h] by y³+ay−h=0. Over K(h), the degree of k is three: the intermediate
field K(h,k)⊂K(y) has index dividing both 3 and 4 and thus index one.
The monic minimal relation F(x,U,V) of k over R[h] therefore has V-degree
three. Its coefficients belong to R[U] by integrality and normality, and
F is also minus Res_y(y³+ay−U,y⁴+by²+cy−V).

Every summand A_i(x,h)k^i, i=0,1,2, has y-degree congruent respectively
to 0,1,2 modulo three. Their leading terms cannot cancel, over R or after
any specialization x=ξ. Consequently

    E/(x−ξ)E = C[h_ξ,k_ξ] ⊂ C[y]                         (1)

injectively; the specialized relation is exactly the monic minimal
relation, and every fiber has degree semigroup <3,4> and delta three.
The same facts hold over K and its algebraic closure. In particular,
whole-algebra generation specializes to generation of the actual fiber
algebra, not a thicker quotient or just a set-theoretic image.

The trace-dual conductor calculation in the independently audited
[low-delta packet](../low_delta_source_conductor_exclusion_20261020/LOW_DELTA_QUADRATIC_AND_RELATIVE_GLOBALIZATION.md) applies over each field and gives

    C(x,y)=F_V(x,h,k)/h_y
          =y⁶+(a+b)y⁴+2cy³+(a²−ab+b²)y²
             +(a+b)cy+a(a−b)²+c².                         (2)

Thus C is monic of degree six in R[y], and C(ξ,y)C[y] is the exact
normalization conductor of every fiber. Formula (2) is a polynomial
identity, so specializations where h_y has roots cause no division issue.
The conductor is intrinsic and does not depend on which plane generating
pair is subsequently used. Its roots are exactly normalization points
over singular points. The global containment is also witnessed directly:

    C=h²+(b−a)k+2ch+a(a−b)²+c²,
    Cy=hk+ck+(a−b)²h,
    Cy²=k²+(a−b)h²−a(a−b)k+c(a−b)h.

Since R[y] is free over R[h] on 1,y,y², these identities give CR[y]⊂E.
Conversely any element of the global conductor is divisible by C over
K[y], and monic division makes it divisible over R[y]. Thus CR[y] is
also the exact global conductor, although only its fiberwise equality
is needed below.

## 2. Immersion and the exceptional generic singularity types

If h_y and k_y vanish simultaneously anywhere, then P_y and Q_y vanish
there for every P,Q∈E, contradicting J(P,Q)=j≠0. Hence the exact resultant

    Res_y(h_y,k_y)=16a³−48a²b+36ab²+27c²=9I,
    I=3c²+4a(b−2a/3)²                                  (3)

has no complex zero as a polynomial in x. Therefore I∈C*.
The factor nine in (3) is retained: I is a normalized immersion invariant,
not the raw resultant.

For distinct parameters r,t with h(r)=h(t), put s=r+t and p=rt.
Dividing the two coordinate equalities by r−t gives

    p=s²+a,
    H(s):=s³−(b−2a)s−c=0,
    (r−t)²=−3s²−4a,
    h(r)=h(t)=−c−(b−a)s.                                (4)

Conversely the roots of T²−sT+(s²+a) for a root of H give a collision
whenever distinct. The resultant of H and −3s²−4a is −9I, so under (3)
all these pairs are distinct-parameter pairs.

Two exceptional regimes force all coefficients to be constant.

(a) If b=a in R, then 9I=4a³+27c²∈C*. If a or c is constant, so is the
other. Otherwise gcd(a,c)=1, and differentiating gives
12a²a'+54cc'=0. Hence c divides a' and a² divides c'. Their degree
inequalities are simultaneously deg c≤deg a−1 and 2deg a≤deg c−1,
an impossibility. No Mason–Stothers or elliptic-curve theorem is needed.

(b) If H has a repeated root over an algebraic closure of K, either
H=s³, or its unique repeated root t is rational over K (explicitly
 t=−3c/[2(b−2a)]). The equations H(t)=H'(t)=0 give

    b−2a=3t²,   c=−2t³.

Since t is integral over the integrally closed ring R, t∈R. This includes
t=0 in the triple-root case. Substitution in (3) gives

    I=(4/9)(a+3t²)(4a+3t²)².                            (5)

Both factors are polynomial units. Their difference gives constant a;
then t², and therefore t, are constant. Thus a,b,c are constant.

The fixed-coefficient case is already excluded by the independently
audited [fixed-cylinder theorem](../fixed_cubic_quartic_embedding_classes_20261020/FINITE_EMBEDDING_CLASSES_AND_FIXED_CYLINDERS.md), which handles all ramified, nodal,
tacnodal and triple-point cases, not only squarefree conductors.

We may henceforth assume b−a≠0 in K and H squarefree over K. Its three
roots give three different image h-values by (4), and the six parameters
are disjoint. Thus the three equal-value conditions are independent on
Kbar[y], by interpolation at six distinct points. Their equalizer has
codimension three. Since Kbar[h,k] also has codimension three and lies in
that equalizer, the two algebras coincide. In particular all three
singularities are ordinary nodes: their delta contributions exhaust the
fiber's total delta three. Direct transversality can also be checked:

    h_y(r)k_y(t)−h_y(t)k_y(r)=(t−r)(3s²+4a)H'(s).

Every factor on the right is nonzero in the present regime. The
conductor (2) is squarefree generically,
and equals the product of the three pair quadratics:

    C(y)=Res_s(H(s),y²−sy+s²+a).                         (6)

This codimension argument also proves that equality at the three pairs
is sufficient for membership in E_K after descent from Kbar. No arbitrary
abstract gluing algebra is substituted for the actual plane algebra.

## 3. The generic embedding label is algebraic, not chosen continuously

We use the exact reduced-pair classification established in the separate
[fixed cubic/quartic continuation](../fixed_cubic_quartic_embedding_classes_20261020/FINITE_EMBEDDING_CLASSES_AND_FIXED_CYLINDERS.md). Here is its extension to any
characteristic-zero field L, and why the resulting automorphism is over L.
For A=L[y³+ay,y⁴+by²+cy], every ordered generating pair equals

    Φ(h,k),  or  Φ(k,h²+(3c/2)h),   Φ∈Aut_L L[U,V].      (7)

The second alternative is available precisely when c≠0 and 2a=3b.

Perform elementary triangular degree reductions over L: whenever one
coordinate degree divides the other, subtract the suitable power with
leading coefficient ratio in L. The process terminates and all inverse
operations remain over L. Apply the classical one-place characteristic
semigroup bound after extending constants to Lbar: a nonprincipal pair
for conductor six has initial degrees at most seven. Its semigroup is
still <3,4>. The characteristic gcd and primitive-generator properties
leave only degree pairs (3,4) and (4,6). All subsequent normalizations
use coefficients in L, not roots or algebraic extensions.

The (3,4) pair affine-normalizes to (h,k). The (4,6) pair normalizes to
(k+βh,h²+γh). The residue-separated free bases in the fixed packet force
β=0 and γ=3c/2: otherwise their deltas are respectively seven and six,
not three. For the remaining pair, c=0 makes both entries even. If
c≠0 but 2a≠3b, its delta is five. Conversely, at b=2a/3 and q=h²+3ch/2,

    q²−k³−2akq+(2a²/3)k²
      +(2a³/27+3c²/4)q−(a⁴/9+ac²)k=(c³/8)h.             (8)

Thus (7) follows over L exactly as parametrized, without precomposition
by an intrinsic automorphism. No assertion that Φ extends over R is made.

Apply (7) to L=K and the actual generic pair P,Q. This selects one of the
two reference pairs once over K. The construction, not an unexplained
continuity argument, supplies that label. In the second case, (3) and
b=2a/3 give I=3c², so c∈C*. Identity (8) now proves that its reference
pair generates E over R globally as well. In either case write the
chosen global reference pair as (u,v), with R[u,v]=E.

The only external classification ingredient here is the classical
nonprincipal degree bound and characteristic-sequence properties, pinned
in PRIMARY_SOURCES.md. The concrete reduced-pair calculations and their
K-descent are explicit above and in the audited dependency.

## 4. Monic conductor descent forces a constant ambient determinant

Choose global representatives P=G(x,u,v), Q=L(x,u,v) with G,L∈R[U,V].
Let

    D(x,y)=(G_U L_V−G_V L_U)(x,u,v)∈R[y].

Over K, (7) gives (P,Q)=Φ(u,v), with det DΦ=λ∈K*. At a singular point
of the generic reference plane curve, the differences G−Φ_1,L−Φ_2 are
multiples of its irreducible equation. The equation and both its partial
derivatives vanish there. Hence D has the value λ at all six generic
normalization points over the three nodes. By squarefreeness,

    C divides D−λ in K[y].                              (9)

Divide D by the monic C in R[y]. Its remainder is in R[y] and equals λ
over K, so λ∈R and (9) holds already in R[y]. This is why conductor
monicity is essential; no denominator-removal assumption is hidden.

Fix any ξ∈C, including fibers at which (7)'s generic automorphism has
poles or the nodes coalesce. The conductor C(ξ,y) has degree six, so its
fiber has a singular point. By (1), (u_ξ,v_ξ) and (P_ξ,Q_ξ) both generate
the same actual fiber algebra. At a singular point of a reduced plane
curve, the classes of U−U₀,V−V₀ form a basis of the two-dimensional
cotangent space m/m². The classes of P−P₀,Q−Q₀ are another generating
pair, and hence another basis. Their transition determinant is exactly
D(ξ,y₀), which is therefore nonzero at every preimage y₀ of that singular
point. From (9), it equals λ(ξ). Thus λ has no complex zero, and

    λ∈C*.                                               (10)

This handles every special-fiber conormal without assuming an ambient
lift there, without assuming constant singularity type, and without
assuming the generic automorphism is defined there.

## 5. Period transport gives an actual algebra-membership condition

For a polynomial automorphism Φ over K of determinant λ, the one-form

    Φ_1 dΦ_2−λU dV

is closed with respect to U,V, and is the differential of a polynomial
in K[U,V]. Polynomial exactness follows by integrating coefficients in
characteristic zero. Along the normalization, its primitive takes equal
values at the two parameters of each node. Consequently, for every
generic node pair r(x),t(x),

    ∫_r^t P Q_y dy = λ ∫_r^t u v_y dy.                   (11)

One can work on a small analytic base open set where all six roots and
the rational coefficients of Φ are regular. The equality is independent
of the integration path because the integrands are polynomials in y.
For any pair A,B whose endpoint values agree, direct differentiation and
integration by parts give

    (d/dx)∫_r^t A B_y dy = ∫_r^t J(A,B)dy.              (12)

The boundary term is the common value of A times the difference between
the total derivatives of the equal endpoint values of B, hence zero.
Since λ is constant, equations (11)–(12) imply that

    T(y)=j y−λ∫J(u,v)dy                                 (13)

has equal values at all three node pairs. The equalizer argument in
Section 2 gives T∈Kbar[h,k], and thus T∈K[h,k] by coefficient-field descent
(or by the residue-separated basis). This is the exact membership
consequence used below. It is stronger than one numerical period identity
but uses all three independent pairs, not just one.

## 6. Standard generic embedding: the first two gaps freeze the family

For (u,v)=(h,k), direct differentiation gives

    J(h,k)=(4a'−3b')y⁴−3c'y³
             +(2ba'−ab')y²+(ca'−ac')y.                  (14)

The primitive (13) has degree at most five. Degree five is a gap of
K[h,k], so membership first forces 4a'−3b'=0, whence

    b=4a/3+B,   B∈C.                                   (15)

Now T has degree at most four, so it lies in the span of 1,h,k. Its y⁴
and y³ coefficients force the respective k- and h-coefficients to be
3λc'/4 and −λ(2ba'−ab')/3. Comparing y² gives

    (3b−2a)c'+2ca'=0.

Using (15), this is ((2a+3B)c)'=0. Write (2a+3B)c=D₀∈C.
If D₀≠0, both polynomial factors are units, so a,c,b are constant. If
D₀=0, either 2a+3B=0, giving a,b constant and I=3c² constant, or c=0,
in which case

    I=(4a/9)(2a+3B)²∈C*

forces a and 2a+3B to be units. Thus all coefficients are constant in
every case. Equation (13) becomes T=j y, impossible in an algebra whose
degree semigroup is <3,4>. This contradicts j≠0.

For clarity, the unused final y coefficient would read
j=−λa(2ba'−ab')/3+3λcc'/4; it gives the same contradiction once frozen.

## 7. Second generic embedding: an unavoidable degree-two remainder

Here b=2a/3 and c∈C*. Put q=h²+(3c/2)h, so (u,v)=(k,q).
If a is constant the fixed theorem already applies. Otherwise the
exceptions in Section 2 do not hold: b−a=−a/3 is generically nonzero,
and the discriminant of H=s³+(4a/3)s−c cannot vanish identically unless
a is constant. Thus the common scalar descent and period argument apply.

Let

    V(y)=∫−(2h+3c/2)(2y⁴+(2a/3)y²+cy)dy.

Since c is constant,

    J(k,q)=−(2h+3c/2)J(h,k)=a' V_y.

The exact reduction identity is

    V=−k²/2−(2a/9)h²+(a²/3)k−(ac/3)h−(c²/4)y².       (16)

In particular (13), modulo K[h,k], has remainder

    j y+(λa'c²/4)y².                                   (17)

There are no nonconstant elements of degree one or two in K[h,k].
Membership forces both coefficients in (17) to vanish, in particular
j=0. This is the desired contradiction. The second embedding's period
is not identified with the standard period; its distinct remainder is
essential.

## 8. Exact controls, hostiles, and actual-source consequence

Run:

    python3 verify_moving_cubic.py
    python3 -O verify_moving_cubic.py

The script checks the raw/normalized immersion resultant, relation and
conductor identities, collision formulas, repeated-root factorization,
primitive gap equations, second-generator recovery, and its exact
primitive remainder. Controls include a generic three-node fiber, the
triple-point fiber, a repeated-H immersed fixed fiber, a ramified fiber,
and the nonconstant immersed moving family

    a=x, b=2x/3, c=1.

The last has I=3 and globally admits both generating reference pairs.
It decisively refutes any inference that immersion alone freezes the
coefficients. A standard-class-only hostile is c=x, a=3(1−x²)/4,
b=(3−x²)/2: again I=3, but 2a−3b=−3 excludes the second class.
Only whole generation plus the period argument excludes a
Keller pair. The previously audited globalization hostile remains a
separate boundary: a generic cubic/quartic presentation, even with
saturation and immersed noninjective fibers, is not a global one.

For an actual noninvertible source with generic conductor degree six,
the low-delta theorem gives generic semigroup <3,4> but no automatic
global monic presentation. Combining the present result with the separate
[saturation globalization dichotomy](../cubic_quartic_globalization_embedding_boundaries_20261020/CUBIC_QUARTIC_GLOBALIZATION_AND_EMBEDDING_BOUNDARIES.md) yields only the following:

    such a source cannot admit global monic cubic/quartic generators;
    under that dichotomy's hypotheses it must have a special fiber
    containing a quadratic (but no linear element).

No full degree-six exclusion or improved universal conductor bound is
claimed. The missing global generators remain a genuine boundary.
