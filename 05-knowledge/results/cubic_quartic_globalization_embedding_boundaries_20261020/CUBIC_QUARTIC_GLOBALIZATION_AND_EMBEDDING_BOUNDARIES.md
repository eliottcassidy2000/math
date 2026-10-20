# Cubic/quartic source algebras: two sharp failure boundaries

20 October 2026. **PROVED below; independently reviewed PASS.**

See `review/AUDIT.md` for the proof and primary-source audit, including
74 producer checks and 56 independent exact controls.

## Outcome and scope

The audited low-delta packet leaves a hypothetical noninvertible planar Keller
source of generic conductor degree six in the following sector:

    S=C[x,P,Q] ⊂ B=C[x,y],   K=C(x),
    Γ(S_K)=<3,4>.

It does not give cubic/quartic generators over C[x]. This note identifies two
genuine obstacles to extending the quadratic-source argument:

1. **Base saturation does not globalize a cubic**, even with a monic original
   generator, immersed primitive slices, a principal conductor, and a collision
   on every slice. An explicit family has generic delta three and special
   delta one, with no global monic cubic.
2. **Unique plane embedding fails inside the immersed rational <3,4> sector.**
   The ordinary-triple-point algebra C[y³,y⁴+y] has two inequivalent plane
   embeddings. An explicit period calculation proves inequivalence.
3. **A precise positive survivor remains:** the primitive cubic is canonical up
   to scalar and constants. Its nonunit leading coefficient forces a
   lower-degree special-fiber element. If no fiber contains a nonconstant
   polynomial of degree at most two, a global monic cubic/quartic presentation
   of the entire algebra does follow.
4. In the actual noninvertible Keller setting, the lower-degree element cannot
   be linear, by the cited one-line theorem. Thus failure of this globalization
   route is localized at genuine quadratic special fibers.

The examples are **not Keller maps**. The note neither excludes all generic
cubic/quartic Keller source algebras nor improves the audited universal
generic-conductor lower bound from six to eight. No historical novelty is
claimed for nonunique plane embeddings or for the elementary period invariant.

## Inheritance and concept board

Closest proved mechanism: Sections 3 and 6 of the frozen
[low-delta packet](../low_delta_source_conductor_exclusion_20261020/LOW_DELTA_QUADRATIC_AND_RELATIVE_GLOBALIZATION.md),
including its [independent review](../low_delta_source_conductor_exclusion_20261020/review/AUDIT.md),
give the generic <3,4> alternative and base-saturation argument. The
[moving quadratic-source packet](../cylinder_generator_period_obstructions_20261020/MOVING_HYPERELLIPTIC_GENERATOR_OBSTRUCTION.md)
gives entire generator obstructions only in its stated quadratic models.

Canonical hostile: the previously audited family
C[x,y³,xy²+y] is saturated but loses its only collision on x=0.
The new hostile retains a collision there. Corrected near miss: one must not
confuse the Jacobian of selected normal-form generators with the Jacobian of
an arbitrary entire generating pair. Least-used sidecar: the ratio between a
collision period and the cotangent determinant at a singular point.

Concept board:

- Generic degree filtration versus its integral C[x]-lattice.
- Degree-drop fibers versus injective fibers.
- Entire generation versus ambient equivalence of embeddings.
- Collision periods versus local cotangent determinants.
- Characteristic degree sequences versus field/sheet degree.

Portfolio: the anchor is relative cubic/quartic globalization, the niche is
inequivalent plane presentations, and the wildcard is the local-determinant
versus collision-period test.

The first bridge below preserves the generic curve and saturation but loses
global monicity. The second preserves the abstract algebra and its
normalization, conductor and local singularity, but changes its plane
embedding. Neither loss is repaired by retaining conductor degree alone.

## 1. A saturation lemma, with its full hypotheses

Let R=C[x], B=R[y], and let T=R[f,g] with f monic in y of positive degree n.
Assume

    C(f_a,g_a)=C(y) for every a∈C,

where f_a=f(a,y), and likewise for g_a. Then:

    T/(x−a)T → C[y]

is injective for every a. In particular, T is R-saturated in B:

    0≠d(x)∈R, b∈B, d b∈T  ⇒  b∈T.               (1)

**Proof.** B is free over D=R[f] with basis 1,y,...,y^(n−1).
Multiplication by g has a monic characteristic polynomial H(x,u,v) of
v-degree n. Generic primitivity makes it the minimal polynomial. Because g
is integral over the normal domain R[u], its minimal polynomial has
coefficients there. Thus

    T≅R[u,v]/(H),

and 1,g,...,g^(n−1) is a D-basis. Specialization of the multiplication matrix
gives H(a,u,v). Slice primitivity and deg f_a=n make this the minimal
polynomial of g_a over C(f_a), again of degree n. This proves the claimed
injection. If (x−a)b∈T, the injection says (x−a)b∈(x−a)T, so cancellation
gives b∈T. Factoring d into linear factors proves (1). QED.

The same setup also permits the trace-dual conductor formula over the
two-dimensional base D, not just over a PID. The two power bases have
perfect trace pairings with their respective derivative denominators:

    B†=(1/f_y)B,     T†=(1/H_v(f,g))T.

Their anti-triangular pairing matrices have unit anti-diagonal. The
double-dual argument therefore gives

    c=H_v(f,g)/f_y ∈ B,    conductor(T,B)=cB.       (2)

No Smith normal form over D is asserted. On each fiber the usual one-variable
formula supplies deg c_a=2δ_a when c_a is its conductor. The explicit
calculations below verify that specialization here has exactly that property.

## 2. A saturated immersed family with no global monic cubic

Set

    P=y⁴,     Q=xy³+y²+y,     S=R[P,Q].             (3)

### 2.1 Immersion, primitivity, and noninjectivity on every slice

The derivative pair is

    P_y=4y³,     Q_y=3xy²+2y+1.

They never vanish together, because the first only vanishes at y=0, where
the second is 1. Every slice is immersed.

For a≠0 the degrees of P_a,Q_a are 4,3. The field extension
[C(y):C(P_a,Q_a)] divides both degrees and is one. For a=0, put

    v=y²+y,     u=y⁴.

The identity

    u=v²+v−(2v+1)y

gives y=(v²+v−u)/(2v+1), so this slice is primitive as well. Section 1 applies:
S is base-saturated, its quotient specialization is its polynomial image, and
B is its finite normalization.

For a≠0, the monic cubic Q_a/a and monic quartic P_a generate an algebra
with degree semigroup <3,4>. Indeed the basis 1,P_a,P_a² over C[Q_a]
has leading degrees 0,4,8 modulo 3. Therefore δ_a=3.

At a=0, make the constant source translation t=y+1/2. Then

    y²+y=t²−1/4,
    y⁴=t⁴−2t³+(3/2)t²−(1/2)t+1/16,

and hence

    S_0=C[t², t(t²+1/4)].                         (4)

Its conductor is (t²+1/4)C[t], its delta is one, and it has an ordinary
node, with normalization preimages t=±i/2. In the original coordinate the
collision pair is y=(-1±i)/2. Thus the x=0 slice is not injective.

For a≠0, choose r with r²=−1/a. The distinct points r,−r have the
same P_a-value and the same Q_a-value, since

    Q_a(r)−Q_a(−r)=2r(ar²+1)=0.

Consequently **no slice is injective**. The globalization failure is not
explained by an exceptional affine-line fiber.

### 2.2 Exact relation and entire conductor

The monic defining relation, with u=P and v=Q, is

    H(x,u,v)
      =v⁴−(4x+2)uv²−4uv−4x²u²v−x⁴u³
        +(2x²−4x+1)u²−u.

It is the characteristic polynomial of multiplication by Q in the basis
1,y,y²,y³ over R[P]. Generic primitivity makes it minimal. Formula (2)
gives the entire conductor:

    c=(xy²+1)(x²y⁴+2xy³+2y²+2y+1).               (5)

The exact identity is H_v(P,Q)=4y³c. Thus the conductor has generic
y-degree six and specializes at x=0 to

    c_0=2y²+2y+1=2(t²+1/4).                       (6)

For a≠0, c_a has degree six and is the slice conductor; for a=0 it is
the degree-two conductor from (4). There is no assumption of constant
fiber delta or a monic conductor in this calculation.

### 2.3 Why no global monic cubic exists

Over K, let

    h=Q/x=y³+(y²+y)/x,      k=P=y⁴.

The degree filtration of S_K has gaps 1,2,5. In particular,

    (S_K)_{≤3}=K⊕K h.

Every monic cubic in S_K is h+r(x), r∈K. Its coefficients of y² and y
are both 1/x, so it cannot belong to R[y]. Therefore S contains no monic
cubic in y. No change of its selected generators repairs this fact.
Nor can an R-defined affine source change y↦u y+b(x), u∈C*, repair it:
pulling back a monic cubic under such a change would give a monic cubic
in the original coordinate after constant scaling.

For completeness,

    S∩R[y]_{≤3}=R⊕R(xy³+y²+y).                    (7)

If rQ has polynomial coefficients for r∈K, its y-coefficient r is
polynomial. Conversely the displayed R-module is clearly contained in S.

### 2.4 Exact failure boundary

The ordinary planar Jacobian of the displayed pair is

    J(P,Q)=−4y⁶,

not a nonzero constant. Thus (3) is not a Keller counterexample.
It refutes the implication from the following collection of properties
to a global monic cubic/quartic presentation:

    monic original coordinate + primitive immersed slices
    + base saturation + principal conductor
    + generic semigroup <3,4> + no injective slice.

The missing ingredient, if a stronger globalization theorem holds for actual
Keller sources, must use more of constant Jacobian than slice immersion and
the prohibition on injective source lines.

## 3. What saturation really does globalize

### 3.1 The canonical primitive cubic and its degree-drop divisor

Let S be any R-saturated R-subalgebra of B with
Γ(S_K)=<3,4>. There is a polynomial

    H=α(x)y³+β(x)y²+γ(x)y∈S,
    α≠0,   gcd(α,β,γ)=1,                         (8)

and

    S∩R[y]_{≤3}=R⊕R H.                          (9)

H is unique up to multiplication by a nonzero complex scalar; allowing a
constant term adds an arbitrary element of R.

**Proof.** Choose a generic cubic, remove its constant coefficient, and clear
denominators. Saturation permits division by the gcd of the remaining three
coefficients, giving (8). The generic degree-≤3 space has basis 1,H.
For any polynomial rH+s with r,s∈K, evaluation at y=0 gives s∈R.
Its other coefficients are rα,rβ,rγ. Their polynomiality and the primitive
coefficient condition force r∈R, by a Bezout identity in R. This proves (9)
and uniqueness. QED.

In particular:

- A global monic cubic exists if and only if α is a unit.
- At each zero a of α, H_a is nonconstant of degree at most two.
- If no slice contains a linear polynomial, then β(a)≠0 at each zero of
  α. Consequently gcd(α,β)=1, and every such degree drop is exactly
  quadratic.

The last statement applies to an actual noninvertible Keller source:
a linear element in its slice algebra would make that slice C[y], hence
would make the Keller map injective on a source line. The cited
Gwoździewicz theorem would then make it invertible.

The leading-coefficient ideal of the cubic filtration is precisely (α).
For the hostile (3), it is (x), and H_0=y²+y. This identifies the missing
integral coordinate, rather than merely recording a failed normal form.

### 3.2 A sufficient relative globalization theorem

Assume in addition that, for every a∈C, the polynomial image S_a contains
no nonconstant element of degree at most two. Then S has monic generators
H,K of respective degrees three and four over R:

    S=R[H,K].                                    (10)

An R-defined source translation and changes of these chosen generators put
the entire algebra in the form

    R[t³+a(x)t, t⁴+b(x)t²+c(x)t],
    a,b,c∈R.                                     (11)

**Proof.** Section 3.1 makes α a unit, so H can be chosen monic. Subtract
its constant term. Choose a generic quartic in S_K, subtract its constant
term and a K-multiple of H to remove its y³ coefficient, and clear
denominators. Divide the gcd of the remaining coefficients using
saturation. The resulting element has the form

    K₀=d(x)y⁴+e(x)y²+f(x)y,
    d≠0,   gcd(d,e,f)=1.

A zero of d would give a prohibited nonconstant degree-at-most-two
element in that special fiber. Thus d is a unit and a monic quartic K
exists in S.

Let T=R[H,K]. Every slice has degrees three and four, hence is primitive.
Section 1 makes T saturated. Over K the algebra T_K has semigroup
<3,4>, so T_K=S_K: both have codimension three in K[y].
For s∈S, some nonzero base polynomial d satisfies ds∈T. Saturation of
T gives s∈T, proving (10).

Finally translate y by one third of the quadratic coefficient of H,
subtract the generators' constant terms, and subtract a polynomial multiple
of the resulting monic cubic from the quartic to remove its cubic term.
All coefficients remain in R. This gives (11). QED.

For an actual noninvertible Keller source of generic conductor degree six,
this proves the following **disjunction**, not an exclusion:

    either a global entire monic cubic/quartic presentation exists,
    or some special fiber contains a quadratic.

The hypothesis excluding special-fiber quadratics is sufficient, not
claimed necessary. It has not been established for all Keller sources.
The generic absence of a quadratic does not exclude a quadratic special
fiber, as (3) demonstrates.

## 4. Two inequivalent embeddings of an immersed <3,4> curve

Let

    h=y³,     k=y⁴+y,     A=C[h,k],
    z=8h²+12h+3=8y⁶+12y³+3.                      (12)

The standard implicit equation is

    F(U,V)=V³−U(U+1)³.

Its normalization conductor is (y³+1)²C[y]. The three roots of y³=−1
map to the ordinary triple point (U,V)=(-1,0). The derivative pair
(3y²,4y³+1) has no common zero, so the normalization is immersed.

The exact polynomial identity

    z²−64k³−9=8h                                 (13)

proves C[k,z]=A. Thus these pairs define closed embeddings of Spec A
into the plane, with normalization parametrizations

    i₁(y)=(h(y),k(y)),      i₂(y)=(k(y),z(y))

of degrees (3,4) and (4,6), respectively. The maps i₁ and i₂ from the
normalization A¹ are not themselves embeddings: each identifies the three
triple-point preimages. In the latter presentation, the cancellation (13)
recovers degree three.

### 4.1 A hypothetical exact ambient lift has determinant four

Suppose a polynomial automorphism Φ of the plane satisfied

    Φ(h(y),k(y))=(k(y),z(y)).                       (14)

Let T(U,V)=(V,8U²+12U+3). The two components of Φ−T vanish on F=0,
so each is a multiple of F. At the singular point p=(-1,0) both F and
its gradient vanish; hence

    DΦ(p)=DT(p),     det DΦ(p)=−(16(−1)+12)=4.

The Jacobian determinant of a polynomial automorphism is a nonzero
constant, so det DΦ is identically four.

### 4.2 Collision periods force a different determinant

For any polynomial automorphism with constant determinant λ, the
polynomial one-form

    Φ*(U dV)−λ U dV

is closed, hence polynomial-exact on C². Pullback to the normalization
and integration between any two points with equal (h,k)-image kills
the exact term. Applied to (14), this requires

    ∫_r^s k dz = 4 ∫_r^s h dk                    (15)

for distinct roots r,s of y³=−1.

The two polynomial primitives are

    L₁(y)=(4/7)y⁷+(1/4)y⁴,
    L₂(y)=(24/5)y¹⁰+12y⁷+9y⁴.

At y³=−1 they reduce to

    L₁(y)=(9/28)y,     L₂(y)=−(9/5)y.

Since r≠s, the first period is nonzero, and

    (L₂(s)−L₂(r))/(L₁(s)−L₁(r))=−28/5≠4.

This contradicts (15). No exact ambient lift (14) exists.

### 4.3 Why this also rules out unparametrized ambient equivalence

Every abstract automorphism of A extends to its normalization C[y], hence
is y↦αy+β, with α≠0. It preserves the conductor divisor, whose support
consists of the three roots of y³+1. Their average is zero, so β=0.
Preserving the root set then gives α³=1.

Conversely each y↦αy with α³=1 preserves A and lifts in the standard
embedding by the ambient linear map

    (U,V)↦(U,αV).

If any ambient automorphism carried the image of i₁ to the image of i₂,
its induced abstract isomorphism would differ from the fixed identification
by one of these automorphisms. Precomposing by the inverse linear lift
would give (14), already contradicted. Therefore the two embedded curves
are ambiently inequivalent.

An exact implicit equation for the second image, with u=k,v=z, is

    G(u,v)=4096u⁶−128u³v²+384u³+v⁴−6v²−8v−3.

It is included as a verification control, not needed for the period proof.

### 4.4 Consequence for the next Keller step

This is an explicit rational, immersed, delta-three obstruction to replacing
the distinct-prime embedding input in the quadratic-source proof by an
unqualified <3,4> uniqueness assertion. It also shows why the whole-algebra
identity alone does not license transport of the standard collision period
using the cotangent determinant.

It does **not** exhibit a generating Keller pair, settle how many embedding
classes an arbitrary member has, or imply that each member has two classes.
The failure of a uniqueness-based proof is not a counterexample to a
Keller-sector exclusion.

## 5. Bounded primary-literature pass

The source record `PRIMARY_SOURCES.md` gives the actual retrieval routes.
The mathematical arguments in Sections 1–4 are self-contained apart from
the already pinned one-line theorem in the Keller corollaries.

Abhyankar–Sathaye (1996), Theorems (1.3) and (1.4), separates prime-exponent
uniqueness from composite-exponent nonuniqueness. Its printed p.1063 also
discusses reduced characteristic degree sequences. Sathaye–Stenerson's
Section 3.1 gives the numerical bounds and primitive-generator property.
Gurycz (2001), Section 7, distinguishes equivalence of images from extension
of a specified isomorphism. This is why Section 4.3 is essential.

For clarity, the following is a **numerical consequence of the cited
characteristic-sequence axioms and bounds**, not a new embedding theorem.
A nonprincipal characteristic degree sequence generating <3,4>, ordered so
that δ₀>δ₁, has one of the two patterns

    (4,3),     (6,4,3).

Indeed its conductor is six. The cited bound gives δ₀,δ₁≤7, while every
later entry is at most five; all entries lie in <3,4>, and the primitive
entries 3 and 4 must appear. The possible initial nondividing pairs in
{3,4,6,7} are (4,3),(6,4),(7,3),(7,4),(7,6). Each pair involving 7
already has gcd one and cannot acquire the missing primitive generators.
The pair (6,4) has gcd two and forces a final entry three.

The two exact embeddings in Section 4 realize both patterns for the same
abstract curve. No uniform relative classification over C[x], nor
parameterized all-embedding theorem over C(x), has been imported from
this numerical observation.

## 6. Reproducibility, status, and next precise question

Run from this packet's directory:

    python3 verify_cubic_quartic_boundaries.py
    python3 -O verify_cubic_quartic_boundaries.py

The checker uses explicit exceptions rather than Python assert statements.
It checks the characteristic matrix and relation, the trace-dual conductor
identity, reduced divided-difference resultant, specialized node,
finite degree filtrations, primitive-cubic coefficient module, both
embedding equations, the generation identity, local determinant, period
ratio, normalization symmetries, and the bounded characteristic-sequence
enumeration. These are **FINITE-EXACT controls**; the all-fiber and
globalization proofs are the arguments above.

The next rigorous question is now narrower:

> Does constant Jacobian forbid the quadratic degree-drop fibers identified
> by the primitive cubic/quartic lattices, or supply a replacement for
> global monicity that controls all relevant plane-embedding classes?

Neither saturation alone nor <3,4> embedding uniqueness can answer it.
The generic conductor-six Keller sector remains open in this packet.
