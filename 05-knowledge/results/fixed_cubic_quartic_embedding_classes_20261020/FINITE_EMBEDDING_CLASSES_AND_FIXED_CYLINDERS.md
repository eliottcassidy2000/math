# At most two plane embeddings for fixed <3,4> polynomial curves

20 October 2026. **PROVED; independently reviewed PASS.**

See `review/AUDIT.md` for the independent proof and primary-source audit.
The separately reviewed `EXACT_EMBEDDING_COUNT_APPENDIX.md` determines
when the two possible classes are distinct, without changing this proof.

This is a separate continuation of the cubic/quartic boundary packet.
It does not alter the frozen globalization hostile or inequivalence proof.

## Outcome and exact scope

Let

    A=C[h(y),k(y)] ⊂ C[y],
    h=y³+a y,   k=y⁴+b y²+c y,   a,b,c∈C.        (1)

Then every generating pair of A is carried, by an ambient polynomial
automorphism, to one of at most two explicit pairs:

    (h,k),
    (k,h²+(3c/2)h),                              (2)

where the second pair is included only if it generates A. In fact, it
generates exactly when c≠0 and 2a=3b. This determines candidate generation,
not whether its ambient class is distinct from the first in every case.

In particular A has at most two plane embeddings up to ambient equivalence.
The bound is sharp: A=C[y³,y⁴+y] has exactly two, by the independently
proved inequivalence witness in the
[boundary packet](../cubic_quartic_globalization_embedding_boundaries_20261020/CUBIC_QUARTIC_GLOBALIZATION_AND_EMBEDDING_BOUNDARIES.md).

Consequently the fixed cylinder

    E=C[x]⊗C A ⊂ C[x,y]

has no entire generating pair with nonzero constant (x,y)-Jacobian.
More strongly, a Jacobian lying in C[x] must vanish identically. If the
normalization has a collision, this follows from the finite-class period
lemma. If it has a ramification point, it follows directly from the chain
rule. Since δ=3, at least one of these two cases occurs.

This is a **fixed-algebra** theorem. The coefficients a,b,c in (1) are
complex constants. No uniform classification of moving embeddings over C[x]
or C(x), global cubic/quartic presentation of a general actual source,
or exclusion of the full generic conductor-six Keller sector is asserted.

Every fixed subalgebra of C[y] with degree semigroup <3,4> admits (1):
choose monic elements of degrees three and four. They generate a subalgebra
of the same degree semigroup and codimension three, hence the whole algebra.
Translate y to remove the cubic's quadratic term, remove the quartic's cubic
term using the cubic, and subtract constants.

## 1. The one classical reduction input

We use the classical characteristic-semigroup theorem for a plane curve
with one place at infinity, in the following precise form:

> For a plane presentation with coordinate pole orders d,e, its degree
> semigroup has a characteristic sequence beginning with d,e. After
> elementary triangular degree reductions, either the curve is a line or
> the first two orders are nonprincipal. For a nonprincipal sequence with
> semigroup conductor C, both initial orders are at most C+1.

This is **CITED**, not newly proved. The initial-order identification and
degree-reduction procedure were directly checked in Abhyankar–Sathaye
(1996), printed p.1062, Section 2. The numerical bound was directly
checked in Sathaye–Stenerson, Section 3.1(4), printed p.9 of the
author-hosted 2014 text. `PRIMARY_SOURCES.md` records both sources.
We do not use the general Abhyankar–Singh equivalence classification.

For the concrete algebra (1), gcd(3,4)=1 gives C(h,k)=C(y). The standard
power basis 1,k,k² over C[h] has degrees 0,4,8 modulo three. Therefore

    Γ(A)=<3,4>,   δ(A)=3,   C(Γ)=6.               (3)

Every other plane embedding of Spec A gives a generating pair p,q∈A
after choosing an isomorphism with A. Its normalization is the same
C[y], with one point at infinity. The cited result applies.

If one coordinate degree divides the other, subtract the corresponding
power of the smaller coordinate with the unique leading coefficient
that cancels. This is an ambient triangular automorphism and decreases
the sum of degrees. Iteration terminates. A constant coordinate at
termination would make A singly generated, hence a polynomial ring in
one variable. Since A has normalization C[y] and proper codimension
three, that is impossible. Thus the terminal pair is nonprincipal.

By the cited bound, both degrees lie in

    Γ(A)∩[1,7]={3,4,6,7}.

The only nonprincipal characteristic patterns generating <3,4>, ordered
by their first two entries, are (4,3) and (6,4,3). Here is the short
numerical argument. The nondividing initial pairs are
(4,3),(6,4),(7,3),(7,4),(7,6). A pair involving 7 has gcd one, so the
characteristic sequence cannot have a later entry; it would fail to
contain both primitive semigroup generators 3 and 4. The pair (6,4)
has gcd two and forces the remaining generator three. Hence the
terminal coordinate degrees, up to swapping, are

    (3,4) or (4,6).                              (4)

Only this bounded degree reduction is imported from the classical
semigroup theory. The classification of all pairs with degrees (4,6)
below is a direct calculation with free modules.

## 2. A degree-separated module lemma

Suppose p,q∈C[y], p monic of degree n, and C(p,q)=C(y). Then
C[p,q] is free over C[p] with basis

    1,q,...,q^(n−1).                             (5)

Indeed q is integral over the normal ring C[p], and its minimal polynomial
over C(p) is monic of degree n with coefficients in C[p].

Suppose a C[p]-basis b₀,...,b_(n−1) has degrees d_i in distinct residue
classes modulo n. Then no leading cancellation is possible between
different summands of

    Σ_i A_i(p)b_i.

Its degree is the maximum of the nonzero summands' degrees. Thus the degree
set is exactly the disjoint union of d_i+nN₀. If r_i is the least
nonnegative residue of d_i, then

    dim_C(C[y]/C[p,q]) = Σ_i (d_i−r_i)/n.         (6)

This formula is exact, not an asymptotic bound: it counts the missing
nonnegative degrees in each residue class.

## 3. The reduced degree-(3,4) case

The degree filtration of A has

    A_≤3=C⊕C h,
    A_≤4=C⊕C h⊕C k.

Consequently every pair of degrees (3,4) is obtained from (h,k)
by an invertible affine transformation. This gives one ambient class.

## 4. The reduced degree-(4,6) case

The next relevant filtration step is

    A_≤6=C⊕C h⊕C k⊕C h².

Any generating pair of degrees (4,6) therefore affine-normalizes to

    p=k+βh,     q=h²+γh,     β,γ∈C.             (7)

To see this explicitly, first scale p's k-coefficient to one and remove
its constant. Scale q's h²-coefficient to one, subtract its k-coefficient
times p, and remove its constant. The remaining h-coefficient is γ.
These are ambient affine changes and preserve entire generation.

Assume throughout this section that C[p,q]=A. In particular p is monic
of degree four, the field is C(y), and (5) applies with n=4.

### 4.1 A nonzero β would force delta seven

Put r=q²−p³. Its y^11 coefficient is −3β. If β≠0, then deg r=11.
Replacing q² by r and q³ by qr changes the basis (5) by a unimodular
matrix over C[p]. Therefore

    1,q,r,qr

is also a C[p]-basis. Its degrees are 0,6,11,17, distinct modulo four.
Formula (6) gives

    δ(C[p,q])=(0+6+11+17−(0+1+2+3))/4=7.

This contradicts δ(A)=3. Thus β=0.

### 4.2 The wrong γ would force delta six

Now p=k. Define

    ρ=4a−3b,
    r=q²−p³−ρ p q.                              (8)

The degree-12 terms cancel. There is no degree-11 term; ρpq cancels
the degree-10 term. The degree-nine coefficient is

    [y⁹]r=2γ−3c,                               (9)

because pq has no degree-nine term.

If γ≠3c/2, then deg r=9. Again 1,q,r,qr is obtained unimodularly
from 1,q,q²,q³: explicitly r=q²−ρpq−p³ and
qr=q³−ρpq²−p³q. Its degrees are 0,6,9,15, distinct modulo four.
Formula (6) gives

    δ(C[p,q])=(0+6+9+15−(0+1+2+3))/4=6,

again a contradiction. Therefore γ=3c/2.

This proves that **every** reduced degree-(4,6) generating pair lies
in the single possible class represented by

    (k,h²+(3c/2)h).

No inference from a lone leading degree to an entire semigroup was made:
the displayed free bases prove that there are no hidden later cancellations.
The at-most-two bound is proved without deciding whether the remaining
candidate generates. The following optional refinement decides that issue.

### 4.3 Exact criterion for generation by the second candidate

Put q=h²+(3c/2)h. Then

    C[k,q]=A  if and only if  c≠0 and 2a=3b.    (10)

If c=0, both k and q are even polynomials in y, whereas h is odd and
nonconstant. Thus they cannot generate A.

Suppose c≠0 and 2a≠3b. In (8), additionally subtract ηk², where

    η=−2a(a−b).

The resulting element r₁ has degree seven and leading coefficient
−c(2a−3b)/2. Since k has degree four, gcd(4,7)=1 gives
C(k,r₁)=C(y), hence C(k,q)=C(y). The same unimodular basis change
now gives the C[k]-basis

    1,q,r₁,qr₁

of degrees 0,6,7,13. Formula (6) gives δ(C[k,q])=5, not three.
This excludes generation in the remaining necessary case.

Conversely, if b=2a/3, direct expansion gives

    q²−k³−2a kq+(2a²/3)k²
      +(2a³/27+3c²/4)q−(a⁴/9+a c²)k=(c³/8)h.  (11)

For c≠0 this recovers h from k,q, proving generation. This refinement
does not claim that the second generating pair is inequivalent to the
first for every a,c. That stronger distinction is only needed, and only
asserted here, for the sharp example below.

## 5. Sharpness

For a=b=0,c=1, put h=y³ and k=y⁴+y. If

    q=h²+(3/2)h,

then

    8(q²−k³)+6q=h.                              (12)

Hence C[k,q]=A. The affine reparametrization z=8q+3 gives the second
embedding in the preceding packet, where a direct cotangent/period
argument proves it inequivalent to (h,k). Thus the upper bound two
is attained.

The exact scalar distinction there is:

    required ambient determinant at the triple point = 4,
    collision-period ratio for (k,z) versus (h,k) = −28/5.

Its intrinsic normalization automorphisms lift in the standard embedding,
so the proof excludes unparametrized ambient equivalence, not merely a
specified coordinate map.

## 6. Fixed-cylinder consequence

Let A be as in (1) and suppose distinct r,s∈C have

    h(r)=h(s),     k(r)=k(s).

Set E=C[x]⊗C A. The independently reviewed
[fixed-cylinder finite-period lemma](../cylinder_generator_period_obstructions_20261020/FIXED_CYLINDER_FINITE_PERIOD_LEMMA.md)
applies to finitely many embedding classes, not only to a unique one.
The theorem just proved supplies that finiteness with at most two classes.

For every entire generating pair

    C[x,P,Q]=E,

it follows that

    I(x)=∫_r^s P(x,y)Q_y(x,y)dy

is constant. If J(P,Q)=j(x)∈C[x], differentiation and integration by
parts give

    I'(x)=(s−r)j(x)+[P Q_x]_r^s=(s−r)j(x).

The boundary term is zero because every element of E takes equal values
at r and s. Since r≠s, j=0.

The finite-period mechanism retains the missing embedding label: the
local conormal determinant is a nonzero constant in x, while each
embedding class supplies only finitely many reference period values.
A polynomial in x with finite image is constant. It never equates the
two inequivalent embeddings' standard periods.

The remaining ramified case is elementary. If h'(y₀)=k'(y₀)=0,
every element of E has zero y-derivative along y=y₀ by the chain rule.
Thus P_y(x,y₀)=Q_y(x,y₀)=0 and J(P,Q)(x,y₀)=0. If J=j(x), this
forces j=0 without any period argument.

These cases exhaust the fixed algebras. A proper finite birational
normalization map with no collisions and no ramification would be an
isomorphism, contrary to δ=3. Equivalently a nonnormal complex plane curve
has either multiple branches over a point or a unibranch singularity,
and the latter forces both derivatives to vanish. Thus the fixed-cylinder
Jacobian conclusion holds with no extra immersion or collision hypothesis.
This does not falsely classify the cusp C[y³,y⁴] as having a collision.

## 7. Boundaries and reproducibility

- The algebra must be fixed over C, not merely generically isomorphic
  to a fixed algebra after an unproved rational change of parameters.
- Entire generation over C[x] is required for the cylinder conclusion;
  generic generation alone is insufficient.
- The second candidate in (2) generates exactly when c≠0 and 2a=3b.
  It is not asserted to give a class different from the first in every
  such case. Only an upper bound is asserted generally; sharpness is
  certified for the displayed triple curve.
- The classical one-place semigroup reduction is a cited all-embedding
  input. Finite symbolic checks do not prove it.
- The constant-Jacobian conclusion concerns fixed cylinders: collision
  cylinders use periods, ramified cylinders use the chain rule. It does
  not cover all moving cubic/quartic actual sources.

Run the exact controls from this directory:

    python3 verify_finite_embedding_classes.py
    python3 -O verify_finite_embedding_classes.py

The controls check the coefficient identities, unimodular basis changes,
residue-separated gap counts, the exact second-candidate generation
criterion, the sharp generating identity, the reduced
characteristic patterns, and selected positive and negative generating
examples. The universal proofs are Sections 1–6.
