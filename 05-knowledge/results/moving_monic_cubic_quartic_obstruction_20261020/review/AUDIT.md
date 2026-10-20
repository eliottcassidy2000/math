# Independent audit: entire moving cubic/quartic generator obstruction

20 October 2026.

## Verdict and precise scope

**PASS within the stated entire, globally monic hypotheses.** For

    R=C[x], h=y³+a(x)y, k=y⁴+b(x)y²+c(x)y, a,b,c∈R,

no pair P,Q generating R[h,k] as an entire R-algebra has nonzero constant
(x,y)-Jacobian. The normalization paragraph correctly extends this to
any global monic degree-three/degree-four presentation by a polynomial
source translation of Jacobian one and changes of the chosen generators.

The proof does not assume that an arbitrary actual generic conductor-six
source has these global generators. It does not exclude the full
conductor-six sector, prove a new universal conductor lower bound, or
settle the planar Jacobian conjecture. It also does not prove the stronger
moving-family assertion with an arbitrary nonconstant Jacobian in C[x].
The constant, nonzero Jacobian assumption is used to force uniform immersion.

No mathematical correction to the audited manuscript was required.
During the outline audit the reviewer supplied explicit global conductor
membership identities and the transverse-node identity, now included in
the producer proof and checked independently below. The second-class
primitive sign and its factor a′ were checked before the written proof
was finalized. Parent-owned status promotion, verified dependency links,
and the documented extra hostile control do not alter the argument.
One summary-only edit retained the phrase “under that dichotomy's
hypotheses” in the README's actual-source consequence. Removing the
documented status/link additions and extra hostile paragraph from the
final proof reproduces the exact prepromotion proof hash; the analogous
checker comparison also matches. See promotion_check.json.

All reviewer writes are confined to this review directory. No mathematical
repository, producer source, or external service was modified by the
reviewer. Publication belongs to the parent.

## 1. Flat degree filtration and exact conductor specialization

The intermediate-field degree divides both three and four, hence is one.
Therefore k has degree three over C(x,h). Since k is integral over the
normal domain R[h], its monic minimal polynomial has coefficients there.
The power basis 1,k,k² is valid over R[h], not merely over its fraction
field. Its distinct degree residues prevent cancellation after every base
specialization as well as generically. This proves the claimed injection
of the actual fiber algebra and delta three on every fiber.

The conductor calculation is not a generic formula incorrectly extended
across a bad fiber. The monic implicit equation specializes to the actual
minimal equation on each fiber, since the degrees remain three and four.
The inherited trace-dual identity consequently computes the exact
conductor on each fiber. Independently constructing the implicit
resultant and the collision resultant gives the same monic polynomial

    C=y⁶+(a+b)y⁴+2cy³+(a²−ab+b²)y²
        +(a+b)cy+a(a−b)²+c².

The direct identities for C, Cy, and Cy² show C R[y]⊂E over R itself.
The normalization is free over R[h] on 1,y,y², so those three identities
suffice for the whole containment. Conversely, an element in the global
conductor is divisible by C generically. Polynomial division by the monic
C takes place over R, and its remainder vanishes generically, hence
identically. Thus the global conductor equality is also justified.

This is sufficient to identify every root of every specialized C with a
normalization point above a singular point. Indeed the curve is reduced,
its normalization is C[y], and regularity and normality coincide for
one-dimensional local rings of complex algebraic curves.

## 2. Immersion, exceptional loci, and independent node constraints

The chain rule makes every element of E have zero y-derivative wherever
h_y and k_y vanish together. A nonzero constant Jacobian forbids such a
point. Both derivative leading coefficients are nonzero constants, so
the resultant has no root in the affine x-line and is a polynomial unit.
The exact normalization is Res(h_y,k_y)=9I, not I.

Both exceptional regimes are completely handled.

- If b=a, the unit equation makes a and c coprime. Differentiation gives
  c|a′ and a²|c′. If both are nonconstant, their two degree inequalities
  contradict one another. If either is constant, the unit equation makes
  the other constant too.
- A repeated root of H is rational unless H=s³, where the root is zero.
  The rational root is integral over R because H is monic, so it is a
  polynomial. The displayed factorization of I into
  (a+3t²)(4a+3t²)² times 4/9 makes both factors constant units; it then
  makes a,t,b,c constant. No polynomial Pell or unproved elliptic
  assertion is being used.

The inherited fixed-cylinder result genuinely covers all fixed
coefficient cases. It distinguishes collisions from ramification and
retains the possible two embedding classes. In the present proof
uniform immersion already excludes the ramified fixed models, but
invoking that broader reviewed theorem is harmless.

Outside the exceptional regimes, the separation resultant is −9I and
the common h-values are −c−(b−a)s. Thus the three roots of the squarefree
H give three different image points and six different normalization
points. The independently verified tangent identity is

    h′(r)k′(t)−h′(t)k′(r)=(t−r)(3s²+4a)H′(s).

Every factor is nonzero. More importantly for the subsequent argument,
the three equal-value conditions are independent by polynomial
interpolation at six points. Their equalizer has codimension three,
the same as the actual algebra. This proves equality of the two spaces
and provides the sufficient algebra-membership criterion the period
argument needs. It is not merely a necessary gluing condition.

## 3. The classification really descends to C(x)

All triangular cancellations are performed over the given coefficient
field: the needed leading coefficient ratio belongs to that field.
The terminal degrees do not change after extending constants. One may
therefore use the classical characteristic-sequence theorem over an
algebraic closure solely to bound and enumerate those terminal degrees.
There is no descent of an arbitrarily chosen algebraic-closure
automorphism.

The inherited degree-separated free-module calculation was reread and
its decisive coefficients were independently recomputed. The only
terminal degrees are (3,4) and (4,6). The latter affine-normalizes over
the original field to (k+βh,h²+γh). Nonzero β gives delta seven; after
β=0, an incorrect γ gives delta six. With γ=3c/2, the off-locus
candidate gives delta five. The residue-separated bases establish these
actual codimensions, not just the existence of a leading cancellation.
When c=0 the candidate is even and cannot generate the odd cubic.
The universal recovery identity proves generation on c≠0,2a=3b.

Reversing the explicit elementary and affine operations supplies an
exactly parametrized ambient automorphism over C(x). No intrinsic
normalization automorphism or continuity choice is missing. Uniqueness
of the chosen class is unnecessary.

For the second reference, the immersion unit gives I=3c² and hence
c∈C*. The recovery identity then has a constant invertible denominator,
so that reference generates over the entire R, not only generically.
This point is essential to the special-fiber cotangent argument.

## 4. Scalar descent across all special fibers

For global representatives of P,Q in the selected reference coordinates,
the difference from the generic ambient automorphism is a multiple of
the generic implicit equation. At a singular point the equation and its
first partial derivatives vanish. Consequently the pulled determinant
D equals the generic constant-in-(U,V) determinant λ at each of the six
generic conductor roots.

Generic squarefreeness implies C|(D−λ) over C(x)[y]. Monic division over
R[y] proves at once that λ∈R and that this same congruence holds globally.
This is a valid argument about remainders; it does not assume a global
ambient lift or discard denominator poles.

At every complex fiber, whole generation gives both pairs as generators
of the same reduced plane-curve algebra. At a singular point its
cotangent space has dimension two. Each pair supplies a basis, so its
transition determinant D is nonzero. The specialized conductor has
degree six and detects a singular point on every fiber, including fibers
where generic nodes coalesce. The global congruence thus gives
λ(ξ)≠0 for every ξ∈C and proves λ∈C*.

The generic-only control (xh,k) has determinant x and loses generation
at x=0. It confirms why generic generation cannot replace this argument.

## 5. Moving endpoints, exact primitives, and both gap obstructions

For a polynomial ambient automorphism Φ of determinant λ, the form
Φ₁dΦ₂−λU dV is polynomial-exact in the reference variables. Evaluation
of its primitive at equal reference endpoints therefore gives the
claimed period identity even when Φ depends rationally on x.

Differentiating a moving-endpoint period does not discard boundary terms.
After integration by parts they are the common value of A multiplied by
the difference of the total x-derivatives of the two equal B-values.
They vanish exactly. A local analytic base neighborhood with regular
roots justifies the calculation; it is also an identity in a finite
algebraic extension with the extended derivation. Since λ is constant,
the required primitive

    j y−λ∫J(reference)dy

takes equal values at all three pairs. The codimension argument proves
membership after scalar extension; the degree-separated basis, or
faithful extension of vector spaces, descends membership to C(x).

For the standard reference, the degree-five coefficient vanishes first.
The degree-two remainder is exactly

    −λ((3b−2a)c′+2ca′)y²/4.

Together these give b=4a/3+B and (2a+3B)c constant. If that constant
is nonzero, both polynomial factors are units. If it is zero, the two
possibilities are 2a+3B=0 or c=0; in each case the unit equation settles
the remaining coefficients. There is no omitted zero-factor case.
Once the family is fixed the primitive is jy, contradicting the gap one.

For the second reference, writing T=−J(k,q)/a′ formally, the integral of
T has remainder +c²y²/4. Accordingly the integral of J(k,q) has remainder
−a′c²y²/4, and the required primitive has remainder

    j y+λa′c²y²/4.

This checks the signs and factors in the manuscript's V convention.
The independent degree-one and degree-two gaps force j=0. The proof
does not equate the periods of the two inequivalent reference embeddings.

## 6. Actual-source consequence and excluded overclaims

The preceding globalization packet's Sections 3.1–3.2 were checked in
full. For an R-saturated algebra with generic semigroup <3,4>, excluding
all special-fiber elements of degree one or two first makes the primitive
cubic monic, then makes a quartic monic. The resulting subalgebra is
saturated and generically equal, hence entirely equal to the source.
The present theorem eliminates that branch.

In the actual noninvertible Keller setting covered by the cited
hypotheses, a linear element in a fiber makes that source line injective
and is ruled out by the inherited one-line theorem. The remaining
consequence is therefore existence of a genuine quadratic special fiber.
It is not absence of all special fibers or an exclusion of generic
conductor degree six. The final manuscript preserves this distinction.

## 7. Primary-source check and exact reproducibility

The author-hosted Sathaye–Stenerson PDF was independently opened anew at
https://www.ms.uky.edu/~sohum/ma561/affgeom/papnew.pdf and the inherited
local copy was read. Section 2.2 fixes the ground-field convention;
Section 3.1(4) gives the conductor-plus-one bound, with the neighboring
gcd and primitive-generator properties used in the finite enumeration.
The passage to an algebraic closure in the present proof satisfies that
convention. The local extraction is a reading aid, not a publication file.

The Abhyankar–Sathaye initial-coordinate and reduction input is retained
from the independently checked fixed packet. A new DOI request returned
an internal retrieval error; indexed original text from the cited paper
was available, but no fresh successful full-PDF retrieval is claimed.
The reviewer does not import the broader general equivalence theorem
mentioned there.

The independent checker imports no producer code and uses explicit
exceptions, so its checks survive Python optimization. Its 64 checks
include two independent resultant constructions of the conductor,
global conductor membership, the transverse-node identity, both terminal
degree exclusions, both gap reductions, exactness of a nontrivial
ambient-automorphism form, and moving/fixed/ramified hostile controls.
Both independent normal and optimized runs have identical outputs.
The producer's normal and optimized scripts were separately replayed;
all 64 producer checks passed and both outputs match the frozen receipts.
All seven pinned dependency files independently match both their recorded
hashes and the published repository HEAD versions. These comparisons are
recorded in REVIEW_SUMMARY.json and dependency_verification.json.

Reproduce from the packet directory:

    python3 review/independent_checks.py
    python3 -O review/independent_checks.py
    python3 verify_moving_cubic.py
    python3 -O verify_moving_cubic.py

The all-pair theorem rests on the audited proof and the precisely scoped
classical input, not on the finite controls. Source hashes and final
replay receipts are pinned in the adjacent JSON records.
