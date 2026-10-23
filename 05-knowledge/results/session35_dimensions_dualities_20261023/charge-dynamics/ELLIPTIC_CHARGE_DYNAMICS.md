# An exact rational charge orbit with a transcendental rotation angle

Status: PROVED elementary identities and limiting distribution; the angle-transcendence step uses the explicitly CITED Gelfond–Schneider theorem. FINITE-EXACT checks are separate from the all-index proofs. This is a new-to-session connection, not a literature-wide novelty claim.

## 1. Input, indexing and branch

The inherited Coulomb coordinate is C(t)=13−6(t+t^−1). The session's formal intermediate-E7.5 level-one charge is c*=38/5. We use this number as an arithmetic input only: none of the following proves existence, rationality, or a representation category of an intermediate-E7.5 VOA.

Set

    u=(13−c*)/6=9/10,
    t=(9+i sqrt(319))/20,
    θ=arccos(9/20) in (0,π),
    t=e^(iθ),
    S_m=t^m+t^−m,   c_m=13−6S_m,   m≥0.

Thus c_0=1, c_1=c*=38/5. The subscript zero here means exponent zero, not the seed; it avoids an indexing collision with descriptions calling the seed c0. The opposite square-root branch replaces t by t^−1 and θ by −θ modulo 2π, and leaves every c_m unchanged. The minimal polynomial is 10t²−9t+10; |t|=1 exactly. The field is Q(sqrt(−319)).

## 2. Infinite order and transcendence are different statements

If t were a root of unity, both t and t^−1 would be algebraic integers. Their sum u would therefore be a rational algebraic integer, hence an integer. But u=9/10. Consequently t has infinite multiplicative order and θ/π is irrational. In particular t^m is never ±1 for m≥1, so

    1<c_m<25  for every m≥1.

For the stronger assertion, suppose β=θ/π were algebraic. It is already irrational. Use the logarithm log(−1)=iπ, so the corresponding value of (−1)^β is exp(iπβ)=t. Gelfond–Schneider says this value is transcendental, contradicting the displayed quadratic equation. Hence θ/π is transcendental. Every argument branch divided by π, β+2k, and the reciprocal branches −β+2k, remains transcendental.

This does not say that the algebraic number t is transcendental. It also does not deduce transcendence of θ just by multiplying β by π, since products of transcendental numbers need not be transcendental. No such inference is used.

CITED theorem: Gelfond–Schneider, in its complex-power form (every value for algebraic base other than 0 or 1 and algebraic irrational exponent). The inherited original research record is A. O. Gelfond (1934), https://www.mathnet.ru/rus/im4924 . Retrieval of the English record https://www.mathnet.ru/eng/im4924 timed out during this pass; the original paper was not freshly audited. An independently accessible author-hosted proof is Richard Elman, *Lectures on Abstract Algebra*, §74, Theorem 74.29 and Corollary 74.30, printed p.442 / PDF p.450: https://www.math.ucla.edu/~rse/algebra_book.pdf . Its statement and arbitrary-logarithm proof were checked. This standard theorem, not computation, supplies transcendence.

## 3. Exact arithmetic and a growing-denominator certificate

The reciprocal trace recurrence gives

    S_0=2, S_1=9/10,
    S_(m+1)=(9/10)S_m−S_(m−1),
    c_(m+1)=(9/10)c_m−c_(m−1)+143/10.

In particular every c_m is rational. Define A_m=10^m S_m. Then

    A_0=2, A_1=9,
    A_(m+1)=9A_m−100A_(m−1),
    c_m=(13·10^m−6A_m)/10^m.

For m≥1, A_m is odd and not divisible by 5, by induction modulo 2 and 5. For m≥2 the numerator of c_m is 2 modulo 4 and nonzero modulo 5; for m=1 the fraction is 38/5 directly. Thus the reduced denominator of c_m is exactly

    10^m/2,  m≥1.

These unbounded, strictly increasing denominators prove that the charges are pairwise distinct. This is an exact arithmetic certificate, not a numerical test of irrationality.

The inherited powering identity gives another exact invariant:

    Δ(c_m)=Δ(c*) U_(m−1)(9/20)^2,
    Δ(c*)=(c*−1)(c*−25)=−9·319/25,
    Δ(c)=(c−1)(c−25).

The U-factor is rational and nonzero, since its vanishing would imply t^(2m)=1. Every positive power therefore retains discriminant square class −319 and quadratic field Q(sqrt(−319)). This extends the inherited c<1 square-class preservation to this infinite-order elliptic seed. It does not extend to finite-order seeds: c=13 has t=i and its second power hits branch value 25.

## 4. All-integer equidistribution, with a direct proof

For every nonzero integer k, t^k≠1, so the geometric sum yields

    (1/N) sum_(m=1)^N exp(ikmθ)
      = t^k(1−t^(kN))/(N(1−t^k)) → 0.

For k=0 the average is 1. Therefore averaging any trigonometric polynomial along mθ gives its circle integral. Uniform approximation of continuous periodic functions by trigonometric polynomials extends this to every continuous function on the circle. This proves the required irrational-rotation equidistribution directly; no quantitative rate uniform in k is claimed.

Push normalized circle measure forward through c=13−12cosφ. On (0,π) the derivative is 12sinφ=sqrt((c−1)(25−c)); the other half-circle gives the same contribution. Hence for every continuous f on [1,25],

    lim_(N→∞) (1/N) sum_(m=1)^N f(c_m)
      = ∫_1^25 f(c) dc / [π sqrt((c−1)(25−c))].

For 1≤a<b≤25 the limiting proportion in [a,b] is

    [arccos((13−b)/12)−arccos((13−a)/12)]/π.

Interval convergence follows by continuous upper/lower approximations because the limiting measure has no atoms, including at the integrable endpoint singularities. Every nonempty interior interval has positive measure, proving density in [1,25] and infinitely many visits to each such interval. Actual positive-index values never equal the endpoints.

Consequences include exact limiting mean 13, variance 72, and centered moments

    mean((c_m−13)^(2k+1)) → 0,
    mean((c_m−13)^(2k)) → 36^k binomial(2k,k).

These follow by expanding (t^m+t^−m)^j and keeping the constant Fourier coefficient. For any fixed lag r≥0 the same argument gives

    lim (1/N) sum_(m=1)^N (c_m−13)(c_(m+r)−13)=36S_r.

Thus the marginal arcsine law does not make the sequence independent or random. The lag-one correlation coefficient is S_1/2=9/20 exactly.

## 5. Why doubling is not the theorem

The inherited quadratic polynomial P_2(c)=25−(13−c)^2/6 sends c_m to c_(2m). Iteration starting at c* samples only m=1,2,4,8,.... The geometric-sum proof above applies to all consecutive integers m, not this lacunary subsequence. Irrationality or even transcendence of θ/π alone does not prove density or equidistribution of the specific doubling orbit. No such claim is made.

The arcsine measure is also invariant under P_2, because circle doubling preserves normalized circle measure. This is a statement about a measure, not about genericity of this particular initial point. It must not be silently substituted for a pointwise orbit theorem.

Nor do the real positive mirror formulas from the inherited c<1 branch apply here: the present reciprocal parameters lie on the unit circle, not on t>1. Complex exponentiation requires its own branch choices and supplies no positive mirror-base interpretation in this note.

## 6. Recovered connections, scope and reproduction

Load-bearing inherited notes: `session35-lie/COULOMB_MIRROR_DUALITY.md` (quotient and real locus) and `session35-lie/COULOMB_POWERING_ADDENDUM.md` (Chebyshev powering and discriminant identity). The same S-recurrence has characteristic trace 3 in the golden-ratio/path-forest example, but trace 9/10 here. The former is hyperbolic and escapes; this one is elliptic and has an arcsine-distributed rational charge sequence. This is a concrete change of dynamical regime of the same quotient, not an assertion that forest counts are the present rational charges.

A bounded prior-repository search recovered the two-step-walk arcsine occurrence in `05-knowledge/results/unit_distance_polygon_ladder_20261017.md`, lines 84–87. It concerns a radial random-walk convolution, not this deterministic orbit. The only transported fact here is the explicit cosine pushforward geometry; no elliptic-curve or graph-count conclusion is imported. A search hit in THM-3012 was a *lemniscatic* arcsine and is not the probability law above.

Run `python3 session35-charge-dynamics/checks.py`. The JSON receipt records exact arithmetic checks, quadratic-field comparisons, Chebyshev/discriminant checks and hostile finite-order controls. The proofs of all-index denominator growth and equidistribution stand independently of that finite run. No VOA-realization claim, no dense-doubling claim, no numerical transcendence inference, and no Git operation is part of this deliverable.
