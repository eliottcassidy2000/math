# Independent Coulomb/mirror audit

Status: PASS, major errors0. Scope: COULOMB_MIRROR_DUALITY.md and coulomb_checks.py/json. Powering addenda are outside this frozen report until explicitly added separately. No Git mutation.

## Quotient and arithmetic

The quadratic inverse has discriminant(c−1)(c−25), with reciprocal product1. Algebraicity equivalence holds in both directions; for rational c the inverse roots are rational exactly when the discriminant is a rational square, otherwise quadratic irrational. The field extension bound is valid even for transcendental c, interpreted as Q(c)⊂Q(t).

The imaginary-part formula verifies the complete real locus. Positive real t gives c≤1, negative real t gives c≥25, and the unit circle gives [1,25]. Overlap only occurs at t=±1. Both branch values and the unramified poles0,∞ are correct. The Laurent invariant-ring proof by symmetric-power recurrence is valid.

## Mirror domain and transcendence

Recovered the Session31 mirror README by Git object read. A distinct positive mirror pair necessarily has1<x<e<y: log(x)/x is injective for x≤e and negative below1, then strictly decreasing above e. Setting t=y/x>1 gives the stated formulas. Conversely the elementary logarithm inequalities give exactly these inequalities and the mirror identity for all t>1. Reciprocal t swaps the two bases. Charge decreases strictly from1 to−∞, so the unordered-pair bijection is exact. The entire diagonal x=y is not parametrized; only its distinguished coalescence(e,e) is the limiting branch point.

For rational c<1 and nonsquare discriminant, t is positive algebraic irrational. The exponent1/(t−1) is algebraic irrational (if it were rational, t would be rational); t is neither0 nor1. Gelfond–Schneider applies to the positive real power and makes X transcendental. Since t is nonzero algebraic, Y=tX must also be transcendental. In the square case both are rational powers of positive rational numbers and algebraic. The exception at c1 is correctly excluded, not patched by substitution into a singular exponent.

The prior integer-mirror theorem was independently checked: if n≥2 is integer and r=a/b rational in lowest terms, n^a b^(nb)=a^(nb) forces b1. The smaller member of a distinct integer mirror pair must be2, and the unique partner is4. Any other algebraic r would be irrational, making n^r transcendental while r^n is algebraic, contradiction. Thus the present warning about transcendental t is justified by inherited examples, not speculation.

## Exact VOA hypotheses and sources

Verified Huang's primary paper https://arxiv.org/pdf/q-alg/9505020 , Theorem2.1, printed pp15–16: simple L(c,0) is rational precisely for coprime minimal-model parameters p,q>1, including c0. Its next paragraph explicitly states rationality of finite tensor products of these minimal models. Therefore the note's criterion on the reduced numerator and denominator of positive rational t is exact, including exclusion of integer ratios and t1. It is not a criterion for every VOA of central charge c.

The large tensor controls are correct: five Lee–Yang factors and34 Ising factors total−5;114 Ising factors total57. An independently verified smaller control is L(c_(2,5),0)⊗L(c_(3,5),0), whose charges−22/5 and−3/5 total−5. Regularity also follows from Dong–Li–Mason https://arxiv.org/pdf/q-alg/9508018 , Theorem3.13 and Proposition3.3. Thus there is no confusion between tensor-product rationality and simple Virasoro rationality.

Integer minimal-model charge forces pq|6 because gcd(pq,(p−q)^2)=1, hence only p,q2,3 and c0. The Ising fusion relations give the stated positive dimensions1 and sqrt2. The positivity/unitary qualifier matters and is present.

The original Gelfond article link was not re-fetched here. Instead verified a complete author-written proof source: Richard Elman, https://www.math.ucla.edu/~rse/algebra_book.pdf , Theorem74.29 (printed p442; PDF page450). It states the exact algebraic-base/nonrational-algebraic-exponent hypotheses used. Wang's original publisher PDF could not be retrieved; Huang's explicit primary-paper restatement supplies the needed iff without that retrieval.

## Independent finite checks

The independent executable traverses rational central charges rather than rational t, checks inverse roots and reciprocal identities, establishes the integer discriminant list through factorization144, and verifies rational mirror and both tensor controls using exact fractions. Analytic and transcendence conclusions are proof-level checks, not claims inferred from floating-point tests.

No Langlands equivalence, categorical functor, arbitrary-transcendental-parameter conclusion, or major open-problem closure is justified or claimed.
