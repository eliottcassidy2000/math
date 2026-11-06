# Source context: Soule's two-section height strategy

## Pinned source and status

Christophe Soule, *On the abc conjecture*, arXiv:2611.05217v1,
3 November 2026, 39 PDF pages.

- Public source: https://arxiv.org/abs/2611.05217v1
- SHA256 of the supplied PDF:
  f85bf570334601930ae862615c4f13365d2cd45370bb0b853252bf15445d9abf.
- Source status here: **PREPRINT CLAIMS, NOT INDEPENDENTLY AUDITED**.

The PDF's first page states a genus-two effective Mordell-type bound
with discriminant coefficient 12+epsilon, and an abc-type upper bound
c<=A(epsilon) rad(abc)^(24+epsilon) for coprime positive a+b=c.
This is not the usual 1+epsilon abc exponent. Section 7 also claims a
uniform number-field statement and a consequence excluding Siegel zeros
for negative quadratic discriminants. None is a premise of our source
theorems or promoted to repository canon here.

## The method and its load-bearing hypotheses

The introduction, pp. 1-3, describes the Parshin/Moret-Bailly strategy:
pass from a genus-two curve to an etale cover carrying a non-isotrivial
smooth proper genus-six family. A relative dualizing self-intersection
becomes a height on the base.

1. Section 1.7, Proposition 2, and Section 1.8, Corollary 1, compare
   arithmetic heights using the ordinary degree of line bundles on a
   projective curve. Degree-zero error terms are bounded by epsilon times
   an ample reference height plus O([F:Q]); the constant depends on epsilon.
2. Section 1.9, equation (2), invokes Minkowski successive minima for a
   rank-two adelic vector space with an integral model. It chooses TWO
   linearly independent sections with a bound on the SUM of their heights.
3. Section 2.2, Proposition 3, assumes ample line bundles with admissible
   metrics and sections whose divisors do not meet. These hypotheses
   permit an arithmetic-intersection bound from the two section heights.
4. Proposition 9, pp. 22-23, specifies an exceptional set Phi and the
   relevant dimensions and degree calculations. In Section 6.3.2 the
   section space has dimension two away from Phi. Section 6.3.3 separately
   invokes basepoint freeness to turn independence into no common zero.
   Independence alone is not substituted for this condition.
5. Sections 6.3.6-6.3.7 compare two estimates on the SAME fiber quantity:

       -L_P^2 <= (3/2+epsilon) omega_P^2+6 delta(F)+O([F:Q]),
       -L_P^2 >= (5/3-epsilon) omega_P^2-O([F:Q]).

   The strict coefficient gap 5/3-3/2=1/6 is the formal source of the
   discriminant coefficient 36. This algebraic gap only has content if
   all preceding geometric, metric, and base-change inputs hold uniformly.
6. Section 6.3.8 handles the finite exceptional fibers separately by a
   bounded normalized height. Finite exception handling is not a license
   to ignore an infinite degeneration locus in another problem.

## The exact transfer made here

Source: the need to retain two independent sections and distinguish
exceptional degeneration before using a two-measurement estimate.

Target: the joint physical source pencil with columns (w,q), transported
by one literal source word. These are two response boundaries, not two
line bundles, two lattice vectors chosen by Minkowski, or two independent
physical sources. The second boundary q is generally nonphysical.

Map: evaluate the two boundaries on four fixed physical sources. The
result is multiplication by the invertible matrix P Q, determinant four.
Therefore rank is preserved exactly. The simultaneous G/H measurements
and common-word certificate survive. No modular reduction, gcd-source
normalization, or independent choice of the two boundaries is made.

New consequences proved without the preprint's major claims:

- Rank-one degeneration forces an explicit nontrivial integer factor of
  every target in the synchronized compiler, not just repeated-return
  targets.
- Rank two is detected by three exact determinants from four tiny
  physical probes, with a quantitative separation bound.
- Rank-one branches may be removed from the exact prime-target inverse
  without losing witnesses or minimum root-change costs.

Destroyed information: these ranks and determinants do not establish
integer target occupancy, physicality of a primitive boundary quotient,
or existence for every prime. The full synchronized columns and a
predecessor recipe remain necessary for construction.

Decisive hostiles: a non-return rank-one state (OJOOJ); a rank-two state
with composite target 51; equal physical w with different auxiliary q;
and the inherited failure of primitive source normalization. The finite
prime comparison in verify.py is a complete numerical-height test, not
an extrapolation of all-order prime coverage.

The preprint's strict height-gap mechanism suggests a separate future
question: find two genuinely independent physical estimates for the same
target fiber, with a positive uniform gap. This note supplies the exact
independence guard but does not manufacture the missing estimates.
