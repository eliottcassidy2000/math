# A marked cubic closes the Keller Fourier remainder

2026-10-17, clock_duality_oct15. **PROVED, scoped and self-audited:** the
identities and all-prime bounds below for the explicit three-dimensional
map. **FINITE-EXACT:** independent character-coefficient verification at
31,412 nonzero prime/prime-power frequencies. No independent referee or novelty claim.
**OPEN:** the planar Jacobian conjecture and any transport to Collatz or
six-colour forcing.

## Inheritance and the decisive retained coordinate

The incoming pre-audit checkpoint c3fc4332ed contains
[THM-4664, Fourier concentration on the dual Jelonek cubic](../../01-canon/theorems/THM-4664-fourier-transform-of-the-thm1300-keller-map-mod-p-concentrates-on-the-dual-of-its-asymptotic-jelonek-cubic.md).
Its exact decomposition is the closest mechanism; its unproved uniform
remainder estimates are the target. The hostile is its own warning that
counting roots while forgetting phases gives the wrong fibre statistic.
The underused coordinate is **the coefficient of the reciprocal source
variable**, rather than a second-moment or genus label.

Write u=1+xy and retain the literal map

\[
 F=(u^3z+y^2u(4+3xy),\quad
 y+3xu^2z+3xy^2(4+3xy),\quad 2x-3x^2y-x^3z).
\]

Its Jacobian is -2. For an odd prime p, put
S(a,b,c)=sum_(x,y,z in F_p) psi(aF_1+bF_2+cF_3), where
psi(t)=exp(2pi i t/p). Roots below are **distinct field elements**, never
counted with algebraic multiplicity.

## K1. Exact Fourier formula, all odd primes

For a!=0, set R(t)=a t^3+3b t^2-c. Then

\[
 \boxed{S(a,b,c)=p^2\psi(-4b^2/a)[a^2c=4b^3]
       -p\sum_{R(t)=0}\psi(a t^2+4bt).}                 \tag{1}
\]

For a=0,b!=0 and for a=b=0,c!=0 respectively,

\[
 S(0,b,c)=-p\sum_{3bt^2=c}\psi(4bt),\qquad
 S(0,0,c)=p^2.                                        \tag{2}
\]

**Direct proof.** For x!=0 use the bijective change t=y+1/x, retaining
x and z. The coefficient of z in the phase is x^3 R(t).
Sum over z first. On R(t)=0, the remaining phase simplifies exactly to

\[
             a t^2+4bt+(at+2b)/x.                     \tag{3}
\]

Since sum_(x!=0) psi(d/x)=p[d=0]-1, this part of S is
p sum_(R(t)=0) psi(at^2+4bt)(p[at+2b=0]-1).
For a!=0 the distinguished root t=-2b/a exists exactly when
a^2c=4b^3; its phase is -4b^2/a. The plane x=0 contributes zero
unless a=b=0, when it contributes p^2. This proves all cases, including
characteristic three. At the zero frequency the original sum is p^3.
No asymptotic equidistribution theorem is used. The same argument proves
(1)--(2) over every finite field of odd cardinality q, with p replaced by q
and any nontrivial additive character. This is an orthogonality identity,
not a prime-field-only exponential parametrization.

There is also a direct reconciliation with THM-4664's remainders. For a!=0,
with G(a)=sum_y psi(ay^2), they are

\[
 E=p\sum_{R(t)=0}\psi(at^2+4bt)-\psi(-4b^2/a)G(a),
 \qquad S_\Gamma=\psi(-4b^2/a)(p-G(a)).                \tag{4}
\]

For E substitute w=ct in that theorem's sum (its summation variable c
is nonzero). For S_Gamma use r=1/x and
F_1=y^2+3ry+2r^2, F_2=4y+6r, then complete the square.
The two Gauss terms cancel in -E+S_Gamma. In particular the stated uniform
remainder obligation is discharged, with |E|<=3p+sqrt(p) and
|S_Gamma|<=p+sqrt(p). Formula (3) is a shorter proof of the final answer.

## K2. The dual cone is exactly the large spectrum

For every prime p>=7,

\[
 \{|S(a,b,c)|>p^2/2:(a,b,c)\ne0\}
   =\{a^2c=4b^3\}\setminus\{0\}.                     \tag{5}
\]

Off the cone (1)--(2) give |S|<=3p. On the cone:

* a!=0,b!=0: R(t)=a(t-b/a)(t+2b/a)^2, so
  S=p psi(-4b^2/a)((p-1)-psi(9b^2/a)), of magnitude at least p(p-2).
* a!=0,b=0,c=0: S=p(p-1).
* a=b=0,c!=0: S=p^2.

These inequalities prove (5). There are exactly p^2-1 nonzero cone
frequencies. The geometry of the cone as the projective dual of the
asymptotic Jelonek cubic is inherited from THM-4664, not re-claimed here.
The same statement holds for q>=7 of characteristic at least5. In
characteristic3 the two displayed roots coalesce; that case is covered
by (1), but is not used in the factorization argument just given.

The second component of the cubic discriminant is a useful hostile.
For p>3, disc_t R=27c(4b^3-a^2c). At c=0 with ab!=0 the root t=0
is double, yet (1) gives

\[
                  S(a,b,0)=-p(1+\psi(-3b^2/a)).        \tag{6}
\]

Thus a double root alone does not explain the large coefficient. The
distinguished condition at+2b=0 is essential. This is a precise instance
of preserving a marked phase through a quotient.

## K3. Exact energy concentration, every prime p>=5

Let E_all=sum_(frequency!=0)|S|^2 and E_cone be its sum on the cone.
Then

\[
 \begin{aligned}
 E_{all}&=p^3(p-1)(p^2+2),\\
 E_{cone}&=p^2(p-1)(p^3-p^2+4p-3),\\
 \frac{E_{cone}}{E_{all}}
   &=1-\frac{p^2-2p+3}{p(p^2+2)}.                    \tag{7}
 \end{aligned}
\]

For the cone energy, square the three explicit cases above and use
sum_(a,b!=0) psi(9b^2/a)=-(p-1), valid for p>3.

For completeness the total energy has an elementary independent count.
Use target coordinates (A,B,C). When C!=0, the fibre recipe inherited
from THM-4664 identifies preimages with the **simple** roots of
w^3-2w^2+BCw-2AC^2; the deleted locus is precisely its derivative-zero
roots. Three preimages therefore correspond to unordered triples of
distinct roots summing to2. There are (p-1)(p-2)/6 such triples, and
p-1 choices of C. When C=0, x=0 gives one preimage for every (A,B),
and x!=0 gives two more exactly when B^2-16A is a nonzero square:
p(p-1)/2 targets. Hence the number of three-point fibres is

\[
 N_3=(p-1)^2(p-2)/6+p(p-1)/2=(p-1)(p^2+2)/6.
\]

The only possible fibre sizes are0,1,3, so the number of ordered
distinct-source collisions is6N_3. Parseval gives E_all=6p^3 N_3.
This proves (7), rather than extrapolating its finite numerical pattern.
Every step also proves (7) for prime powers q of characteristic at least5,
with p replaced by q.

## Exact checks, controls, and non-transfer

[Checker](../../04-computation/experiments/keller_cubic_fourier_exact_20261017d.py),
[receipt](../../04-computation/experiments/escape_renewal_20261017d/keller.json).
The explicit universe is every nonzero frequency over p=3,5,7,11,13,17,19:
15,788 cases. Direct evaluation of F produces an integer vector of phase
counts. The proposed expression is a second integer vector. Their
difference must be constant across all p entries, equivalently zero
modulo the prime cyclotomic polynomial. No floating-point FFT is used.
Symbolic checks verify the Jacobian, both substitutions, and (3).
Direct fibre-count Parseval agrees with (7). The characteristic-three
control is included for (1)--(2), outside the scope of (7). An additional
15,624 nonzero frequencies over F25=F5[t]/(t^2-2) are checked using the
trace character; this independently tests the prime-power extension.

The reusable move is **sum the linear variable before bounding the
remaining sum; preserve the coefficient that can vanish**. Its source is
the explicit polynomial phase, target the marked cubic root set, preserved
predicate exact Fourier value, and lost information the x-coordinate
whose orthogonality is retained in the factor p[at+2b=0]-1.
The hostile (6) distinguishes this from unmarked singularity counting.

For Collatz, the source-dependent legal parity guard and integer clock
must remain in any analogous summation; a finite-field affine word does
not provide such a guard. For six-colour forcing, the actual boundary
extension relation must survive; a polygon moment or point count does
not supply it. This result closes an explicit three-dimensional Fourier
calculation, not the planar Jacobian conjecture or either universal target.
