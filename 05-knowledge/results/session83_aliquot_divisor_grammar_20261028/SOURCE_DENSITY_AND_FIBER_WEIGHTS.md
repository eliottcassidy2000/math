# The three classes under two different sampling rules

Status: PROVED elementary identities and density statement. These concern source integers, not the density of untouchable targets.

Let I be the prime powers p^a (a>=1), II the products pq of distinct primes, and III all other integers n>=2. Put1 in a separate boundary class. This is the priority partition from the clarified three-proper-divisor rule.

## 1. Exact source generating functions

For real u>1 let P(u)=sum_p p^(-u). Absolute convergence justifies

\[
D_I(u)=\sum_{a\ge1}P(au),\qquad
D_{II}(u)=\tfrac12(P(u)^2-P(2u)),\qquad
D_{III}(u)=\zeta(u)-1-D_I(u)-D_{II}(u).
\]

The factor1/2 forgets the order of the two distinct primes; subtracting P(2u) excludes their diagonal. Prime squares remain in I. This is an exact decomposition of Dirichlet series, not a claim of multiplicativity of any individual class indicator.

## 2. Class III has natural density one

Here is an elementary proof that avoids using a numerical census as asymptotic evidence. Fix the first K primes and let X_K(n) count how many divide n. Over their common period M=product p, the divisibility indicators are independent Bernoulli random variables, with

\[
\mu_K=\sum_{j\le K}1/p_j,\quad
\operatorname{Var}(X_K)=\sum_{j\le K}(1/p_j)(1-1/p_j)\le\mu_K.
\]

The Chinese remainder theorem proves the asserted finite-period distribution. Every integer having at most two distinct prime factors satisfies X_K<=2. For mu_K>2, Chebyshev gives its upper natural density at most mu_K/(mu_K-2)^2. The reciprocal-prime sum diverges, so this bound tends to zero as K tends to infinity. Thus integers with at most two distinct prime factors have density zero. In particular I and II have density zero, and III has density one.

For completeness, reciprocal-prime divergence follows from Euler's product: if sum1/p converged, the products product_(p<=N)(1-1/p)^(-1) would be bounded, since -log(1-x)=x+O(x²). But each product is at least the harmonic sum through N, a contradiction.

This proof also shows that every fixed bound on the number of distinct prime factors gives a density-zero set. It does not assert that every III source has three distinct prime factors:12,18,20,28 and36 are immediate counterexamples.

## 3. Regrouping by aliquot target changes the measure

Define f_C(t)=#{n in C:s(n)=t}. For t>=2 this is finite. Then

\[
\sum_{2\le t\le T}f_C(t)
=\#\{n\in C:2\le s(n)\le T\}.
\]

This is not #{n in C:n<=T}. For instance the prime-square source (t-1)^2 lies over t whenever t-1 is prime. Small target rows can include sources quadratic in the cutoff. The infinitely many primes all lie in F(1), so even that boundary row cannot be assigned a finite uniform source weight.

Consequently class III having source density one says nothing by itself about the fraction of target rows that are empty, or about class III's share among sources whose aliquot values are bounded. Losing this distinction would manufacture an unjustified conclusion about untouchable-number density.

For an exact finite-period version of the argument, the polynomial

\[
\prod_{p\le p_K}((p-1)+z)
\]

counts residues modulo M by how many of those primes divide them. The sum of its coefficients of z^0,z^1,z^2 divided by M is exactly Pr(X_K<=2). This polynomial retains the support multiplicity that a scalar average would discard.

## 4. Checks

Run `python3 check_source_density.py` (and with `-O`). It compares the finite-period polynomial to brute CRT residue enumeration for the first five primes, checks coefficient counts at source cutoff100000 using both divisor-factorization classes and independently generated prime powers/distinct semiprimes, and tests class III examples with only two prime factors. The asymptotic proof is the argument above, not extrapolation from those checks.
