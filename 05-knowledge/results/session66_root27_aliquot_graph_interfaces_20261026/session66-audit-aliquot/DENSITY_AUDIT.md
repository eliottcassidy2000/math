# Independent exact lower-density certificate audit

2026-10-26. **PASS for both numerical enclosures and their stated consequence.**
This is a cited theorem instantiated by a new session calculation, not a
proof of density existence, an exact density, or a claim of a record bound.

## Source theorem and coefficients

[Chen–Zhao, Nonaliquot numbers](https://publi.math.unideb.hu/paper/1551/download/10_5486_PMD_2011_4820.pdf),
Publ. Math. Debrecen78/2 (2011),439–442, Theorem1 (PDF page2), states for
each fixed positive M that `U(X)>=g_M X+o_M(X)`, where

    g_M = sum_(d|M) phi(M/d)/(M/d) max(0,1/(2d)-1/(sigma(2d)-2d)).

Thus `liminf U(X)/X>=g_M`. The normalization is among **all** integers,
not a proportion conditional on even targets; no extra factor2 applies.
The theorem's displayed M has exponents6,5,4,3,2 on2,3,5,7,11 and exponent1
on13,17,19,23,29,31,37,41. The sign and each denominator were checked against
the primary theorem and its proof. The result is not a statement that U/X
converges. The o-term is for fixed M, not uniform in a growing M(X).

## Independent calculation

The producer tracks the rational weight phi(M/d)/(M/d). The independent
checker instead streams divisors and integer complement totients, using
the algebraically equivalent positive term

    phi(M/d)*(sigma(2d)-4d)/(2*M*(sigma(2d)-2d)).

It handles the two-adic factor separately from odd divisors and checks
small M using literal divisors, gcd-counted totients and exact Fractions.
Every positive term is rounded down at denominator D=10^24. If L is the
sum of integer floors and t>0 is the number of positive terms, then

    L/D <= g_M < (L+t)/D.

The strict upper bound follows because each individual rounding error is
strictly less than1/D; no floating-point operation is used.

| M | Divisors | Positive terms | L | L+t |
|---|---:|---:|---:|---:|
| 53133049021337845560000 |645120|644034|60275780518518619919468|60275780518518620563502|
| 43 times preceding M |1290240|1288370|60452936134095005620815|60452936134095006909185|

Both endpoints match the producer exactly, in normal and optimized Python.
In particular

    liminf U(X)/X > 0.060452936134095.

The certified enlarged lower endpoint is strictly greater than the old
upper endpoint. This establishes the computed improvement without a
literature-priority claim. Source theorem plus integer-certified evaluation
is the dependency chain; finite census fractions are not used.

## Monotonicity proof audit

For a new odd prime q, writing `h(d)=max(0,1/(2d)-1/s(2d))` gives the
claimed difference `sum w_M(d)[h(qd)-h(d)/q]`. When h(d)>0 its bracket is
`sigma(2d)/(q*s(2d)*(q*s(2d)+sigma(2d)))`, strictly positive. When h(d)=0,
the bracket is nonnegative directly. Hence adding43 increases the bound;
d6 is an available strictly-positive witness since s(12)=16>12.

The stronger root-authored identity also checks:
`g_M=E[H(gcd(R,M))]`, `H(d)=d*h(d)`, R uniform modulo M. The exact gcd
probability is phi(M/d)/M. If M divides M', a uniform residue modulo M'
projects uniformly modulo M and its two gcds are divisibility-ordered.
The ratio sigma(2d)/(2d), and hence H(d), is nondecreasing under divisibility.
Therefore `g_M<=g_M'`. This does not allow the fixed-M o-term to be used
uniformly with M depending on X. The checker tests all gcd-expectation
identities and divisibility refinements for moduli at most240.

## Repairs found and final boundary

An initial generic helper omitted the extra factor3 in sigma(2d) when M
was odd. The hostile M9 has exact g9=1/126, whereas that version returned0.
The two published even-M computations were unaffected. The producer repaired
the missing factor and added the M9 regression control.

A second generic boundary was that t=0 gives g_M=0, so L+t=0 is not a strict
upper bound. The producer repaired this by using `L+max(1,t)`; the supplied
two certificates both have t>0 and are unchanged. The final generic helper
is independently checked on odd M9 and exact-zero M1,2,3, among other cases.
Both repairs are included in final script SHA256
`d187b2e9f574191b3bef82c25c8e85ba89a79e88c79aebc35e3c531c2c07e5a2`.
No unresolved numerical, source, or proof blocker remains.

Reproduce the independent audit with `python3 audit_density.py` and
`python3 -O audit_density.py`; compare their redirected JSON outputs. The
frozen source inputs and both results are retained in this audit directory.
