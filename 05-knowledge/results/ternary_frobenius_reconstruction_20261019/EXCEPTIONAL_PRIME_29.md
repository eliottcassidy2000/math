# The degree-42 reconstruction defect in characteristic 29

Status: explicit finite coefficient certificate, 2026-10-19.

Let C f(x,y,z)=f(y,z,x), A f(x,y,z)=f(-y,x-y,z), and N=I+A+A². Write Alt(a,b,c) for the determinant with rows (x^a,x^b,x^c), (y^a,y^b,y^c), (z^a,z^b,z^c).

The following table specifies a polynomial E over the integers: for each row, add coefficient times Alt(a,b,c).

| a | b | c | coefficient |
|---:|---:|---:|---:|
| 2 | 8 | 32 | 4 |
| 2 | 12 | 28 | 1 |
| 2 | 14 | 26 | 2 |
| 2 | 16 | 24 | -13 |
| 2 | 18 | 22 | -8 |
| 4 | 6 | 32 | -12 |
| 4 | 10 | 28 | 11 |
| 4 | 12 | 26 | 1 |
| 4 | 14 | 24 | -5 |
| 4 | 16 | 22 | -10 |
| 4 | 18 | 20 | -11 |
| 6 | 8 | 28 | -9 |
| 6 | 10 | 26 | -8 |
| 6 | 12 | 24 | 7 |
| 6 | 14 | 22 | -5 |
| 6 | 16 | 20 | 11 |
| 8 | 10 | 24 | 8 |
| 8 | 12 | 22 | 14 |
| 8 | 14 | 20 | 12 |
| 8 | 16 | 18 | 2 |
| 10 | 12 | 20 | -14 |
| 10 | 14 | 18 | -1 |
| 12 | 14 | 16 | -1 |

FINITE-EXACT. Every coefficient of NE is divisible by 29, and their greatest common divisor is exactly 29. Hence the reduction of E is a nonzero homogeneous element of degree 42 satisfying CE=E and NE=0 over F_29. Nonzero is immediate from the coefficient one of x²y¹²z²⁸. Every term is divisible by x²y²z², so restriction to every coordinate plane is zero.

The 23 alternants have disjoint exponent multisets and expand to 138 terms. Each exponent is even. Consequently E is divisible over the integers by

    x²y²z² (x²-y²)(y²-z²)(z²-x²),

and the quotient is a symmetric degree-30 polynomial in x,y,z, equivalently a symmetric degree-15 polynomial in x²,y²,z². This divisibility follows from alternating polynomials in the three squared variables; no numerical factorization is needed.

The certificate proves delta_29 <= 42 < min(3*29+3,2*29+12)=70. It therefore disproves the universal-equality conjecture suggested by the first eight prime thresholds. The two all-prime construction bounds remain valid. A smaller modular class can arise from additional prime-specific torsion, so rational-family searches alone cannot prove modular minimality.

## Verification

Run `python3 -B verify_exception_29.py` or `python3 -B -O verify_exception_29.py` from this directory. The verifier reconstructs E directly from the table, literally substitutes the three linear maps, and checks all integer norm coefficients before reducing them. It also checks the degree, support, cyclicity, alternation, parity of exponents, nonzero coefficient, and integer norm coefficient gcd. It uses only standard-library routines from `verify_periods_independent.py`, with no production search or matrix code.

`next_thresholds.json` records a separate production all-degree search: first interior kernels at (p,n)=(23,58),(29,42),(31,74), each of dimension one. Minimality is a separate rank claim and is not needed for the counterexample. The review receipt identifies which of these extended rank claims have been independently reproduced.

No representation-theoretic explanation or all-prime exceptional-degree law is claimed. In particular, this is an even-degree interior class, distinct from the odd-degree supported Frobenius classes and their all-precision theorem.

## Finite 29-adic lift obstruction

A second FINITE-EXACT calculation concerns only cyclic all-interior inputs of degree 42. Let M be their integral norm matrix in the distinct cyclic-orbit basis; it has 274 columns and 946 monomial rows. The displayed E gives a primitive mod-29 kernel vector. To lift it to a cyclic all-interior solution modulo 29², one must solve

    M g = -NE/29  (mod 29).

The matrix M has rank 273 modulo 29, whereas adjoining this right-hand side raises the rank to 274. Therefore no such lift exists. The exact ranks are in `exception_29_lift.json`; the production check is `probe_exception_29_lift.py`, with independent reproduction in the review artifacts.

By the characteristic-zero restriction theorem, M has full column rank over Q. Its mod-29 nullity one says exactly one nonzero Smith factor is divisible by 29. The failure to lift a generator of the mod-29 kernel forces that factor's 29-adic valuation to be exactly one. Thus the 29-primary torsion of coker M is precisely Z/29. This is an even-degree, interior, finite-prime statement, separate from the all-degree odd supported theorem.
