# Independent audit of rank 30/31/32 certificates

**VERIFIED:** the supplied points on ICARM curves 273, 302 and 4239 are
respectively 30, 31 and 32 independent rational points. This upgrades the
initial point-membership-only audit; it does not certify an upper rank bound.
No BSD, GRH, floating height determinant or numerical regulator is needed.

## Exact mechanism

For an odd good-reduction prime p, complete the square in the Weierstrass
model to obtain

    Y² = g(x) = x³ + b2*x²/4 + b4*x/2 + b6/4.

If r is an Fp-root of g, the rational 2-torsion point T=(r,0) gives the
usual 2-descent character

    chi_r(O)=1,
    chi_r(T)=g'(r),
    chi_r(P)=x(P)-r otherwise,

with values in Fp*/Fp*². The derivative is nonzero because reduction is
good. This is a group homomorphism. In the ordinary chart, intersect a line
Y=mx+b with the monic cubic: at x=r its three roots give

    (x(P)-r)(x(Q)-r)(x(R)-r)=(mr+b)².

Since P+Q+R=O and chi_r(-P)=chi_r(P), this is the multiplicative relation.
The derivative value at T is exactly the extension across the omitted
zero-coordinate case; O has neutral value. Equivalently this is the
standard function-divisor descent map for div(x-r)=2(T)-2(O).
Compose with the group-homomorphic good reduction E(Q)->E(Fp) and the
quadratic character to obtain a homomorphism E(Q)->F2 annihilating 2E(Q).
All supplied row evaluations avoid T and O, and the checker verifies that.

The rows across different primes form a matrix on the named rational points.
Full column rank proves their classes independent modulo 2E(Q).
To pass to Z-independence one must eliminate rational 2-torsion: a
rootless reduction of 4g at an odd good prime does so. Indeed a rational
root of the monic g is p-integral and would reduce to a root modulo p.

Suppose now a nonzero integer relation exists. Its coefficient vector is
zero modulo 2 by full column rank. Divide all coefficients by 2; the new
combination is 2-torsion and hence zero. Repeating contradicts the finite
2-adic valuation of the nonzero coefficient vector. Thus no integer
relation exists.

## Independent computational reconstruction

`check_independence_independent.py` imports no producer modules. It uses:

- exact point membership and independently recomputed discriminants;
- trial-division primality and good-reduction tests;
- direct finite-field square sets, rather than producer Euler-power signs;
- direct-search modular inverses;
- dense, column-wise F2 elimination, rather than producer integer pivot rows;
- exhaustive rootlessness at the named 2-torsion exclusion primes.

The result is ranks **30, 31, 32**, with exclusion primes **59, 59, 7**.
Normal and `python -O` outputs agree byte-for-byte. These are unconditional
lower bounds. Exact ranks, analytic bounds, family generic ranks and any
record-holder maximality remain distinct claims requiring their own evidence.

The method is a finite decision certificate despite the enormous point
coordinates. Its success depends on retaining the named curve, points, good
prime, cubic root and character row, not merely a scalar numerical regulator.
