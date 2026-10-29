# Five Ogg values as cycle lengths: explicit arithmetic models

**Status:** FINITE-EXACT multiplicative-order and polynomial certificates;
PROVED consequences of the cyclic-kernel theorem. The Ogg-prime list itself
is CITED, as documented in [the geometry note](RICHELOT_FROBENIUS_AND_FANO.md).

The parameters must be typed carefully. In this note p is an Ogg **value used
as a branch-cycle length**. Every example has ground-field characteristic5,
not characteristic p. There is no asserted equivalence between Ogg rationality
in characteristic p and the binary multiplicative-order criterion.

| Ogg value p | ord_p(2) | `(1,p)` stable-kernel count |
|---:|---:|---:|
| 2 | not applicable | odd-cycle theorem does not apply |
| 3 | 2 | 0 |
| 5 | 4 | 0 |
| 7 | 3 | 2 |
| 11 | 10 | 0 |
| 13 | 12 | 0 |
| 17 | 8 | 0 |
| 19 | 18 | 0 |
| 23 | 11 | 2 |
| 29 | 28 | 0 |
| 31 | 5 | 8 |
| 41 | 20 | 0 |
| 47 | 23 | 2 |
| 59 | 58 | 0 |
| 71 | 35 | 2 |

The intersection is therefore exactly `{7,23,31,47,71}`, with genera
`{3,11,15,23,35}`. This intersection is reproducible arithmetic, but it does
not by itself explain why the Ogg list arises in moonshine. Non-Ogg values
73 and89 also pass the binary criterion.

## Five concrete curves over F5

For the following irreducible polynomials, take `C_p: y²=f_p(x)`:

| p | Polynomial f_p over F5 | Genus | Rational Lagrangian kernels |
|---:|---|---:|---:|
| 7 | x^7+x+1 | 3 | 2 |
| 23 | x^23+x²+1 | 11 | 2 |
| 31 | x^31+x+1 | 15 | 8 |
| 47 | x^47+x^4+x^3+4x+1 | 23 | 2 |
| 71 | x^71+x^4+2x+3 | 35 | 2 |

The odd-degree model has one rational branch point at infinity and a single
Frobenius cycle on its p finite roots. The kernel counts follow from the
all-order theorem, not by exhaustive enumeration of high-dimensional vector
spaces. In particular no enumeration of 2^35 or 2^36 codewords is performed.

For each displayed prime degree p, the exact irreducibility certificate is

`x^(5^p)=x mod f_p` and `gcd(f_p,x^5-x)=1`.

These are complete Rabin conditions because the degree is prime. The
polynomials are squarefree, so the curves are smooth hyperelliptic models
with the stated genera. The checker searches low coefficient vectors in a
fixed base5 order and retains the first certificate in its declared finite
search family. Candidate masks are6,26,6,771,638 respectively, with bit/base5
conventions printed explicitly. The earlier search restricted to correction
terms of degree at most3 found no p47 example; allowing degree4 supplies the
displayed one. This bounded search failure is not a nonexistence claim about
genus23 curves or sparse polynomials in general.

These are examples of rational **polarized 2-power isogeny kernels**, not
examples of high elliptic Mordell–Weil rank, finite Sha, or superspecial
Jacobians. Codomain equations have not been computed. The comparison with
Ogg primes remains an intersection of two separately defined predicates.

## Reproduction

Run `python check_ogg_cycle_models.py` and `python -O check_ogg_cycle_models.py`.
All arithmetic is exact over F5; the output includes all fifteen order-table
rows, the five low-coefficient vectors and21 explicit acceptance checks.
