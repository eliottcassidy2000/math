# Sylvester, Mersenne, and Lossless Encoding

## The Sylvester sequence IS the Phi_6 orbit from 2

Sylvester: a_1=2, a_{n+1} = a_1·a_2·...·a_n + 1.
= 2, 3, 7, 43, 1807, 3263443, ...

This is ALSO: Phi_6(2)=3, Phi_6(3)=7, Phi_6(7)=43, Phi_6(43)=1807.
Same recurrence: a_{n+1} = a_n^2 - a_n + 1 = a_n(a_n-1) + 1.

## The Egyptian fraction decomposition of unity

1 = 1/2 + 1/3 + 1/7 + 1/43 + 1/1807 + ...

The UNIT decomposes into reciprocals of {doubler, curvature, forbidden, bridge, ...}.
Each partial sum + remainder = 1 exactly. PERFECTLY LOSSLESS.

Remainders: 1/2, 1/6, 1/42, 1/1806, ...
= reciprocals of the running products 2, 6, 42, 1806, ...
= reciprocals of {doubler, flat, Hurwitz, Hurwitz·43, ...}

The remainder after 1/2 + 1/3 + 1/7 = 1/42.
The Hurwitz number IS the denominator of the remainder after exhausting the Platonic primes.

## Mersenne and Sylvester share EXACTLY {3, 7}

Mersenne primes: 3, 7, 31, 127, 8191, ...
Sylvester terms: 2, 3, 7, 43, 1807, ...
Intersection: {3, 7}. The curvature and the forbidden.

After {3, 7}: Mersenne goes to 31, 127 (binary structure).
Sylvester goes to 43, 1807 (cyclotomic structure).
They DIVERGE. Binary and cyclotomic capture different aspects.
The AGREEMENT at {3, 7} = the lossless zone where both encodings coincide.

## 43 as the bridge

43 = Phi_6(7) = the 4th Sylvester term.
43 = 1 mod 3 = SPLITS in Eisenstein integers.
43 appears in BOTH the Sylvester chain and the Eisenstein splitting primes.
43 = Hurwitz + 1 = 42 + 1. The first number past the Hurwitz boundary.

The Sylvester and Lucas-Lehmer chains both contain 7:
Sylvester: 2, 3, 7, 43, 1807 (via x -> x^2 - x + 1)
Lucas-Lehmer: 4, 7, 47, 2207 (via x -> x^2 - 2)
Difference at x=7: 47-43 = 4 = the Cayley period.

## 6174 in the Sylvester basis

6174 = 2^1 · 3^2 · 7^3 = a_1^1 · a_2^2 · a_3^3.
Staircase exponents in the Sylvester basis. The Kaprekar constant
IS the Sylvester product with ascending powers.

## Lossless encoding theorem

Every product of Sylvester terms can be recovered from its value:
if P = a_1^{e_1} · a_2^{e_2} · ... · a_k^{e_k}, the factorization is unique
because the Sylvester terms are pairwise coprime (any common factor d
of a_i and a_j with i<j would divide a_j - (product-of-earlier) = 1).

The Sylvester sequence gives a UNIQUE FACTORIZATION for any number
expressible as a product of its terms. This is LOSSLESS.
