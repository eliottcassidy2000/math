# A local permutive gate in complementary checksum words

For an r-bit word v, write its bits b_0,...,b_(r-1) from the right
and its complement as vbar=2^r-1-v. The two-tape surcharge satisfies

$$
J(v):=\frac{p(v)+p(\overline v)}2
=\sum_{d=1}^{r-1} w_d\sum_{i=0}^{r-1-d}(b_i\oplus b_{i+d}),
\qquad w_1=1,\quad w_d=2^{d-2}\ (d\ge2).
$$

Indeed each zero bit in p adds its strictly higher suffix divided by
the appropriate power of two, plus the next digit. For complementary
words the ordered-bit event is exactly XOR, and the extra weight at
distance one comes from the next-digit term.

Modulo two only lags one and two survive. Their interior bits cancel
pairwise, giving the boundary formula J(v)≡b_1 XOR b_(r-2) for r≥2.
For r=3 the two positions coincide, so this parity is forced zero.
This makes concrete why demanding independent endpoint and checksum
parities can fail.

With all sites except j fixed, the exact change on toggling j from
0 to 1 is the sum over lags d of w_d times the sum of (1-2b)
at each available neighbor j-d,j+d. For any t≥1, ensure at least
t+2 neighbors on both sides; set opposite guards in the two sites
at every lag 1,...,t and equal guards at lag t+1. Lags through t
cancel. The contribution at lag t+1 is ±2^t. At lag t+2 the equal
weight on the two present sites makes divisibility by 2^(t+1)
automatic, and all later lags are divisible anyway. Thus toggling
j changes J by exactly 2^t modulo 2^(t+1).

Given correct lower bits, precisely one value at j supplies the
requested bit t without disturbing those lower bits. This is a local,
conditional analogue of Rule 30's left permutivity, not a conjugacy
or a solution of a Rule 30 problem. It works on literal complementary
addresses. Sharing guards across many such gates without changing
previous bits or losing absolute placement remains open.
