# Every finite bank of dyadic windows has a synchronized composite subgrid

Let f_1,...,f_b be finitely many integer polynomials in d variables.
Consider their counts on the dyadic grid
N_i(k)=f_i(2^{k_1},...,2^{k_d}). Suppose at some common allowed corner
a all N_i(a) are odd integers larger than 7. Let L be the lcm of those
integers and M the multiplicative order of 2 modulo L. Then all points
with k_j = a_j mod M (in any allowed direction) satisfy

    N_i(a) divides N_i(k), for every i.

This is direct polynomial evaluation modulo L. Hence, wherever each
N_i(k)>7 N_i(a), no arm gives a prime or 7 times a prime. Indeed a
fixed divisor d>7 of a prime would have to be the prime, and a
divisor d>7 of 7p contains p, so 7p≤7d. The same proof with L a
multiple of the desired base moduli works even when some arms'
base values are themselves already composite.

In particular, for the finite bank consisting of the reported
one-axis and two-axis windows in this session, choose a large
common corner (treat a missing axis as a dummy variable), then
take the forward subgrid far enough from the corner. Every arm
is blocked simultaneously. The base argument works without any
need for distinct factors, prime values at the corner, CRT
compatibility of different phase systems, or factorization of L.
It also works for phase-fixed branches with polynomial output
in the free powers of 2, by first fixing their legitimate phases.
Changing between two tracked readers can widen a window but
cannot remove this obstruction if the resulting polynomial bank
is finite.

The conclusion concerns the *parameters* of these grids. It
does not say that an integer produced on a blocked grid cannot
be produced in another way, or that an infinite varying-slope
source grammar misses an integer. Nor does it imply the primes
are sparse in the rest of a given family's grid. What it rules
out is the intuition that exact contiguous windows alone make
prime and 7p witnesses uniformly available at all sufficiently
large choices of their height parameters.
