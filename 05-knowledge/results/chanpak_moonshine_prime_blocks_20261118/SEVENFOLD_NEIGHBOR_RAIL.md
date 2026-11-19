# A two-parameter sevenfold rail inside the neighbor triple

This is a follow-up to TWO_PORT_NEIGHBOR_FAMILIES.md, Family III.
It classifies divisibility on the whole family, not just selected
heights. It does not prove that infinitely many quotients are prime
or that all quotients are tournament counts in their own right.

Let N(m,k) be the central count of the triple, so actual tournaments
give N-2, N and N+2:

    N(m,k) = A_m 2^k + B_m,
    A_m = 16·2^m + 24·3^m - 16,
    B_m = 3(2^m-1).

For m,k>=1, reducing modulo 7 gives the complete table:

    m mod 6       0  1  2  3  4  5
    (A,B) mod 7  (3,0) (4,3) (5,2) (4,0) (0,3) (0,2)
    k mod 3       —  0  0  —  —  —

Here 2 has order 3 and 3 has order 6 modulo 7. The table follows
by substituting each m class and trying 2^k=1,2,4. Thus every
(m,k) in the two stated windows is a physical sevenfold target.

For m≡2 (mod 6) the quotient N/7 is always divisible by 3,
because N≡2^k(2^m-1) (mod 3). For m≡1 (mod 6) the quotient
is coprime to 3, and of the two neighbors exactly one is
divisible by 3: the lower for odd k and the upper for even k.
In that m class A_m is a unit modulo 7. Write k=3+3j.
Since 2^3=1+7, modulo 49 we have
2^k≡8(1+7j). Consequently N/7 modulo 7 is an affine function
of j with nonzero slope. Exactly one j class modulo 7 must be
discarded to get v_7(N)=1. This is a modular eligibility theorem,
not a primality theorem.

This rail contains on the order of (log X)^2 distinct central
counts below X even when m≡1 mod6 and v_7(N)=1 are imposed.
For 1<=m<=k, v_2(N+3)=m, because N+3=3·2^m+A_m2^k and
8 divides A_m. Recover m from the valuation and then k from
the linear read in 2^k. Taking 1<=m<=k up to a small constant
times log X gives the lower bound: the m and k congruences
and one discarded j class cost only a constant fraction.
Both m and k are O(log X) when N<=X, giving the matching
upper bound.

## Paired prime spot checks

At m=1 and k=3+6j the quotient is q_j=(704·64^j+3)/7,
with q_0=101 and q_(j+1)=64q_j-27. It is divisible by 7
exactly at j≡4 (mod 7), as also follows by reading the
recurrence modulo 7. At j=0 the literal upper outputs are
705, 707=7·101, and prime 709. A separate count-101 tournament
has outgoing masks [8,13,49,84,99,11,39]. The finite five-template
atlas lists 101 as a bank miss.

At j=2 we get prime q_2=411941 and N=2883587. Its lower witness
is the inherited F1 one-arc join at addresses b=3,e=84871
(see POLYGON_COMPLETION_AND_PRIME_COUNTS.md); the upper witness
comes from the present rail. These are independent witnesses,
not an assertion that one is obtained from the other by deleting
a vertex. The code checks the first case by literal subset DP
and verifies both primes and the second F1 address exactly.

Run node sevenfold_neighbor_rail_replay.js.
