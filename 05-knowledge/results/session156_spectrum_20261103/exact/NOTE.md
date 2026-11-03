# Exact 7p counts from a one-apex family

## Status

A rigorous, explicit infinite family of tournament Hamilton-path counts divisible by 7, plus two finite verified examples whose quotients are prime. This does NOT prove infinitely many primes in the family, nor the proposed global {7,21} spectrum. No claim of literature novelty or of non-overlap with earlier exact constructions is made.

Let d>=1. Take the transitive tournament on vertices 1,...,5+d (i->j iff i<j), and add one vertex r. Put r->i precisely for i in {1,4,5}, and i->r for all other base vertices. Equivalently the orientation word (1 means r beats the base vertex) is 1 00 11 0^d. The tournament is strong: r reaches every base vertex through vertex 1, and every base vertex reaches r through the last base vertex. Its Hamilton-path count is exactly

    H(d) = 19*2^d - 13.

### Direct proof of count

Delete r from a Hamilton path. The two remaining blocks are individually increasing sequences. If S is the set before r, its maximum i, when present, must beat r, while the minimum j of the complement, when present, must be beaten by r. The empty-S and full-S cases each give one path. For nonempty blocks and i<j, all integers must belong to one of the blocks, so i=j-1; here the only allowed such pair is (i,j)=(3,4), contributing one. For i>j, all integers less than j belong to S, all greater than i belong to the complement, while the i-j-1 intermediate vertices are independently assigned to either block. Thus that pair contributes 2^(i-j-1).

The allowed pairs with i in {2,3} and j=1 contribute 1+2=3. For each tail index i=6+t, 0<=t<d, the allowed successors j=1,4,5 contribute (16+2+1)*2^t=19*2^t. Total is 1+1+1+3+19*(2^d-1)=19*2^d-13. All assignments are valid and give distinct Hamilton paths.

Since powers of 2 have order 3 modulo 7, H(d) is divisible by 7 exactly for d=2+3k, k>=0. Put

    p_k = H(2+3k)/7 = (76*8^k-13)/7.

Then p_k grows strictly; p_3=5557 is prime and the 17-vertex tournament with d=11 has H=38,899=7*5557. Also p_11=93,262,146,997 is prime, giving d=35, order 41, and H=652,835,028,979. These finite primality claims are certified by trial division through the square root in the accompanying Python commands below.

## Structural sieves, not prime generation

For even k, 76*8^k-13 is divisible by 3; for k congruent to 1 mod 4 it is divisible by 5. As p_k>=9 and p_1=85, neither gives a prime quotient. Thus primality forces k congruent to 3 mod 4. Also 8^k=(1+7)^k=1+7k (mod 49), so

    p_k = 2-k (mod 7).

Since all p_k>7, primality also requires k not congruent to 2 mod 7. Combined, k must be one of 3,7,11,15,19,27 modulo 28. These necessary conditions are not sufficient (e.g. p_7=22,769,077=103*221,059). It remains open here whether infinitely many p_k are prime.

In fact the family has arbitrarily high exact 7-adic valuations: for each t>=1 there is a unique k mod 7^(t-1) for which 7^t divides 76*8^k-13. To prove, for t=1 every k works. Inductively vary k to k+a*7^(t-1), a=0,...,6. Since (1+7)^(7^(t-1)) = 1+7^t (mod 7^(t+1)), the seven values of the expression divided by 7^t represent every residue mod 7, as 76*8^k is a 7-adic unit. Exactly one lifts. Choosing any of the other six lifts gives valuation exactly t. This says nothing about prime quotients for high valuations.

## Finite independent verification

`check.cpp` encodes the tournament orientations explicitly and computes all Hamilton-path counts with endpoint/subset DP. Its output is:

d=2 H=63
d=5 H=595
d=8 H=4851
d=11 H=38899

The direct combinatorial proof above is independent of this implementation and covers every d. Python exact trial division verifies the finite primality claims:

    from math import isqrt
    for p in (5557, 93262146997):
        print(p, all(p%d for d in range(2,isqrt(p)+1)))

The first output count also agrees with the previously reported one-apex recurrence (S155 primes/README.md). Prior sparse-bank theorem applies to its separate named mixed-prefix bank, and does not obstruct the present one-apex family.
