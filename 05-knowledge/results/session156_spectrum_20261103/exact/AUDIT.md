Independent audit of the one-apex family

Verdict: All checked mathematical claims are correct. In particular, the vertex counts 17 and 41 include r; d+5 is only the base-tournament order. The construction does not establish infinitely many prime quotients, and the source correctly avoids that claim.

Strongness: Let N=5+d. Since r->1, 1->v for every base v>1, and v->N->r for every base v<N (while N->r directly), r reaches and is reachable from every vertex. Hence every vertex reaches every other via r.

Count: Every Hamilton path is uniquely determined by the set S of base vertices before r, since the induced transitive blocks must each be in increasing order. Its validity is equivalent to (S empty or max(S)->r) and (complement empty or r->min(complement)). The source's i<j boundary pair is uniquely (3,4), and its i>j weighted sum is 3+19 sum_{t=0}^{d-1}2^t; the two empty-block paths give H(d)=19*2^d-13. As an independently implemented cross-check, enumerate the 2^(d+5) cut sets S, not endpoint/subset dynamic programming. This produces d=2: 63; d=5: 595; d=8: 4851; d=11: 38899. Direct exhaustive permutation enumeration independently gives d=1: 25, d=2: 63, d=3: 139.

Divisibility and examples: 19*2^d-13 is 5*2^d+1 modulo 7. As ord_7(2)=3, it vanishes exactly for d congruent to 2 mod 3. Then p_k=(76*8^k-13)/7. Trial division by every integer from 2 through floor(sqrt(p)) shows p_3=5,557 and p_11=93,262,146,997 prime. Their respective square-root bounds are 74 and 305388. Direct arithmetic verifies H(11)=38,899=7*5,557 and H(35)=652,835,028,979=7*93,262,146,997. Also p_7=22,769,077=103*221,059.

Sieve: For even k, p_k is divisible by 3; for k=1 mod 4, p_k is divisible by 5. As p_k>=9 and p_1=85, no exception is a prime quotient. This restricts prime-producing k to 3 mod 4. Modulo 49, 8^k=1+7k, so p_k=2-k mod 7; as p_k>7, exclude k=2 mod 7. The combined residues are exactly 3,7,11,15,19,27 mod 28. These are necessary only.

Valuations: Write F(k)=76*8^k-13. For each t>=1, k mod 7^(t-1) determines F(k) mod 7^t. Every k works at t=1. If 7^t divides F(k), then for 0<=a<7,
 F(k+a*7^(t-1))/7^t = F(k)/7^t + a*76*8^k (mod 7).
The coefficient of a is a unit mod 7, so exactly one of the seven lifts is a root at level t+1 and the other six have valuation exactly t. Thus there is one root class at each level, and for each t there are infinitely many nonnegative k with v_7(F(k))=t. Direct checks give root classes k=0 mod 1, 2 mod 7, 9 mod 49, 303 mod 343, 1675 mod 2401 for t=1,...,5.

Reproducible independently written cut-set and trial-division checks:

from math import isqrt

def cut_count(d):
    n = 5 + d
    def edge(a,b):
        if a == 0: return b in (1,4,5)
        if b == 0: return a not in (1,4,5)
        return a < b
    total = 0
    for mask in range(1 << n):
        a = [i for i in range(1,n+1) if mask >> (i-1) & 1]
        b = [i for i in range(1,n+1) if not (mask >> (i-1) & 1)]
        total += (not a or edge(a[-1],0)) and (not b or edge(0,b[0]))
    return total

for d in (2,5,8,11):
    print(d,cut_count(d),19*2**d-13)
for p in (5557,93262146997):
    print(p,isqrt(p),all(p % q for q in range(2,isqrt(p)+1)))
for t in range(1,6):
    print(t,[k for k in range(7**(t-1)) if (76*pow(8,k,7**t)-13) % (7**t) == 0])
