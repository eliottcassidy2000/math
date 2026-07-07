#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S29 (HYP-4592) -- REFUTE the "prime q=3N+2" criterion.

S28 (HYP-4582) claimed the canonical bordered-AP S_N={1..N-2}u{N,3(N-1)} achieves the
mediant 3/(3N+2) ONLY when q=3N+2 is prime.  Extending the range REFUTES this: N=25
(q=77=7*11, COMPOSITE) IS a gap-witness.  The governing law is N==1 mod 6 (concurrent
HYP-4572 trichotomy), not primality.  Uses the S20-FIXED fast-M (checks ALL a).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

def Mfast(S):
    S=sorted(set(S)); Q=set()
    for v in S: Q.add(2*v)
    for a,b in combinations(S,2): Q.add(a+b); Q.add(abs(a-b))
    Q.discard(0); best=F(0)
    for q in Q:
        for a in range(1,q):
            mn=min(min((v*a)%q,q-((v*a)%q)) for v in S)
            val=F(mn,q)
            if val>best: best=val
    return best

def isprime(m): return m>1 and all(m%d for d in range(2,int(m**.5)+1))

if __name__=="__main__":
    print("S_N={1..N-2}u{N,3(N-1)}: mediant 3/(3N+2) achieved <=> 3N+2 prime?\n")
    print(f"  {'N':>2} {'3N+2':>5} {'prime?':>7} {'M(S_N)':>8} {'mediant':>8} {'ACHIEVED':>9}")
    ach=[]
    for N in range(4,27):
        S=sorted(set(list(range(1,N-1))+[N,3*(N-1)]))
        Sp=[x//reduce(gcd,S) for x in S]
        M=Mfast(Sp); med=F(3,3*N+2); q=3*N+2; p=isprime(q)
        got = (M==med and F(1,N+1)<M<F(2,2*N+1))
        if got: ach.append(N)
        print(f"  {N:>2} {q:>5} {str(p):>7} {str(M):>8} {str(med):>8} {str(got):>9}")
    print(f"\n  mediant ACHIEVED (via S_N) at N = {ach}")
    print(f"  3N+2 PRIME at N = {[N for N in range(4,27) if isprime(3*N+2)]}")
    print("  => sets DIFFER: N=25 (q=77=7*11 composite) achieves; N=5,9,15,17,23 prime FAIL.")
    print("     'prime q' REFUTED; law is N==1 mod 6 (HYP-4572), with a sporadic N=31 exception.")
