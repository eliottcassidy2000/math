#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S29 (HYP-4592) -- the BINDING WITNESS per N (argmax provenance).

For the canonical S_N={1..N-2}u{N,3(N-1)}, report the ACTUAL argmax witness that sets
M: its denominator q and the pair/doubling that produced it.  Reveals WHY N=25 (q=77
composite) still witnesses (intended sum(5,3(N-1))=3N+2 wins) and WHY N=31 misses (a
doubling dbl(16,16)=32 caps M at the floor 1/32).  Confirms the N==1 mod 6 law and its
sporadic N=31 exception (concurrent HYP-4572 trichotomy), independently.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations

def Mfull(S):
    """Return (M, best_q, binding_provenance): the argmax witness."""
    S=sorted(set(S)); Q={}
    for v in S: Q.setdefault(2*v,("dbl",v,v))
    for a,b in combinations(S,2):
        Q.setdefault(a+b,("sum",a,b)); Q.setdefault(abs(a-b),("dif",a,b))
    Q.pop(0,None); best=F(0); bq=None; bp=None
    for q,prov in Q.items():
        for a in range(1,q):
            mn=min(min((v*a)%q,q-((v*a)%q)) for v in S)
            val=F(mn,q)
            if val>best: best=val; bq=q; bp=prov
    return best,bq,bp

def isprime(m): return m>1 and all(m%d for d in range(2,int(m**.5)+1))

if __name__=="__main__":
    print("S_N={1..N-2}u{N,3(N-1)}: is it a GAP-WITNESS (M strictly inside (1/(N+1),2/(2N+1)))?")
    print("  binding = the actual argmax pair/denominator that sets M.\n")
    print(f"  {'N':>2} {'N%6':>3} {'3N+2':>5} {'pr?':>4} {'M(S_N)':>8} {'in-gap':>6} {'binding q':>9} {'binding pair':>16}")
    witness=[]
    for N in range(6,38):
        S=sorted(set(list(range(1,N-1))+[N,3*(N-1)]))
        Sp=[x//reduce(gcd,S) for x in S]
        M,bq,bp=Mfull(Sp); q=3*N+2
        ing = F(1,N+1)<M<F(2,2*N+1)
        if ing: witness.append(N)
        pair=f"{bp[0]}({bp[1]},{bp[2]})"
        print(f"  {N:>2} {N%6:>3} {q:>5} {str(isprime(q))[0]:>4} {str(M):>8} {str(ing):>6} {bq:>9} {pair:>16}")
    print(f"\n  GAP-WITNESS at N = {witness}")
    print(f"  == exactly 6k+1?  {witness == [N for N in range(6,38) if N%6==1]}  (False => sporadic exception, e.g. N=31)")
    print("  => N=25 (q=77 composite) witnesses via sum(5,72)=77; N=31 MISSES via dbl(16,16)=32 => M=1/32=floor.")
