import sys
from fractions import Fraction
from functools import reduce
from math import gcd
def covered(S, r):
    arcs=[]
    for v in S:
        for m in range(0, v+1):
            lo=Fraction(m,v)-Fraction(r,v); hi=Fraction(m,v)+Fraction(r,v)
            if hi<0 or lo>1: continue
            arcs.append((lo if lo>0 else Fraction(0), hi if hi<1 else Fraction(1)))
    arcs.sort(); cur=Fraction(0)
    for lo,hi in arcs:
        if lo>cur: return False
        if hi>cur: cur=hi
        if cur>=1: return True
    return cur>=1
def is_tight_or_level(S,k):
    floor=Fraction(1,k+1)
    # M <= floor ? (would mean tight, M=floor exactly since LRC says M>=floor)
    if covered(S, floor): 
        # M <= floor and M>=floor => M=floor (TIGHT)
        return "TIGHT (M=floor=1/%d)"%(k+1), floor
    # else M>floor: find level. test mediants a/(a(k+1)-1)
    for a in range(2, 40):
        q=a*(k+1)-1; r=Fraction(a,q)
        if covered(S,r):
            # M<=a/q; check not covered just below to confirm exactly this mediant
            below=r-Fraction(1,(4*max(S)**2*len(S)))
            if not covered(S,below):
                return "level a=%d (M=%d/%d)"%(a,a,q), r
            else:
                return "BELOW a=%d mediant (deeper, M<%d/%d)"%(a,a,q), r
    return "above all small mediants", None
print("=== verify F(211,5)={1..209,211,1050} : tight or level-5? ===")
for (mult,k) in [(5,211),(4,31),(3,13)]:
    S=sorted(set(list(range(1,k-1))+[k,mult*(k-1)]))
    desc,val=is_tight_or_level(S,k)
    print(f"  F({k},{mult}) (k-1={k-1}): {desc}")
print()
print("=== is a=5 reachable by F(k,5) at OTHER k? scan k where k-1 div by 2,3,5,7 ===")
# a=4 worked at k=31,61,91 (k-1 div 30). a=5 candidate: k-1 div 210? or finer.
for k in [211,421,631,841, 71,141, 106, 176, 246]:  # various
    S=sorted(set(list(range(1,k-1))+[k,5*(k-1)]))
    if len(S)!=k: continue
    desc,val=is_tight_or_level(S,k)
    fac=[p for p in [2,3,5,7] if (k-1)%p==0]
    print(f"  F({k},5) k-1={k-1} (primes {fac}): {desc}")
