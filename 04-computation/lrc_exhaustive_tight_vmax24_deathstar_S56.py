import sys
from math import gcd
from functools import reduce
from itertools import combinations
def is_tight(fam):
    Q=2*fam[-1]+2; hit=False
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if 14*m>q: return False
            if 14*m==q: hit=True
    return hit
found=[]
for w in range(13,25):                         # Vmax = w up to 24 (GW has 24)
    cnt=0
    for rest in combinations(range(1,w),12):
        fam=list(rest)+[w]
        if reduce(gcd,fam)!=1: continue
        if is_tight(fam): found.append(tuple(fam)); cnt+=1
    print(f"Vmax={w}: new tight={cnt}  total={len(found)}",flush=True)
print("ALL primitive tight families Vmax<=24:",found,flush=True)
