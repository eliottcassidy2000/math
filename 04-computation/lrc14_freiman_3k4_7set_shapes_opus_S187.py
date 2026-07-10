from itertools import combinations
from math import gcd
from functools import reduce
def rss(A): return len(set(A[i]+A[j] for i in range(len(A)) for j in range(i+1,len(A))))
B=18
by={}
for combo in combinations(range(1,B+1),6):
    A=(0,)+combo
    if reduce(gcd,[A[i+1]-A[i] for i in range(6)])!=1: continue
    r=rss(A)
    if r<=13: by.setdefault(r,[]).append(A)
for r in [11,12,13]:
    sets=by.get(r,[]); spreads=[a[-1] for a in sets]
    print(f"|A+^A|={r}: {len(sets)} primitive shapes, spread range [{min(spreads)},{max(spreads)}]")
    # show the max-spread ones (potential stability violations)
    mx=max(spreads)
    for A in sets:
        if A[-1]>=9:
            print(f"   spread {A[-1]:>2}: {A}  (verify |A+^A|={rss(A)})")
print("\n=> the near-minimal FINITE FAMILY (affine-normalized primitive shapes):")
print(f"   |A+^A|=11: {len(by.get(11,[]))} (the AP), 12: {len(by.get(12,[]))}, 13: {len(by.get(13,[]))}")
