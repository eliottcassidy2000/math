from itertools import combinations
from math import comb
def E3(S):
    Sset=set(S); return sum(1 for a in S for b in S if a+b in Sset)  # ordered pairs (a,b) in S^2, a+b in S
# Exhaustive: among all k-subsets of {1..N}, find max E3 and ALL maximizers; check they are exactly dilated intervals c*{1..k}
for k in [3,4,5]:
    N = 6*k
    best=-1; maxers=[]
    for S in combinations(range(1,N+1), k):
        e=E3(S)
        if e>best: best=e; maxers=[S]
        elif e==best: maxers.append(S)
    # dilated intervals c*{1..k} within range
    dil=[tuple(c*i for i in range(1,k+1)) for c in range(1,N//k+1) if c*k<=N]
    dilset=set(dil)
    all_dil = all(m in dilset for m in maxers)
    print(f"k={k}: max E3={best} (C(k,2)={comb(k,2)}), #maximizers={len(maxers)}, all dilated-intervals? {all_dil}")
    if not all_dil:
        print("   NON-dilate maximizers:", [m for m in maxers if m not in dilset][:5])
    else:
        print(f"   maximizers = {maxers[:6]}{'...' if len(maxers)>6 else ''}")
# also confirm r(s_j) <= j-1 bound is tight only for intervals: print r-profile
def rprofile(S):
    S=sorted(S); Sset=set(S)
    return [sum(1 for a in S if (S[j]-a) in Sset and a < S[j]) for j in range(len(S))]
print("r-profile {1..5}:", rprofile([1,2,3,4,5]), "(should be 0,1,2,3,4)")
print("r-profile {2,4,6,8,10}:", rprofile([2,4,6,8,10]))
print("r-profile {1,2,3,4,6}:", rprofile([1,2,3,4,6]))
