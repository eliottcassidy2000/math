"""
The razor-thin line: where is M(S) relative to 1/14? Disproof = M < 1/14 strictly somewhere.
LRC(14) is TIGHT at {1..13} (non-covering, M=1/14 exactly). Map the boundary.
"""
from fractions import Fraction as F
def norm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
def M_exact(S):
    C=set()
    for i in range(len(S)):
        for j in range(len(S)):
            for d in (S[i]+S[j],abs(S[i]-S[j])):
                if d:
                    for m in range(d+1): C.add(F(m,d))
        for m in range(2*S[i]+1): C.add(F(m,2*S[i]))
    best=F(0);arg=F(0)
    for t in C:
        if 0<t<1:
            v=min(norm(s*t) for s in S)
            if v>best: best,arg=v,t
    return best,arg
def is_cov(S): return all(any(s%q==0 for s in S) for q in range(2,15))

thr=F(1,14)
print(f"threshold 1/14 = {float(thr):.5f}\n")
sets={
 "{1..13} (TIGHT, non-covering)": list(range(1,14)),
 "{1..11,13} (12-speed core)": [1,2,3,4,5,6,7,8,9,10,11,13],
 "{2..14} (covering)": list(range(2,15)),
 "{1..12,14} (covering)": list(range(1,13))+[14],
 "{1..11,13,14} (covering)": [1,2,3,4,5,6,7,8,9,10,11,13,14],
 "{1..11,13,84} (THM-523 min~7/89?)":[1,2,3,4,5,6,7,8,9,10,11,13,84],
}
for nm,S in sets.items():
    M,t=M_exact(S); cov=is_cov(S)
    margin=M-thr
    print(f"  {nm:36s}: M={M}={float(M):.5f}  margin={float(margin):+.5f}  covering={cov}  witness t={t}")
print()
# search for the LOWEST M covering set with bounded speeds (near-tight covering)
print("searching low-M covering 13-sets (speeds<=18, sample near {1..14}):")
import itertools, random
random.seed(0)
best=(F(1),None)
cands=set()
# structured: {1..13} with one swapped to make covering (add 14, drop k)
for drop in range(1,14):
    S=sorted(set(range(1,15))-{drop})
    if len(S)==13 and is_cov(S): cands.add(tuple(S))
# random covering near small speeds
tries=0
while len(cands)<60 and tries<200000:
    tries+=1; S=tuple(sorted(random.sample(range(1,19),13)))
    if is_cov(S): cands.add(S)
for S in cands:
    M,_=M_exact(list(S))
    if M<best[0]: best=(M,S)
print(f"  lowest M found: {best[0]}={float(best[0]):.5f} (margin {float(best[0]-thr):+.5f}) at {best[1]}")
print(f"  1/14={float(thr):.5f}, 7/89={float(F(7,89)):.5f}")
