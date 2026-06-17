from fractions import Fraction as Fr
from math import gcd
from functools import reduce
from itertools import combinations
import sys

def danger_arcs(v):
    w=Fr(1,14*v); A=[]
    for k in range(v+1):
        c=Fr(k,v); lo=c-w; hi=c+w
        if lo<0: A+=[(Fr(0),hi),(1+lo,Fr(1))]
        elif hi>1: A+=[(lo,Fr(1)),(Fr(0),hi-1)]
        else: A.append((lo,hi))
    return A
_cache={}
def darcs(v):
    if v not in _cache: _cache[v]=danger_arcs(v)
    return _cache[v]
def L_exact(S):
    A=[]
    for v in S: A.extend(darcs(v))
    A=sorted((a,b) for a,b in A if b>a)
    tot=Fr(0); cl=ch=None
    for a,b in A:
        if ch is None: cl,ch=a,b
        elif a<=ch:
            if b>ch: ch=b
        else: tot+=ch-cl; cl,ch=a,b
    if ch is not None: tot+=ch-cl
    return 1-tot
def lcm(a,b): return a*b//gcd(a,b)
def lcmset(S): return reduce(lcm,S,1)

# drop-6 family: {1..5,7..13} + w. Disprove claim: min positive L = 19/10626 at w=69.
core6=[1,2,3,4,5,7,8,9,10,11,12,13]
best=None; tightw=[]
for w in range(14,3001):
    if w in core6: continue
    S=core6+[w]
    L=L_exact(S)
    if L==0: tightw.append(w)
    if L>0 and (best is None or L<best[0]): best=(L,w)
print(f"drop-6 family {{1..5,7..13,w}}: tight w={tightw}; min positive L={best[0]}={float(best[0]):.6f} at w={best[1]}")
print(f"  19/10626 = {float(Fr(19,10626)):.6f}; match? {best[0]==Fr(19,10626) and best[1]==69}")
print(f"  L still ABOVE 1/1260={float(Fr(1,1260)):.6f}? {best[0]>Fr(1,1260)}")
sys.stdout.flush()

# Below-1/1260 hunt: ALL single perturbations of AP, w up to 2000 (the provable region)
print()
print("Single-perturbation global min positive L (all drops e, all w<=2000):")
gmin=None
APset=set(range(1,14))
for e in range(1,14):
    base=sorted(APset-{e})
    for w in range(14,2001):
        if w in APset: continue
        S=base+[w]
        L=L_exact(S)
        if L>0 and (gmin is None or L<gmin[0]): gmin=(L,e,w)
print(f"  min positive single-pert L = {gmin[0]}={float(gmin[0]):.6f} (drop {gmin[1]}, add {gmin[2]})")
print(f"  == 1/1260? {gmin[0]==Fr(1,1260)}; nothing below 1/1260 found: {gmin[0]>=Fr(1,1260)}")
sys.stdout.flush()

# 2-drop scan for anything below 1/1260 (entries to 60, this is the uncontrolled regime)
print()
print("2-drop scan (drop 2 from AP, add 2 new entries <=60) hunting for L < 1/1260:")
below=[]
allnew=[w for w in range(14,61)]
drops=list(combinations(range(1,14),2))
cnt=0
target=Fr(1,1260)
for (e1,e2) in drops:
    base=sorted(APset-{e1,e2})
    for i in range(len(allnew)):
        for j in range(i+1,len(allnew)):
            S=base+[allnew[i],allnew[j]]
            cnt+=1
            L=L_exact(S)
            if 0<L<target:
                below.append((L,S))
print(f"  scanned {cnt} configs; # with 0<L<1/1260: {len(below)}")
if below:
    below.sort()
    for L,S in below[:5]: print("   ", L, float(L), S)
else:
    print("   NONE below 1/1260 (champion 1/1260 survives 2-drop scan)")
