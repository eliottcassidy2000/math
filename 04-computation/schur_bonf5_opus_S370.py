# opus-2026-07-17-S370 -- THE SCHUR DIRECTION: is additive richness the
# DILATION-INVARIANT predictor that min-speed failed to be?
#
# CRUCIAL PROPERTY: the Schur count (# of triples a+b=c inside V) is
# DILATION-INVARIANT -- scaling every speed by k preserves a+b=c exactly.
# So unlike min speed (MISTAKE-154/156, THM-1055), it is an ADMISSIBLE
# predictor of a dilation-invariant quantity like BONF5.  That is the whole
# reason this direction can succeed where the threshold direction could not.
from fractions import Fraction as F
from itertools import combinations
import random
LAM=F(1,14)
def teeth01(x):
    w=LAM/x; out=[]
    for j in range(0,x+1):
        a,b=max(F(j,x)-w,F(0)), min(F(j,x)+w,F(1))
        if a<b: out.append((a,b))
    return out
def inter(u,v):
    out,i,j=[],0,0
    while i<len(u) and j<len(v):
        a,b=max(u[i][0],v[j][0]), min(u[i][1],v[j][1])
        if a<b: out.append((a,b))
        if u[i][1]<v[j][1]: i+=1
        else: j+=1
    return out
def meas(iv): return sum(b-a for a,b in iv)
def union(ivs):
    ivs=sorted(ivs); out=[]
    for a,b in ivs:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def uncovered(V):
    allv=[]
    for x in V: allv.extend(teeth01(x))
    return 1-sum(b-a for a,b in union(allv))
def bonf5(V):
    T=[teeth01(v) for v in V]; n=len(V)
    S=[F(0)]*6
    def rec(i,cur,k):
        if k==5: return
        for j in range(i,n):
            nxt = inter(cur,T[j]) if k>0 else T[j]
            if not nxt: continue
            S[k+1]+=meas(nxt)
            rec(j+1,nxt,k+1)
    rec(0,None,0)
    return 1-S[1]+S[2]-S[3]+S[4]-S[5]
def schur_count(V):
    s=set(V); return sum(1 for a,b in combinations(sorted(V),2) if a+b in s)

print("(5) SCHUR COUNT vs BONF5 -- is additive richness the real predictor?")
print("    (Schur count is DILATION-INVARIANT, so unlike min-speed it is admissible.)")
random.seed(370)
rows=[]
for _ in range(26):
    V=sorted(random.sample(range(3,140),13))
    rows.append((schur_count(V), float(bonf5(V)), float(uncovered(V)), min(V)))
rows.sort()
print("    schur  BONF5      uncovered   minspeed")
for sc,b,u,ms in rows:
    print(f"    {sc:5d}  {b:+.5f}   {u:.5f}     {ms:4d}")
lo=[r for r in rows if r[0]<=1]; hi=[r for r in rows if r[0]>=3]
if lo and hi:
    print()
    print(f"    schur<=1 (n={len(lo):2d}): median BONF5 = {sorted(x[1] for x in lo)[len(lo)//2]:+.5f}")
    print(f"    schur>=3 (n={len(hi):2d}): median BONF5 = {sorted(x[1] for x in hi)[len(hi)//2]:+.5f}")
    print(f"    BONF5>0 rate: schur<=1 -> {sum(1 for x in lo if x[1]>0)}/{len(lo)},"
          f"  schur>=3 -> {sum(1 for x in hi if x[1]>0)}/{len(hi)}")

print()
print("(6) THE KNOWN FAILURE FAMILIES -- are they additively rich?")
for name,V in [("THM-1055 primitive failure",[27,36,46,70,101,114,117,121,140,160,194,277,293]),
               ("AP {1..13} (tight)",list(range(1,14))),
               ("AP d=8",[1+8*i for i in range(13)]),
               ("odd sum-free-ish",[2*i+1 for i in range(13)])]:
    print(f"    {name:28s} schur={schur_count(V):3d}  BONF5={float(bonf5(V)):+.5f}  unc={float(uncovered(V)):.5f}")
