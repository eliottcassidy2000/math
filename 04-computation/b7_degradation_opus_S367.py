from fractions import Fraction as F
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
def fold_inter(S):
    cur=teeth01(S[0])
    for x in S[1:]:
        cur=inter(cur,teeth01(x))
        if not cur: break
    return cur
def lowerk(S):
    k=len(S); p=(2*LAM)**k
    for i in range(k-1):
        f=1-F(7*S[i],S[i+1])
        if f<=0: return F(0)
        p*=f
    return p
random.seed(9)
print("DEGRADATION OF THE CONTAINMENT FLOOR WITH k (separated regime, ratios 8-10):")
print("  k :  exact/floor  (1.0 = tight)")
prev=None
for k in (2,3,4,5,6):
    rat=[]
    for _ in range(12):
        S=[random.randint(1,3)]
        for _ in range(k-1): S.append(S[-1]*random.randint(8,10))
        ex=meas(fold_inter(S)); lb=lowerk(S)
        if lb>0: rat.append(float(ex/lb))
    rat.sort(); med=rat[len(rat)//2]
    g = f"   x{med/prev:.1f} vs previous" if prev else ""
    print(f"  {k} :  {med:8.2f}{g}")
    prev=med
print()
print("  => the floor is TIGHT only at k=2 (where THM-1012/1025 proved it sharp)")
print("     and degrades by a roughly constant factor per level thereafter.")
