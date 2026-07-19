"""Is inf M(W) over SPREAD (rho>=6.5) far-from-AP covering-2..12 12-cores bounded away from 1/13?
If yes -> a crude gap c>0 exists (provable target). If it approaches 1/13 -> sharp Freiman needed.
death-star-S58i."""
from fractions import Fraction as F
from math import gcd
import random,time
def Mexact(V):
    mx=max(V);Qs=set()
    for i in range(len(V)):
        for j in range(i,len(V)):
            s=V[i]+V[j]
            for d in range(2,2*mx+1):
                if s%d==0:Qs.add(d)
    best=F(0)
    for q in Qs:
        for a in range(1,q):
            if gcd(a,q)!=1:continue
            m=min(min((v*a)%q,q-((v*a)%q)) for v in V)
            if F(m,q)>best:best=F(m,q)
    return best
def covers212(W): return all(any(w%k==0 for w in W) for k in range(2,13))
def ham(W):
    best=12
    for d in range(1,max(W)//12+2):
        ap=set(d*i for i in range(1,13)); best=min(best,12-len(set(W)&ap))
    return best
random.seed(17); THR=F(1,13); mn=(None,F(1)); n=0; below13=0
t0=time.time()
while time.time()-t0<170:
    hi=random.choice([26,30,40,50,60])
    W=sorted(random.sample(range(2,hi+1),12))
    if not covers212(W): continue
    rho=max(W)/min(W)
    if rho<6.5: continue
    h=ham(W)
    if h<7: continue
    M=Mexact(W); n+=1
    if M<THR: below13+=1; print("!!! far covering-2..12 12-core with M<1/13:",W,M,flush=True)
    if M<mn[1]: mn=(W,M); print("new min M(W)=%s=%.6f margin=%+.6f at %s (rho=%.1f H=%d)"%(M,float(M),float(M-THR),W,rho,h),flush=True)
print("searched %d spread far covering-2..12 cores; M<1/13 count=%d"%(n,below13))
print("INFIMUM M(W)=%s=%.6f, margin over 1/13 = %.6f"%(mn[1],float(mn[1]),float(mn[1]-THR)))
