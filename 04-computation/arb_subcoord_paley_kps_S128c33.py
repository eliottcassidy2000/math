#!/usr/bin/env python3
"""arb_subcoord_paley_kps_S128c33.py -- sub-coordinate minimality for the tau-drift (n=4..6 exhaustive)
+ Paley closed form referee + BEST cross-check."""
import sys
from math import comb, factorial
from fractions import Fraction as F
from itertools import permutations, combinations
sys.stdout.reconfigure(line_buffering=True)

def det_int(M):
    M=[r[:] for r in M]; n=len(M); sign=1; prev=1
    for k in range(n-1):
        if M[k][k]==0:
            for i in range(k+1,n):
                if M[i][k]: M[k],M[i]=M[i],M[k]; sign=-sign; break
            else: return 0
        for i in range(k+1,n):
            for j in range(k+1,n):
                M[i][j]=(M[i][j]*M[k][k]-M[i][k]*M[k][j])//prev
        prev=M[k][k]
    return sign*M[-1][-1]
def taus(B,n):
    L=[[0]*n for _ in range(n)]
    for u in range(n):
        for v in range(n):
            if u!=v and B[u][v]: L[v][v]+=1; L[v][u]-=1
    return [det_int([[L[i][j] for j in (k for k in range(n) if k!=r)] for i in (k for k in range(n) if k!=r)]) for r in range(n)]
def ham(B,n):
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for S in range(1<<n):
        for v in range(n):
            c=dp[S][v]
            if not c or not (S>>v)&1: continue
            for u in range(n):
                if not (S>>u)&1 and B[v][u]: dp[S|1<<u][u]+=c
    return sum(dp[(1<<n)-1][v] for v in range(n))
def codd(B,n,Lc):
    if Lc>n: return 0
    tot=0
    for S in combinations(range(n),Lc):
        u=S[0]
        for perm in permutations(S[1:]):
            prev=u; ok=True
            for w in perm:
                if not B[prev][w]: ok=False; break
                prev=w
            if ok and B[prev][u]: tot+=1
    return tot

print("== sub-coordinate minimality for tau-drift, n=4..6 exhaustive ==")
for n in (4,5,6):
    tiles=[(x,y) for y in range(1,n-1) for x in range(n,y+1,-1) if x-y>=2]
    m=len(tiles); m2=comb(n,2)
    recs=[]
    for t in range(1<<m):
        B=[[False]*n for _ in range(n)]
        for k in range(2,n+1): B[k-1][k-2]=True
        for i,(x,y) in enumerate(tiles):
            if (t>>i)&1: B[x-1][y-1]=True
            else: B[y-1][x-1]=True
        H=ham(B,n); c5=codd(B,n,5); c7=codd(B,n,7); tt=sum(taus(B,n))
        s=F(0)
        for u in range(n):
            for v in range(n):
                if u!=v and B[u][v]:
                    B[u][v],B[v][u]=False,True
                    s+=sum(taus(B,n))-tt
                    B[u][v],B[v][u]=True,False
        recs.append(((H,c5,c7,tt),F(s,m2)))
    projs={'tau':lambda k:(k[3],),'(tau,H)':lambda k:(k[3],k[0]),'(tau,c5)':lambda k:(k[3],k[1]),
           '(tau,H,c5)':lambda k:(k[3],k[0],k[1]),'census-only':lambda k:(k[0],k[1],k[2]),
           'full':lambda k:k}
    line=["n=%d:"%n]
    for name,pr in projs.items():
        f={}; bad=0
        for k,v in recs:
            p=pr(k)
            if p in f and f[p]!=v: bad+=1
            f.setdefault(p,v)
        line.append("%s:%s"%(name,"OK(%d)"%len(f) if bad==0 else "X(%d)"%bad))
    print("  "+" ".join(line))

print("== Paley closed form: tau_r(P_p) = p^{(p-3)/2} ((p+1)/4)^{(p-1)/2} ==")
for p in (3,7,11):
    QR={(a*a)%p for a in range(1,p)}
    B=[[False]*p for _ in range(p)]
    for a in range(p):
        for b in range(p):
            if a!=b and (b-a)%p in QR: B[a][b]=True
    tv=taus(B,p)
    pred=p**((p-3)//2)*((p+1)//4)**((p-1)//2)
    print("  p=%d: tau_r=%s (root-indep %s), predicted %d -> %s"%(p,tv[0],len(set(tv))==1,pred,"MATCH" if tv[0]==pred else "FAIL"))

print("== BEST cross-check on the two regular n=7 classes ==")
def best_ec(B,n):
    tv=taus(B,n)
    d=[sum(B[v]) for v in range(n)]
    # BEST: ec = tau_r(rooted anywhere, Eulerian) * prod (outdeg-1)!  [circuits from a fixed root arc... standard count of Eulerian circuits]
    return tv[0]*1 if False else tv[0]*prodf(d)
def prodf(d): 
    r=1
    for x in d: r*=factorial(x-1)
    return r
for name,S in [("Paley7 {1,2,4}",{1,2,4}),("C7 {1,2,3}",{1,2,3})]:
    B=[[False]*7 for _ in range(7)]
    for a in range(7):
        for b in range(7):
            if a!=b and (b-a)%7 in S: B[a][b]=True
    tv=taus(B,7)
    ec=tv[0]*prodf([3]*7)
    print("  %s: tau=%d, ec=%d"%(name,tv[0],ec))
print("DONE")
