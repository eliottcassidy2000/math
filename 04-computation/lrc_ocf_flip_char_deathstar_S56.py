from fractions import Fraction as F
from math import gcd, comb
import sys
def M_arg(fam):
    Q=2*max(fam)+2; best=F(0); arg=None
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q); arg=(a,q)
    return best,arg
def resid_adj(fam):
    (M,(a,q))=M_arg(fam); R=sorted((v*a)%q for v in fam); n=len(R)
    return [[i!=j and 0<2*((R[i]-R[j])%q)<q for j in range(n)] for i in range(n)],n,M
def H_count(adj,n):
    full=(1<<n)-1; dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        row=dp[mask]
        for v in range(n):
            d=row[v]
            if not d or not (mask>>v)&1: continue
            av=adj[v]
            for w in range(n):
                if (mask>>w)&1 or not av[w]: continue
                dp[mask|1<<w][w]+=d
    return sum(dp[full][v] for v in range(n))
def R13():
    n=13; return [[any((k:=(j-i)%n)==kk for kk in range(1,7)) for j in range(n)] for i in range(n)],n
PIN=3551083
print("Robust pinning: H for more M<1/13 families (expect 3551083):")
for name,fam in [("ladder m=7",list(range(1,13))+[1274]),
                 ("ladder m=14",list(range(1,13))+[2548]),
                 ("dilate 11*",[11*v for v in range(1,13)]+[2002]),
                 ("dilate 7* m=2",[7*v for v in range(1,13)]+[7*364])]:
    adj,n,M=resid_adj(fam); H=H_count(adj,n)
    print(f"  {name}: M={float(M):.4f} H={H}  ==PIN: {H==PIN}"); sys.stdout.flush()
# Characterize which flips of R_13 give the pinned H
a13,n=R13()
print("\nWhich single flips (i,0) of R_13 give H=PIN? (by chord distance d=(0-i)%13):")
byd={}
for i in range(1,13):
    A=[row[:] for row in a13]; A[0][i],A[i][0]=A[i][0],A[0][i]
    h=H_count(A,13); d=min(i,13-i); byd.setdefault(d,[]).append(h)
for d in sorted(byd): print(f"  chord distance {d}: H={byd[d][0]}  (=PIN: {byd[d][0]==PIN})")
print(f"\n=> M<1/13 tournament = R_13 with ONE flip of a distance-6 chord (the antipodal edge), H=3551083 PINNED.")
