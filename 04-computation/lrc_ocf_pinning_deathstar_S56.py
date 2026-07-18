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
def resid_tourney(fam):
    (M,(a,q))=M_arg(fam); R=sorted((v*a)%q for v in fam); n=len(R)
    adj=[[False]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i!=j and 0<2*((R[i]-R[j])%q)<q: adj[i][j]=True
    return adj,n,M
def H_count(adj,n):
    # number of Hamiltonian PATHS (Redei, odd). bitmask DP: dp[mask][v]=#paths on mask ending at v
    full=(1<<n)-1
    # incoming lists: preds[v] = [u : u->v]
    preds=[[u for u in range(n) if adj[u][v]] for v in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for v in range(n):
            d=dp[mask][v]
            if not d or not (mask>>v)&1: continue
            for w in range(n):
                if (mask>>w)&1: continue
                if adj[v][w]:      # extend v->w
                    dp[mask|1<<w][w]+=d
    return sum(dp[full][v] for v in range(n))
def c3(adj,n):
    s=[sum(adj[i]) for i in range(n)]
    return comb(n,3)-sum(comb(x,2) for x in s)
# rotational R_13
def R13():
    n=13; adj=[[False]*n for _ in range(n)]
    for i in range(n):
        for k in range(1,7): adj[i][(i+k)%n]=True
    return adj,n
def flip(adj,i,j):  # reverse edge between i,j
    A=[row[:] for row in adj]; A[i][j],A[j][i]=A[j][i],A[i][j]; return A

print("H (Hamiltonian-path count, Redei-odd) of residue tournaments:")
for name,fam in [("deep well {1..12,182}",list(range(1,13))+[182]),
                 ("ladder m=2",list(range(1,13))+[364]),
                 ("ladder m=3",list(range(1,13))+[546]),
                 ("dilate 3*",[3*v for v in range(1,13)]+[546]),
                 ("dilate 5*",[5*v for v in range(1,13)]+[910])]:
    adj,n,M=resid_tourney(fam); H=H_count(adj,n)
    print(f"  {name}: M={float(M):.4f}  H={H}  c3={c3(adj,n)}"); sys.stdout.flush()
a13,n=R13(); print(f"\n  pure R_13 (rotational): H={H_count(a13,n)}  c3={c3(a13,n)}")
# R_13 + one flip: try a few flips
for (i,j) in [(0,1),(0,6),(0,7)]:
    A=flip(a13,i,j); print(f"  R_13 + flip({i},{j}): H={H_count(A,13)}  c3={c3(A,13)}")
print("\nContrast M>1/13:")
for name,fam in [("{1..11,13,84}",list(range(1,12))+[13,84]),("{1..13} tight AP",list(range(1,14)))]:
    adj,n,M=resid_tourney(fam); print(f"  {name}: M={float(M):.4f} H={H_count(adj,n)} c3={c3(adj,n)}"); sys.stdout.flush()
