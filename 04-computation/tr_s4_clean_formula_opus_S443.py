#!/usr/bin/env python3
"""
Clean global tr(S^4) formula + the SC4 mechanism (opus-S443 follow-up).
Confirm: (i) tr(S^4) = clean poly in (n, Sum s_i^2, SC4) with n-INDEPENDENT coeffs;
(ii) SC4 (# strongly-connected 4-subtournaments) decouples from (score,c3) at EXACTLY n=5
     -- the mechanism behind var(lambda^2)'s n=5 wall.
"""
import itertools, numpy as np
from math import comb
import numpy.linalg as la

def edges_iter(n):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(1<<len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            adj[i][j]=1 if (bits>>k)&1 else 0; adj[j][i]=1-adj[i][j]
        yield adj
def tr4(adj,n):
    S=np.array([[0.0 if i==j else (1.0 if adj[i][j] else -1.0) for j in range(n)] for i in range(n)])
    S2=S@S; return int(round(np.trace(S2@S2)))
def scores(adj,n): return [sum(adj[i]) for i in range(n)]
def c3(adj,n):
    return sum(1 for i,j,k in itertools.combinations(range(n),3)
               if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]))
def sc4(adj,n):
    cnt=0
    for quad in itertools.combinations(range(n),4):
        c=sum(1 for i,j,k in itertools.combinations(quad,3)
              if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]))
        if c==2: cnt+=1
    return cnt

# global fit with n-independent features: try tr4 = A*n^2 + B*n + C*sum_s2 + D*SC4 + E
print("Global tr(S^4) fit, features {n^2, n, sum s_i^2, SC4, 1}, pooled n=3..6:")
rows=[]; ys=[]
for n in range(3,7):
    for adj in edges_iter(n):
        s=scores(adj,n); ss2=sum(v*v for v in s)
        rows.append([n*n, n, ss2, sc4(adj,n), 1]); ys.append(tr4(adj,n))
A=np.array(rows,float); y=np.array(ys,float)
coef,_,_,_=la.lstsq(A,y,rcond=None); pred=A@coef
print(f"   coeffs [n^2,n,sum_s2,SC4,1] = {[round(c,3) for c in coef]}  exact={np.allclose(pred,y,atol=1e-6)}")

# the score part alone should equal 2*sum s_i(s_i-1)+... : try tr4 = 2*sum s_i^2 - stuff.
# known: tr(S^2) = -n(n-1). Guess tr(S^4)=2*(sum s_i^2)+... let's fit {sum_s2, SC4, n, n^2,1} already did.
# cleaner: express via c3.  ss2 = 2*sumC(s,2)+ (n)(n-1)/1 ... test tr4 = a*n(n-1) + 8*sum s_i^2 + 64*SC4 + ...
for guess_name, feats in [("{n(n-1), sum_s2, SC4}", lambda n,ss2,sc: [n*(n-1), ss2, sc, 1])]:
    rows=[]; ys=[]
    for n in range(3,7):
        for adj in edges_iter(n):
            s=scores(adj,n); ss2=sum(v*v for v in s)
            rows.append(feats(n,ss2,sc4(adj,n))); ys.append(tr4(adj,n))
    A=np.array(rows,float); y=np.array(ys,float); coef,_,_,_=la.lstsq(A,y,rcond=None); pred=A@coef
    print(f"   fit {guess_name}: {[round(c,3) for c in coef]}  exact={np.allclose(pred,y,atol=1e-6)}")

# SC4 mechanism: threshold of SC4 | (score,c3)
print("\nSC4 decoupling threshold (the mechanism behind var's n=5 wall):")
def sig(a,n): return (tuple(sorted(scores(a,n))), c3(a,n))
for n in range(4,7):
    byY={}
    for adj in edges_iter(n):
        byY.setdefault(sig(adj,n),set()).add(sc4(adj,n))
    split=any(len(v)>1 for v in byY.values())
    print(f"   n={n}: SC4 determined by (score,c3)? {'NO (SPLIT)' if split else 'yes'}")

# confirm var(=tr4) fully determined by (score, SC4) [c3 is a score function, so (score,SC4)]
print("\nvar(=tr4) determined by (score, SC4) alone?")
for n in range(3,7):
    byY={}
    for adj in edges_iter(n):
        byY.setdefault((tuple(sorted(scores(adj,n))), sc4(adj,n)),set()).add(tr4(adj,n))
    print(f"   n={n}: {'yes' if all(len(v)==1 for v in byY.values()) else 'NO'}")
print("\n => var(lambda^2) is a 4-SUBTOURNAMENT-CENSUS invariant (score + SC4); it decouples from")
print("    (score,c3) exactly when SC4 does = n=5. The 'quaternion wall' = the degree-4 invariant")
print("    first escaping the degree-<=2 data.")
