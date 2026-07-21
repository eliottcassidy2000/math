#!/usr/bin/env python3
"""
CONCRETE NEXT STEPS (opus-S443):
 (1) Resolve tr(S^4) as an EXACT formula in the 4-subtournament census -> the finer invariant that
     determines var(lambda^2), and the 32-step insertion index (completes THM-1930).
 (2) Test the OCTONION-WALL direction: does joint degree-2 data (score,c3) push the decoupling
     threshold past n=5? (THM-1935 second wall.)

tr(S^4)=Sum lambda^4. Closed 4-walks in the +-1 skew matrix decompose by #distinct vertices:
2-vertex (const in n), 3-vertex (score/c3 terms), 4-vertex (the genuine 4-cycles = the strongly-
connected-4-subtournament count SC4). We fit tr(S^4) = a*n(n-1) + b*Sum C(s_i,2) + c*SC4 + d and
find exact integer coeffs, then re-run the decoupling matrix with the coarser Y determined.
"""
import itertools, numpy as np
from fractions import Fraction as F

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
    t=0
    for i,j,k in itertools.combinations(range(n),3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]): t+=1
    return t
def sc4(adj,n):
    """# of induced 4-subtournaments that are strongly connected (the (1,1,2,2) type)."""
    cnt=0
    for quad in itertools.combinations(range(n),4):
        # strongly connected iff every vertex has in- and out-degree >=1 within the quad AND connected;
        # for 4-tournaments: SC <=> score seq is (1,1,2,2) <=> #3-cycles among the 4 == 2
        c=0
        for i,j,k in itertools.combinations(quad,3):
            if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]): c+=1
        if c==2: cnt+=1     # (1,1,2,2) strongly-connected 4-tournament has exactly 2 cyclic triples
    return cnt

from math import comb
print("(1) tr(S^4) EXACT FORMULA fit against {n(n-1), sum C(s_i,2), SC4}")
print("="*70)
import numpy.linalg as la
for n in range(3,7):
    rows=[]; ys=[]
    for adj in edges_iter(n):
        s=scores(adj,n); sumCs2=sum(comb(v,2) for v in s)
        feat=[n*(n-1), sumCs2, sc4(adj,n), 1]
        rows.append(feat); ys.append(tr4(adj,n))
    A=np.array(rows,dtype=float); y=np.array(ys,dtype=float)
    # exact integer solve on a few independent rows, then verify on ALL
    coef,res,rk,sv=la.lstsq(A,y,rcond=None)
    coef_r=[round(c,4) for c in coef]
    pred=A@coef
    exact = np.allclose(pred,y,atol=1e-6)
    print(f" n={n}: tr(S^4) = {coef_r[0]}*n(n-1) + {coef_r[1]}*SumC(s_i,2) + {coef_r[2]}*SC4 + {coef_r[3]}   exact={exact}")

# global fit (all n pooled) to get n-independent coefficients
print("\n GLOBAL fit (all n=3..6 pooled):")
rows=[]; ys=[]
for n in range(3,7):
    for adj in edges_iter(n):
        s=scores(adj,n); sumCs2=sum(comb(v,2) for v in s)
        rows.append([n*(n-1), sumCs2, sc4(adj,n), 1]); ys.append(tr4(adj,n))
A=np.array(rows,dtype=float); y=np.array(ys,dtype=float)
coef,_,_,_=la.lstsq(A,y,rcond=None)
pred=A@coef
print(f"   tr(S^4) = {round(coef[0],3)}*n(n-1) + {round(coef[1],3)}*SumC(s_i,2) + {round(coef[2],3)}*SC4 + {round(coef[3],3)}")
print(f"   exact over ALL tournaments n=3..6: {np.allclose(pred,y,atol=1e-6)}")

# (2) octonion-wall: threshold of var (=tr4) determined by degree-1 vs degree-2 Y
print("\n" + "="*70)
print("(2) DECOUPLING THRESHOLD with coarser Y (does degree-2 push the wall past 5?)")
def sig_score(adj,n): return tuple(sorted(scores(adj,n)))
def sig_sc(a,n): return (sig_score(a,n), c3(a,n))              # (score,c3)
def sig_sc4(a,n): return (sig_score(a,n), c3(a,n), sc4(a,n))   # (score,c3,SC4)
Ys={'score':sig_score,'(score,c3)':sig_sc,'(score,c3,SC4)':sig_sc4}
data={}
for n in range(3,7):
    data[n]=[(tr4(a,n), {k:f(a,n) for k,f in Ys.items()}) for a in edges_iter(n)]
for Yname in Ys:
    row=[]; thr=None
    for n in range(3,7):
        byY={}
        for tr4v,sigs in data[n]:
            byY.setdefault(sigs[Yname],set()).add(tr4v)
        split=any(len(s)>1 for s in byY.values())
        row.append('SPLIT' if split else 'det')
        if split and thr is None: thr=n
    print(f"  var(=tr4) | {Yname:<16}: {row}   threshold={thr if thr else '>6'}")
print("\n READING: if tr(S^4)=affine in SC4, then var is a 4-vertex invariant; it decouples from Y")
print(" exactly when SC4 does. The threshold under (score,c3) vs (score,c3,SC4) locates the wall.")
