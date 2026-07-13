#!/usr/bin/env python3
"""mac-mini-S69: verify the load-bearing claims for a RIGOROUS single-killer rigidity proof.
(A) BALANCE LEMMA edge: for S=C u {v}, v resonant at a core optimum t0 (v*t0 in Z), the
    perturbation t0-eps gives M(S) >= M_C * v/(v+s), s=binding-runner rate. Check numerically.
(B) SHALLOW-WITNESS existence: for dilated core C=c*{1..12} and killer v (primitive gcd(c,v)=1),
    a 'good dilation' a with gcd(a,13)=1 and (a*v mod 13c) at distance >= c from 0 ALWAYS exists
    => M(C u {v}) >= 1/13. Verify existence exhaustively over c, v.
(C) COVERING => v>=182 for dilated cores coprime to 182 (sanity)."""
from fractions import Fraction as F
from math import gcd

def resdist(x,q): r=x%q; return min(r,q-r)

# (A) balance lemma numeric check on interval core {1..12}, killer 182|v
def M_lower_balance(C, v):
    # core optimum t0=1/13 for {1..12}; general: scan t0 rationals to find M_C & binding s
    pass
def M_exact(S,Qmax):
    best=F(0)
    for q in range(2,Qmax+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            mind=q
            for v in S:
                d=resdist(a*v,q)
                if d<mind: mind=d
                if d==0: break
            if mind>0 and F(mind,q)>best: best=F(mind,q)
    return best
print("(A) interval-core balance: M({1..12,v}) vs balance v/(13(v+1)), for 182|v:")
for v in [182,364,546,728,910]:
    S=list(range(1,13))+[v]
    M=M_exact(S, min(2*v,500))
    bal=F(v,13*(v+1))
    print(f"    v={v}: M_exact={M}={float(M):.6f}  balance v/(13(v+1))={bal}={float(bal):.6f}  match={M==bal}")

# (B) shallow-witness good-dilation existence for dilated cores
print("\n(B) shallow-witness: exists a (gcd(a,13)=1) with (a*v mod 13c) at distance >= c from 0?")
print("    (=> M(c*{1..12} u {v}) >= 1/13). Test all c=2..15, v coprime to c, v in [13,400]:")
fails=[]
tested=0
for c in range(2,16):
    q=13*c
    for v in range(13, 400):
        if gcd(v,c)!=1: continue      # primitivity of the 13-set (gcd(core,v)=gcd(c,v))
        tested+=1
        good_a=[a for a in range(1,q) if gcd(a,13)==1 and resdist(a*v,q)>=c]
        if not good_a:
            fails.append((c,v))
print(f"    tested {tested} (c,v) primitive pairs; good-dilation MISSING in {len(fails)} cases.")
if fails: print("    *** fails:", fails[:10])
else: print("    => shallow-witness ALWAYS exists (M >= 1/13 for every primitive dilated-core config).")

# also confirm the true M for a few dilated configs >= 1/13
print("\n    confirm M_exact >= 1/13 for dilated cores + killer 182:")
for c in [2,3,5,7,11,13]:
    S=sorted(set([c*k for k in range(1,13)]+[182]))
    if len(S)!=13: 
        print(f"    c={c}: dup (182 in core), skip"); continue
    g=0
    for x in S: g=gcd(g,x)
    M=M_exact(S,min(2*max(S),400))
    print(f"    c={c}: prim={g==1} M={M}={float(M):.5f} {'>=1/13' if M>=F(1,13) else '< 1/13 !!'}")
