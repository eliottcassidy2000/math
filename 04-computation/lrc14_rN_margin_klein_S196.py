#!/usr/bin/env python3
"""
klein-2026-07-09-S196: characterize opus-S165's r_N < 1 (the LRC(14) capstone).

opus-S165: j* <= N  <=>  S_N := sum_{j=1}^N W(j/Vmax) > 0,  N = ceil(7(k-1)/6).
S_N = N(6/7)^k + Corr_N,  r_N = |Corr_N|/(N(6/7)^k).  r_N < 1 => S_N>0 => j*<=N.
Exact complete-residue AP (E={0..k-1}, Vmax=k) is the r_N=1 boundary (S_N=0, tight LRC).

Questions:
 (1) Is j* <= N ALWAYS (for Vmax>N)?  Is the margin 1-r_N bounded below?
 (2) What is the worst (max r_N / min margin) hard cluster, and is it the near-AP?
 (3) RIGIDITY: is S_N=0 (no good period in 1..N) ONLY the tight AP (Vmax<=N)?
"""
import numpy as np
from math import ceil, gcd
rng=np.random.default_rng(196)
INV7=1/7

def W(x,E,Vmax):
    ph=np.sort((np.array(E)*x)%Vmax)/Vmax
    g=np.diff(ph); g=np.append(g,1-ph[-1]+ph[0])
    return np.maximum(g-INV7,0).sum()
def jstar(E,Vmax,Jmax=None):
    Jmax=Jmax or Vmax
    for j in range(1,Jmax+1):
        if W(j,E,Vmax)>1e-12: return j
    return None
def SN(E,Vmax,N): return sum(W(j,E,Vmax) for j in range(1,N+1))
def rN(E,Vmax,N,k):
    main=N*(6/7)**k; return abs(SN(E,Vmax,N)-main)/main
def sample_hard(k,Vmax,n,tries=50000):
    out=[];t=0
    while len(out)<n and t<tries:
        t+=1
        rest=rng.choice(np.arange(1,Vmax),k-1,replace=False)
        E=tuple(sorted([0]+[int(x) for x in rest]))
        if max(E)<6*Vmax/7: continue
        if jstar(E,Vmax,1) is not None: continue  # hard = j=1 fails
        out.append(E)
    return out

print("(1)+(2) max r_N (min margin) over hard clusters; is j*<=N always?")
print(f"{'k':>3} {'N':>3} {'Vmax':>6} {'#hard':>6} {'max j*':>7} {'j*<=N?':>7} {'max r_N':>8} {'min margin':>10} {'worst near-AP?':>13}")
for k in (8,9,10,11,12,13):
    N=ceil(7*(k-1)/6)
    for Vmax in (2*k, 3*k, 91, 300):
        if Vmax<=N: continue
        hs=sample_hard(k,Vmax,200)
        if not hs: continue
        js=[jstar(E,Vmax) for E in hs]
        mxj=max(j for j in js if j)
        rs=[(rN(E,Vmax,N,k),E) for E in hs]
        mr,wE=max(rs)
        # evenness of worst
        ph=np.sort(np.array(wE)%Vmax); g=np.diff(ph); g=np.append(g,Vmax-ph[-1]+ph[0])
        ev=g.max()/g.mean()
        print(f"{k:>3} {N:>3} {Vmax:>6} {len(hs):>6} {mxj:>7} {str(mxj<=N):>7} {mr:>8.3f} {1-mr:>10.3f} {ev:>13.3f}")

print("\n(3) RIGIDITY: the exact complete-residue AP E={0..k-1} at small Vmax -- does S_N=0 (no good period)?")
for k in (8,11,13):
    N=ceil(7*(k-1)/6)
    E=tuple(range(k))
    for Vmax in range(k, N+3):
        js=jstar(E,Vmax)
        s=SN(E,Vmax,N)
        tag="TIGHT (no good j in 1..N)" if (js is None or js>N) else f"good at j*={js}"
        print(f"  k={k} N={N} E={{0..{k-1}}} Vmax={Vmax}: S_N={s:.4f}  {tag}")
