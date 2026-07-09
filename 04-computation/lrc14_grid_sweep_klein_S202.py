#!/usr/bin/env python3
"""
klein-2026-07-09-S202: adversarial sweep -- where does (GP) close the good period, and is it SOUND?
(GP): R_grid_abs < E[W] - (6/7)/Vmax  =>  good period exists (Sum_{j>=1}W(j/V)>0).
Keeps R0 (n.e=0, density floor) SIGNED; only R_grid absolute. Q: over HARD clusters (7*spread in
(6Vmax, 7Vmax], i.e. j=1 fails and non-degenerate), for which Vmax/spread does (GP) hold uniformly,
and does it EVER over-claim (fire when no good period)? Find Vmax0.
"""
import numpy as np
from math import gcd
from functools import reduce
rng=np.random.default_rng(202)
INV7=1/7
def Phi(E,N):
    E=np.array(sorted(set(E)),float); y=np.arange(N)/N
    ph=np.mod(np.outer(y,E),1.0); ph.sort(axis=1)
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    return np.maximum(g-INV7,0).sum(axis=1)
def stats(E,Vmax,L=64):
    N=L*Vmax; P=Phi(E,N); EW=P.mean(); F=np.fft.rfft(P)/N
    mmax=(N//2)//Vmax; Rga=2*sum(abs(F[m*Vmax]) for m in range(1,mmax+1))
    return EW,Rga
def realGP(E,Vmax):
    E=np.array(sorted(set(E)))
    for j in range(1,Vmax):
        p=np.unique((E*j)%Vmax); g=np.diff(p); g=np.append(g,Vmax-p[-1]+p[0])
        if g.max()>Vmax/7+1e-12: return True
    return False
def hard(E,Vmax):
    E=np.array(sorted(set(E))); p=E%Vmax; p.sort(); g=np.diff(p); g=np.append(g,Vmax-p[-1]+p[0])
    return g.max()<=Vmax/7+1e-9   # j=1 fails (maxgap<=V/7)
def genhard(mode):
    if mode==0: d=int(rng.integers(2,9)); E=sorted(set([d*i+int(rng.integers(-1,2)) for i in range(13)]))  # near-AP
    elif mode==1: E=sorted(set([7*int(rng.integers(0,14)) for _ in range(7)]+[int(rng.integers(1,98)) for _ in range(6)]))  # 7-struct
    elif mode==2: E=sorted(set(rng.choice(range(1,160),12,replace=False).tolist()+[0]))  # random/dissoc
    else: E=list(range(13))  # tight AP
    E=[e-min(E) for e in E]
    return E if len(E)==13 else None

print("SOUNDNESS + Vmax0 sweep (k=13 hard clusters). over_claim MUST be 0.\n")
buckets={}  # ratio-band -> [n, n_holds, n_realGP, n_overclaim]
overclaims=[]; fails_with_GP=[]
for _ in range(3000):
    E=genhard(int(rng.integers(0,4)))
    if E is None: continue
    sp=max(E)
    if reduce(gcd,E)!=1: continue
    for c in [1.02,1.1,1.3,1.6,2.0,3.0]:
        Vmax=int(sp*c)+1
        if Vmax<=sp or not hard(E,Vmax): continue
        EW,Rga=stats(E,Vmax); rhs=EW-(6/7)/Vmax; holds=Rga<rhs; real=realGP(E,Vmax)
        band=round(Vmax/sp,1)
        b=buckets.setdefault(band,[0,0,0,0]); b[0]+=1; b[1]+=holds; b[2]+=real; b[3]+=(holds and not real)
        if holds and not real: overclaims.append((E,Vmax,Rga,rhs))
        break
print(f"{'Vmax/sp':>7} {'#hard':>6} {'#(GP)holds':>10} {'#realGP':>8} {'#overclaim':>10} {'GP-cover%':>9}")
for band in sorted(buckets):
    n,nh,nr,noc=buckets[band]
    print(f"{band:>7} {n:>6} {nh:>10} {nr:>8} {noc:>10} {100*nh/max(n,1):>8.1f}%")
print(f"\nTOTAL over-claims (GP fires, NO good period) = {len(overclaims)}  (MUST be 0 for soundness)")
for E,V,r,rhs in overclaims[:4]: print(f"  OVERCLAIM E={E} V={V} Rga={r:.4f} rhs={rhs:.4f}")
print("\nReading: GP-cover% -> 100 as Vmax/sp grows = Vmax0 boundary. Below it, no-GP clusters (tight-AP-like)")
print("correctly excluded (exact-check territory). Soundness = 0 over-claims (the tight AP is respected).")
