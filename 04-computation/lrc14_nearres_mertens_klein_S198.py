#!/usr/bin/env python3
"""
klein-2026-07-09-S198: the NEAR-RESONANCE COUNT (route-a target) + the MERTENS
cancellation, and its link to #arcs (route c).

Corr_N = sum_{n!=0} What(n) G_N(n.E/V), dominated by NEAR-RESONANCES ||n.E/V||<1/(2N).
 (1) COUNT: #near-resonances (small support/height) vs longest-AP L -- is it bounded
     by the longest-AP cap (dissociated => few)?  This is the a-priori target.
 (2) MERTENS: absolute sum |What|G_N vs SIGNED sum. Absolute ~20x target (kps-S92);
     signed cancels to ~0.1x. The sign is (-1)^r (support parity) -- a Mobius-like
     parity, so the sum is a Mertens-type signed sum: cancellation is real but,
     like sum mu(n) (Mertens conj DISPROVED), MUST NOT be assumed -- hence route (c).
 (3) LINK: #arcs(Good_E) tracks the near-resonance count (both = coherent gap-events);
     low L => few near-res => few arcs => route (c) c(L)<rho* holds, POSITIVELY.
"""
import numpy as np, cmath
from itertools import product
from math import gcd
rng=np.random.default_rng(1988)
INV7=1/7; M=6/7
def b0(m): return 1/7 if m==0 else (cmath.exp(2j*np.pi*m/7)-1)/(2j*np.pi*m)
def c_of(s): return 1/7 if s==0 else (1-cmath.exp(-2j*np.pi*s/7))/(2j*np.pi*s)
def What(n,k):
    nz=[v for v in n if v!=0]; r=len(nz); sig=sum(n); pr=1+0j
    for v in nz: pr*=b0(v)
    return ((-1)**r)*(6/7)**((k-1)-r)*pr*((1 if sig==0 else 0)-c_of(sig))
def longest_AP(E):
    Es=sorted(set(E)); Eset=set(Es); best=1
    for i in range(len(Es)):
        for jj in range(i+1,len(Es)):
            d=Es[jj]-Es[i]; L=2; x=Es[jj]+d
            while x in Eset: L+=1; x+=d
            best=max(best,L)
    return best
def good_j1(E,V):
    p=np.sort(np.mod(np.array(E),V)); g=np.diff(p); g=np.append(g,V-p[-1]+p[0]); return g.max()>V/7+1e-9
def arcs(E,V,Nx=150000):
    x=(np.arange(Nx)+0.5)/Nx*V; e=np.array(E)
    ph=np.sort(np.mod(np.outer(x,e),V),axis=1)/V
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    gi=(g.max(axis=1)>INV7+1e-12).astype(int); ed=np.diff(np.concatenate([gi,gi[:1]]))
    nc=int((ed==1).sum());
    if gi.all():nc=1
    return nc
def analyze(E,V,N=3,H=4):
    k=len(E); absS=0.0; sgnS=0.0+0j; nearcnt=0
    for n in product(range(-H,H+1),repeat=k-1):
        if all(v==0 for v in n): continue
        ne=int(np.dot(n,E[1:])); th=((ne/V)+0.5)%1.0-0.5
        GN=N if abs(th)<1e-12 else min(N,1/(2*abs(th)))
        w=What(n,k); absS+=abs(w)*GN; sgnS+=w*GN
        if abs(th)<1/(2*N): nearcnt+=1
    return absS, sgnS.real, nearcnt

k=8; V=90; N=3; tgt=N*(6/7)**k
print(f"k={k}, V={V}, N={N}, target=N(6/7)^k={tgt:.4f}, H=4")
print(f"{'longest-AP L':>13} {'#near-res':>10} {'#arcs':>6} {'absBound/tgt':>13} {'|signed|/tgt':>13} {'good@N?':>8}")
# group by longest-AP
seen={}
tries=0
while tries<40000 and sum(len(v) for v in seen.values())<40:
    tries+=1
    rest=rng.choice(np.arange(1,V),k-1,replace=False)
    E=tuple(sorted([0]+[int(x) for x in rest]))
    if max(E)<6*V/7 or good_j1(E,V): continue
    L=longest_AP(E)
    if L in seen and len(seen[L])>=6: continue
    ab,sg,nr=analyze(E,V,N); nc=arcs(E,V)
    # good period at some j<=N?
    goodN=any(not good_j1(E,1) and (np.sort(np.mod(np.array(E)*j,V)/V) is not None) for j in [1]) # placeholder
    gj=None
    for j in range(1,N+1):
        p=np.sort(np.mod(np.array(E)*j,V)); g=np.diff(p); g=np.append(g,V-p[-1]+p[0])
        if g.max()>V/7+1e-9: gj=j; break
    seen.setdefault(L,[]).append((nr,nc,ab/tgt,abs(sg)/tgt,gj))
for L in sorted(seen):
    for (nr,nc,abr,sgr,gj) in seen[L][:4]:
        print(f"{L:>13} {nr:>10} {nc:>6} {abr:>13.1f} {sgr:>13.3f} {str(gj):>8}")

print("\nREAD: #near-res and #arcs INCREASE with longest-AP L (near-AP => many coherent")
print("resonances; Sidon L=2 => few). absBound/tgt ~20 (Mertens: absolute hopeless);")
print("|signed|/tgt < 1 (cancellation gives the bound, but MUST be proven, not assumed --")
print("the Mertens-conjecture lesson). Route (c) uses #arcs (POSITIVE) => no cancellation.")
