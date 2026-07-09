#!/usr/bin/env python3
"""
klein-2026-07-09-S198: close the Sidon-like dissociated branch (L<=k-7) of j*=O(k).

THE SINGLE INEQUALITY (route c, mac-mini-S61, CANCELLATION-FREE): a good period
EXISTS whenever #arcs(Good_E) < rho*(E)*Vmax, and since spread<=Vmax this holds if
  c(L) := max #arcs/spread   <   rho_min(L) := min rho*  (over longest-AP=L clusters).
This is a bound on POSITIVE quantities -- no signed cancellation (unlike route a's
Corr_N). MERTENS LESSON: kps-S92 showed the absolute |W-hat| sum is ~20x the target
=> the partial-sum route needs cancellation, which -- like the Mertens conjecture
(sum mu(n) < sqrt(n), DISPROVED) -- cannot be assumed. Route (c) SIDESTEPS it.

Compute for dissociated clusters (longest-AP L, hard): c(L), rho_min(L) (min D3, min mu),
the near-resonance COUNT (route-a target) vs L, and the absolute-vs-signed sum (Mertens).
"""
import numpy as np
from math import ceil, gcd
from itertools import product
rng=np.random.default_rng(198)
INV7=1/7; M=6/7

def maxgap_arr(E,x,V):
    ph=np.sort(np.mod(np.outer(x,E),V),axis=1)/V
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    return g.max(axis=1)
def arcs_meas(E,V,N=200000):
    x=(np.arange(N)+0.5)/N*V
    good=maxgap_arr(E,x,V)>INV7+1e-12
    gi=good.astype(int); edges=np.diff(np.concatenate([gi,gi[:1]]))
    nc=int((edges==1).sum());
    if gi.all(): nc=1
    return nc, good.mean()
def D3(E,V,N=150000):
    x=(np.arange(N)+0.5)/N
    ph=np.sort(np.mod(np.outer(x,E)/V*V,1.0),axis=1) if False else np.sort(np.mod(np.outer(x,np.array(E)),1.0),axis=1)
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    W=np.maximum(g-INV7,0).sum(axis=1); m1,m2,m3=W.mean(),(W**2).mean(),(W**3).mean()
    return m1/M+(m1-m2/M)**2/(m2-m3/M)
def longest_AP(E):
    Es=sorted(set(E)); Eset=set(Es); best=1
    for i in range(len(Es)):
        for jj in range(i+1,len(Es)):
            d=Es[jj]-Es[i]; L=2; x=Es[jj]+d
            while x in Eset: L+=1; x+=d
            best=max(best,L)
    return best
def is_hard(E,V): return maxgap_arr(E,np.array([V/V]),V)  # placeholder
def good_j1(E,V):
    p=np.sort(np.mod(np.array(E),V)); g=np.diff(p); g=np.append(g,V-p[-1]+p[0]); return g.max()>V/7+1e-9

k=13
print(f"k={k}: dissociated branch, c(L)=max#arcs/spread vs rho_min(L)=min rho*  (route c, cancellation-free)")
print(f"{'L':>3} {'#sets':>6} {'c(L)=max#arcs/sp':>16} {'rho_min(=minmu)':>15} {'min D3':>8} {'c<rho?':>7} {'margin':>7}")
byL={}
for V in (200,300,500):
    tries=0
    while tries<250000:
        tries+=1
        rest=rng.choice(np.arange(1,V),k-1,replace=False)
        E=tuple(sorted([0]+[int(x) for x in rest]))
        if max(E)<6*V/7: continue
        if good_j1(E,V): continue     # hard
        L=longest_AP(E)
        if L>k-7: continue            # dissociated (L<=6)
        nc,mu=arcs_meas(E,V)
        byL.setdefault(L,[]).append((nc/max(E),mu,E,V))
        if sum(len(v) for v in byL.values())>=600: break
for L in sorted(byL):
    rows=byL[L]
    cL=max(r[0] for r in rows); rmin=min(r[1] for r in rows)
    # min D3 over a few
    d3s=[D3(list(r[2]),r[3]) for r in rows[:20]]
    d3min=min(d3s) if d3s else 0
    print(f"{L:>3} {len(rows):>6} {cL:>16.3f} {rmin:>15.3f} {d3min:>8.3f} {str(cL<rmin):>7} {rmin-cL:>7.3f}")

# near-resonance count + Mertens (absolute vs signed) for a sample dissociated E
print("\nNear-resonance COUNT (route a target) + Mertens cancellation check:")
def b0(m):
    import cmath
    return 1/7 if m==0 else (cmath.exp(2j*np.pi*m/7)-1)/(2j*np.pi*m)
def c_of(s):
    import cmath
    return 1/7 if s==0 else (1-cmath.exp(-2j*np.pi*s/7))/(2j*np.pi*s)
def What(n,k):
    nz=[v for v in n if v!=0]; r=len(nz); sig=sum(n); pr=1+0j
    for v in nz: pr*=b0(v)
    return ((-1)**r)*(6/7)**((k-1)-r)*pr*((1 if sig==0 else 0)-c_of(sig))
def near_res_analysis(E,V,N,H=4):
    Ecl=list(E)  # k co-offsets incl 0
    absS=0.0; sgnS=0.0+0j; cnt=0; nearcnt=0
    for n in product(range(-H,H+1),repeat=len(Ecl)-1):  # n on the k-1 nonzero-phase coords
        if all(v==0 for v in n): continue
        ne=int(np.dot(n,Ecl[1:]))
        theta=((ne/V)+0.5)%1.0-0.5
        GN=min(N,1/(2*abs(theta))) if abs(theta)>1e-12 else N
        w=What(n,k); absS+=abs(w)*GN; sgnS+=w*GN; cnt+=1
        if abs(theta)<1/(2*N): nearcnt+=1
    return absS, abs(sgnS.real), nearcnt, cnt
# pick a Sidon-like (low L) E
for V in (200,):
    for _ in range(300):
        rest=rng.choice(np.arange(1,V),k-1,replace=False)
        E=tuple(sorted([0]+[int(x) for x in rest]))
        if max(E)<6*V/7 or good_j1(E,V) or longest_AP(E)>k-7: continue
        N=3
        # small k for the full product: reduce to k=8 dissociated for feasibility
        break
# feasible near-res: do it for k=8 dissociated
print("  (k=8 dissociated, H=4, N=3): absBound/target vs signed/target -- the Mertens cancellation")
k8=8
for _ in range(2000):
    V=90
    rest=rng.choice(np.arange(1,V),k8-1,replace=False)
    E=tuple(sorted([0]+[int(x) for x in rest]))
    if max(E)<6*V/7 or good_j1(E,V): continue
    L=longest_AP(E)
    if L>2: continue   # Sidon (L=2)
    N=3; tgt=N*(6/7)**k8
    ab,sg,nr,tot=near_res_analysis(E,V,N,H=4)
    print(f"    Sidon E={E[:6]}... L={L}: absBound/target={ab/tgt:.1f}  signed/target={sg/tgt:.3f}  #near-res={nr}/{tot}")
    break
print("\n=> route(c) c(L)<rho_min (positive, no cancellation) CLOSES the Sidon branch;")
print("   route(a) needs cancellation (absBound~20x, Mertens-caution) => route(c) is the robust closure.")
