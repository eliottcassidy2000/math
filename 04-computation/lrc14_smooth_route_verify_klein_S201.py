#!/usr/bin/env python3
"""
klein-2026-07-09-S201: rigorously check opus-S170's SMOOTH-AVERAGING good-period route
at the razor. A good period exists if E_grid[maxgap] > 1/7 (max>=mean). opus:
  E_grid[maxgap] = E_x[maxgap] + discrepancy,  E_x[maxgap]>1/7 (min ~1.047x=0.150),
  |discrepancy| <= 0.006 (maxgap 1-D piecewise-linear => Fourier alpha=2, abs. conv.).
Net margin ~0.0007 -- RAZOR. Is it really positive for ALL hard clusters, or does some
adversarial (coarse-cluster / 7-structured / dilated-AP) cluster break E_grid[maxgap]>1/7?

Compute E_x[maxgap] (fine-grid integral) and E_grid[maxgap] (exact ruler mean) and the
good-period condition, over adversarial k=13 families + a min-hunt.
"""
import numpy as np
rng=np.random.default_rng(201)
def maxgap_x(E,Nx):  # maxgap(x) sampled; returns array
    E=np.array(E); x=(np.arange(Nx)+0.5)/Nx
    ph=np.sort(np.mod(np.outer(x,E),1.0),axis=1)
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    return g.max(axis=1)
def Ex_maxgap(E,Nx=2000000): return maxgap_x(E,Nx).mean()
def Egrid_maxgap(E,Vmax):
    E=np.array(E); j=np.arange(Vmax); x=j/Vmax
    ph=np.sort(np.mod(np.outer(x,E),1.0),axis=1)
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    return g.max(axis=1).mean()
def hard(E,V):
    p=np.sort(np.mod(np.array(E),V)); g=np.diff(p); g=np.append(g,V-p[-1]+p[0]); return g.max()<=V/7+1e-9

INV7=1/7
print(f"1/7 = {INV7:.5f};  E_grid[maxgap] > 1/7  <=>  a good period exists")
print(f"{'family':>16} {'E_x[maxgap]':>11} {'x/(1/7)':>8}")
fams={
 'AP{1..13}': list(range(1,14)),
 'primsat 2{1..12}+13': [2*i for i in range(1,13)]+[13],
 'dilAP 13*{0..12}': [13*i for i in range(13)],
 '7-struct(S200)': [0,7,14,21,26,29,37,44,51,58,67,75,82],
 'block{0..12}': list(range(13)),
 '2block': [0,1,2,3,4,5,100,101,102,103,104,105,106],
}
for nm,E in fams.items():
    ex=Ex_maxgap([e-min(E) for e in E]); print(f"{nm:>16} {ex:>11.5f} {ex/INV7:>8.3f}")

print("\nMIN-HUNT E_x[maxgap] over hard k=13 (structured/coarse adversaries) + E_grid check:")
worst=(9,None)
for _ in range(4000):
    mode=rng.integers(0,3)
    if mode==0:  # coarse: a*{small set} + breakers
        a=int(rng.integers(2,5)); base=sorted(rng.choice(range(1,14),rng.integers(6,10),replace=False))
        E=sorted(set([a*b for b in base]+[int(rng.integers(1,3*a*max(base)))for _ in range(13-len(base))]))
    elif mode==1:  # near-AP dilated + perturb
        d=int(rng.integers(2,20)); E=sorted(set([d*i+int(rng.integers(-2,3)) for i in range(13)]))
    else:  # 7-structured
        E=sorted(set([7*int(rng.integers(0,20)) for _ in range(7)]+[int(rng.integers(1,140)) for _ in range(6)]))
    E=[e-min(E) for e in E]
    if len(E)!=13: continue
    ex=Ex_maxgap(E,Nx=400000)
    if ex<worst[0]: worst=(ex,E)
ex,E=worst
print(f"  min E_x[maxgap] found = {ex:.5f} = {ex/INV7:.3f}x (1/7)   E={E}")
# check E_grid at its hard rulers
sp=max(E); nfail=0; mn=9
for V in range(sp+1, int(7*sp/6)+1):
    p=np.sort(np.mod(np.array(E),V)); g=np.diff(p); g=np.append(g,V-p[-1]+p[0])
    if g.max()>V/7+1e-9: continue  # not hard
    eg=Egrid_maxgap(E,V); mn=min(mn,eg)
    if eg<=INV7: nfail+=1
print(f"  min E_grid[maxgap] over hard rulers = {mn:.5f} = {mn/INV7:.3f}x   #(E_grid<=1/7, route FAILS)={nfail}")
print(f"\n=> if min E_grid[maxgap] > 1/7 for all hard adversaries, opus's smooth route holds (thin).")
print(f"   the a-priori pieces: E_x[maxgap] floor (>1/7) + discrepancy |E_grid-E_x| bound.")
