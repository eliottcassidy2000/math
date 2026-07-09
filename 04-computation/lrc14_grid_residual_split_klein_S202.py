#!/usr/bin/env python3
"""
klein-2026-07-09-S202: transport the OCF/decay truncation constant to an a-priori good-period bound.

opus-S171 claimed R_abs = Sum_{Vmax|n.e} 2|What(n)| < (6/7)^k for ALL clusters. kps-S98 REFUTED
uniformity: the TOTAL absolute sum exceeds (6/7)^k at small/mid spread (1.55@s50) -- essential
cancellation. kps reverted to LEM-013.

THE REFINEMENT (this script): kps took EVERYTHING absolute. Split instead:
  R = E_grid[W]-(6/7)^k = R0_signed + R_grid,
  R0_signed = E[W]-(6/7)^k   (the n.e=0 EXACT relations = the DENSITY FLOOR; keep SIGNED, it's bounded
                              because E[W]>=0, in fact >= a positive floor)
  R_grid    = E_grid[W]-E[W]  (the PURE GRID resonances n.e=mVmax, m!=0; bound ABSOLUTELY)
Then E_grid[W] = E[W] + R_grid >= E[W] - R_grid_abs. So a good period exists (E_grid[W]>0) once
  ** R_grid_abs < E[W] **  -- keeping R0 signed SIDESTEPS mac-mini's/kps's cancellation wall (which
was in the n.e=0 exact-relation part, the AP decorrelation). R_grid_abs -> 0 as Vmax grows regardless
of cluster (grid resonances pushed to high height => LEM-011 geometric decay). Q: does it, how fast,
and is R_grid_abs < E[W] for Vmax>Vmax0 including the TIGHT AP (opus/kps's blind spot)?

Phi(y) = uncovered measure of {frac(e_i y)}, W = Sum(gap-1/7)_+. E[W]=Phi_hat(0); E_grid = mean over
y=j/Vmax; R_grid_abs = Sum_{m>=1} 2|Phi_hat(m Vmax)| via FFT of Phi at resolution N=L*Vmax.
"""
import numpy as np
INV7=1/7
def Phi_samples(E,N):
    E=np.array(sorted(set(E)),dtype=float); y=np.arange(N)/N
    ph=np.mod(np.outer(y,E),1.0); ph.sort(axis=1)
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    return np.maximum(g-INV7,0).sum(axis=1)   # Phi(y) = W at shift y
def analyze(E,Vmax,L=96):
    k=len(set(E)); N=L*Vmax
    Phi=Phi_samples(E,N)
    EW=Phi.mean()                              # continuum E[W] (approx, fine grid)
    Egrid=Phi[::L].mean()                      # mean over y=j/Vmax  (j=0..Vmax-1)
    F=np.fft.rfft(Phi)/N                        # Phi_hat(M)=F[M]
    mmax=(N//2)//Vmax
    Rgrid_abs=2*sum(abs(F[m*Vmax]) for m in range(1,mmax+1))
    main=(6/7)**k
    return dict(EW=EW,Egrid=Egrid,Rgrid_signed=Egrid-EW,Rgrid_abs=Rgrid_abs,
                R0_signed=EW-main,R_total=Egrid-main,main=main)

clusters={
 'tightAP {0..12}':      list(range(13)),
 'near-AP 3*{0..9}+int': [0,3,6,8,9,12,15,18,21,24,27][:11]+[30,31],  # k=13 near-AP
 '7-struct(MISTAKE128)': [0,7,14,21,26,29,37,44,51,58,67,75,82],
 'dissociated B2-ish':   [0,1,3,7,12,20,30,44,65,80,96,112,129],
}
print("R_grid_abs < E[W] ?  (good period exists if yes)   main=(6/7)^13=%.5f\n"%((6/7)**13))
for nm,E in clusters.items():
    sp=max(E)-min(E)
    print(f"[{nm}] spread={sp}")
    print(f"  {'Vmax':>5} {'E[W]':>8} {'R_grid_abs':>11} {'ratio g/EW':>10} {'R_grid*Vmax':>11} {'Rtot/main':>9} {'Egrid>0?':>8}")
    for Vmax in [sp+1, int(sp*1.07)+1, 2*sp, 4*sp, 8*sp]:
        if Vmax<=sp: continue
        d=analyze([e-min(E) for e in E],Vmax)
        rat=d['Rgrid_abs']/d['EW'] if d['EW']>1e-12 else float('inf')
        print(f"  {Vmax:>5} {d['EW']:>8.5f} {d['Rgrid_abs']:>11.5f} {rat:>10.3f} {d['Rgrid_abs']*Vmax:>11.3f} {d['R_total']/d['main']:>9.3f} {str(d['Egrid']>1e-9):>8}")
    print()
print("KEY: if R_grid_abs < E[W] for Vmax>Vmax0 (ratio<1), keeping R0 signed CLOSES the good period")
print("a-priori where kps's total-absolute failed. R_grid*Vmax ~ const => R_grid_abs = O(1/Vmax) transport.")
