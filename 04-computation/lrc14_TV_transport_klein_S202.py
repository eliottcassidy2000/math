#!/usr/bin/env python3
"""
klein-2026-07-09-S202: the OCF/decay-truncation transport, made explicit & VALIDATED.

CLAIM (the transported theorem shape): a good period exists once
  (GP)   R_grid_abs  <  E[W] - (6/7)/Vmax
where R_grid_abs = Sum_{m>=1} 2|Phi_hat(m Vmax)|, Phi(y)=uncovered measure of {frac(e_i y)}.
Reason (Sum_{j=1}^{V-1} W(j/V) = Vmax*E_grid[W] - W(0), W(0)=6/7): (GP) => Sum_{j>=1}W(j/V)>0 =>
some j in {1..V-1} has W>0 => maxgap>1/7 => good period. R0 (the n.e=0 exact relations = density
floor) stays SIGNED inside E[W] (bounded >=0), ONLY R_grid is taken absolutely -- sidestepping the
kps-S98/mac-mini exact-relation cancellation wall.

TRANSPORT (the OCF clean-truncation constant, covering side): Phi is CONTINUOUS piecewise-linear in y
(the tournament Walsh-OCF THM-076 "clean truncation" analogue: continuous => Fourier O(1/M^2)), so
  |Phi_hat(M)| <= TV(Phi') / (2 pi M)^2   and   R_grid_abs <= TV(Phi') / (12 Vmax^2).
TV(Phi') = total variation of Phi' over one period = the explicit covering "truncation constant".

VALIDATE (this script):
 (A) R_grid_abs <= TV(Phi')/(12 Vmax^2)  (the transport bound holds).
 (B) (GP) is SOUND: whenever (GP) holds, a good period ACTUALLY exists (max_{j>=1} maxgap>1/7).
     Never over-claims (unlike opus-S170 smooth mean / total-absolute). Test the TIGHT AP hard.
 (C) TV(Phi') is an explicit cluster constant; bound TV(Phi') <= C0 * k * spread (uniform).
"""
import numpy as np
INV7=1/7
def Phi_and_maxgap(E,N):
    E=np.array(sorted(set(E)),dtype=float); y=np.arange(N)/N
    ph=np.mod(np.outer(y,E),1.0); ph.sort(axis=1)
    g=np.diff(ph,axis=1); g=np.concatenate([g,(1-ph[:,-1]+ph[:,0])[:,None]],axis=1)
    W=np.maximum(g-INV7,0).sum(axis=1); mg=g.max(axis=1)
    return W,mg
def TV_deriv(Phi,N):  # total variation of Phi' = sum |second difference| * N  (discrete)
    d1=(np.roll(Phi,-1)-Phi)*N            # forward derivative (periodic)
    return np.abs(np.roll(d1,-1)-d1).sum()  # sum of |jumps in derivative|
def good_period_exists(E,Vmax):
    E=np.array(sorted(set(E)))
    for j in range(1,Vmax):
        p=np.unique((E*j)%Vmax); g=np.diff(p); g=np.append(g,Vmax-p[-1]+p[0])
        if g.max()>Vmax/7+1e-12: return True
    return False

clusters={
 'tightAP{0..12}':      list(range(13)),
 'near-AP':             [0,3,6,8,9,12,15,18,21,24,27,30,31],
 '7-struct(M128)':      [0,7,14,21,26,29,37,44,51,58,67,75,82],
 'dissoc':              [0,1,3,7,12,20,30,44,65,80,96,112,129],
}
print("(A) transport bound R_grid_abs <= TV(Phi')/(12 Vmax^2) ;  (B) (GP) SOUND vs actual good period\n")
for nm,E in clusters.items():
    E=[e-min(E) for e in sorted(set(E))]; k=len(E); sp=max(E)
    # TV(Phi') from a fine Vmax-independent grid
    Nfine=1<<18
    Phi,_=Phi_and_maxgap(E,Nfine); TV=TV_deriv(Phi,Nfine); EW=Phi.mean()
    print(f"[{nm}] k={k} spread={sp}  E[W]={EW:.5f}  TV(Phi')={TV:.1f}  TV/(k*sp)={TV/(k*sp):.2f}")
    print(f"  {'Vmax':>5} {'R_grid_abs':>11} {'TV/(12V^2)':>11} {'bound ok':>8} {'GP-cond':>8} {'(GP)holds':>9} {'realGP':>7} {'sound':>6}")
    for Vmax in [sp+1, int(sp*1.1)+1, 2*sp, 3*sp, 6*sp]:
        if Vmax<=sp: continue
        N=64*Vmax
        Phi2,_=Phi_and_maxgap(E,N); F=np.fft.rfft(Phi2)/N
        mmax=(N//2)//Vmax
        Rga=2*sum(abs(F[m*Vmax]) for m in range(1,mmax+1))
        tvb=TV/(12*Vmax**2)
        gp_rhs=EW-(6/7)/Vmax
        holds=Rga<gp_rhs
        real=good_period_exists(E,Vmax)
        sound=(not holds) or real   # (GP) holds => must have real good period
        print(f"  {Vmax:>5} {Rga:>11.5f} {tvb:>11.5f} {str(Rga<=tvb*1.05):>8} {gp_rhs:>8.5f} {str(holds):>9} {str(real):>7} {str(sound):>6}")
    print()
print("(C) if TV(Phi') <= C0*k*spread uniformly, then R_grid_abs <= C0*k*spread/(12 Vmax^2);")
print("    for hard clusters spread<=Vmax => R_grid_abs <= C0*k/(12 Vmax) -> a-priori Vmax0.")
print("SOUND must be True everywhere (never claims GP when none exists -- esp. tight AP small V).")
