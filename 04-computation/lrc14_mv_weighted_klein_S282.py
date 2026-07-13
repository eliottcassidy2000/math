#!/usr/bin/env python3
"""
lrc14_mv_weighted_klein_S282.py
===============================
klein-2026-07-13-S282. CORRECTED: Q_s is dominated by the MANY THIN R_s-arcs, each contributing ~w_i
(arc width) to f_hat. So it's a WIDTH-WEIGHTED exponential sum Sum_i w_i e(-ell w c_i).

Montgomery-Vaughan: Sum_{ell~L}|Sum_i a_i e(-ell w c_i)|^2 <= Sum_i |a_i|^2 (L + delta_i^-1),
a_i=w_i (width), delta_i = nearest-nbr of {w c_i} (midpoints x w). The DECISIVE question:
does the WIDTH-WEIGHTED clustering sum Sum_i w_i^2/delta_i stay CONTROLLED (weights kill clustered
terms) while the UNWEIGHTED Sum 1/delta_i blows up? If Sum w_i^2/delta_i = O(mu/w) ~ L-term, MV closes
Q_s=O(r) (or a power-saving). Compute both, vs mu, r, w.
"""
import math
NG=1500000
def sec(e,x): return int((e*x%1.0)*7.0)%7
def occ(E,x):
    o=0
    for e in E: o|=1<<sec(e,x)
    return o
def Rs_arcs(E,s):
    arcs=[]; inR=False; start=0.0
    for k in range(1,NG):
        x=k/NG; o=occ(E,x); cur=(7-bin(o).count("1")==1) and not((o>>s)&1)
        if cur and not inR: start=(k-0.5)/NG; inR=True
        if (not cur) and inR: arcs.append((start,(k-0.5)/NG)); inR=False
    if inR: arcs.append((start,1.0))
    return arcs
def nn_sep(betas):
    P=sorted(betas); r=len(P); dd=[1.0]*r
    order=sorted(range(r), key=lambda i:P[i])
    B=[P[i] for i in order]
    for a in range(r):
        for b in range(a+1,r):
            v=abs(B[b]-B[a]); v=min(v,1-v)
            if v>0:
                if v<dd[a]: dd[a]=v
                if v<dd[b]: dd[b]=v
            if B[b]-B[a]>0.5: break
    return dd

print("width-weighted vs unweighted clustering sums for {w c_i} (midpoints x w), w=997")
print("does Sum w_i^2/delta_i stay ~O(mu/w) (MV closes) while Sum 1/delta_i blows up?")
print("="*80)
w=997
print("  {:26s} {:>4} {:>4} {:>7} {:>9} {:>11} {:>11} {:>11}".format(
      "E'","diam","r","mu","mu/w","Sum wi^2/di","Sum 1/di","L-term(=r? )"))
for E in [[0,1,2,3,4,5,6],[0,3,7,15,30,55,90],[0,10,27,55,99,150,199],[0,20,50,90,140,170,199]]:
    for s in [1]:
        arcs=Rs_arcs(E,s); r=len(arcs)
        if r<2: continue
        wid=[b-a for a,b in arcs]; mu=sum(wid)
        cmid=[((a+b)/2*w)%1.0 for a,b in arcs]
        dd=nn_sep(cmid)
        Sw=sum(wid[i]**2/dd[i] for i in range(r))
        Su=sum(1.0/dd[i] for i in range(r))
        print("  {:26s} {:>4} {:>4} {:>7.4f} {:>9.2e} {:>11.3e} {:>11.2e} {:>11}".format(
              str(E),max(E),r,mu,mu/w,Sw,Su,r))
print("-"*80)
print("  MV: Q_s/(2pi)^2 ~ Sum_ell |f_hat(ell w)|^2 <~ (over dyadic) Sum_i w_i^2(L+1/d_i)/L^2-terms.")
print("  If Sum w_i^2/d_i ~ mu/w (comparable to the L=1 diagonal Sum w_i^2 ~ mu/w): weighted MV gives")
print("  Q_s = (2pi w)^2 * O(mu/w) = O(w*mu) ... [check w-dependence]; vs unweighted Sum 1/d_i ~ r^2 (blows up).")
print("\ndone.")
