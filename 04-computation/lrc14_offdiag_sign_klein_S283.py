#!/usr/bin/env python3
"""
lrc14_offdiag_sign_klein_S283.py
================================
klein-2026-07-13-S283 (owner: prove the Weyl cancellation for arc midpoints).

Reframing: Q_s = diag + offdiag, where
  diag  = Sum_i 4 pi^2 {w w_i}(1-{w w_i})   [CLOSED FORM, arc widths w_i; RIGOROUS; = O(w mu)=O(r)]
  offdiag = Q_s - diag = (2pi w)^2 Sum_{i!=j} [cross-correlation of arc i,j at xw].
If offdiag <= 0 (arcs ANTI-correlate under xw), then Q_s <= diag = O(r) RIGOROUSLY => density row closes.

TEST: is offdiag <= 0 ROBUSTLY (all clusters, all sectors, clean AND resonant w)? If yes, "prove
offdiag<=0" is the clean concrete target replacing the vague Weyl cancellation.
"""
import math
NG=1500000
LMAX=800
def sec(e,x): return int((e*x%1.0)*7.0)%7
def occ(E,x):
    o=0
    for e in E: o|=1<<sec(e,x)
    return o
def Rs_arcs(E,s):
    arcs=[]; inR=False; st=0.0
    for k in range(1,NG):
        x=k/NG; o=occ(E,x); cur=(7-bin(o).count("1")==1) and not((o>>s)&1)
        if cur and not inR: st=(k-0.5)/NG; inR=True
        if (not cur) and inR: arcs.append((st,(k-0.5)/NG)); inR=False
    if inR: arcs.append((st,1.0))
    return arcs
def Qs_endpoints(arcs,w):
    # Q_s = Sum_{ell=1}^{LMAX} 2*|U_s(ell w)|^2/ell^2, U_s(N)=Sum_arcs (e(-N a)-e(-N b))
    tot=0.0
    for ell in range(1,LMAX+1):
        N=ell*w; re=im=0.0
        for (a,b) in arcs:
            re+= math.cos(-2*math.pi*N*a)-math.cos(-2*math.pi*N*b)
            im+= math.sin(-2*math.pi*N*a)-math.sin(-2*math.pi*N*b)
        tot+= 2*(re*re+im*im)/(ell*ell)
    return tot
def frac(x): return x-math.floor(x)
def diag(arcs,w):
    d=0.0
    for (a,b) in arcs:
        wi=b-a; f=frac(w*wi); d+= 4*math.pi*math.pi*f*(1-f)
    return d

print("Q_s vs diag (closed form) vs offdiag=Q_s-diag; is offdiag <= 0 robustly?")
print("="*76)
print("  {:26s} {:>2} {:>4} {:>9} {:>9} {:>9} {:>7}".format("E'","s","r","Q_s","diag","offdiag","off<=0?"))
allneg=True
for E in [[0,1,2,3,4,5,6],[0,3,7,15,30,55,90],[0,10,27,55,99,150,199],[0,1,2,28,29,30,15]]:
    for w in [997, None]:  # None => resonant (lcm-ish); use a resonant multiple
        ww = w if w else 60  # for {0,1,2,28,29,30,15} lcm huge; just test w=997 and a small resonant 60
        for s in [1,3]:
            arcs=Rs_arcs(E,s); r=len(arcs)
            if r<2: continue
            q=Qs_endpoints(arcs,ww); dg=diag(arcs,ww); off=q-dg
            if off>1e-6: allneg=False
            print("  {:26s} {:>2} {:>4} {:>9.3f} {:>9.3f} {:>9.3f} {:>7}".format(str(E),s,r,q,dg,off,"YES" if off<=1e-6 else "NO!"))
        break  # one w per cluster to keep it short
print("-"*76)
print("  offdiag <= 0 for all above:", allneg)
print("  If robust: Q_s <= diag = Sum_i 4pi^2{w w_i}(1-{w w_i}) <= 4pi^2 w mu = O(r) RIGOROUS => row closes.")
print("  => the clean concrete target is: PROVE offdiag<=0 (R_s-arcs anti-correlate under xw).")
print("\ndone.")
