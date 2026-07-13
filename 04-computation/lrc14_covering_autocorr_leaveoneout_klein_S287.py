#!/usr/bin/env python3
"""
lrc14_covering_autocorr_leaveoneout_klein_S287.py
=================================================
klein-2026-07-13-S287 (owner: build the middle-order x-integral for the covering side).

FIX to the first run: opus's eps_v is the LEAVE-ONE-OUT covariance
    eps_v = Cov(1_{D_v}, 1_{G'_{~v}}),   G'_{~v} = intersection of good arcs of speeds w != v,
so D_v CAN overlap G'_{~v} (unlike the full G', which excludes D_v and gives the trivial eps=-|G'|/7).
Peeling identity:  L_full = (6/7)|G'_{~v}| - eps_v.

THE POSITIVE-DEFINITE x-INTEGRAL (covering THM-729 analogue), now on the leave-one-out set:
    |eps_v|^2 <= (6/49) * [ (1/v) sum_{j=0}^{v-1} A_{~v}(j/v) - |G'_{~v}|^2 ],           (*)
A_{~v}(tau)=meas(G'_{~v} ^ (G'_{~v}-tau))  (>=0, positive-definite). Bracket = v-grid discrepancy (>=0).

Reports: (1) does (*) hold; (2) the ALIGNMENT RATIO rho_v = |eps_v|^2 / [(6/49)disc_v] in [0,1]
(how tight is Cauchy-Schwarz -- rho near 1 = sharp, rho near 0 = loose); (3) the CERTIFICATE margin
L_cert(v) = (6/7)|G'_{~v}| - sqrt((6/49)disc_v)  vs the true L_full -- does the positive-definite bound
certify L_full>0? (4) ACID TEST: deep well vs near-AP residue.
"""
import numpy as np
NG = 1 << 21
THR = 1.0/14.0
t = np.arange(NG, dtype=np.float64)/NG

def good_indicator(S):
    g = np.ones(NG, dtype=np.float64)
    for w in S:
        frac=(w*t)%1.0; dist=np.minimum(frac,1.0-frac)
        g *= (dist>=THR)
    return g

def autocorr(g):
    G=np.fft.rfft(g); return np.fft.irfft(G*np.conj(G),n=NG)/NG

def Dv_ind(v):
    frac=(v*t)%1.0; dist=np.minimum(frac,1.0-frac); return (dist<THR).astype(np.float64)

def analyze(S,name,peel=None):
    Lfull=good_indicator(S).mean()
    cores = sorted(set(S)) if peel is None else peel
    rows=[]
    for v in cores:
        Snv=[w for w in S if w!=v]
        gnv=good_indicator(Snv); Lnv=gnv.mean()
        A=autocorr(gnv); Abar=A.mean()   # =Lnv^2
        Dv=Dv_ind(v)
        eps=(Dv*gnv).mean()-(1.0/7.0)*Lnv
        idx=np.round((np.arange(v)/v)*NG).astype(np.int64)%NG
        disc=A[idx].mean()-Abar
        rhs=(6.0/49.0)*disc
        rho=(eps*eps)/rhs if rhs>0 else float('nan')
        Lcert=(6.0/7.0)*Lnv-np.sqrt(max(rhs,0.0))
        rows.append((v,Lnv,eps,eps*eps,rhs,disc,rho,Lcert))
    return Lfull,rows

FAM=[([1,2,3,4,5,6,7,8,9,10,11,12,182],"deep well {1..12,182}"),
     ([1,2,3,4,5,6,7,8,9,10,11,13,84], "near-AP residue {1..11,13,84}"),
     ([2,3,4,5,6,7,8,9,10,11,12,13,14],"small-far {2..14} (far elt=14)"),
     ([1,3,4,5,6,7,8,9,10,11,12,13,182],"variant {1,3..13,182}")]

print("="*104)
print("COVERING MIDDLE-ORDER x-INTEGRAL, leave-one-out.  |eps_v|^2 <= (6/49)[v-grid disc of A_{~v}]   (NG=%d)"%NG)
print("  peeling identity L_full=(6/7)|G'_{~v}|-eps_v.  rho=tightness in [0,1]. L_cert=(6/7)|G'_{~v}|-sqrt(rhs).")
print("="*104)
acid=[]
for S,name in FAM:
    Lfull,rows=analyze(S,name)
    print("\n%s   L_full=|G'|=%.6f"%(name,Lfull))
    print("   %4s %10s %11s %11s %11s %8s %10s  %s"%("v","|G'_~v|","eps_v","eps_v^2","(6/49)disc","rho_v","L_cert","(*)"))
    allok=True; rhos=[]; best_cert=-1e9
    for v,Lnv,eps,eps2,rhs,disc,rho,Lcert in rows:
        ok=eps2<=rhs+1e-12; allok&=ok; rhos.append(rho); best_cert=max(best_cert,Lcert)
        if v in (1,) or v>=12 or not ok:
            print("   %4d %10.6f %11.6f %11.3e %11.3e %8.4f %10.6f  %s"%(v,Lnv,eps,eps2,rhs,rho,Lcert,"OK" if ok else "FAIL"))
    print("   -> (*) all cores: %s ; median rho=%.4f (tightness) ; BEST L_cert=%.6f (true L=%.6f)"%(
          allok,float(np.median(rhos)),best_cert,Lfull))
    acid.append((name,Lfull,float(np.median(rhos)),best_cert,max(r[5] for r in rows)))

print("\n"+"="*104)
print("ACID TEST (mac-mini-S83: metric object must rank residue MORE stuck than deep well).")
print("   %-34s %10s %10s %12s %12s"%("family","L_full","med_rho","best_Lcert","max_disc"))
for name,L,rho,cert,md in acid:
    print("   %-34s %10.6f %10.4f %12.6f %12.4e"%(name,L,rho,cert,md))
print("""
Reading:
 - (*) holding + rho<1 quantifies HOW LOOSE Cauchy-Schwarz is (tight covering bound needs rho->1).
 - L_cert>0 would CERTIFY lonely times from the positive-definite x-integral alone; L_cert<0 = too loose.
 - If best_Lcert is MORE negative (or disc/L ratio larger) for the residue, the autocorrelation
   discrepancy correctly flags it as the binding (more-stuck) family -- a faithful METRIC surrogate.
""")
print("done.")
