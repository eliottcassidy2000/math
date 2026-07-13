#!/usr/bin/env python3
"""
lrc14_chi_dft_klein_S278.py
===========================
klein-2026-07-13-S278 (owner: carry out the k-dim Weyl estimate on the endpoint sum).

KEY (S277): the coupled endpoint inner sum Sum_{j<e'} chi(j) e(-Nj/e') = e'*chihat(N mod e') EXACTLY
(chi periodic mod e' in j => 1-D DFT over Z/e'). So U_s^{e'}(N) = e' Sum_sigma e(-N sigma/7e') chihat_{s,sigma}(N mod e').
The bound then needs: (0) net imbalance chihat(0)=mean chi = O(1)? (resonance term), and
(k) off-DC |chihat(kappa)| decay.

chi_{s}(j): at boundary s (e' ENTERS sector s, winding j): -1 if others cover exactly {0..6}\{s}
   (R_s-END). At boundary s+1 (e' LEAVES s): +1 if others cover exactly {0..6}\{s} (R_s-START).
   "others" = E'\{e'} (incl 0). cond_s(x) = [occ_others == {0..6}\{s}].

Compute: NET(e',s)=Sum_j chi_s(j) (the imbalance = e'*chihat_s(0)); max|NET|; scaling vs e', k, Sigma-e'';
and the DFT profile max_kappa |chihat(kappa)| (off-DC). Prediction: LARGEST offset => NET=0 exactly.
"""
import math
from math import gcd
from functools import reduce
def sec(e,x): return int((e*x%1.0)*7.0)%7
def occ_others(E, eprime, x):
    o=0
    for e in E:
        if e!=eprime: o|=1<<sec(e,x)
    return o
def cond_s(E,eprime,x,s):
    # others cover exactly {0..6}\{s}: 6 sectors, missing exactly s
    o=occ_others(E,eprime,x)
    return (bin(o).count("1")==6) and not ((o>>s)&1)

def chi_seq(E, eprime, s):
    """chi_s(j) for j=0..e'-1: -cond at enter (boundary s), +cond at leave (boundary s+1)."""
    ch=[]
    for j in range(eprime):
        x_enter=(j + s/7.0)/eprime
        x_leave=(j + ((s+1)%7)/7.0)/eprime
        c = (1.0 if cond_s(E,eprime,x_leave,s) else 0.0) - (1.0 if cond_s(E,eprime,x_enter,s) else 0.0)
        ch.append(c)
    return ch

def dft_maxabs_and_dc(ch):
    e=len(ch);
    # DC = mean* e' = sum(ch); chihat(0)=sum/e'
    dc=sum(ch)
    # max off-DC |e'*chihat(kappa)| = max_{kappa!=0} |sum_j ch[j] e(-kappa j/e')|
    mx=0.0
    for kappa in range(1,e):
        re=sum(ch[j]*math.cos(-2*math.pi*kappa*j/e) for j in range(e))
        im=sum(ch[j]*math.sin(-2*math.pi*kappa*j/e) for j in range(e))
        mx=max(mx, math.hypot(re,im))
    return dc, mx

print("net endpoint imbalance NET(e',s)=Sum_j chi_s(j) [=e'*chihat(0)], and off-DC max, per offset")
print("(prediction: NET=0 for the LARGEST offset; |NET| small overall)")
print("="*76)
clusters=[
  [0,1,2,3,4,5,6],[0,1,2,3,4,5,7],[0,1,2,4,7,10,13],[0,1,2,28,29,30,15],
  [0,1,2,3,5,7,11],[0,1,5,6,10,11,15],[0,1,2,3,4,10,20],
]
globmaxNET=0; globmaxOff=0
for E in clusters:
    Emax=max(E); sumall=sum(E)
    print(f"\n  E'={E} (max={Emax}, sum={sumall})")
    for ep in sorted(E):
        if ep==0: continue
        sig_e=sum(e for e in E if e not in (0,ep))  # sum of OTHER nonzero offsets
        maxNET=0; maxoff=0; sNET=[]
        for s in range(7):
            ch=chi_seq(E,ep,s); dc,offmx=dft_maxabs_and_dc(ch)
            sNET.append(int(round(dc))); maxNET=max(maxNET,abs(dc)); maxoff=max(maxoff,offmx)
        tag=" <-LARGEST" if ep==Emax else ""
        globmaxNET=max(globmaxNET,maxNET); globmaxOff=max(globmaxOff,maxoff)
        print(f"    e'={ep:3d} (sig_others={sig_e:3d}): max_s|NET|={maxNET:.0f}  max_s off-DC={maxoff:.2f}  NET_s={sNET}{tag}")
print("-"*76)
print(f"  GLOBAL max|NET| = {globmaxNET:.0f}   GLOBAL max off-DC = {globmaxOff:.2f}")
print(f"  => if max|NET| small (O(1) or O(k)) and off-DC bounded, U_s^e'(N)=e'*chihat is controlled.")
print("\ndone.")
