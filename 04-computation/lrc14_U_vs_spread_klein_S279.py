#!/usr/bin/env python3
"""
lrc14_U_vs_spread_klein_S279.py
===============================
klein-2026-07-13-S279. #cond_s-arcs GROWS ~diam (rainbow_arcs test). So U_s^{e'}<=#arcs is NOT bounded;
the S276/S278 boundedness must come from DFT cancellation. DECISIVE: how does |U_s^{e'}|=e'|chihat(kappa)|
scale as the OTHERS spread (diam grows)?
  O(1)     => density crux is bounded (maybe provable directly)
  ~sqrt(#arcs) => 2nd-moment/large-sieve bound; SUFFICES for density (slack: err~sqrt(diam)/d->0), tractable
  ~#arcs   => no cancellation; hard multi-linear (Gowers), same as covering
Compute max_{s,kappa} e'|chihat| for e'=large offset D, with 6 SPREAD others of growing diameter.
"""
import math
def sec(e,x): return int((e*x%1.0)*7.0)%7
def occ_others(E,ep,x):
    o=0
    for e in E:
        if e!=ep: o|=1<<sec(e,x)
    return o
def cond_s(E,ep,x,s):
    o=occ_others(E,ep,x); return (bin(o).count("1")==6) and not ((o>>s)&1)
def chi_seq(E,ep,s):
    ch=[]
    for j in range(ep):
        xe=(j+s/7.0)/ep; xl=(j+((s+1)%7)/7.0)/ep
        ch.append((1.0 if cond_s(E,ep,xl,s) else 0.0)-(1.0 if cond_s(E,ep,xe,s) else 0.0))
    return ch
def maxdft_and_T(ch):
    e=len(ch); mx=0.0
    for kappa in range(e):
        re=sum(ch[j]*math.cos(2*math.pi*kappa*j/e) for j in range(e))
        im=sum(ch[j]*math.sin(2*math.pi*kappa*j/e) for j in range(e))
        mx=max(mx,math.hypot(re,im))
    return mx, sum(1 for c in ch if c!=0)

print("max_{s,kappa} e'|chihat| (=max|U_s^{e'}|) for e'=D=997 (prime, large), 6 SPREAD others of growing diam")
print("="*76)
print("  {:34s} {:>6} {:>8} {:>6} {:>10}".format("6 others (spread)","diam","maxDFT","maxT","DFT/sqrtT"))
D=997
fams=[
  ("{0,1,2,3,4,5}",     [0,1,2,3,4,5]),
  ("{0,1,2,3,4,20}",    [0,1,2,3,4,20]),
  ("{0,1,4,9,16,25}",   [0,1,4,9,16,25]),
  ("{0,7,17,31,50,73}", [0,7,17,31,50,73]),
  ("{0,10,23,41,66,99}",[0,10,23,41,66,99]),
  ("{0,20,50,90,140,199}",[0,20,50,90,140,199]),
]
for name,others in fams:
    E=sorted(others+[D])
    mxdft=0;mxT=0
    for s in range(7):
        ch=chi_seq(E,D,s); md,T=maxdft_and_T(ch)
        mxdft=max(mxdft,md); mxT=max(mxT,T)
    r = mxdft/math.sqrt(mxT) if mxT>0 else 0
    print("  {:34s} {:6d} {:8.2f} {:6d} {:10.2f}".format(name,max(others),mxdft,mxT,r))
print("-"*76)
print("  maxT ~#cond_s-arcs grows ~diam. If maxDFT stays ~O(1): bounded (provable). If ~sqrt(maxT):")
print("  2nd-moment bound suffices for density (slack). If ~maxT: hard multilinear (like covering).")
print("\ndone.")
