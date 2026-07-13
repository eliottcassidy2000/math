#!/usr/bin/env python3
"""
lrc14_chi_largeeprime_klein_S278.py
===================================
klein-2026-07-13-S278. Confirm the KEY bound at LARGE e': is e'|chihat(kappa)| = |Sum_{j<e'} ch_j e(-kappa j/e')|
BOUNDED (O(1)) as e' grows, for the swing-difference sequence ch_j = cond_s(leave_j)-cond_s(enter_j)?
If yes, |U_s^{e'}(N)| <= max_kappa e'|chihat| = O(1) => |S| = O(k) (the estimate closes).

Test e'=D as the large offset with fixed small others {0,1,2,3,4,5} (and a couple mixed clusters).
Report max over s,kappa of e'|chihat|, and the # nonzeros T, as e' grows.
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
def maxdft(ch):
    e=len(ch); mx=0.0
    for kappa in range(e):
        re=sum(ch[j]*math.cos(2*math.pi*kappa*j/e) for j in range(e))
        im=sum(ch[j]*math.sin(2*math.pi*kappa*j/e) for j in range(e))
        mx=max(mx,math.hypot(re,im))
    return mx, sum(1 for c in ch if c!=0), int(round(sum(ch)))

print("large-e' test: max_{s,kappa} e'|chihat| (=max|U_s^{e'}| contribution) as e' grows")
print("="*72)
print("  {:38s} {:>4} {:>7} {:>5} {:>6}".format("cluster (eprime=large offset)","ep","maxDFT","maxT","maxNET"))
for others,label in [([0,1,2,3,4,5],"{0..5,D}"), ([0,1,2,3,4,5,6],"{0..6,D}"), ([0,1,7,12,20,33],"spread6+D")]:
    for D in [50,100,200,400]:
        E=sorted(others+[D])
        if len(E)!=len(others)+1: continue
        mxdft=0;mxT=0;mxNET=0
        for s in range(7):
            ch=chi_seq(E,D,s); md,T,net=maxdft(ch)
            mxdft=max(mxdft,md);mxT=max(mxT,T);mxNET=max(mxNET,abs(net))
        print(f"  {label+' D='+str(D):38s} {D:4d} {mxdft:7.2f} {mxT:5d} {mxNET:6d}")
print("-"*72)
print("  => if maxDFT stays bounded (~O(1)) as e'->400, the estimate CLOSES: |U_s^e'|=O(1), |S|=O(k).")
print("     (T may grow, but the DFT/Fourier-linf norm is what bounds U; watch maxDFT vs maxT.)")
print("\ndone.")
