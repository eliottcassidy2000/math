#!/usr/bin/env python3
"""
lrc14_BS_bandlimited_signed_macmini_0618s7.py  (mac-mini-2026-06-18-S7) ANGLE A diagnosis

CRUX DIAGNOSTIC for the Beurling-Selberg route.  The route integrates a degree-N majorant of
1_{S7}; the integral is a FINITE SIGNED lattice sum (cancellation preserved).  Two questions:

 (Q1) Does the SIGNED band-limited approximation of meas(S7) (truncate the exact relation-lattice
      sum to |n_i|<=N, KEEPING signs) converge to the true meas(S7) as N grows?  If yes, a
      band-limited certificate is in principle possible.  We compute this exactly-ish (float, but
      the lattice sum is finite per N) for consec k=7,8.
 (Q2) The Beurling-Selberg MAJORANT replaces hat 1_{B_T}(n) by hat V(n) (|n|<=N) with hat V(0)=
      |B_T|+1/(N+1) and |hat V(n)-hat 1(n)| <= 1/(N+1).  The resulting integral is a GUARANTEED
      upper bound.  The overshoot vs true meas(S7) is the price of band-limiting.  We model the
      BEST CASE (a hypothetical perfect majorant with hat V(n)=hat 1(n) for n!=0, only the n=0
      mode inflated by 1/(N+1) per factor) to see the unavoidable main-term inflation:
        each factor's 0-mode becomes (|B_T|+1/(N+1)) -> M7-like term with inflated base.
      Compute M7_infl(k,N) = sum_T (-1)^|T| (1-|T|/7 + 1/(N+1))^{k-1} and compare to cap_k.

This isolates whether the inflation alone (even with a PERFECT high-mode majorant) already breaks
the bound, independent of the tail-cancellation issue.
"""
import sys, itertools, math, cmath
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def measS7_geom(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E:
            y=(e*xm)%1; secs.add(int(y*7))
        if len(secs)==7: total+=x1-x0
    return total

def M7(k):
    s=F(0)
    for t in range(0,7):
        s += F((-1)**t * math.comb(6,t)) * F(7-t,7)**(k-1)
    return s

cap8=0.38146; cap_float={7:0.0,8:0.38146,9:0.4943,10:0.6044,11:0.7253,12:0.8571,13:1.0}

def hat_BT(T, n):
    if n==0: return 1.0 - len(T)/7.0
    s=0j
    for j in T:
        s += cmath.exp(-2j*math.pi*n*j/7.0)*(1-cmath.exp(-2j*math.pi*n/7.0))/(-2j*math.pi*n)
    return -s

# (Q1) signed band-limited lattice sum of meas(S7), truncating each n_i to |n|<=N.
# DP over partial relation sum to avoid (2N+1)^{k-1} blowup. Complex coeffs.
def measS7_bandlimited(E, N):
    Enz=[e for e in E if e!=0]
    total=0.0
    for r in range(0,7):
        for T in itertools.combinations(range(1,7), r):
            sgn=(-1)**r
            # DP: dict mapping partial sum s -> accumulated complex product, over factors i with e_i.
            # factor i contributes sum_{n=-N..N} hat_BT(T,n) e at shift n*e_i. We need total s=0.
            dp={0:1.0+0j}
            for e in Enz:
                nd={}
                for s,acc in dp.items():
                    for n in range(-N,N+1):
                        c=hat_BT(T,n)
                        if c==0: continue
                        ns=s+n*e
                        nd[ns]=nd.get(ns,0j)+acc*c
                dp=nd
            total += sgn*dp.get(0,0j).real
    return total

# (Q2) best-case main-term inflation with a perfect high-mode majorant.
def M7_inflated(k,N):
    s=0.0; inf=1.0/(N+1)
    for t in range(0,7):
        s += ((-1)**t)*math.comb(6,t)*((7-t)/7.0 + inf)**(k-1)
    return s

if __name__=="__main__":
    print("="*86)
    print("ANGLE A diagnosis: signed band-limited convergence + best-case main-term inflation")
    print("="*86)
    print("\n(Q1) signed band-limited lattice sum -> true meas(S7)?  (keeps cancellation)")
    for name,E in [("consec k=7",list(range(7))),("consec k=8",list(range(8))),
                   ("dissoc k=8 2^i",[0,1,3,7,15,31,63,127]),
                   ("Sidon k=8",[0,1,3,7,12,20,30,44]),
                   ("generic k=8",[0,5,13,27,41,58,79,97])]:
        g=float(measS7_geom(E))
        print(f"\n  {name}: true meas(S7)={g:.5f}")
        for N in [3,5,7,9,11]:
            v=measS7_bandlimited(E,N)
            print(f"     N={N:>2}: signed band sum = {v:.5f}   (err {v-g:+.5f})")
    print("\n(Q2) best-case inflation: even a PERFECT high-mode majorant inflates the n=0 base by 1/(N+1):")
    print(f"  {'k':>3}{'cap_k':>9} | M7_inflated(k,N) for N=:   {'6':>9}{'12':>9}{'24':>9}{'48':>9}{'96':>9}  (M7 exact)")
    for k in [8,9,10,11,12,13]:
        row=f"  {k:>3}{cap_float[k]:>9.4f} | "
        for N in [6,12,24,48,96]:
            row+=f"{M7_inflated(k,N):>9.4f}"
        row+=f"   {float(M7(k)):>8.4f}"
        print(row)
    print("\nIf M7_inflated > cap_k even at large N, the band-limiting inflation alone breaks the bound")
    print("for the dangerous rows -> a pure interval-majorant-per-factor certificate cannot close k=8..11.")
