#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_routeB_corner_telescope_opus_0620s3.py   (opus-2026-06-20-S3)

BREAKTHROUGH probe.  From S2: U(E) = 1 - E[N] + sum_a w_a E[(N-a)_+].
The linear term -E[N] and the a=1 corner are NOT independent:
   E[(N-1)_+] = E[N] - P(N>=1) = E[N] - (1 - p_0).
So  -E[N] + 1*E[(N-1)_+] = -(1-p_0) = p_0 - 1.

More generally REWRITE everything in the SURVIVAL basis G_a = P(N>=a):
   E[(N-a)_+] = sum_{b>a} G_b = sum_{b=a+1}^{6} G_b.
   E[N] = sum_{b=1}^{6} G_b.
Then U = 1 - E[N] + sum_a w_a sum_{b>a} G_b  =  1 + sum_{b=1}^{6} c_b G_b
for some coefficients c_b.  Since p_0 = 1 - G_1, U is an AFFINE functional of the
survival vector (G_1,...,G_6).  The CLEAN question:  what are the c_b, and does
consec MAXIMIZE/EXTREMIZE the survival vector in a coordinatewise or majorization
sense weighted by c_b?

This reframes Route B as: an inequality on the SURVIVAL function G_a(E)=P(N_E>=a),
a=1..6 (a monotone decreasing sequence in a).  G_a(E) = S_?(E)?  Actually
G_a = P(at least a sectors empty) = P(N>=a) = sum_{t>=a} p_t.
Note S_r = E[C(N,r)] = sum_t C(t,r) p_t -- different basis (binomial), but both
linear in p.  We connect them.

GOAL:
  (1) get exact c_b for U4 and L_y, k=8,9,10.  Show U = 1 + sum_b c_b G_b.
  (2) which G_b does consec maximize?  (cut-by-cut on the survival ladder)
  (3) the DECISIVE structural test:  c_b signs.  If c_b for the LARGE b (deep tail,
      b=4,5,6) are POSITIVE and consec maximizes those G_b, while the SMALL-b
      coefficients are NEGATIVE and consec MINIMIZES those G_b, then EVERY term
      points the same way -> a genuine cut-by-cut proof that survives the
      non-separability of the moments.
  (4) verify the telescoped identity p_0-1 = -E[N]+E[(N-1)_+] exactly on consec.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def dist_p(E):
    E = sorted(set(E))
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e+1):
            bps.add(F(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [F(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        hit = set()
        for e in E:
            v = (e*mid) % 1
            hit.add((v.numerator*7)//v.denominator)
        t = sum(1 for j in range(1,7) if j not in hit)
        p[t] += (hi-lo)
    return p

def phi_Ly(k):
    if k==8: return [F((t-1)*(t-2)*(t-4)*(t-5),40) for t in range(7)]
    if k in (9,10): return [F(-(t-2)*(t-3)*(t-6),36) for t in range(7)]
    return [F((t-3)*(t-4),12) for t in range(7)]
PHI_U4 = [F(1),F(0),F(0),F(0),F(0),F(1),F(5)]

def survival(p): return [sum(p[t] for t in range(a,7)) for a in range(7)]  # G_0..G_6

def phi_to_survival_coeffs(phi):
    """phi(t)=phi[0]+sum_{b=1}^6 c_b * 1{t>=b}.  Because 1{t>=b} steps.
       Apply: U=sum_t phi[t] p_t = phi[0]*1 + sum_b c_b G_b  where c_b=phi[b]-phi[b-1]."""
    c={b: phi[b]-phi[b-1] for b in range(1,7)}
    return phi[0], c

def apply_phi(p,phi): return sum(p[t]*phi[t] for t in range(7))
def consec(k): return list(range(k))

if __name__=="__main__":
    for name,getphi in [("U4",lambda k:PHI_U4),("L_y",phi_Ly)]:
        print("="*78); print(f"FUNCTIONAL {name}: U = phi0 + sum_b c_b * G_b   (G_b=P(N>=b))")
        print("="*78)
        for k in (8,9,10):
            phi=getphi(k)
            phi0,c=phi_to_survival_coeffs(phi)
            print(f"\n--- k={k} phi={[str(x) for x in phi]} ---")
            print(f"  survival coeffs: U = {phi0} + " +
                  " + ".join(f"({c[b]})*G_{b}" for b in range(1,7)))
            signs={b:('+' if c[b]>0 else ('-' if c[b]<0 else '0')) for b in range(1,7)}
            print(f"  c_b signs (b=1..6): {[signs[b] for b in range(1,7)]}")

            C=consec(k); pc=dist_p(C); Gc=survival(pc); Uc=apply_phi(pc,phi)
            # verify reconstruction
            Urec=phi0+sum(c[b]*Gc[b] for b in range(1,7))
            print(f"  consec U={Uc}={float(Uc):.6f}  reconstruction matches: {Urec==Uc}")
            print(f"  consec G_b (b=1..6) = {[float(Gc[b]) for b in range(1,7)]}")

            span=k+4 if k==8 else k+3
            bank=[[0]+list(r) for r in itertools.combinations(range(1,span+1),k-1)]
            # cut-by-cut: for each b, does consec EXTREMIZE G_b in the direction c_b wants?
            # c_b>0 wants G_b max; c_b<0 wants G_b min.
            wrongdir={b:0 for b in range(1,7) if c[b]!=0}
            beats=0
            for E in bank:
                p=dist_p(E); G=survival(p); U=apply_phi(p,phi)
                if U>Uc+F(1,10**15): beats+=1
                for b in wrongdir:
                    if c[b]>0 and G[b]>Gc[b]+F(1,10**15): wrongdir[b]+=1   # someone bigger (bad for max)
                    if c[b]<0 and G[b]<Gc[b]-F(1,10**15): wrongdir[b]+=1   # someone smaller (bad for min)
            print(f"  beats consec on U: {beats}/{len(bank)}")
            print(f"  per-cut: # shapes where consec is NOT extremal in c_b's direction: {wrongdir}")
            allgood=all(v==0 for v in wrongdir.values())
            print(f"  => cut-by-cut survival certificate {'HOLDS (each G_b extremal!)' if allgood else 'FAILS at some b'}")

    # telescope identity check
    print("\n"+"="*78); print("TELESCOPE IDENTITY: -E[N]+E[(N-1)_+] = p_0 - 1 ?"); print("="*78)
    for k in (8,9):
        pc=dist_p(consec(k))
        EN=sum(t*pc[t] for t in range(7))
        c1=sum((t-1)*pc[t] for t in range(2,7))  # E[(N-1)_+]
        lhs=-EN+c1; rhs=pc[0]-1
        print(f"  k={k}: -E[N]+E[(N-1)_+]={lhs}  p_0-1={rhs}  EQUAL={lhs==rhs}")
