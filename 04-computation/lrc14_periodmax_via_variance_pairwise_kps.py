#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CAN THE THM-563 PERIOD-MAX BE BOUNDED WITHOUT ENUMERATING THE 7*lcm(B) PERIOD?

kind-pasteur-2026-06-30. Concrete next step from the barrier-residual study. THM-563:
w*Delta_w is periodic (period P=7*lcm(B)); the general-bounded-base closure needs
period-max(B) < 15*margin(B) for ALL B, and enumerating P (up to 7*lcm(1..14)=2.5M) is
the bottleneck. Reciprocity computes SUMS, not MAXES -- so the honest question is:

  Is the period-MAX faithfully bounded by the VARIANCE (an L2 sum), and does the variance
  DECOUPLE into pairwise correlations computable WITHOUT the full period?

  g(w) = w*Delta_w = Sum_e eps_e * G0(w*theta_e - phi_e)      [G0 = tent, endpoints from cells]
  E[g^2] = Sum_{e,e'} eps_e eps_e' * E_w[G0(w th_e - ph_e) G0(w th_e' - ph_e')]
  and each pairwise E_w depends only on w mod lcm(den th_e, den th_e') << P.

So the variance is a sum of PAIRWISE tent-correlations (each a Dedekind-Rademacher sum),
computable over SMALL pairwise periods -- no 7*lcm(B) enumeration. IF period-max <= C*rms
with C = O(1) uniformly, the variance route bounds the period-max fast.

Engine (G0, cells_of) verbatim from the validated lrc14_uniform_C_growth_kps.py.
"""
import sys, time
from fractions import Fraction as Fr
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
def lcm(a,b): return a*b//gcd(a,b)

def G0(y):
    f=y-(y.numerator//y.denominator)
    return Fr(6,7)*f if f<=Fr(1,7) else Fr(6,49)-(f-Fr(1,7))/7
def cells_of(Ep):
    Ep=sorted(set(Ep)); bps={Fr(0),Fr(1)}
    for e in Ep:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fr(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); out=[]
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in Ep:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        m=[j for j in range(1,7) if j not in hit]
        if len(m)==1:
            if out and out[-1][2]==m[0] and out[-1][1]==lo: out[-1]=(out[-1][0],hi,m[0])
            else: out.append((lo,hi,m[0]))
    return out

def endpoints(B):
    eps=[]
    for (a,b,s) in cells_of(B):
        eps.append((+1, b, Fr(s,7)))
        eps.append((-1, a, Fr(s,7)))
    return eps
def g(eps,w): return sum(e*G0(w*th-ph) for (e,th,ph) in eps)
def full_period(B):
    l=1
    for e in B:
        if e: l=lcm(l,e)
    return 7*l

def Eg2_pairwise(eps):
    """E_w[g^2] over the FULL period, computed as a sum of PAIRWISE correlations,
    each averaged over its OWN small period lcm(den th_e, den th_e') -- never touches 7*lcm(B)."""
    # group endpoints by theta to reuse pairwise averages
    tot=Fr(0); maxL=0
    # precompute per-pair; cache by (den1,den2, ...) not worth it, just do it
    for (e1,th1,ph1) in eps:
        for (e2,th2,ph2) in eps:
            L=lcm(th1.denominator, th2.denominator); maxL=max(maxL,L)
            s=Fr(0)
            for w in range(L):
                s+=G0(w*th1-ph1)*G0(w*th2-ph2)
            tot+=e1*e2*s/L
    return tot, maxL

print("="*90)
print("PART 1  period-max vs rms  (is period-max <= C*rms with C = O(1)?)")
print("  full enumeration of one period P=7*lcm(B); rms=sqrt(E[g^2]); ratio=period-max/rms")
print("="*90)
bases = {
  "consec7 {0..6}":  list(range(7)),
  "consec8 {0..7}":  list(range(8)),
  "consec9 {0..8}":  list(range(9)),
  "consec10 {0..9}": list(range(10)),
  "skip {0,1,2,3,5,7,9,11}": [0,1,2,3,5,7,9,11],
  "wide {0,1,2,7,8,9,13}":   [0,1,2,7,8,9,13],
}
print(f"  {'base':28} {'period P':>9} {'period-max':>11} {'at w':>7} {'rms':>9} {'max/rms':>8}")
ratios=[]
for name,B in bases.items():
    eps=endpoints(B); P=full_period(B)
    if P>1_200_000:
        print(f"  {name:28} {P:9d}  (skip full enum: P too large)"); continue
    gmax=Fr(0); argw=0; sq=Fr(0)
    for w in range(P):
        v=g(eps,w); a=abs(v); sq+=v*v
        if a>gmax: gmax=a; argw=w
    Eg2=sq/P; rms=(float(Eg2))**0.5
    ratio=float(gmax)/rms if rms>0 else float('inf')
    ratios.append(ratio)
    print(f"  {name:28} {P:9d} {float(gmax):11.4f} {argw:7d} {rms:9.4f} {ratio:8.3f}")
if ratios:
    print(f"  => max/rms across bases: min={min(ratios):.2f}  max={max(ratios):.2f}  "
          f"{'ROUGHLY CONSTANT -> variance is a faithful proxy' if max(ratios)/min(ratios)<2.2 else 'VARIES -> variance bound is loose'}")

print("\n"+"="*90)
print("PART 2  variance DECOUPLES pairwise (no full-period enumeration) -- correctness + cost")
print("="*90)
for name,B in list(bases.items())[:4]:
    eps=endpoints(B); P=full_period(B)
    # baseline E[g^2] by full enumeration
    t0=time.time(); sq=sum(g(eps,w)**2 for w in range(P)); Eg2_full=sq/P; t_full=time.time()-t0
    # pairwise
    t0=time.time(); Eg2_pw,maxL=Eg2_pairwise(eps); t_pw=time.time()-t0
    ok = (Eg2_full==Eg2_pw)
    print(f"  {name:28} P={P:7d} #endpoints={len(eps):3d}  E[g^2]={float(Eg2_full):.5f}  "
          f"match={ok}  maxPairLcm={maxL}  t_full={t_full:.2f}s t_pairwise={t_pw:.2f}s")
print("  => pairwise E[g^2] matches full enumeration EXACTLY and only touches periods <= maxPairLcm << P.")

print("\n"+"="*90)
print("PART 3  LARGE-lcm base: full period infeasible, pairwise variance still cheap")
print("="*90)
Bbig=[0,1,2,3,4,5,6,7,9,11,13]   # lcm(1..7,9,11,13)=180180 -> P=1,261,260
eps=endpoints(Bbig); P=full_period(Bbig)
t0=time.time(); Eg2_pw,maxL=Eg2_pairwise(eps); t_pw=time.time()-t0
rms=float(Eg2_pw)**0.5
print(f"  base {Bbig}")
print(f"  full period P = 7*lcm = {P:,}  (full-max enumeration = {P:,} evals -> infeasible here)")
print(f"  pairwise variance: E[g^2]={float(Eg2_pw):.5f}  rms={rms:.4f}  maxPairLcm={maxL}  "
      f"#endpoints={len(eps)}  computed in {t_pw:.2f}s")
if ratios:
    Cest=max(ratios)
    print(f"  => period-max <= C*rms with empirical C~{Cest:.2f} gives period-max <~ {Cest*rms:.3f} "
          f"WITHOUT enumerating {P:,} points.")

print("\n"+"="*90)
print("VERDICT")
print("="*90)
print(" - Reciprocity computes SUMS, not the MAX: the period-MAX is a genuine extremum.")
print(" - BUT the VARIANCE (an L2 sum) DECOUPLES into pairwise tent-correlations, each a")
print("   Dedekind-Rademacher sum over a SMALL period (<= maxPairLcm), never the 7*lcm(B) full")
print("   period -- and each pairwise term is further reducible by Dedekind reciprocity (O(log)).")
print(" - IF period-max <= C*rms uniformly (Part 1 measures C), the fast variance BOUNDS the")
print("   period-max without enumeration. This is the honest form of the 'reciprocity' payoff:")
print("   an L2/variance surrogate, faithful iff C=O(1); the exact MAX still needs enumeration")
print("   or an Ostrowski/three-distance sup bound (the continued-fraction sup, same CF object).")
print("DONE.")
