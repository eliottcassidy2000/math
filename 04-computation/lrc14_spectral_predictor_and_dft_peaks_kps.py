#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
(A) Which order parameter predicts p0 best -- additive energy, spectral flatness,
    or spread?  and (B) does the DFT of the ACTUAL extremizer {0,1,2,30,31,32,60,61}
    predict WHERE w*Delta_w peaks?

Context (kind-pasteur-2026-06-30, follows lrc14_riesz_depth_ladder_kps.py):
  The repo found "consec maximizes p0" but additive energy is a NON-MONOTONE proxy
  (1368+ within-span inversions; consec-max via raw additive energy REFUTED). Reframing
  C asked: is a SPECTRAL functional (Wiener entropy / spectral flatness) a cleaner
  monotone?  The Sidon result (depth-ladder script) says the Delta-driver is SPREAD,
  so we also test |E-E| (difference-set size) as a predictor.

  Additive energy is the 4th spectral moment  AE = ||\hat 1_E||_4^4 ; spectral flatness
  is the geometric/arithmetic-mean ratio of the SAME power spectrum -- more robust to a
  few large peaks. This is a clean apples-to-apples "which spectral functional tracks
  p0's ordering" test (Kendall tau over a fixed k-set family).

Engine (G0, cells_of, wD, p0) copied verbatim from the verified
lrc14_uniform_C_growth_kps.py / lrc14_wide_multiscale_p0_kps.py.
"""
import sys, itertools
from fractions import Fraction
from functools import reduce
from math import gcd
import numpy as np
from scipy.stats import kendalltau
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

# ---- VERIFIED ENGINE (verbatim) ----
def G0(y):
    f=y-(y.numerator//y.denominator)
    return Fraction(6,7)*f if f<=Fraction(1,7) else Fraction(6,49)-(f-Fraction(1,7))/7
def cells_of(Ep):
    Ep=sorted(set(Ep)); bps={Fraction(0),Fraction(1)}
    for e in Ep:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
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
def wD(cells,w):
    return sum(G0(w*b-Fraction(s,7))-G0(w*a-Fraction(s,7)) for (a,b,s) in cells)
def p0(E):
    E=sorted(set(E)); bps={Fraction(0),Fraction(1)}
    for e in E:
        if e==0: continue
        for a in range(0,7*e+1): bps.add(Fraction(a,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); tot=Fraction(0)
    for i in range(len(bps)-1):
        lo,hi=bps[i],bps[i+1]
        if hi==lo: continue
        mid=(lo+hi)/2; hit=set()
        for e in E:
            v=e*mid; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        if len(hit)==7: tot+=hi-lo
    return tot
# ------------------------------------

def add_energy(E):  # difference additive energy = #{(a,b,c,d): a-b=c-d} = sum_delta d(delta)^2
    from collections import Counter
    d=Counter()
    for a in E:
        for b in E: d[a-b]+=1
    return sum(v*v for v in d.values())
def diffset_size(E):
    return len({a-b for a in E for b in E})
def spectral_flatness(E, M=211):
    # geomean/mean of |\hat 1_E(j/M)|^2 over j=1..M-1  (Wiener entropy in [0,1]; small=peaky/structured)
    E=np.array(sorted(E)); j=np.arange(1,M)
    ph=np.exp(2j*np.pi*np.outer(E, j)/M)        # |E| x (M-1)
    fh=ph.sum(axis=0); power=np.abs(fh)**2
    power=np.maximum(power,1e-12)
    return float(np.exp(np.mean(np.log(power)))/np.mean(power))

print("="*80)
print("(A) PREDICTOR SHOOT-OUT: which order parameter tracks p0's ordering best?")
print("    family: all k=8 subsets E of {0..13} with 0 in E  (C(13,7)=1716 sets)")
print("    predictors: additive energy AE | spectral flatness SF | diff-set size |E-E| | width")
print("    metric: Kendall tau_b vs p0  (|tau|->1 = fewer inversions = cleaner monotone)")
print("="*80)
sets=[(0,)+c for c in itertools.combinations(range(1,14),7)]
rows=[]
for E in sets:
    rows.append((E, float(p0(E)), add_energy(E), spectral_flatness(E), diffset_size(E), max(E)))
P=np.array([r[1] for r in rows]); AE=np.array([r[2] for r in rows], float)
SF=np.array([r[3] for r in rows]); DS=np.array([r[4] for r in rows], float); WD=np.array([r[5] for r in rows], float)
def tau(x): return kendalltau(x, P).correlation
print(f"  #sets = {len(rows)}")
print(f"  tau(additive energy AE, p0) = {tau(AE):+.4f}")
print(f"  tau(spectral flatness SF, p0) = {tau(SF):+.4f}   (expect NEGATIVE: flat spectrum -> low cover)")
print(f"  tau(diff-set size |E-E|, p0) = {tau(DS):+.4f}   (spread)")
print(f"  tau(width max(E), p0)        = {tau(WD):+.4f}")
# argmax/argmin sanity
imaxP=int(np.argmax(P)); print(f"\n  p0-argmax set = {rows[imaxP][0]}  p0={P[imaxP]:.4f}  (is it consecutive {tuple(range(8))}?)")
imaxAE=int(np.argmax(AE)); iminSF=int(np.argmin(SF))
print(f"  AE-argmax set = {rows[imaxAE][0]}   SF-argmin set = {rows[iminSF][0]}")
# per the repo claim, restrict to |tau| and count strict inversions AE vs SF
def inversions(x, sign):
    # count pairs where predictor (scaled by sign) disagrees with p0 strict order
    xs=sign*x; n=len(xs); bad=0; tot=0
    # sample to keep it O(n^2) but cheap: full is 1.47M, fine
    for i in range(n):
        di=xs[i]-xs[i+1:]; dp=P[i]-P[i+1:]
        conc=(di*dp); tot+=len(dp)
        bad+=int(np.sum((conc<0)))
    return bad, tot
biAE,tot=inversions(AE,+1); biSF,_=inversions(SF,-1); biDS,_=inversions(DS,+1)
print(f"\n  strict discordant pairs (out of {tot}):")
print(f"    AE (+): {biAE} ({100*biAE/tot:.1f}%)   SF (-): {biSF} ({100*biSF/tot:.1f}%)   |E-E| (+): {biDS} ({100*biDS/tot:.1f}%)")
best=min([("additive energy",abs(tau(AE))),("spectral flatness",abs(tau(SF))),("diff-set size",abs(tau(DS)))], key=lambda t:-t[1])
print(f"  => BEST single predictor by |tau|: {best[0]} (|tau|={best[1]:.4f})")

print("\n"+"="*80)
print("(B) DOES THE DFT OF THE EXTREMIZER PREDICT WHERE w*Delta_w PEAKS?")
print("    E' = {0,1,2,30,31,32,60,61}  (= clean tower B(+)30B minus the element 62)")
print("="*80)
Ep=[0,1,2,30,31,32,60,61]
# B1: defected-Riesz factorization  \hat 1_{E'} = \hat 1_B(xi)\hat 1_B(30 xi) - e(62 xi)
B=[0,1,2]
def fh(S,xi): return sum(np.exp(2j*np.pi*s*xi) for s in S)
xis=np.array([k/97.0 for k in range(1,60)])
lhs=np.array([fh(Ep,x) for x in xis])
rhs=np.array([fh(B,x)*fh(B,30*x)-np.exp(2j*np.pi*62*x) for x in xis])
print(f"  B1 defected-Riesz identity  max|\\hat1_E' - (\\hat1_B(xi)\\hat1_B(30xi) - e(62xi))| = {np.max(np.abs(lhs-rhs)):.2e}")
print(f"     => the extremizer IS an exact 2-fold Riesz product MINUS one term (a 'defected' Riesz product).")
# B2: power spectrum of E' -- where are its side-peaks?  (predict multiples of 1/30)
N=1800; jj=np.arange(1,N//2)
power=np.abs(np.array([sum(np.exp(2j*np.pi*s*j/N) for s in Ep) for j in jj]))**2
# local maxima
loc=[j for j in range(1,len(power)-1) if power[j]>power[j-1] and power[j]>=power[j+1]]
loc.sort(key=lambda j:-power[j])
print(f"  B2 top side-peaks of |\\hat 1_E'(xi)|^2  (xi = j/{N}); nearest k/30 and k/60:")
for j in loc[:6]:
    xi=jj[j]/N; k30=round(xi*30); k60=round(xi*60)
    print(f"       xi={xi:.4f}  power={power[j]:5.1f}   |xi - {k30}/30|={abs(xi-k30/30):.4f}   |xi - {k60}/60|={abs(xi-k60/60):.4f}")
# B3: w*Delta_w peaks, and the DFT of the discrepancy sequence in w
cells=cells_of(Ep); W=630
wd=np.array([float(wD(cells,w)) for w in range(1,W+1)])
order=np.argsort(-np.abs(wd))
print(f"  B3 top |w*Delta_w| peaks over w=1..{W}:  (w mod 7, w mod 30)")
for idx in order[:8]:
    w=idx+1; print(f"       w={w:4d}  |wD|={abs(wd[idx]):.4f}   (w mod 7 = {w%7}, w mod 30 = {w%30})")
# FFT of the discrepancy sequence -> dominant periods in w
sp=np.abs(np.fft.rfft(wd - wd.mean()))
freqs=np.fft.rfftfreq(len(wd))
top=np.argsort(-sp)[:8]
print(f"  B3' dominant FREQUENCIES of the w->w*Delta_w sequence (period = 1/freq in units of w):")
for t in top:
    if freqs[t]==0: continue
    per=1/freqs[t]; print(f"       period ~ {per:6.2f} in w   (amp {sp[t]:.1f})   [7?{abs(per-7)<0.6}  30?{abs(per-30)<2}  15?{abs(per-15)<1}]")
print("\nDONE.")
