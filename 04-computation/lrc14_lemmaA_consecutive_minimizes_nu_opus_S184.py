"""
lrc14_lemmaA_consecutive_minimizes_nu_opus_S184.py  (opus-2026-07-09-S184)

LEMMA A (the lone open hfloor node): nu(E) >= nuConsec(k) -- the CONSECUTIVE co-offset cluster
{0,1,..,k-1} MINIMIZES nu(E)=meas{x: maxgap(frac(e_i x)) > 1/7} (the good-set measure). This is the
MEASURE-level version of "the AP/interval is the extremal / hardest" (my LEM-015 = discrete Schur-triple
version). VERIFY: (a) consecutive achieves (near) the min nu; (b) nu tracks INVERSELY with the additive
structure E3=#Schur-triples of the cluster -- more resonance (AP) => orbit more equidistributed => smaller
maxgap => smaller good-set => smaller nu. This gives the fleet the conceptual bridge for the crux: the AP
minimizes nu BECAUSE it maximizes E3 (additive coherence = resonance = hardness, per kps-S113/mac-mini-S68).
"""
import numpy as np
from collections import Counter

NG=200003
X=(np.arange(1,NG))/NG   # x in (0,1)
def maxgap_measure(E):
    """nu(E) = fraction of x with maxgap{frac(e*x)} > 1/7."""
    E=np.array(sorted(set(E)))
    good=np.zeros(len(X),dtype=bool)
    # for each x compute maxgap of {frac(e*x)}: sort the k values, circular gaps
    P=(np.outer(E,X))%1.0           # k x |X|
    P=np.sort(P,axis=0)
    gaps=np.diff(P,axis=0)
    wrap=(P[0]+1.0)-P[-1]
    mg=np.maximum(gaps.max(axis=0),wrap)
    return float((mg>1.0/7).mean())
def E3(S):
    Sset=set(S); return sum(1 for a in S for b in S if a+b in Sset)

k=13
clusters=[
 ("consecutive {0..12}", list(range(0,13))),
 ("interval {1..13}", list(range(1,14))),
 ("near-2-AP {0,2,3..12,14}", [0,2,3,4,5,6,7,8,9,10,11,12,14]),
 ("2*AP {0,2,..24}", [2*i for i in range(0,13)]),
 ("Sidon-ish", [0,1,3,7,12,20,30,44,65,80,96,122,147]),
 ("random A", [0,5,9,13,18,26,31,40,55,63,77,91,110]),
 ("GAP 4x4", sorted(set(a+5*b for a in range(4) for b in range(4)))[:13]),
]
print("="*88)
print(f"LEMMA A: does the CONSECUTIVE cluster minimize nu(E)=meas{{maxgap>1/7}}?  (k={k})")
print("="*88)
print(f"  {'cluster':>26} {'nu(E)':>9} {'E3 Schur':>9}")
rows=[]
for name,E in clusters:
    if len(set(E))!=13: 
        print(f"  {name:>26}  (size {len(set(E))}, skip)"); continue
    nu=maxgap_measure(E); e3=E3(E); rows.append((name,nu,e3))
nu_consec=[r[1] for r in rows if "consecutive" in r[0]][0]
for name,nu,e3 in rows:
    flag=" <-- MIN?" if abs(nu-min(r[1] for r in rows))<1e-4 else ""
    print(f"  {name:>26} {nu:>9.4f} {e3:>9}{flag}")
print("-"*88)
print(f"  nu(consecutive) = {nu_consec:.4f}; min over sample = {min(r[1] for r in rows):.4f}")
nus=[r[1] for r in rows]; e3s=[r[2] for r in rows]
print(f"  corr(E3, nu) = {np.corrcoef(e3s,nus)[0,1]:+.3f}  (expect NEGATIVE: more Schur => smaller good-set)")
print("  READING: if consecutive ~ min nu AND corr(E3,nu)<0, then the AP minimizes nu BECAUSE it")
print("  maximizes E3 (LEM-015) -- the measure-level bridge for Lemma A (the lone open hfloor crux).")
print("="*88)
