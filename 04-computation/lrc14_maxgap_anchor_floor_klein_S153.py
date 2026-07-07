#!/usr/bin/env python3
"""
klein-2026-07-07-S153.  CONSTRUCTIVE anchor-floor for  inf_E E[maxgap] > 1/7  (extends kps-S57).

kps-S57 (HYP-4747): Route-1 density floor  mu_{1/7} = P_x(maxgap{frac(e_i x)} > 1/7)
  >= (7/6)(E[maxgap] - 1/7)   [reverse Markov],  reducing to the open target inf_E E[maxgap] > 1/7,
with the partial E[origingap]=1/7 (=> E[maxgap] >= 1/7 NON-strict). Remaining: "max beats typical
with margin."

klein-S153 NEW: a CONSTRUCTIVE, provable lower bound on E[maxgap] via a FINITE ANCHOR SET.
For any finite A = {a_1,...,a_r} subset of the circle,
      maxgap(x)  >=  max_{a in A} gap_a(x)         (the max gap dominates any single anchor's gap)
so   E[maxgap]  >=  E[ max_{a in A} gap_a(x) ].
KEY IDENTITY (verify): E[gap_a(x)] = 2/(k+1) = 1/7 at k=13, for EVERY anchor a and EVERY family
(length-biased origin gap; structure-independent). Then
      E[max_{a in A} gap_a]  =  1/7  +  E[ max_a gap_a  -  gap_{a_1} ]  >  1/7
with a strictly positive margin whenever the anchor gaps are not perfectly correlated. As |A| grows,
E[max_a gap_a] increases to E[maxgap]. So a FINITE anchor set gives a PROVABLE floor > 1/7 -- the
"max beats typical" margin made constructive.

This script: (1) verifies E[gap_a] = 2/(k+1) structure-independent; (2) computes the anchor floors
E[max_{a in A} gap_a] for growing A and shows they exceed 1/7 with margin and rise to E[maxgap];
(3) finds the WORST family (min E[maxgap]) = the AP, and its anchor floor margin.
"""
import math, random

def points(E, x):
    return sorted((e*x) % 1.0 for e in E)

def gaps_of(pts):
    n=len(pts); return [ (pts[(i+1)%n]-pts[i]) % 1.0 if i<n-1 else (pts[0]+1.0-pts[-1]) for i in range(n) ]

def gap_at(pts, a):
    """length of the gap (of the cyclic point set) containing anchor a in [0,1)."""
    n=len(pts)
    a = a % 1.0
    # find the two consecutive points straddling a
    # nearest point <= a (cyclically) and next point
    for i in range(n):
        lo = pts[i]; hi = pts[(i+1)%n] if i<n-1 else pts[0]+1.0
        aa = a if a>=pts[0] else a+1.0
        if lo <= aa < hi:
            return hi-lo
    # a below pts[0]: in wrap gap
    return (pts[0]+1.0-pts[-1])

def analyze(E, NG=100000, anchorsets=None):
    k=len(E)
    if anchorsets is None:
        anchorsets = {
            "{0}":[0.0],
            "{0,1/2}":[0.0,0.5],
            "{0,1/2,1/4,3/4}":[0.0,0.5,0.25,0.75],
            "{j/8}":[j/8 for j in range(8)],
            "{j/16}":[j/16 for j in range(16)],
            "{j/40}":[j/40 for j in range(40)],
        }
    sum_maxgap=0.0
    sum_gap_a={name:0.0 for name in anchorsets}        # E[max over that anchor set]
    sum_ega=0.0; ne=0                                   # E[gap at a single anchor 0], for the identity
    sum_ega_half=0.0
    for j in range(NG):
        x=(j+0.5)/NG
        pts=points(E,x)
        g=gaps_of(pts)
        mg=max(g)
        sum_maxgap+=mg
        sum_ega += gap_at(pts,0.0)
        sum_ega_half += gap_at(pts,0.3178)   # a generic anchor
        for name,A in anchorsets.items():
            sum_gap_a[name]+=max(gap_at(pts,a) for a in A)
    res={"Emaxgap":sum_maxgap/NG, "Ega0":sum_ega/NG, "Ega_gen":sum_ega_half/NG}
    for name in anchorsets:
        res["anchor "+name]=sum_gap_a[name]/NG
    return res

random.seed(7)
THRESH=1.0/7.0
print("="*88)
print(f"ANCHOR FLOOR for E[maxgap]  (threshold 1/7={THRESH:.5f};  2/(k+1) identity)")
print("="*88)

def rand_family(k):
    while True:
        s=sorted(random.sample(range(1,60),k))
        if math.gcd(*s)==1: return s

for k in [8,10,13]:
    print(f"\n----- k={k}  (2/(k+1) = {2/(k+1):.5f}) -----")
    fams=[("AP {1..k}",list(range(1,k+1)))]
    for r in range(3): fams.append((f"spread#{r}",rand_family(k)))
    for name,E in fams:
        R=analyze(E, NG=60000)
        idcheck = f"E[gap@0]={R['Ega0']:.4f} E[gap@gen]={R['Ega_gen']:.4f} (2/(k+1)={2/(k+1):.4f})"
        print(f"  {name:14s} E[maxgap]={R['Emaxgap']:.4f}  | {idcheck}")
        line="      anchor floors: "
        for a in ["{0}","{0,1/2}","{0,1/2,1/4,3/4}","{j/8}","{j/16}","{j/40}"]:
            line+=f"{a}={R['anchor '+a]:.4f}  "
        print(line)

print()
print("="*88)
print("WORST-CASE (min E[maxgap]) hunt at k=13 -- is the AP the minimizer, and does a FINITE")
print("anchor set already exceed 1/7 for it (=> provable floor)?")
print("="*88)
# focus AP at k=13, finer grid, bigger anchor sets
AP=list(range(1,14))
big_anchors={f"{{j/{m}}}":[j/m for j in range(m)] for m in [8,16,32,64,128]}
R=analyze(AP, NG=120000, anchorsets=big_anchors)
print(f"AP k=13: E[maxgap]={R['Emaxgap']:.5f}, E[gap@0]={R['Ega0']:.5f} (1/7={THRESH:.5f})")
for m in [8,16,32,64,128]:
    name=f"{{j/{m}}}"
    v=R['anchor '+name]
    print(f"   anchor set {name:9s} (|A|={m:3d}): E[max_a gap_a]={v:.5f}   margin over 1/7 = {v-THRESH:+.5f}"
          + ("  >1/7 OK" if v>THRESH else "  <=1/7"))
print()
print("READOUT: if E[gap@a]=2/(k+1)=1/7 is structure-independent, and a FINITE anchor set gives")
print("E[max_a gap_a] > 1/7 for the AP (the minimizer), then inf_E E[maxgap] > 1/7 is provable via")
print("that finite anchor floor (each anchor's mean is 1/7; the max over anchors is a computable")
print("order statistic > 1/7). This turns kps-S57's open 'margin' into a finite computation.")
