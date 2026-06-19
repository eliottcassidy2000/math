#!/usr/bin/env python3
"""
lrc14_sector_majorize_macmini_0618s7b.py  (mac-mini-2026-06-18-S7, ANGLE B)

The AP partial-sum-walk picture: for E={0..k-1}, residues hit at theta = {floor(e theta) mod 7}
= partial sums of a mechanical word.  For GENERAL E, residues hit = {floor(e theta) mod 7: e in E}
= a SUBSET of the AP's partial sums (those at indices e in E).  KEY OBSERVATION TO TEST:

  For ANY E with max(E)=k-1 (same span as AP_k) but |E| possibly < k, the residues hit are a
  SUBSET of the AP_k residues at the same theta.  So at every theta, AP_k covers >= any sub-E.
  => meas(S7(E)) <= meas(S7(AP_{max(E)+1})) when E subset {0..max E}.  (POINTWISE domination!)

But the competition is over |E|=k FIXED, E primitive, NOT E subset {0..k-1}. The AP_k has
span k-1 (minimal span for k points). Any other |E|=k primitive has span >= k-1, with EQUALITY
iff E=AP_k (consecutive). Larger span = the offsets are 'stretched', dilating the mechanical
word -- does stretching ALWAYS lose coverage?

TESTS:
 (1) POINTWISE subset domination: for E subset {0..k-1}, meas(S7(E)) <= meas(S7(AP_{k})) at the
     SAME k cells? Confirm meas(S7(E)) <= meas(S7(AP_{maxE+1})) (E with a hole loses).
 (2) The real claim: among |E|=k primitive, span(E) MINIMAL (=k-1, AP) maximizes meas(S7).
     Test: sort competitors by span; is meas(S7) monotone-ish decreasing in span?
 (3) A CLEANER monotonicity: meas(S7) as a function of E under "spreading moves" (e_i -> e_i+1
     for the top element, keeping primitivity). Does each spreading move decrease meas(S7)?
     If meas(S7) is monotone under elementary spreads, AP (the most-compressed) is the max.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)

def measS7(E):
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in Enz:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        res=set(int(7*e*xm)%7 for e in E)
        if len(res)==7: total+=x1-x0
    return total

def span(E): return max(E)-min(E)
def primitive(E):
    g=0
    for e in E: g=gcd(g,e)
    return g==1

# (1) pointwise subset domination
print("="*92)
print("(1) POINTWISE subset domination: E subset {0..k-1} => meas(S7(E)) <= meas(S7(AP_k))")
print("    (residues of E are a subset of AP_k residues at every theta)")
print("="*92)
for k in [8,9,10]:
    AP=tuple(range(k)); apv=measS7(AP)
    ok=True; worst=None
    # all subsets of {0..k-1} containing 0, of size 7..k (need at least 7 to cover)
    for sz in range(7,k+1):
        for sub in itertools.combinations(range(1,k),sz-1):
            E=(0,)+sub
            v=measS7(E)
            if v>apv: ok=False; worst=(E,v)
    print(f"  k={k}: AP_k={float(apv):.6f}; all subsets E of {{0..{k-1}}} have meas(S7(E))<=AP_k? {ok}"
          + ("" if ok else f"  VIOLATION {worst}"))

# (2) span monotonicity among |E|=k primitive
print()
print("="*92)
print("(2) Among |E|=k primitive normalized (min=0): is meas(S7) decreasing in span?")
print("    span=k-1 unique to AP. Show AP is max; show correlation meas(S7) vs span.")
print("="*92)
def gen(k,maxE):
    out=[]
    for rest in itertools.combinations(range(1,maxE+1),k-1):
        E=(0,)+rest
        if not primitive(E): continue
        out.append(E)
    return out
for k in [8,9]:
    maxE=k+3
    shapes=gen(k,maxE)
    data=[(measS7(E),span(E),E) for E in shapes]
    data.sort(reverse=True)
    AP=tuple(range(k)); apv=measS7(AP)
    # bucket max meas by span
    byspan={}
    for v,s,E in data:
        if s not in byspan or v>byspan[s][0]: byspan[s]=(v,E)
    print(f"  k={k}, box maxE<={maxE}: AP meas={float(apv):.6f} (span {k-1})")
    print(f"    top overall: {[(float(v),E) for v,s,E in data[:3]]}")
    print(f"    max meas(S7) per span: " + ", ".join(f"span{s}:{float(byspan[s][0]):.4f}" for s in sorted(byspan)))

# (3) elementary spreading move monotonicity
print()
print("="*92)
print("(3) Elementary spread monotonicity: replace the largest element e_max -> e_max+1.")
print("    Does meas(S7) decrease (weakly) under each such spread?  Test from AP outward.")
print("="*92)
for k in [8,9]:
    E=list(range(k)); v0=measS7(tuple(E))
    print(f"  k={k}: start AP meas={float(v0):.6f}")
    cur=E[:]; prev=v0; chain=[(tuple(cur),v0)]
    mono=True
    for step in range(6):
        cur=cur[:]; cur[-1]+=1  # spread the top
        if not primitive(tuple(cur)):
            # skip to next primitive
            while not primitive(tuple(cur)): cur[-1]+=1
        v=measS7(tuple(cur))
        dec = v<=prev
        if not dec: mono=False
        chain.append((tuple(cur),v))
        print(f"     -> E={tuple(cur)} meas={float(v):.6f}  {'(dec)' if dec else '(*** INCREASE ***)'}")
        prev=v
    print(f"     spreading-top monotone decreasing: {mono}")
print("\nDONE.")
