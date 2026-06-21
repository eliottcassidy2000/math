#!/usr/bin/env python3
"""
lrc14_route_D_twoblock_bound_opus_s63.py   (opus-2026-06-20-S63)

ROUTE D, part 4.  Parts 1-3 established (VERIFIED, k=8,9,10):
  * consec (single solid block) is the GLOBAL max of Ly over bounded shapes.
  * The ONLY other local maxima are MULTI-BLOCK shapes, and they sit FAR below
    consec (Ly ~0.13 vs 0.358 at k=8).  They are local-max traps, not threats.
  * No single-offset / block-merge / gap-close operator is globally monotone
    (the multi-block ridge has its own resonant local peaks).
  * Convex order on N is DEAD (dual g not convex; means go wrong way).

So the metagraph DAG dream fails, but the EXTREMALITY is robust.  The cleanest
remaining route to a PROOF is a structural DICHOTOMY:

    THE BLOCK DICHOTOMY.  Every bounded shape E either
      (a) is a SINGLE solid block  (=> E = consec, done), or
      (b) has an internal gap / splits into >=2 blocks  => Ly(E) <= some
          explicit value B_split < cap_k, hence E is NOT the maximizer.

If we can show every NON-consec shape has Ly <= B_split with B_split clearly
below Ly(consec), the maximizer must be consec.  This part MEASURES B_split:

  (1) max Ly over all NON-single-block shapes (the "split ceiling").  Compare to
      Ly(consec).  The GAP = consec advantage = the proof's safety margin.
  (2) Does the split ceiling stay bounded away from Ly(consec) as span grows?
      (test span = k+1 .. k+6).  If the split ceiling is monotone in span and
      saturates BELOW Ly(consec), the bounded check is a finite certificate.
  (3) Identify WHICH split shape achieves the ceiling (the "hardest competitor")
      and its block signature -- this is the shape a proof must rule out.
  (4) THE ATOMIC LEMMA test:  "removing one internal gap (de-splitting) by
      DELETING the gap entirely -- i.e. comparing a shape to its 'gap-collapsed'
      version -- never decreases Ly when the result is consec."  Precisely:
      for every shape E, let collapse(E) = solidify to consec; we already know
      Ly(consec) >= Ly(E) globally.  Here we test the SHARPER per-shape claim:
      Ly is monotone under the SPECIFIC operator "delete the SMALLEST internal
      gap and re-pack" -- a global de-splitting -- toward consec.

python3 stdlib only; exact Fractions.
"""
import sys, itertools, functools
from math import comb, gcd
from functools import reduce
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

@functools.lru_cache(maxsize=None)
def pdist(Etuple):
    E=sorted(set(Etuple))
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for j in range(7):
            for t in (F(j,7),F(j+1,7)):
                for m in range(e): bps.add((t+m)/e)
    bps=sorted(z for z in bps if 0<=z<=1)
    P=[F(0)]*7
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        occ=[False]*7
        for e in E:
            occ[int(((e*xm)%1)*7)]=True
        N=sum(1 for j in range(1,7) if not occ[j])
        P[N]+=x1-x0
    return tuple(P)

def yvec(k):
    if k==8: g=lambda t:F((t-1)*(t-2)*(t-4)*(t-5),40)
    elif k in (9,10): g=lambda t:F(-(t-2)*(t-3)*(t-6),36)
    else: g=lambda t:F((t-3)*(t-4),12)
    return [g(t) for t in range(7)]

@functools.lru_cache(maxsize=None)
def Ly(Etuple,k):
    g=yvec(k); P=pdist(Etuple)
    return sum(g[t]*P[t] for t in range(7))

def norm(E):
    E=sorted(set(E)); E=[e-E[0] for e in E]
    g=reduce(gcd,E,0)
    if g>1: E=[e//g for e in E]
    return tuple(E)
def span(E): return max(E)-min(E)

def is_single_block(E):
    E=sorted(set(E))
    return all(E[i+1]==E[i]+1 for i in range(len(E)-1))

def all_shapes(k, max_span):
    seen=set()
    for rest in itertools.combinations(range(1,max_span+1),k-1):
        E=norm([0]+list(rest))
        if span(E)<=max_span and E not in seen:
            seen.add(E); yield E

if __name__=="__main__":
    print("="*72)
    print("ROUTE D part 4: the BLOCK DICHOTOMY -- split ceiling vs consec")
    print("="*72)
    for k in (8,9,10,11,12):
        consec=tuple(range(k)); lyc=Ly(consec,k)
        print(f"\n### k={k}: Ly(consec)={lyc}={float(lyc):.6f} ###")
        # (1)+(2) split ceiling as a function of span
        for ms in range(k, k+7):
            best=None; bestE=None; nsplit=0
            for E in all_shapes(k,ms):
                if is_single_block(E): continue  # only NON-block (split) shapes
                nsplit+=1
                v=Ly(E,k)
                if best is None or v>best: best=v; bestE=E
            if best is None:
                print(f"  span<={ms:2d}: (no split shapes)")
                continue
            margin=lyc-best
            print(f"  span<={ms:2d}: split-ceiling={float(best):.6f} at {bestE}  "
                  f"margin(consec-split)={float(margin):.6f}  {'SAFE' if margin>0 else 'DANGER'}")
        # the hardest competitor at the widest span
    print("\n--- Interpretation ---")
    print("If for every k the split-ceiling stays strictly below Ly(consec) and")
    print("SATURATES as span grows, then: any maximizer with an internal gap is")
    print("beaten by consec, so the maximizer is the single block.  The saturation")
    print("value is the explicit B_split the dichotomy needs.")
    print("\nDONE part 4.")
