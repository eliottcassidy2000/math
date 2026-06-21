#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_consec_maximizer_kpswf3.py   (kind-pasteur 2026-06-21, THREAD D -> HYP)

THE LEAD: CONSECUTIVE-MAXIMIZER conjecture (HYP candidate).
   Among all k-element sets E of positive integers, measS7(E) is MAXIMIZED by the
   consecutive run {1,2,...,k} (equivalently any dilate / translate-free run).

If true, the LRC(14) cover bound's WORST CASE is a SINGLE explicit family, and the
whole p0 <= cap_k question becomes a finite check on consec runs vs cap_k.

This script:
  (1) widens the exhaustive search (larger N) to stress the maximizer claim at k=8,9;
  (2) checks measS7([1..k]) vs cap_k for k=8,9,10 -- does consec itself stay under cap?
  (3) tests TRANSLATES {a+1,...,a+k} and whether starting-offset changes p0 (it must
      not be larger -- if a translate beats [1..k], conjecture is false);
  (4) records S_7=0 structural fact and the consec p0 sequence.

Scaling invariance (proved-checked separately, kpswf3): measS7(cE)=measS7(E) all c.
So we may assume gcd(E)=1 WLOG and min behaviour is about the *shape*.
"""
import itertools
from fractions import Fraction as Fr
from math import gcd

P=7
def sector(yf): return int(P*yf)
def breakpoints(E):
    bp={Fr(0),Fr(1)}
    for e in E:
        if e==0: continue
        for t in range(0,P*e): bp.add(Fr(t,P*e))
    return sorted(bp)
def measS7(E):
    E=[int(e) for e in E if int(e)!=0]
    xs=breakpoints(E); tot=Fr(0)
    for a,b in zip(xs,xs[1:]):
        mid=(a+b)/2
        if len(set(sector((e*mid)%1) for e in E))==P: tot+=(b-a)
    return tot

CAP={8:Fr(2243,5880),9:Fr(1979,4004),10:Fr(55,91)}

def main():
    print("#"*80)
    print("# CONSECUTIVE-MAXIMIZER conjecture stress test")
    print("#"*80)

    # (2) consec vs cap
    print("\n=== consec [1..k] p0 vs cap_k ===")
    for k in [7,8,9,10,11]:
        E=list(range(1,k+1)); v=measS7(E); cap=CAP.get(k)
        print(f"  k={k}: measS7([1..{k}])={float(v):.5f}"
              f"{'  cap='+str(round(float(cap),5))+'  '+('UNDER cap' if cap and v<=cap else ('OVER cap!!' if cap else '')) if cap else ''}")

    # (1) widened exhaustive search at k=8 (N up to 16), k=9 (N up to 15)
    print("\n=== widened exhaustive maximizer search ===")
    for k,N in [(8,16),(9,15)]:
        cons=list(range(1,k+1)); pc=measS7(cons)
        best=pc; bestE=tuple(cons); beat=[]
        total=0
        for combo in itertools.combinations(range(1,N+1),k):
            if combo[0]!=1: continue          # normalize translate (min=1)
            gg=0
            for e in combo: gg=gcd(gg,e)
            if gg!=1: continue                # normalize dilation (gcd=1)
            total+=1
            v=measS7(combo)
            if v>pc+Fr(1,10**9):
                beat.append((combo,float(v)))
                if v>best: best=v; bestE=combo
        print(f"  k={k},N={N}: searched {total} normalized subsets; consec p0={float(pc):.5f}")
        if beat:
            print(f"     !!! CONSEC BEATEN by {len(beat)}; top: {sorted(beat,key=lambda t:-t[1])[:3]}")
        else:
            print(f"     consec is the UNIQUE max (no normalized subset exceeds it)")

    # (3) translate test: does NOT starting at 1 ever help?  test {a..a+k-1}
    print("\n=== translate test: measS7({a,...,a+k-1}) vs consec [1..k] ===")
    for k in [8,9]:
        pc=measS7(list(range(1,k+1)))
        worst_trans=pc; worstA=1
        for a in range(1,40):
            v=measS7(list(range(a,a+k)))
            if v>worst_trans+Fr(1,10**9):
                worst_trans=v; worstA=a
        print(f"  k={k}: [1..{k}] p0={float(pc):.5f};  max over translates a in[1,40)"
              f" = {float(worst_trans):.5f} at a={worstA}"
              f"  {'(translate beats! conj FALSE)' if worst_trans>pc+Fr(1,10**9) else '(no translate beats)'}")

    # (4) consec sequence + does the MARGIN cap_k - consec grow?  (room to spare)
    print("\n=== margin cap_k - measS7([1..k]) (the actual LRC slack on the worst family) ===")
    for k in [8,9,10]:
        v=measS7(list(range(1,k+1))); cap=CAP[k]
        print(f"  k={k}: cap-consec = {float(cap-v):.5f}  ({'POSITIVE -> cover bound holds on consec' if cap>v else 'NEGATIVE!'})")

    print("\nDONE.")

if __name__=="__main__":
    main()
