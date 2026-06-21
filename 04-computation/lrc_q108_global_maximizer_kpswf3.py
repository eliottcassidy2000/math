#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc_q108_global_maximizer_kpswf3.py   (kind-pasteur 2026-06-21, THREAD D -> HYP)

GLOBAL CONSECUTIVE-MAXIMIZER (canon convention: 0 may be in E; e=0 always covers
sector 0).  CLAIM:  among all k-subsets E of {0,1,2,...} (primitive: gcd=1),
   measS7(E)  is maximized by the consecutive block  C_k = {0,1,...,k-1}.

This UNIFIES the LRC(14) ledger: the wide-bound max (single-far {0,2,..,14,29}=0.372
at k=9) is BELOW consec_9 (0.438); consec is narrow (span k-1<=14 for k<=15) so it
lives in the L3 finite-check branch.  If consec is the GLOBAL argmax then:
   * the binding constraint of the whole LRC(14) cover bound is the single family C_k;
   * L4-L7 (wide cases) are all DOMINATED -> "everything is <= consec < cap".

Tests (exact rational, canon convention with 0):
 (1) measS7(C_k) and cap_k, margins, for k=7..12.
 (2) Exhaustive: over all primitive k-subsets of {0..N}, is C_k the unique argmax?
     (k=8 N up to 15; k=9 N up to 14).  Report any set exceeding consec.
 (3) Dilation invariance with 0: measS7(cE)=measS7(E)?  (the search-space reducer)
 (4) Spread stress: AP {0,d,2d,...}, two-cluster {0..t}u{f..}, vs consec.
"""
import itertools, random
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
    E=[int(e) for e in E]
    has0 = 0 in E
    Enz=[e for e in E if e!=0]
    xs=breakpoints(Enz) if Enz else [Fr(0),Fr(1)]
    tot=Fr(0)
    for a,b in zip(xs,xs[1:]):
        mid=(a+b)/2
        secs=set()
        if has0: secs.add(0)
        for e in Enz: secs.add(sector((e*mid)%1))
        if len(secs)==P: tot+=(b-a)
    return tot

CAP={8:Fr(2243,5880),9:Fr(1979,4004),10:Fr(55,91)}

def gcdl(L):
    g=0
    for x in L: g=gcd(g,x)
    return g

def main():
    print("#"*80)
    print("# GLOBAL CONSECUTIVE-MAXIMIZER  C_k={0,..,k-1}  (canon convention, 0 in E)")
    print("#"*80)

    # (1) consec vs cap
    print("\n=== (1) measS7(C_k) vs cap_k ===")
    cons_seq={}
    for k in range(7,13):
        Ck=list(range(0,k)); v=measS7(Ck); cons_seq[k]=v; cap=CAP.get(k)
        print(f"  k={k}: measS7({{0..{k-1}}})={float(v):.5f}"
              f"{('  cap='+str(round(float(cap),5))+'  margin='+str(round(float(cap-v),5))+'  '+('UNDER' if v<=cap else 'OVER!')) if cap else ''}")

    # (3) dilation invariance with 0
    print("\n=== (3) dilation invariance measS7(cE)=measS7(E) (canon, with 0) ===")
    ok=True
    for E in [[0,1,2,3,4,5,6,7],[0,1,3,4,9],[0,2,5,11,13]]:
        base=measS7(E)
        for c in [2,3,5,7,11]:
            if measS7([c*e for e in E])!=base: ok=False
    print(f"  measS7(cE)==measS7(E) for c in {{2,3,5,7,11}}: {'HOLDS' if ok else 'FAILS'}")

    # (2) exhaustive argmax (primitive subsets of {0..N})
    print("\n=== (2) EXHAUSTIVE: is C_k the global argmax over primitive k-subsets of {0..N}? ===")
    for k,N in [(8,15),(9,14)]:
        Ck=list(range(0,k)); pc=cons_seq[k]
        nbeat=[]; tot=0; maxv=pc; maxE=tuple(Ck)
        for combo in itertools.combinations(range(0,N+1),k):
            # primitive: gcd of nonzero =1 (dilation-normalize); and to dedupe translates
            nz=[e for e in combo if e!=0]
            if gcdl(nz)!=1: continue
            tot+=1
            v=measS7(combo)
            if v>pc+Fr(1,10**9): nbeat.append((combo,float(v)))
            if v>maxv: maxv=v; maxE=combo
        print(f"  k={k}, N={N}: {tot} primitive subsets; consec={float(pc):.5f}; "
              f"global max={float(maxv):.5f} at {maxE}")
        if nbeat:
            print(f"     !!! {len(nbeat)} BEAT consec; top: {sorted(nbeat,key=lambda t:-t[1])[:3]}")
        else:
            print(f"     CONSEC IS THE UNIQUE GLOBAL ARGMAX (over {{0..{N}}})")

    # (4) spread stress vs consec
    print("\n=== (4) spread stress: nothing wide should exceed consec ===")
    for k in [8,9]:
        pc=cons_seq[k]
        worst=pc; worstdesc="consec"
        # AP with diff d (dilation-equiv to consec only if... actually {0,d,2d,..}=d*{0,1,..} => equiv)
        for d in [1,2,3,5,7,11,13]:
            E=[i*d for i in range(k)]; v=measS7(E)
            if v>worst+Fr(1,10**9): worst=v; worstdesc=f"AP d={d}"
        # two-cluster {0..t-1} u {f..f+(k-t)-1}
        t=k//2
        for f in [10,15,20,30,50,100,300]:
            E=list(range(0,t))+list(range(f,f+(k-t))); v=measS7(E)
            if v>worst+Fr(1,10**9): worst=v; worstdesc=f"2clust f={f}"
        # single far: consec_{k-1} u {f}
        for f in [15,20,29,50,100]:
            E=list(range(0,k-1))+[f]; v=measS7(E)
            if v>worst+Fr(1,10**9): worst=v; worstdesc=f"single-far f={f}"
        print(f"  k={k}: consec={float(pc):.5f}; best spread config={float(worst):.5f} ({worstdesc})"
              f"  {'consec dominates' if worst<=pc+Fr(1,10**9) else 'SPREAD BEATS consec!'}")

    print("\nDONE.")

if __name__=="__main__":
    main()
