#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_verify_reduction_1_DE_kps-S2-wf.py  -- sections (D) decisive hunt + (E) floor audit,
plus a new (F) EXTREME-PUSH: drive M(S) as close to 1/14 from above as possible with
coordinated large speeds, and a (G) k=4 / k=5 coordinated stress test, all bounded so
they terminate.  Exact fractions throughout.
"""
import sys
from math import gcd
from fractions import Fraction as F
from functools import reduce
import itertools as it

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def gmin(S, t): return min(nrm(v * t) for v in S)
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def exact_M(S):
    b = F(0); at = None
    for t in cand(S):
        v = gmin(S, t)
        if v > b: b = v; at = t
    return b, at

def danger_arcs(v):
    w = F(1, 14 * v); out = []
    for k in range(v + 1):
        c = F(k, v); lo, hi = c - w, c + w
        if lo < 0: out.append((F(0), hi)); out.append((1 + lo, F(1)))
        elif hi > 1: out.append((lo, F(1))); out.append((F(0), hi - 1))
        else: out.append((lo, hi))
    return [(a, b) for (a, b) in out if a < b]
def union_arcs(intervals):
    iv = sorted(intervals)
    if not iv: return []
    out = []; cl, ch = iv[0]
    for lo, hi in iv[1:]:
        if lo <= ch:
            if hi > ch: ch = hi
        else: out.append((cl, ch)); cl, ch = lo, hi
    out.append((cl, ch)); return out
def safe_arcs(S):
    danger = union_arcs([iv for v in set(S) for iv in danger_arcs(v)])
    out = []; cur = F(0)
    for lo, hi in danger:
        if lo > cur: out.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1: out.append((cur, F(1)))
    return out
def meas(arcs): return sum((b - a for a, b in arcs), F(0))
def L(S): return meas(safe_arcs(S))
def primitive(S): return reduce(gcd, S) == 1
def is_covering(S):
    Sset = set(S)
    return all(any(v % q == 0 for v in Sset) for q in range(2, 15))
def banner(t): print("\n" + "=" * 78); print(t); print("=" * 78)

def main():
    small = list(range(1, 14))

    banner("(D) DECISIVE HUNT (bounded): any L=0 covering 13-set in coordinated families?")
    def scan(label, gen):
        checked=0; l0=0; worstM=None; worstM_S=None; minposL=None; tight=[]
        for S in gen:
            S = sorted(set(S))
            if len(S)!=13 or not is_covering(S) or not primitive(S): continue
            checked+=1
            Lv = L(S)
            if Lv==0:
                l0+=1
                Mv,_ = exact_M(S)
                if worstM is None or Mv<worstM: worstM=Mv; worstM_S=tuple(S)
                if Mv<=F(1,14): tight.append((Mv,tuple(S)))
            else:
                if minposL is None or Lv<minposL: minposL=Lv
        mp = f"{float(minposL):.6f}" if minposL is not None else "n/a"
        print(f"  {label}: checked={checked} L=0-survivors={l0} min-positive-L={mp}")
        if worstM is not None:
            print(f"     min M among L=0 survivors: {worstM}={float(worstM):.8f} at {worstM_S} "
                  f"(<1/14? {worstM<F(1,14)}; =1/14? {worstM==F(1,14)})")
        for Mv,Sx in tight[:6]:
            tag = "COUNTEREXAMPLE" if Mv<F(1,14) else "TIGHT(=1/14)"
            print(f"     {Sx} M={Mv}={float(Mv):.8f} [{tag}]")
        return worstM, tight

    # (D-i) 3 coordinated mult-of-14, drop 3 smalls, m up to 14
    def gen_d_i():
        for drops in it.combinations(range(1,14),3):
            basec=[v for v in small if v not in drops]
            for ms in it.combinations(range(1,15),3):
                yield basec+[14*x for x in ms]
    wD1,tD1 = scan("(D-i) 3-coord 14m (drop3, m<=14)", gen_d_i())

    # (D-ii) 4 coordinated mult-of-14, drop 4 smalls, m up to 9  (bounded)
    def gen_d_ii():
        for drops in it.combinations(range(1,14),4):
            basec=[v for v in small if v not in drops]
            for ms in it.combinations(range(1,10),4):
                yield basec+[14*x for x in ms]
    wD2,tD2 = scan("(D-ii) 4-coord 14m (drop4, k=4, m<=9)", gen_d_ii())

    # (D-iii) shared-modulus coordinated clusters
    def gen_d_iii():
        for drops in it.combinations(range(1,14),3):
            basec=[v for v in small if v not in drops]
            for modbase in (84,98,168,196,252):
                for start in range(1,6):
                    yield basec+[14*start, 14*start+modbase, 14*start+2*modbase]
    wD3,tD3 = scan("(D-iii) shared-modulus coord clusters", gen_d_iii())

    allM = [w for w in [wD1,wD2,wD3] if w is not None]
    alltight = tD1+tD2+tD3
    realctx = [x for x in alltight if x[0] < F(1,14)]
    print(f"\n  VERDICT (D): L=0 survivors total = {sum(1 for _ in [wD1,wD2,wD3] if _ is not None)} families with survivors;")
    print(f"     min M over all L=0 survivors = {min(allM) if allM else 'NONE (no L=0 survivors at all)'}")
    if allM: print(f"     >= 1/14? {min(allM) >= F(1,14)}")
    print(f"  TIGHT(=1/14): {len([x for x in alltight if x[0]==F(1,14)])}; "
          f"COUNTEREXAMPLES(<1/14): {len(realctx)}")
    if realctx: print(f"  *** {realctx}")

    # ===================================================================
    banner("(F) EXTREME PUSH: drive M(S) -> 1/14+ with coordinated large speeds (exact M)")
    # Search broadly for the covering primitive 13-set with the SMALLEST M, restricting to
    # sets with >=2 coordinated large mult-of-14 (the dangerous regime). Use exact M directly
    # on a curated but broad family (drop combos x small coordinated clusters).
    best = None  # (M, S)
    runner_up = []
    seen = set()
    # 2 coordinated large + parked structure: drop 2 smalls, add 2 mult-of-14 small multiples
    for drops in it.combinations(range(1,14),2):
        basec=[v for v in small if v not in drops]
        for ms in it.combinations(range(1,9),2):
            S = tuple(sorted(set(basec+[14*x for x in ms])))
            if len(S)!=13 or S in seen: continue
            seen.add(S)
            if not is_covering(S) or not primitive(S): continue
            Mv,_ = exact_M(list(S))
            if best is None or Mv<best[0]:
                best=(Mv,S)
            runner_up.append((Mv,S))
    # 3 coordinated large + drop 3
    for drops in it.combinations(range(1,14),3):
        basec=[v for v in small if v not in drops]
        for ms in it.combinations(range(1,8),3):
            S = tuple(sorted(set(basec+[14*x for x in ms])))
            if len(S)!=13 or S in seen: continue
            seen.add(S)
            if not is_covering(S) or not primitive(S): continue
            Mv,_ = exact_M(list(S))
            if best is None or Mv<best[0]:
                best=(Mv,S)
            runner_up.append((Mv,S))
    runner_up.sort(key=lambda x:x[0])
    print(f"  scanned {len(seen)} coordinated covering primitive 13-sets with exact M.")
    print(f"  GLOBAL MIN M = {best[0]}={float(best[0]):.8f} at {best[1]} (>1/14? {best[0]>F(1,14)})")
    print(f"  1/14 = {float(F(1,14)):.8f}.  gap above 1/14 = {float(best[0]-F(1,14)):.8f}")
    print(f"  10 smallest M among coordinated covering sets:")
    for Mv,S in runner_up[:10]:
        print(f"     M={str(Mv):>10s}={float(Mv):.8f}  S={list(S)}")

    # ===================================================================
    banner("(G) k=4,5 COORDINATED STRESS: does growing k drive M toward 1/14?")
    # For increasing k, take {1..(13-k)} u {k coordinated small mult-of-14} (the tightest
    # available covering structure) and report exact M and L.
    for k in range(1, 6):
        basec = list(range(1, 14-k))
        clus = [14*(i+1) for i in range(k)]  # 14,28,...,14k
        S = sorted(set(basec+clus))
        if len(S)!=13:
            # pad/trim
            S = sorted(set(basec+clus))
        cov = is_covering(S); prim = primitive(S)
        Mv,at = exact_M(S) if len(S)==13 else (None,None)
        Lv = L(S) if len(S)==13 else None
        print(f"  k={k}: S={S} |S|={len(S)} cov={cov} prim={prim} "
              f"M={Mv}={float(Mv) if Mv else 0:.8f} L={float(Lv) if Lv else 0:.8f}")
    # also: the MOST coordinated covering set we can build - many mult-of-14
    print("\n  Maximally-coordinated covering 13-sets (as many mult-of-14 as possible):")
    for base_small, clus in [
        ([1,2,3,4,5,6], [14,28,42,56,70,84,98]),   # 6 small + 7 coord = 13
        ([1,2,3,4,5],   [14,28,42,56,70,84,98,112]),  # 5+8=13
        ([1,2,3,4],     [14,28,42,56,70,84,98,112,126]),  # 4+9=13
    ]:
        S = sorted(set(base_small+clus))
        if len(S)!=13: continue
        cov=is_covering(S); prim=primitive(S)
        Mv,at = exact_M(S); Lv=L(S)
        print(f"  S={S} cov={cov} prim={prim} M={Mv}={float(Mv):.8f} L={float(Lv):.8f} "
              f"(>1/14? {Mv>F(1,14)})")

    # ===================================================================
    banner("(E) DECOUPLING-FLOOR AUDIT: floor>0 vs L>0 over coordinated 12-cores + parked w")
    def decoupling_floor(C, w):
        arcs = safe_arcs(C); r = len(arcs)
        return F(6,7)*meas(arcs) - F(r, 7*w)
    audit_cores = [
        ("{1..9}u{14,28,42}",  [1,2,3,4,5,6,7,8,9,14,28,42]),
        ("{1..6}u6coord",      [1,2,3,4,5,6,14,28,42,56,70,84]),
        ("{1..4}u8coord",      [1,2,3,4,14,28,42,56,70,84,98,112]),
    ]
    floor_silent_L_pos = 0; total=0
    for name, C in audit_cores:
        C = sorted(set(C))
        if len(C)!=12: continue
        mC = L(C); arcs=safe_arcs(C); r=len(arcs)
        wstar = F(r,6)/mC if mC>0 else None
        print(f"  {name}: meas={float(mC):.8f} r={r} floor-threshold w*={float(wstar):.1f}")
        tight=False
        for m in range(1,60):
            w=14*m; fl=decoupling_floor(C,w); Lv=L(list(C)+[w]); total+=1
            if fl<=0 and Lv>0: floor_silent_L_pos+=1
            if Lv==0:
                tight=True; Mv,_=exact_M(list(C)+[w])
                print(f"     w={w}: L=0! M={Mv}={float(Mv):.8f}")
        msmall=14
        while decoupling_floor(C,msmall)<=0 and msmall<14*300: msmall+=14
        print(f"     smallest w with proved floor>0: {msmall}; L there={float(L(list(C)+[msmall])):.8f}; "
              f"any tight w in m<=59: {tight}")
    print(f"\n  {floor_silent_L_pos}/{total} (w,core) pairs: proved floor <=0 (SILENT) but true L>0.")
    print("  => Exactly the gap: where the proved floor is silent, L is still positive, but the")
    print("     reduction lacks a proof of that. That silent zone IS GAP A / OPEN-Q-108.")

    banner("DONE")

if __name__ == "__main__":
    main()
