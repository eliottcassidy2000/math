#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""F/G/E sections (FAST): use cheap exact L(S) to screen, exact M only on smallest-L sets.
L(S)>0 <=> M(S)>1/14 is exact, so min-L sets are the closest-to-tight; exact M confirms."""
import sys
from math import gcd
from fractions import Fraction as F
from functools import reduce
import itertools as it

OUT = open(sys.argv[1], "w", encoding="utf-8") if len(sys.argv) > 1 else sys.stdout
def p(*a):
    print(*a, file=OUT); OUT.flush()

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
def banner(t): p("\n" + "=" * 78); p(t); p("=" * 78)

small = list(range(1, 14))

banner("(F) EXTREME PUSH: smallest L (=> closest-to-tight) over coordinated covering 13-sets")
p("  Strategy: L(S) is cheap and L(S)>0 <=> M(S)>1/14 (exact). Find the coordinated covering")
p("  primitive 13-sets with SMALLEST L; those are the closest to violating 1/14. Then run")
p("  exact M on the 15 smallest-L sets to confirm their true gap.")
cand_sets = []
seen = set()
# 2 coordinated large mult-of-14, drop 2 smalls
for drops in it.combinations(range(1,14),2):
    basec=[v for v in small if v not in drops]
    for ms in it.combinations(range(1,9),2):
        S = tuple(sorted(set(basec+[14*x for x in ms])))
        if len(S)!=13 or S in seen: continue
        seen.add(S)
        if not is_covering(S) or not primitive(S): continue
        cand_sets.append((L(list(S)), S))
# 3 coordinated large mult-of-14, drop 3 smalls
for drops in it.combinations(range(1,14),3):
    basec=[v for v in small if v not in drops]
    for ms in it.combinations(range(1,8),3):
        S = tuple(sorted(set(basec+[14*x for x in ms])))
        if len(S)!=13 or S in seen: continue
        seen.add(S)
        if not is_covering(S) or not primitive(S): continue
        cand_sets.append((L(list(S)), S))
cand_sets.sort(key=lambda x:x[0])
p(f"  scanned {len(seen)} coordinated covering primitive 13-sets (cheap L screen).")
nzero = sum(1 for Lv,_ in cand_sets if Lv==0)
p(f"  L=0 (=> M<=1/14 tight/counterexample) count: {nzero}")
p(f"  smallest positive L = {cand_sets[0][0]}={float(cand_sets[0][0]):.8f} at {list(cand_sets[0][1])}")
p(f"  15 smallest-L coordinated covering 13-sets, with EXACT M:")
worstM = None; worstM_S = None
for Lv, S in cand_sets[:15]:
    Mv,_ = exact_M(list(S))
    if worstM is None or Mv < worstM: worstM=Mv; worstM_S=S
    p(f"     L={float(Lv):.8f}  M={str(Mv):>10s}={float(Mv):.8f}  S={list(S)}")
p(f"  --> SMALLEST M among the 15 closest-to-tight coordinated sets: {worstM}={float(worstM):.8f}")
p(f"      1/14={float(F(1,14)):.8f}; gap above = {float(worstM-F(1,14)):.8f}; >1/14? {worstM>F(1,14)}")

banner("(G) k=1..5 COORDINATED STRESS: {1..(13-k)} u {14,28,...,14k}")
for k in range(1, 6):
    basec = list(range(1, 14-k)); clus = [14*(i+1) for i in range(k)]
    S = sorted(set(basec+clus))
    cov = is_covering(S); prim = primitive(S)
    if len(S)==13:
        Lv = L(S); Mv,at = exact_M(S)
        p(f"  k={k}: S={S} cov={cov} prim={prim} M={Mv}={float(Mv):.8f} L={float(Lv):.8f} (>1/14? {Mv>F(1,14)})")
    else:
        p(f"  k={k}: S={S} |S|={len(S)} != 13 (skipped)")
p("\n  Maximally-coordinated covering 13-sets (max # mult-of-14):")
for base_small, clus in [
    ([1,2,3,4,5,6], [14,28,42,56,70,84,98]),
    ([1,2,3,4,5],   [14,28,42,56,70,84,98,112]),
    ([1,2,3,4],     [14,28,42,56,70,84,98,112,126]),
]:
    S = sorted(set(base_small+clus))
    if len(S)!=13: continue
    cov=is_covering(S); prim=primitive(S); Lv=L(S); Mv,_=exact_M(S)
    p(f"  S={S} cov={cov} prim={prim} M={Mv}={float(Mv):.8f} L={float(Lv):.8f} (>1/14? {Mv>F(1,14)})")

banner("(E) DECOUPLING-FLOOR AUDIT: proved floor>0 vs true L>0 (silent-zone = GAP A)")
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
    p(f"  {name}: meas={float(mC):.8f} r={r} floor-threshold w*=r/(6 meas)={float(wstar):.1f}")
    tight=False
    for m in range(1,60):
        w=14*m; fl=decoupling_floor(C,w); Lv=L(list(C)+[w]); total+=1
        if fl<=0 and Lv>0: floor_silent_L_pos+=1
        if Lv==0:
            tight=True; Mv,_=exact_M(list(C)+[w])
            p(f"     w={w}: L=0! M={Mv}={float(Mv):.8f}")
    msmall=14
    while decoupling_floor(C,msmall)<=0 and msmall<14*400: msmall+=14
    p(f"     smallest w with proved floor>0: {msmall}; L there={float(L(list(C)+[msmall])):.8f}; "
      f"any tight w (m<=59): {tight}")
p(f"\n  {floor_silent_L_pos}/{total} (w,core) pairs: proved floor <=0 (SILENT) but true L>0.")
p("  => Where the proved floor is silent, the TRUE L is positive, but the reduction has no")
p("     proof of it. That silent zone is precisely GAP A / OPEN-Q-108.")
banner("DONE")
if OUT is not sys.stdout: OUT.close()
