#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
FAST focused counterexample hunt + resonance/peel probes for LRC(14) OPEN-Q-108.
Companion to lrc_q108_verify_boundarywind_kps-Sx-wf.py (Probe 1/1b already passed:
MU_8=1087/5880 confirmed exactly, Fatou-decorrelation gap confirmed POSITIVE up to
+0.059 -- the decorrelated value is a FLOOR not a majorant, per MISTAKE-082).

Here: PROBE 2 (hunt, span<=120), PROBE 3 (peel V-bound on WIDE bases), PROBE 4
(resonant directions), PROBE 5 (exhaustive small-span window k=8).
EXACT arithmetic.
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd, comb

try:
    sys.stdout.reconfigure(encoding='utf-8', line_buffering=True)
except Exception:
    pass

CAPS = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91), 11: F(66,91), 12: F(6,7)}

def bp(E):
    s = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for j in range(7):
            m = 0
            while True:
                xv = (F(j,7)+m)/e
                if xv >= 1: break
                if xv >= 0: s.add(xv)
                m += 1
    return sorted(b for b in s if 0 <= b < 1)

def p0p1(E):
    E = sorted(set(E)); B = bp(E); a = F(0); b = F(0)
    for lo, hi in zip(B, B[1:]+[F(1)]):
        if hi <= lo: continue
        mid = (lo+hi)/2
        miss = set(range(1,7)) - set(int((e*mid)%1*7) for e in E)
        if len(miss) == 0: a += hi-lo
        elif len(miss) == 1: b += hi-lo
    return a, b

def p0(E): return p0p1(E)[0]

def miss_profile(B):
    Bb = bp(B); prof = {t: F(0) for t in range(7)}
    for lo, hi in zip(Bb, Bb[1:]+[F(1)]):
        if hi <= lo: continue
        mid = (lo+hi)/2
        miss = set(range(1,7)) - set(int((e*mid)%1*7) for e in B)
        prof[len(miss)] += hi-lo
    return prof

def c_t(t, r): return sum((-1)**i*comb(t,i)*(F(7-i,7))**r for i in range(t+1))
def Pr_fatou(B, r):
    prof = miss_profile(B); return sum(prof[t]*c_t(t,r) for t in range(7))

def normalize(E):
    E = sorted(set(E)); E = [e-E[0] for e in E]
    g = 0
    for e in E: g = gcd(g,e)
    if g > 1: E = [e//g for e in E]
    return tuple(E)

def gcd_tuple(E):
    g = 0
    for e in E: g = gcd(g,e)
    return g

def is_wide(E):
    E = sorted(set(E)); return E[0]==0 and E[-2]-E[0] > 14

def trueV(B):
    Bb = bp(B); pts = Bb+[F(1)]; ivs = []
    for lo, hi in zip(pts, pts[1:]):
        if hi <= lo: continue
        mid = (lo+hi)/2
        miss = frozenset(set(range(1,7)) - set(int((e*mid)%1*7) for e in B))
        ivs.append(miss)
    Vtot = 0
    for j in range(1,7):
        prev = False
        for miss in ivs:
            cur = (miss == frozenset({j}))
            if cur and not prev: Vtot += 1
            prev = cur
    return Vtot

# ================= PROBE 2 (fast) =================
def probe2():
    print("="*78); print("PROBE 2 (fast).  COUNTEREXAMPLE HUNT span<=120: wide E with p0>cap_k."); print("="*78)
    worst = {k:(F(0),None) for k in range(8,13)}; viol = []; seen = set(); n = 0
    def consider(E, fam):
        nonlocal n
        E = sorted(set(E))
        if len(E) < 8 or len(E) > 12: return
        if E[0] != 0: E = [e-E[0] for e in E]
        if not is_wide(E): return
        key = normalize(E)
        if key[-2] <= 14 or key[-1] > 120: return
        kk = len(key); tag = (kk,key)
        if tag in seen: return
        seen.add(tag); n += 1
        v = p0(list(key)); cap = CAPS[kk]
        if v > worst[kk][0]: worst[kk] = (v,key)
        if v > cap: viol.append((kk,key,v,cap,fam))
    # (a) consec(k-1)+1far up to 120
    for k in range(8,13):
        base = list(range(k-1))
        for g in range(15,121): consider(base+[g], "c+1far")
    # (b) consec(k-2)+2far boundary 15-60
    for k in range(8,13):
        base = list(range(k-2))
        for g1 in range(15,61):
            for g2 in range(g1+1,81): consider(base+[g1,g2], "c+2far")
    # (c) resonant multiples
    for k in range(8,13):
        base = list(range(k-2))
        for g in range(8,55):
            for mlt in (2,3,4,5,7): consider(base+[g,mlt*g], "res-mult")
    # (d) AP cores + d=7
    for k in range(8,13):
        for d in range(1,9): consider([d*i for i in range(k)], f"AP{d}")
        # near-pinned stretched consec (the MISTAKE-082 family)
        for kbase in range(2,k):
            base = list(range(kbase)); extra = k-kbase
            for step in range(15,41):
                consider(base+[base[-1]+step*(j+1) for j in range(extra)], "stretch")
    # (e) two-block all separations 15-70, all size splits
    for k in range(8,13):
        for s1 in range(1,k):
            s2 = k-s1; A = list(range(s1))
            for sep in range(15,71):
                consider(A+[sep+i for i in range(s2)], "2blk")
    # (e2) shifted single block {0}+{w..w+k-2} -- the 0.202 MISTAKE-080 family
    for k in range(8,13):
        for w in range(15,80):
            consider([0]+[w+i for i in range(k-1)], "0+block")
    # (f) random multi-cluster
    random.seed(424242)
    for k in range(8,13):
        for _ in range(3000):
            ncl = random.randint(2, min(4,k-1))
            cuts = sorted(random.sample(range(1,k), ncl-1)) if ncl>1 else []
            sizes = [b-a for a,b in zip([0]+cuts, cuts+[k])]
            if any(s<1 for s in sizes): continue
            E = []; base = 0; ok = True
            for idx,s in enumerate(sizes):
                start = 0 if idx==0 else base+random.randint(15,40)
                if s==1: clu=[start]
                else:
                    diam = random.randint(s-1, min(13,4*s))
                    pts = sorted(random.sample(range(1,diam+1), s-1)); clu=[start]+[start+p for p in pts]
                E += clu; base = max(E)+1
            consider(E, "rand")
    print(f"  distinct wide primitive sets checked: {n}")
    print(f"  {'k':>3} {'max p0':>10} {'cap_k':>10} {'margin':>10}   argmax")
    for k in range(8,13):
        v,E = worst[k]; cap = CAPS[k]
        print(f"  {k:>3} {float(v):>10.5f} {float(cap):>10.5f} {float(cap-v):>10.5f}   {E}")
    if viol:
        print(f"\n  *** {len(viol)} VIOLATIONS ***")
        for (k,E,v,cap,fam) in viol[:20]:
            print(f"    k={k} p0={float(v):.6f} > cap={float(cap):.6f} fam={fam} E={E}")
    else:
        print("\n  NO violation: every wide set tested has p0 <= cap_k.")
    print()
    return viol, worst

# ================= PROBE 3 peel on wide bases =================
def probe3():
    print("="*78); print("PROBE 3.  Peel (6/49)V/w bound on WIDE / multi-cluster bases."); print("="*78)
    print("  THM-546 proved for residual against a FIXED base; the iterated peel uses")
    print("  bases that are themselves wide. Does |Dw-Phi1| <= (6/49)V(B)/w survive?")
    print(f"  {'B':>28} {'|B|':>3} {'V':>4} {'w':>4} {'|Dw-Phi1|':>11} {'(6/49)V/w':>11} {'ratio':>7}")
    worst = F(0); wr = None
    bases = [list(range(7)), [0,1,2,20,21,22], [0,1,2,3,30,31,32],
             [0,5,10,15,20,25], [0,7,14,21,28], [0,1,2,50,51,100,101]]
    for B in bases:
        p0B,p1B = p0p1(B); Phi = p1B/7; Vb = trueV(B)
        for w in [15,17,23,31,53,101]:
            if w <= max(B): continue
            dw = p0(B+[w]) - p0B; rw = abs(dw-Phi)
            bnd = F(6,49)*Vb/w; ratio = rw/bnd if bnd>0 else F(0)
            if ratio > worst: worst = ratio; wr = (B,w)
            flag = "  VIOL" if ratio>1 else ""
            print(f"  {str(B):>28} {len(B):>3} {Vb:>4} {w:>4} {float(rw):>11.6f} {float(bnd):>11.6f} {float(ratio):>7.3f}{flag}")
    print(f"  worst ratio = {float(worst):.4f}  (>1 => bound BROKEN); worst at {wr}")
    print()
    return worst

# ================= PROBE 4 resonance =================
def probe4():
    print("="*78); print("PROBE 4.  Resonant far scales: p0 - P_r(B) (excess over Fatou FLOOR)."); print("="*78)
    print(f"  {'B':>14} {'far':>16} {'p0':>9} {'P_r':>9} {'p0-P_r':>9}  note")
    B = list(range(6)); maxexc = F(0)
    for far,note in [([7,14],"7,14"),([7,49],"7,49"),([15,30],"g,2g"),([21,42],"21,42"),
                     ([15,45],"g,3g"),([35,70],"5*7,10*7"),([16,24],"sh8"),([15,16],"adj"),
                     ([100,200],"far g,2g")]:
        E = B+far; pv = p0(E); pf = Pr_fatou(B,len(far)); exc = pv-pf
        if exc > maxexc: maxexc = exc
        print(f"  {str(B):>14} {str(far):>16} {float(pv):>9.5f} {float(pf):>9.5f} {float(exc):>+9.5f}  {note}")
    print(f"  MAX excess p0-P_r(B) = {float(maxexc):+.6f}  (this is what R must absorb; FLOOR not majorant)")
    print()
    return maxexc

# ================= PROBE 5 exhaustive small window =================
def probe5():
    print("="*78); print("PROBE 5.  TRULY exhaustive window k=8, span<=24 (test family-completeness)."); print("="*78)
    k = 8; Smax = 24; worst = F(0); wE = None; nviol = 0; ntot = 0
    for S in range(15, Smax+1):
        for mid in itertools.combinations(range(1,S), k-2):
            E = (0,)+mid+(S,)
            if gcd_tuple(E) != 1: continue
            if E[-2] <= 14: continue
            ntot += 1
            v = p0(list(E))
            if v > worst: worst = v; wE = E
            if v > CAPS[k]:
                nviol += 1
                if nviol <= 5: print(f"    VIOL E={E} p0={float(v):.6f}")
    print(f"  k=8 EXHAUSTIVE span<=24: {ntot:,} wide sets, max p0={float(worst):.6f} at {wE},")
    print(f"     cap={float(CAPS[8]):.5f}, margin={float(CAPS[8]-worst):.5f}, violations={nviol}")
    print()
    return nviol, worst

def main():
    print("#"*78); print("# LRC OPEN-Q-108 FAST hunt + peel/resonance/window probes"); print("#"*78); print()
    viol,_ = probe2()
    probe3()
    probe4()
    nv5,_ = probe5()
    print("="*78); print("VERDICT (this script)"); print("="*78)
    print(f"  Probe2 wide-set violations: {len(viol)}")
    print(f"  Probe5 exhaustive k=8 span<=24 violations: {nv5}")
    print("  => objective inequality p0<=cap held everywhere tested." if not viol and nv5==0
          else "  => COUNTEREXAMPLE FOUND.")

if __name__ == "__main__":
    main()
