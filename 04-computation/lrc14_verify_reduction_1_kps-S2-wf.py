#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_verify_reduction_1_kps-S2-wf.py   (kind-pasteur-2026-06-17-S2 adversarial verify)

ADVERSARIAL VERIFICATION of the EASY-DOMINATES-HARD reduction for LRC(14),
focused ONLY on the COORDINATED REGIME: do k>=3 arithmetically-coordinated LARGE
speeds break the peel?

The reduction's claim chain (as stated):
  STEP1 (THM-523, proved): a counterexample must be COVERING (mult of every q in 2..14),
         in particular contains w with 14|w.
  STEP2 (proved): S = C u {w}, |C|=12; LRC(12) => meas(G_C) > 0;
         L(S) = meas(G_C) - meas(G_C cap D_w).
  STEP3 (proved decoupling floor, ONE large w): L(S) >= (6/7) meas(G_C) - r/(7w) > 0
         once w > r/(6 meas(G_C)).
  GAP A (OPEN-Q-108): with k>=3 coordinated large speeds, the one-at-a-time iteration is
         not obviously valid; need meas(G_C) bounded below uniformly in v_max.

WHAT THIS SCRIPT DOES (all EXACT, fractions.Fraction):
  A. Stress the ONE-AT-A-TIME peel chain: append coordinated large speeds w1<...<wk and
     watch L after each append.  Does each step really retain a 6/7 factor, or can a
     coordinated cluster collapse L far below the product of 6/7 factors?
  B. The CRUX: take a 12-CORE C that itself already contains k-1 coordinated large speeds.
     Is meas(G_C) still bounded below?  If a 12-core can have meas(G_C) -> 0 as we add more
     coordinated large speeds, the peel's STEP2 premise (meas(G_C) bounded below) FAILS and
     the reduction does NOT certify M(S) >= 1/14.
  C. Push M(S) toward 1/14 from ABOVE using coordinated large speeds.  Find the infimum of M
     over coordinated covering families and report the minimizer exactly.
  D. Directly hunt for an L(S)=0 covering 13-set among coordinated families (would be a tight
     set or counterexample).  Exact M on every L=0 survivor.
  E. Honest decoupling-floor audit: where (if anywhere) does the proved floor go NEGATIVE
     while L is still positive (floor too weak) vs where L itself goes to 0 (real break)?
"""
import sys
from math import gcd
from fractions import Fraction as F
from functools import reduce
import itertools as it

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

# ---------- exact M tool (validated; verbatim from task spec) ----------
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

# ---------- exact lonely measure G_S at gap 1/14 (arc unions) ----------
def danger_arcs(v):
    w = F(1, 14 * v); out = []
    for k in range(v + 1):
        c = F(k, v); lo, hi = c - w, c + w
        if lo < 0:
            out.append((F(0), hi)); out.append((1 + lo, F(1)))
        elif hi > 1:
            out.append((lo, F(1))); out.append((F(0), hi - 1))
        else:
            out.append((lo, hi))
    return [(a, b) for (a, b) in out if a < b]

def union_arcs(intervals):
    iv = sorted(intervals)
    if not iv: return []
    out = []; cl, ch = iv[0]
    for lo, hi in iv[1:]:
        if lo <= ch:
            if hi > ch: ch = hi
        else:
            out.append((cl, ch)); cl, ch = lo, hi
    out.append((cl, ch))
    return out

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

def banner(t):
    print("\n" + "=" * 78); print(t); print("=" * 78)

# cross-check: L(S)>0 <=> M(S)>1/14 (sanity on a sample)
def cross_check(S):
    Lv = L(S); Mv, at = exact_M(S)
    consistent = ((Lv > 0) == (Mv > F(1, 14)))
    return Lv, Mv, at, consistent


def main():
    banner("(0) SANITY: L(S)>0 <=> M(S)>1/14 on the 84m hard family + drop cores")
    base = [1,2,3,4,5,6,7,8,9,10,11,13]
    ok = True
    for m in [1,2,3]:
        S = sorted(base + [84*m])
        Lv, Mv, at, c = cross_check(S)
        ok &= c
        print(f"  84m m={m}: S has |S|={len(S)} cov={is_covering(S)} prim={primitive(S)} "
              f"L={float(Lv):.8f} M={Mv}={float(Mv):.8f} consistent={c}")
    print(f"  --> L<->M certificate consistent on sample: {ok}")

    # ===================================================================
    banner("(A) ONE-AT-A-TIME PEEL CHAIN: does each coordinated append keep a 6/7 factor?")
    # Reduction STEP3 says appending ONE large w multiplies L by >= ~6/7 (minus r/(7w)).
    # If k>=3 coordinated speeds each kept 6/7, final L >= (6/7)^k * meas(core) > 0.
    # ADVERSARIAL Q: can a coordinated cluster collapse the chain BELOW (6/7)^k?
    print("  Start from a 10-speed easy sub-core; append 3 coordinated large mult-of-14.")
    print("  Report L after each step and the step ratio vs 6/7=0.857143.")
    core0 = [1,2,3,4,5,6,7,8,9,10]
    L0 = L(core0)
    print(f"  core0={core0}  meas(G_core0)={L0}={float(L0):.8f}\n")
    worst_ratio = None; worst_chain = None
    # many coordinated triples: APs, geometric-ish, shared-modulus clusters
    triple_families = [
        ("AP small step 14,28,42",        (14,28,42)),
        ("AP large 98,112,126",           (98,112,126)),
        ("AP huge 1400,2800,4200",        (1400,2800,4200)),
        ("tight cluster 14,16x... no(use 14|w) 98,112,126", (98,112,126)),
        ("7-coprime spread 14*11,14*13,14*17", (14*11,14*13,14*17)),
        ("same residue mod 196: 14,210,406", (14,210,406)),  # all ≡14 mod 196
        ("multiples of 84: 84,168,252",   (84,168,252)),
        ("near-equal big 1414,1428,1442", (1414,1428,1442)),
        ("powers-ish 14,196,2744",        (14,196,2744)),
    ]
    for name, ws in triple_families:
        cur = list(core0); chain = [L(cur)]
        for w in ws:
            cur = cur + [w]; chain.append(L(cur))
        ratios = []
        for i in range(len(chain)-1):
            if chain[i] == 0:
                ratios.append(None)
            else:
                ratios.append(chain[i+1]/chain[i])
        finalS = sorted(set(cur)); fl = L(finalS)
        # track the single worst step ratio (most collapse) over all families
        for rr in ratios:
            if rr is not None and (worst_ratio is None or rr < worst_ratio):
                worst_ratio = rr; worst_chain = (name, ws, ratios)
        print(f"  {name}")
        print(f"     L-chain = {[float(x) for x in chain]}")
        print(f"     ratios  = {[float(r) if r is not None else None for r in ratios]}")
        print(f"     final |S|={len(finalS)} L={float(fl):.8f}  product-of-ratios floor check: "
              f"L>= (6/7)^3*L0 = {float(F(6,7)**3*L0):.8f}? {fl >= F(6,7)**3*L0}")
    print(f"\n  WORST single step-ratio observed: {float(worst_ratio):.6f} "
          f"(< 6/7={6/7:.6f}? {worst_ratio < F(6,7)})  in {worst_chain[0]}")
    print("  INTERPRETATION: a step ratio < 6/7 means that one coordinated append removed MORE")
    print("  than the decoupling-floor 1/7 share -> the one-at-a-time (6/7)^k product is NOT a")
    print("  valid uniform floor; the chain can dip below it.  (Whether L stays >0 is question B/D.)")

    # ===================================================================
    banner("(B) THE CRUX: can a 12-CORE C have meas(G_C) -> 0 with coordinated large speeds?")
    # STEP2's premise: meas(G_C) bounded below for ALL 12-cores.  If a 12-core built from
    # coordinated large speeds drives meas(G_C) toward 0, the reduction's premise FAILS.
    # Build 12-cores = (small base) u (coordinated large speeds), shrink the base, grow the cluster.
    print("  12-cores = small base u coordinated large speeds.  meas(G_C) as cluster grows:")
    print("  (these are 12-speed cores; appending a single parked w gives the 13-set.)")
    crux_cores = [
        ("{1..9} u {14,28,42}",            [1,2,3,4,5,6,7,8,9,14,28,42]),
        ("{1..6} u 6 coord {14k}",         [1,2,3,4,5,6,14,28,42,56,70,84]),
        ("{1..4} u 8 coord {14k}",         [1,2,3,4,14,28,42,56,70,84,98,112]),
        ("{1..3} u 9 coord {14k}",         [1,2,3,14,28,42,56,70,84,98,112,126]),
        ("{1,2} u 10 coord {14k}",         [1,2,14,28,42,56,70,84,98,112,126,140]),
        ("{1..9} u {98,112,126}(big AP)",  [1,2,3,4,5,6,7,8,9,98,112,126]),
        ("{1..6} u 6 coord big {14*20k}",  [1,2,3,4,5,6,280,560,840,1120,1400,1680]),
        ("{1..4} u 8 coord 7-coprime 14p", [1,2,3,4,14*11,14*13,14*17,14*19,14*23,14*29,14*31,14*37]),
    ]
    minmeas = None; minmeas_name = None
    for name, C in crux_cores:
        C = sorted(set(C))
        if len(C) != 12:
            print(f"  {name}: SKIP |C|={len(C)} != 12"); continue
        mC = L(C); arcs = safe_arcs(C); r = len(arcs)
        wid = max((b-a for a,b in arcs), default=F(0))
        prim = primitive(C)
        print(f"  {name:40s}: meas(G_C)={str(mC):>16s}={float(mC):.8f} r={r} widest={float(wid):.6f} prim={prim}")
        if minmeas is None or mC < minmeas:
            minmeas = mC; minmeas_name = name
    print(f"\n  MIN meas(G_C) over these coordinated 12-cores: {minmeas}={float(minmeas):.8f}  ({minmeas_name})")
    print(f"  compare extremal single-drop 7/858={float(F(7,858)):.8f}.")
    print(f"  Is the coordinated min BELOW 7/858? {minmeas < F(7,858)}")
    print("  KEY: if a primitive 12-core has meas(G_C) < 7/858, then 7/858 is NOT the universal")
    print("  floor and OPEN-Q-108's conjectured constant is WRONG (the reduction would need a")
    print("  smaller, possibly zero, uniform bound).")

    # systematic: among ALL primitive 12-cores of form {1..b} u {coordinated 14k}, find min meas
    banner("(B2) SYSTEMATIC MIN meas(G_C) over {1..b} u {AP of mult-of-14}, exact")
    gmin_meas = None; gmin_C = None
    for b in range(2, 11):
        smallbase = list(range(1, b+1))
        need = 12 - b
        if need < 1: continue
        # coordinated AP of mult-of-14: start a, step d
        for a in range(14, 14*8+1, 14):
            for d in range(14, 14*8+1, 14):
                clus = [a + i*d for i in range(need)]
                C = sorted(set(smallbase + clus))
                if len(C) != 12: continue
                if not primitive(C): continue
                mC = L(C)
                if gmin_meas is None or mC < gmin_meas:
                    gmin_meas = mC; gmin_C = C
    print(f"  GLOBAL min meas(G_C) over this scan: {gmin_meas}={float(gmin_meas):.8f}")
    print(f"     at C={gmin_C}")
    print(f"     >= 7/858={float(F(7,858)):.8f}? {gmin_meas >= F(7,858)}   ==0? {gmin_meas==0}")

    # ===================================================================
    banner("(C) PUSH M(S) TOWARD 1/14 FROM ABOVE with coordinated large speeds")
    # Among covering primitive 13-sets with k>=3 coordinated large speeds, minimize M(S).
    # The reduction is SAFE iff the infimum is > 1/14.  Report the minimizer EXACTLY.
    print("  Minimize M(S) over covering primitive 13-sets = (small drops) u (>=3 coord 14k).")
    best = None  # (M, S)
    small = list(range(1, 14))
    cnt = 0
    # 3 coordinated large speeds, drop 3 smalls
    for drops in it.combinations(range(1, 14), 3):
        basec = [v for v in small if v not in drops]
        if len(basec) != 10: continue
        for m1 in range(1, 9):
            for m2 in range(m1+1, 10):
                for m3 in range(m2+1, 11):
                    S = sorted(set(basec + [14*m1, 14*m2, 14*m3]))
                    if len(S) != 13: continue
                    if not is_covering(S) or not primitive(S): continue
                    # fast screen: only compute M if L is small
                    Lv = L(S)
                    cnt += 1
                    if Lv == 0:
                        Mv, _ = exact_M(S)
                        if best is None or Mv < best[0]:
                            best = (Mv, S)
                    else:
                        # L>0 => M>1/14; track the min-L (closest-to-tight) too
                        pass
    # Always also compute exact min-M directly on a curated dangerous subset (small mult-of-14)
    curated = []
    for drops in it.combinations(range(1, 14), 3):
        basec = [v for v in small if v not in drops]
        for trip in [(14,28,42),(14,28,56),(14,42,70),(28,42,56),(14,56,84),(42,84,126)]:
            S = sorted(set(basec + list(trip)))
            if len(S)==13 and is_covering(S) and primitive(S):
                curated.append(S)
    cur_best = None
    for S in curated:
        Mv,_ = exact_M(S)
        if cur_best is None or Mv < cur_best[0]:
            cur_best = (Mv, S)
    print(f"  scanned {cnt} covering+prim 3-coord sets via L-screen.")
    if best is not None:
        print(f"  L=0 survivors among 3-coord: min M={best[0]}={float(best[0]):.8f} at {best[1]} "
              f"(<1/14? {best[0]<F(1,14)})")
    else:
        print(f"  NO L=0 survivors among 3-coord coordinated sets (all have L>0 => all M>1/14).")
    print(f"  curated exact min M among hand-picked dangerous 3-coord sets: "
          f"{cur_best[0]}={float(cur_best[0]):.8f} at {cur_best[1]}  (>1/14? {cur_best[0]>F(1,14)})")

    # ===================================================================
    banner("(D) DECISIVE HUNT: any L(S)=0 covering 13-set in coordinated families?")
    # L=0 is NECESSARY for M<=1/14.  Sweep coordinated families; exact M on every L=0 survivor.
    def scan(label, gen, mmax=None):
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
        if tight:
            for Mv,Sx in tight[:6]:
                tag = "COUNTEREXAMPLE" if Mv<F(1,14) else "TIGHT(=1/14)"
                print(f"     {Sx} M={Mv}={float(Mv):.8f} [{tag}]")
        return worstM, tight

    # (D-i) three coordinated mult-of-14, dropping 3 smalls, wider m-range
    def gen_d_i():
        for drops in it.combinations(range(1,14),3):
            basec=[v for v in small if v not in drops]
            for m1 in range(1,13):
                for m2 in range(m1+1,14):
                    for m3 in range(m2+1,15):
                        yield basec+[14*m1,14*m2,14*m3]
    wD1,tD1 = scan("(D-i) 3-coord 14m (drop3)", gen_d_i())

    # (D-ii) FOUR coordinated mult-of-14, dropping 4 smalls (k=4 regime)
    def gen_d_ii():
        for drops in it.combinations(range(1,14),4):
            basec=[v for v in small if v not in drops]
            for ms in it.combinations(range(1,10),4):
                yield basec+[14*x for x in ms]
    wD2,tD2 = scan("(D-ii) 4-coord 14m (drop4, k=4)", gen_d_ii())

    # (D-iii) coordinated SHARED-MODULUS clusters: all ≡ same residue mod some M (84,98,...)
    def gen_d_iii():
        small = list(range(1,14))
        for drops in it.combinations(range(1,14),3):
            basec=[v for v in small if v not in drops]
            for modbase in (84,98,168,196,252):
                for start in range(1,6):
                    clus=[14*start, 14*start+modbase, 14*start+2*modbase]
                    yield basec+clus
    wD3,tD3 = scan("(D-iii) shared-modulus coord clusters", gen_d_iii())

    # (D-iv) MIXED coordination: large speeds that are NOT all mult of 14 but jointly cover.
    # e.g. one mult of 14 plus coordinated multiples of 84, 70, 90 (covering odd q via different speeds)
    def gen_d_iv():
        # base covers small q; add coordinated large speeds hitting q=11,13,14 etc.
        templates = [
            [1,2,3,4,5,6,7,8,9,10,11,13,14],     # AP+14 (the |S|=13? no, 13 elems already... check)
        ]
        # systematically: keep {1..9}, add 4 coordinated large covering 10..14
        for a in range(1,6):
            for d in range(1,8):
                clus = [14*a, 14*(a+d), 14*(a+2*d), 14*(a+3*d)]
                yield [1,2,3,4,5,6,7,8,9] + clus
    wD4,tD4 = scan("(D-iv) {1..9}+4 coord 14k", gen_d_iv())

    allM = [w for w in [wD1,wD2,wD3,wD4] if w is not None]
    alltight = tD1+tD2+tD3+tD4
    realctx = [x for x in alltight if x[0] < F(1,14)]
    print(f"\n  VERDICT (D): min M over ALL L=0 survivors in coordinated families = "
          f"{min(allM) if allM else 'n/a (no L=0 survivors)'}")
    if allM:
        print(f"     >= 1/14? {min(allM) >= F(1,14)}")
    print(f"  TIGHT (M=1/14): {len([x for x in alltight if x[0]==F(1,14)])}; "
          f"COUNTEREXAMPLES (M<1/14): {len(realctx)}")
    if realctx:
        print(f"  *** COUNTEREXAMPLES: {realctx}")
    else:
        print(f"  No counterexample, no tight covering set found in coordinated families "
              f"(consistent with LRC(14)).")

    # ===================================================================
    banner("(E) DECOUPLING-FLOOR AUDIT: floor>0 vs L>0 over coordinated 12-cores + parked w")
    # The PROVED floor is (6/7)meas(G_C) - r/(7w).  Does it certify the coordinated cases,
    # or does it go NEGATIVE (uninformative) while L stays positive?  That gap is exactly what
    # the reduction LEAVES OPEN.  We quantify: smallest w for which floor>0, vs whether L>0 anyway.
    def decoupling_floor(C, w):
        arcs = safe_arcs(C); r = len(arcs)
        return F(6,7)*meas(arcs) - F(r, 7*w)
    print("  For coordinated 12-cores, the floor threshold w* = r/(6 meas(G_C)) and whether L>0")
    print("  for ALL parked w (14|w) below w* (where the floor is silent):")
    audit_cores = [
        ("{1..9}u{14,28,42}",  [1,2,3,4,5,6,7,8,9,14,28,42]),
        ("{1..6}u6coord",      [1,2,3,4,5,6,14,28,42,56,70,84]),
        ("{1..4}u8coord",      [1,2,3,4,14,28,42,56,70,84,98,112]),
    ]
    floor_fails_but_L_pos = 0; total_checks = 0
    for name, C in audit_cores:
        C = sorted(set(C))
        if len(C)!=12: continue
        mC = L(C); arcs = safe_arcs(C); r=len(arcs)
        if mC==0:
            print(f"  {name}: meas(G_C)=0 -> reduction premise already fails for this core"); continue
        wstar = F(r,6)/mC
        print(f"  {name}: meas={float(mC):.8f} r={r} floor-threshold w*={float(wstar):.2f}")
        # check parked w below and above w*
        any_tight=False
        for m in range(1, 40):
            w = 14*m
            fl = decoupling_floor(C, w)
            Lv = L(list(C)+[w])
            total_checks += 1
            if fl <= 0 and Lv > 0:
                floor_fails_but_L_pos += 1
            if Lv == 0:
                any_tight=True
                Mv,_ = exact_M(list(C)+[w])
                print(f"     w={w}: L=0! M={Mv}={float(Mv):.8f} (tight/ctx)")
        # report the smallest w where floor>0 and confirm L>0 there
        msmall = 14
        while decoupling_floor(C, msmall) <= 0 and msmall < 14*200:
            msmall += 14
        print(f"     smallest w with floor>0: {msmall}; L there = {float(L(list(C)+[msmall])):.8f}; "
              f"any tight w found below: {any_tight}")
    print(f"\n  Across audit: {floor_fails_but_L_pos}/{total_checks} cases where the PROVED floor is "
          f"<=0 (silent) but L>0 anyway.")
    print("  => The floor is a LOSE sufficient condition: where it is silent, L is still positive,")
    print("     but the reduction has NO PROOF of that -- that silent zone IS OPEN-Q-108 / GAP A.")

    banner("DONE")

if __name__ == "__main__":
    main()
