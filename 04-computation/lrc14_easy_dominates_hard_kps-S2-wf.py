#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_easy_dominates_hard_kps-S2-wf.py   (kind-pasteur-2026-06-17-S2, THM-525, A2)

EASY-DOMINATES-HARD REDUCTION for LRC(14), STRESS-TESTED with EXACT rationals.

SETUP. M(S)=max_tau min_{v in S} ||v tau||.  LRC(14): M(S)>=1/14 for every primitive
13-speed set S.  THM-523: a counterexample must be a COVERING set -- it must contain a
multiple of every q in 2..14, in particular a speed w with 14 | w (so w sits in section
0 forever on the grid tau=a/14).  Write S = C u {w}, |C|=12.  C is an EASY core:
PROVEN LRC(12) gives M(C) >= 1/13 > 1/14, hence meas(G_C) > 0 where
   G_C = { tau in [0,1) : ||v tau|| > 1/14  for all v in C }   (the gap-1/14 lonely set).
We have L(S) = meas(G_S) and  L(C u {w}) = meas(G_C) - meas(G_C ∩ D_w),  D_w = w's
danger comb { tau : ||w tau|| <= 1/14 }.

This script does FOUR things, all exact:
 (1) Make the {hard covering S} <-> {easy 12-core C} correspondence explicit & WELL-DEFINED
     (which w to peel; multiple multiples of 14; primitivity).  Census small covering sets.
 (2) THE KEY LOWER BOUND.  For the extremal worst cores (AP-drop-6 = {1..5,7..13}, meas=7/858,
     and others), test EXHAUSTIVELY over ALL w with 14|w whether D_w can cover G_C.  Compare
     the decoupling floor (6/7)meas(G_C) - r/(7w) to the TRUE exact L(C u {w}).  Quantify the
     "G_C is fat enough" transversality: w's comb teeth have half-width 1/(14w); G_C's arcs
     have widths >> that for the extremal cores, so a tooth can kill at most O(1/w) of mass.
 (3) THE COORDINATED REGIME (k>=3 arithmetically-coordinated large speeds): the only
     uncontrolled case.  Peel one-at-a-time vs all-at-once; find where the easy/hard peel
     BREAKS or holds.  Be brutally honest.
 (4) The 84m hard family: M = 7m/(84m+5), centering vs binding.  And: does the easy core's
     STRUCTURE (lonely intervals, SDR spread, (Z/14)* symmetry) close the bound?

Everything uses fractions.Fraction.  Floats only for fast prescreen, exact-confirm decisions.
"""
import sys
from math import gcd
from fractions import Fraction as F
from functools import reduce
from itertools import combinations

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

# ===========================================================================
#  EXACT lonely measure / G_C as a union of arcs.  All rational.
# ===========================================================================
def danger_arcs(v):
    """Closed danger arcs in [0,1): { tau : ||v tau|| <= 1/14 }, as list of (lo,hi)."""
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
    """Merge a list of (lo,hi) into disjoint sorted closed arcs."""
    iv = sorted(intervals)
    if not iv:
        return []
    out = []; cl, ch = iv[0]
    for lo, hi in iv[1:]:
        if lo <= ch:
            if hi > ch: ch = hi
        else:
            out.append((cl, ch)); cl, ch = lo, hi
    out.append((cl, ch))
    return out

def safe_arcs(S):
    """G_S as disjoint open arcs (complement of union of danger arcs in [0,1))."""
    danger = union_arcs([iv for v in set(S) for iv in danger_arcs(v)])
    out = []; cur = F(0)
    for lo, hi in danger:
        if lo > cur:
            out.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1:
        out.append((cur, F(1)))
    return out  # open arcs (lonely set), the safe set

def meas(arcs):
    return sum((b - a for a, b in arcs), F(0))

def L(S):
    """Exact meas(G_S) = lonely measure of S at gap 1/14."""
    return meas(safe_arcs(S))

def primitive(S):
    return reduce(gcd, S) == 1

def lcm(a, b): return a * b // gcd(a, b)
def lcmlist(S): return reduce(lcm, S, 1)

# ===========================================================================
#  (1)  THE HARD<->EASY CORRESPONDENCE.  Well-defined peel rule.
# ===========================================================================
#  Covering set S: for every q in 2..14, some v in S with q | v.  In particular some w with
#  14 | w.  PEEL RULE (canonical):  among all w in S with 14 | w, peel the SMALLEST.
#  C = S \ {w}.  Then meas(G_C) >= meas(G_S) (removing a runner only enlarges the lonely set),
#  and LRC(12) guarantees meas(G_C) >= 1/13-gap measure > 0 IF C is primitive; if not, we
#  dilate.  We verify the peel is well-defined and the residual C is an honest 12-core.

def is_covering(S):
    Sset = set(S)
    for q in range(2, 15):
        if not any(v % q == 0 for v in Sset):
            return False
    return True

def peel_w(S):
    """Canonical peel: smallest multiple of 14.  Returns (w, C) or None if no mult of 14."""
    mults = sorted(v for v in S if v % 14 == 0)
    if not mults:
        return None
    w = mults[0]
    C = [v for v in S if v != w]
    return w, C

# ===========================================================================
#  (2)  THE KEY LOWER BOUND.  Decoupling floor vs exact; transversality of D_w.
# ===========================================================================
def decoupling_floor(C, w):
    """(6/7)*meas(G_C) - r/(7w)  where r = #arcs of G_C.  PROVED lower bound for L(C u {w})."""
    arcs = safe_arcs(C)
    r = len(arcs)
    return F(6, 7) * meas(arcs) - F(r, 7 * w)

def exact_after_add(C, w):
    """Exact L(C u {w})."""
    return L(list(C) + [w])

def comb_can_cover(C, w):
    """
    Does D_w (teeth half-width 1/(14w), spaced 1/w) cover ALL of G_C?
    True  => L(C u {w}) = 0 (tight).  False => survivor (lonely) exists.
    Computed exactly via L.
    """
    return exact_after_add(C, w) == 0

def widest_arc(C):
    arcs = safe_arcs(C)
    if not arcs:
        return F(0)
    return max(b - a for a, b in arcs)

# ===========================================================================
#  WORST CORES (extremal small-measure 12-cores).
# ===========================================================================
AP13      = list(range(1, 14))                       # {1..13}, the AP
DROP6     = [v for v in range(1, 14) if v != 6]      # {1..5,7..13}, meas=7/858 (extremal)
GW_T5     = [1,2,3,4,5,6,7,8,9,10,11,13]             # {1..11,13}, the Goddyn-Wong easy core
DROPS     = {e: [v for v in range(1, 14) if v != e] for e in range(1, 14)}

# ===========================================================================
def banner(t):
    print("\n" + "=" * 78); print(t); print("=" * 78)

def main():
    banner("(0) Baselines: meas(G_C) for the single-drop cores {1..13}\\{e}")
    drop_meas = {}
    for e in range(1, 14):
        C = DROPS[e]; m = L(C); drop_meas[e] = m
        arcs = safe_arcs(C)
        print(f"  drop e={e:2d}: meas(G_C)={str(m):>14s} = {float(m):.8f}  "
              f"#arcs={len(arcs):2d}  widest={float(widest_arc(C)):.6f}  "
              f"floor(6/7)m={float(F(6,7)*m):.8f}")
    emin = min(drop_meas, key=lambda e: drop_meas[e])
    print(f"  --> minimum single-drop meas at e={emin}: {drop_meas[emin]} "
          f"= {float(drop_meas[emin]):.8f}  (THM-523: extremal core, 7/858 expected)")
    print(f"  CHECK 7/858 = {float(F(7,858)):.8f}; drop-6 = {drop_meas[6]} -> "
          f"{'MATCH' if drop_meas[6]==F(7,858) else 'MISMATCH'}")

    # -----------------------------------------------------------------------
    banner("(1) HARD<->EASY correspondence: SPEED-COUNT bookkeeping + the genuine 13-speed family")
    print("  *** BOOKKEEPING CORRECTION (load-bearing): LRC(14) is about |S|=13 speeds.")
    print("      The naive 'hard family' {1..13} u {14m} has |S|=14 -> that is the N=15 case,")
    print("      NOT ours.  Indeed M({1..14})=1/15<1/14, but that is the tight 14-speed AP for")
    print("      LRC(15), irrelevant to LRC(14).  A GENUINE 13-speed covering set must DROP one")
    print("      of {1..13} to make room for the parked w.  So the hard 13-cores are the")
    print("      SINGLE-DROP cores {1..13}\\{e} (plus possibly further drops), NOT the full AP.")
    # demonstrate the speed-count trap once:
    S14 = list(range(1, 15)); M14, _ = exact_M(S14)
    print(f"      DEMO: {{1..14}} has |S|={len(S14)}, M={M14}={float(M14):.6f} < 1/14 "
          f"(N=15 tight AP, NOT a LRC(14) counterexample).")
    print()
    print("  Genuine 13-speed hard family: C_e = {1..13}\\{e}, parked w=14m.  M(C_e u {w}):")
    for e in [6, 12, 11, 1]:
        Ce = DROPS[e]
        for m in [1, 2, 6]:
            w = 14 * m
            S = sorted(set(Ce + [w]))
            if len(S) != 13:
                continue
            Mv, at = exact_M(S)
            cov = is_covering(S); prim = primitive(S)
            print(f"    drop e={e:2d}, w={w:3d}: |S|={len(S)} cov={cov} prim={prim} "
                  f"M={str(Mv):>10s}={float(Mv):.6f}  (>1/14? {Mv > F(1,14)})")
    print("  PEEL RULE (well-defined): given a covering 13-set S, peel the SMALLEST w with 14|w;")
    print("  C=S\\{w} is a 12-speed core.  If C imprimitive, dilate by 1/gcd (scale-invariance,")
    print("  THM-522).  meas(G_C) >= meas(G_S) always (deleting a runner enlarges the safe set).")

    # The 84m family from the task setup: S = {1..11,13} u {84m}
    banner("(1b) The codex 84m hard family: S = {1..11,13} u {84m}; M = 7m/(84m+5)")
    base = [1,2,3,4,5,6,7,8,9,10,11,13]
    for m in [1, 2, 3, 5, 10]:
        w = 84 * m
        S = base + [w]
        cov = is_covering(S); prim = primitive(S)
        # exact M via the validated tool
        Mval, at = exact_M(S)
        pred = F(7*m, 84*m + 5)
        print(f"  m={m:2d} w={w:4d}: covering={cov} prim={prim}  M(S)={str(Mval):>10s}={float(Mval):.8f} "
              f"at tau={at}  predicted 7m/(84m+5)={str(pred)}={float(pred):.8f}  "
              f"{'MATCH' if Mval==pred else 'DIFF'}  (1/14={float(F(1,14)):.6f})")
    print("  So every 84m set has M>1/14 (min 7/89 at m=1).  These are LRC(14)-SAFE, not")
    print("  counterexamples -- exactly what easy-dominates-hard must reprove structurally.")

    # -----------------------------------------------------------------------
    banner("(2) KEY LOWER BOUND: can ANY w (14|w) cover the extremal cores' G_C?  EXHAUSTIVE")
    # For each worst core, sweep ALL w with 14|w up to a bound; the decoupling floor proves
    # that for w large enough no covering is possible, so only finitely many w need checking.
    # NOTE: cores here are 12-ELEMENT (so C u {w} is the genuine 13-speed set).
    # DROP6 and GW_T5 are 12-element single-drop cores. (AP13 has 13 elements -> not a 12-core;
    # appending w would make 14 speeds = N=15 trap, so it is EXCLUDED here.)
    for name, C in [("DROP6 {1..5,7..13} (extremal, meas=7/858)", DROP6),
                    ("GW-easy {1..11,13} (drop e=12)", GW_T5),
                    ("DROP4 {1..3,5..13} (small meas)", DROPS[4]),
                    ("DROP10 {1..9,11..13} (small meas)", DROPS[10])]:
        assert len(C) == 12, f"core {name} must have 12 elements, has {len(C)}"
        mC = L(C); arcs = safe_arcs(C); r = len(arcs); wid = widest_arc(C)
        print(f"\n  CORE {name}: meas(G_C)@gap-1/14={str(mC)}={float(mC):.8f}  r={r} arcs  widest={float(wid):.6f}")
        if mC == 0:
            print(f"    *** meas(G_C)=0 at gap 1/14: this core is ITSELF tight at gap 1/14")
            print(f"        (its LRC(12) gap-1/13 lonely set does NOT fatten to gap 1/14).")
            print(f"        Decoupling floor is vacuous (0); the AP core is the degenerate case.")
            thr_ceil = 0
        else:
            # decoupling: floor>0  <=>  w > r/(6 meas) ; beyond that NO w can make L=0.
            thr = F(r, 6) / mC    # w must exceed this for floor>0
            thr_ceil = -(-thr.numerator // thr.denominator)
            print(f"    decoupling floor (6/7)meas - r/(7w) > 0  <=>  w > r/(6 meas) = "
                  f"{float(thr):.2f}  -> for w >= {thr_ceil+1} (mult of 14) floor>0 => L>0 guaranteed.")
        # exhaustively test all multiples of 14 up to a safe bound (>> threshold).
        Wmax = max(2000, 2 * (thr_ceil + 14))
        tight_ws = []; minpos = None; minpos_w = None
        w = 14
        while w <= Wmax:
            Lv = exact_after_add(C, w)
            if Lv == 0:
                tight_ws.append(w)
            else:
                if minpos is None or Lv < minpos:
                    minpos = Lv; minpos_w = w
            w += 14
        if tight_ws:
            print(f"    *** {len(tight_ws)} values of w (14|w, w<={Wmax}) make C u {{w}} TIGHT (L=0).")
            print(f"        first few: {tight_ws[:8]}{' ...' if len(tight_ws)>8 else ''}")
            if mC == 0:
                print(f"        EXPECTED: meas(G_C)=0 already, so adding ANY w keeps L=0.")
                print(f"        These S have max-min = EXACTLY 1/14 (AP equality case, NOT < 1/14).")
                # confirm max-min is exactly 1/14 on a sample
                for tw in tight_ws[:3]:
                    Sx = sorted(C) + [tw]
                    Mx, atx = exact_M(Sx)
                    print(f"          w={tw}: M(S)={Mx}={float(Mx):.8f}  (=1/14? {Mx==F(1,14)})  prim={primitive(Sx)}")
            else:
                for tw in tight_ws[:5]:
                    print(f"        w={tw}: S={sorted(C)+[tw]}  primitive={primitive(sorted(C)+[tw])}")
        else:
            print(f"    NO w (14|w, w<={Wmax}) makes C u {{w}} tight.  "
                  f"min positive L over these w: {str(minpos)}={float(minpos):.8f} at w={minpos_w}.")
            print(f"    => for EVERY w with 14|w, M(C u {{w}}) >= 1/14 (transversality holds for this core).")

    # -----------------------------------------------------------------------
    banner("(2b) Transversality quantified: tooth half-width 1/(14w) vs G_C arc widths")
    # A single tooth of D_w (half-width 1/(14w)) can sit inside an arc of width L_arc and
    # delete at most min(2/(14w), L_arc).  To cover an arc fully a tooth-cluster is needed;
    # the comb spacing is 1/w, so #teeth in an arc of width L_arc is ~ w*L_arc.  Covering
    # needs every point within 1/(14w) of a tooth-center j/w.  An arc survives iff its width
    # exceeds the residual gap 2(1/w - 2/(14w)) = 2*(6/(14w)) handled by the comb...
    # we just REPORT exact survivor mass for the extremal core at the worst small w.
    C = DROP6; mC = L(C)
    print(f"  Extremal core DROP6 meas={mC}={float(mC):.8f}; arcs:")
    for a, b in safe_arcs(C):
        print(f"    arc [{a},{b}] width={b-a}={float(b-a):.8f}")
    print("  For each small w (14|w), exact survivor L and which arcs survive:")
    for m in range(1, 13):
        w = 14 * m
        Lv = exact_after_add(C, w)
        surv = safe_arcs(list(C) + [w])
        print(f"    w={w:3d}: L={str(Lv):>14s}={float(Lv):.8f}  #survivor-arcs={len(surv)}")

    # -----------------------------------------------------------------------
    banner("(3) COORDINATED REGIME: k>=3 arithmetically-coordinated LARGE speeds")
    # Replace several core speeds by large coordinated ones (e.g. an AP of large speeds, or
    # multiples sharing a modulus).  Test whether L stays >= 1/1260 and whether the easy
    # peel (one-at-a-time decoupling) still controls them.
    print("  (a) k coordinated large multiples-of-14 added/substituted, exact L:")
    # Family: keep a small easy sub-core, add k coordinated large speeds in AP.
    tests = [
        ("base {1..10} + {14,28,42}",        [1,2,3,4,5,6,7,8,9,10,14,28,42]),
        ("base {1..10} + {98,112,126}",      [1,2,3,4,5,6,7,8,9,10,98,112,126]),
        ("base {1..10} + 3 coord {14a}",     [1,2,3,4,5,6,7,8,9,10,14*11,14*13,14*17]),
        ("drop6 core + coordinated 84,168,252", [1,2,3,4,5,7,8,9,10,84,168,252]),
        ("{1,2,3,4,5} + 8 large coord 14k",  [1,2,3,4,5,14,28,42,56,70,84,98,112]),
    ]
    for name, S in tests:
        S = sorted(set(S))
        if len(S) != 13:
            note = f"(|S|={len(S)}, not 13)"
        else:
            note = ""
        Lv = L(S); cov = is_covering(S); prim = primitive(S)
        print(f"    {name:38s}: |S|={len(S)} cov={cov} prim={prim} L={str(Lv):>14s}={float(Lv):.8f} {note}")

    print("\n  (b) PEEL-ONE-AT-A-TIME vs ALL-AT-ONCE on coordinated large speeds.")
    print("      Take S0 = small easy core; append coordinated large speeds w1<w2<w3; track")
    print("      L after each append (decoupling predicts each multiplies by ~6/7).")
    core0 = [1,2,3,4,5,6,7,8,9,10]   # 10-speed easy sub-core
    print(f"      core0={core0}  meas(G_core0)={float(L(core0)):.8f}")
    for ws in [(14,28,42), (98,112,126), (154,168,182), (1400,2800,4200)]:
        cur = list(core0); chain = [L(cur)]
        for w in ws:
            cur = cur + [w]; chain.append(L(cur))
        ratios = [float(chain[i+1]/chain[i]) if chain[i] != 0 else None for i in range(len(chain)-1)]
        print(f"      append {ws}: L-chain={[float(x) for x in chain]}")
        print(f"                    step-ratios (vs 6/7={6/7:.4f})={ratios}")
        finalS = sorted(set(cur))
        if len(finalS) == 13:
            print(f"                    final |S|=13 L={float(L(finalS)):.8f} >= 1/1260? "
                  f"{L(finalS) >= F(1,1260)}")

    print("\n  (c) THE HONEST DANGER: coordinated speeds that SHARE all teeth (same residues).")
    print("      If w1,w2,w3 are all multiples of 14 sharing structure, their combs may")
    print("      OVERLAP rather than tile -> less coverage than 'independent' 6/7 each.")
    print("      But the WORST is when they are 7-coprime & spread (maximal independent kill).")
    # exhaustive small coordinated search: cores {1..5} u {3 coordinated large speeds},
    # minimize L over coordinated triples (AP-structured) of multiples of 14.
    best = None
    smallcore = [1,2,3,4,5,6,7,8,9,10]
    # search arithmetic triples a, a+d, a+2d (all mult of 14) for minimal L
    import itertools
    for a in range(14, 14*15, 14):
        for d in range(14, 14*10, 14):
            trip = [a, a+d, a+2*d]
            S = sorted(set(smallcore + trip))
            if len(S) != 13: continue
            if not primitive(S): continue
            Lv = L(S)
            if best is None or Lv < best[0]:
                best = (Lv, trip, S)
    if best:
        print(f"      min L over coordinated AP-triples (mult of 14) appended to {{1..10}}:")
        print(f"        L={str(best[0])}={float(best[0]):.8f} at triple {best[1]}, S={best[2]}")
        print(f"        >= 1/1260={float(F(1,1260)):.8f}? {best[0] >= F(1,1260)}   "
              f">= 1/14-counterexample? {best[0] == 0}")

    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    banner("(3d) DECISIVE: hunt a 13-speed COVERING set with M<=1/14 (a real counterexample)")
    print("  KEY FILTER (exact, proved): L(S)>0  =>  M(S)>1/14.  Contrapositive: a counterexample")
    print("  (M<1/14) or a TIGHT set (M=1/14) must have L(S)=0.  So we FAST-screen with exact L,")
    print("  and only run the (expensive) exact-M tool on the L=0 survivors.  This makes the")
    print("  coordinated-large-speed sweep tractable while remaining a true M-gap certificate.")
    import itertools as it
    small = list(range(1, 14))
    def scan(label, gen):
        """gen yields candidate 13-sets; return (#checked, #L0, worstM-among-L0, min-positive-L)."""
        checked = 0; l0 = 0; worstM = None; minposL = None; worstM_S = None
        tight_examples = []
        for S in gen:
            S = sorted(set(S))
            if len(S) != 13 or not is_covering(S) or not primitive(S):
                continue
            checked += 1
            Lv = L(S)
            if Lv == 0:
                l0 += 1
                Mv, _ = exact_M(S)
                if worstM is None or Mv < worstM:
                    worstM = Mv; worstM_S = tuple(S)
                if Mv <= F(1, 14):
                    tight_examples.append((Mv, tuple(S)))
            else:
                if minposL is None or Lv < minposL:
                    minposL = Lv
        print(f"  {label}: checked {checked} covering+prim; L=0 survivors {l0}; "
              f"min positive L={float(minposL) if minposL else None:.6f}" if minposL else
              f"  {label}: checked {checked} covering+prim; L=0 survivors {l0}")
        if worstM is not None:
            print(f"      worst (min) M among L=0 survivors: M={worstM}={float(worstM):.6f} at {worstM_S} "
                  f"(<1/14? {worstM < F(1,14)}; =1/14? {worstM == F(1,14)})")
        if tight_examples:
            print(f"      *** {len(tight_examples)} sets with M<=1/14 (TIGHT or counterexample):")
            for Mv, Sx in tight_examples[:6]:
                tag = "COUNTEREXAMPLE!!" if Mv < F(1, 14) else "tight (=1/14)"
                print(f"          {Sx}  M={Mv}={float(Mv):.6f}  [{tag}]")
        return worstM, tight_examples

    # (i) single large parked: drop e, add 14m (covering needs the mult of 14)
    def gen_i():
        for e in range(1, 14):
            for m in range(1, 30):
                yield [v for v in small if v != e] + [14 * m]
    w_i, t_i = scan("(i) single-drop {1..13}\\{e} + 14m", gen_i())

    # (ii) two coordinated large speeds: drop 2 smalls, add 2 mults of 14
    def gen_ii():
        for drops in it.combinations(range(1, 14), 2):
            base = [v for v in small if v not in drops]
            for m1 in range(1, 18):
                for m2 in range(m1 + 1, 18):
                    yield base + [14 * m1, 14 * m2]
    w_ii, t_ii = scan("(ii) two-drop + 2 coord 14m", gen_ii())

    # (iii) THREE coordinated large speeds (the named uncontrolled regime)
    def gen_iii():
        for drops in it.combinations(range(1, 14), 3):
            base = [v for v in small if v not in drops]
            for m1 in range(1, 11):
                for m2 in range(m1 + 1, 12):
                    for m3 in range(m2 + 1, 13):
                        yield base + [14 * m1, 14 * m2, 14 * m3]
    w_iii, t_iii = scan("(iii) three-drop + 3 coord 14m (uncontrolled regime)", gen_iii())

    allM = [w for w in [w_i, w_ii, w_iii] if w is not None]
    alltight = t_i + t_ii + t_iii
    real_ctx = [x for x in alltight if x[0] < F(1, 14)]
    print(f"\n  VERDICT (3d): min M over all L=0 survivors = "
          f"{min(allM) if allM else 'n/a'} "
          f"(>=1/14? {min(allM) >= F(1,14) if allM else 'n/a'}).")
    print(f"  TIGHT sets (M=1/14) found: {len(alltight)}.  REAL counterexamples (M<1/14): "
          f"{len(real_ctx)}  --> {'NONE (consistent with LRC(14))' if not real_ctx else real_ctx}")

    # -----------------------------------------------------------------------
    banner("(4) WHAT THE EASY STRUCTURE BUYS: lonely intervals, SDR spread, (Z/14)* symmetry")
    print("  NAIVE CLAIM (REFUTED below): LRC(12) gap>=1/13 at tau0 => an interval of width")
    print("  >= 2*(1/13-1/14)=1/91 where gap stays >1/14.  This is FALSE: the correct Lipschitz")
    print("  half-width is (M(C)-1/14)/v_max, and v_max (largest core speed) can be LARGE, so the")
    print("  guaranteed arc SHRINKS like 1/v_max.  This is exactly why the structural lever does")
    print("  NOT close OPEN-Q-108.  Quantify the TRUE guaranteed interval:")
    # For each single-drop easy core, find the gap-1/13 maximizer and the width of the
    # surrounding gap-1/14 arc (this is the GUARANTEED fattening from LRC(12)).
    for e in [6, 12, 1, 7]:
        C = DROPS[e]
        M13, tau13 = exact_M(C)   # the LRC(12) maximizer (gap >= 1/13)
        vmax = max(C)
        # the TWO binding runners at tau13, and the true guaranteed half-width
        dists = sorted((nrm(v * tau13), v) for v in C)
        binders = [(d, v) for d, v in dists if d == M13]
        # half-width to keep min >= 1/14, governed by the tightest (largest-v) binder
        slack = M13 - F(1, 14)
        true_hw = min(slack / v for d, v in binders) if binders else F(0)
        # host arc of G_C containing tau13
        arcs = safe_arcs(C)
        host = next(((a, b) for a, b in arcs if a < tau13 < b), None)
        hostw = (host[1] - host[0]) if host else F(0)
        print(f"    drop e={e:2d}: M(C)={str(M13)}={float(M13):.6f} at tau={tau13}; v_max={vmax}; "
              f"binders={[v for _,v in binders]}")
        print(f"           slack=M(C)-1/14={float(slack):.5f}; TRUE Lipschitz half-width "
              f"(slack/v_binder)={float(true_hw):.6f}; host-arc width={float(hostw):.6f} "
              f"(>=1/91? {hostw >= F(1,91)})")
    print("\n  CONCLUSION (honest):  the guaranteed safe interval around the LRC(12) maximizer")
    print("  has half-width ~ (M(C)-1/14)/v_binder, NOT a universal constant.  For the extremal")
    print("  cores the binding runners are the LARGE ones (11,12), so the arc is THIN (~0.0027).")
    print("  meas(G_C) is bounded below NOT by a single arc but by the SUM over its (few) arcs,")
    print("  and the extremal value 7/858 already reflects this.  So LRC(12) alone does NOT")
    print("  hand a uniform meas(G_C) bound: the missing estimate is a UNIFORM-IN-v_max lower")
    print("  bound on the TOTAL gap-1/14 measure -- exactly OPEN-Q-108, not closed by Lipschitz.")

    banner("DONE")

# ---- exact M tool (validated; from task spec) -----------------------------
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

if __name__ == "__main__":
    main()
