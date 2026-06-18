#!/usr/bin/env python3
"""
LRC(14) -- A1: the "perfect middle of section 0" centering geometry.
kps-S2-wf

For covering sets S = C u {w}, w == 0 mod 14, the runner w sits in SECTION 0
on the grid tau=a/14 FOREVER (gridM=0). The true lonely time is OFF the grid.

This script does, all with EXACT rationals (fractions.Fraction):

 (1) THE OPTIMUM GEOMETRY. For S=C u {w}, compute exact M(S) and argmax tau*.
     Then evaluate frac(w*tau*) (the SIGNED position of w inside its 0-band) and
     ||w*tau*||. Question: at the optimum, is w BINDING (||w tau*|| == M, hugging
     the EDGE of section 0) or CENTERED (||w tau*|| ~ 1/2, perfect middle)?
     Sweep the family w=84m, m=1..8, plus w in {14,28,...,168} with the AP core,
     plus alternative cores.

 (2) THE CONSTRUCTIVE CENTERING DEVICE. Define, for threshold h,
        T_w(h) = { tau in [0,1) : ||w*tau|| >= h }   (the SAFE band of w)
     and G_C = { tau : ||v*tau|| >= 1/14 for all v in C }   (easy core gap set).
     M(S) >= 1/14  <==  T_w(1/14) INTERSECTS G_C (then any tau in the intersection
     is a survivor). Compute meas(T_w(1/14) cap G_C) EXACTLY and the SLACK
     = max over tau in G_C of ||w tau|| - 1/14 (how far we can center w while
     staying lonely for the core). Does centering ever fail to clear 1/14?

 (3) THE CANONICAL CENTERING TAU. Is there tau = (2k+1)/(2w) (exact center of a
     w-band, ||w tau||=1/2) that is FORCED into G_C by the core structure?
     Relate to dihedral clock events (14k +/- 1)/(14 v).

All measures are EXACT Fractions via interval union on the circle.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce

# ----------------------------------------------------------------------------
# EXACT M TOOL (verbatim from task spec; THM-524 complete candidate set)
# ----------------------------------------------------------------------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r

def gmin(S, t):
    return min(nrm(v * t) for v in S)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2 * k + 1, 2 * v) <= F(1, 2):
            C.add(F(2 * k + 1, 2 * v)); k += 1
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            for d in (S[i] + S[j], S[j] - S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2):
                        C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C

def M(S):
    b = F(0); at = None
    for t in cand(S):
        v = gmin(S, t)
        if v > b:
            b = v; at = t
    return b, at

def fracpart(x):
    r = x - int(x)
    return r + 1 if r < 0 else r

# ----------------------------------------------------------------------------
# Exact circle-interval union measure.
# Arcs given as (lo, hi) with 0<=lo<hi<=1 (already split across wraparound).
# ----------------------------------------------------------------------------
def union_measure(arcs):
    arcs = sorted([(a, b) for (a, b) in arcs if b > a])
    tot = F(0); cl = ch = None
    for a, b in arcs:
        if ch is None:
            cl, ch = a, b
        elif a <= ch:
            if b > ch: ch = b
        else:
            tot += ch - cl; cl, ch = a, b
    if ch is not None:
        tot += ch - cl
    return tot

def danger_arcs(v, h=F(1, 14)):
    """{ tau : ||v tau|| < h } as a list of (lo,hi) sub-arcs of [0,1)."""
    w = h / v
    A = []
    for k in range(v + 1):
        c = F(k, v); lo = c - w; hi = c + w
        if lo < 0:
            A += [(F(0), hi), (1 + lo, F(1))]
        elif hi > 1:
            A += [(lo, F(1)), (F(0), hi - 1)]
        else:
            A.append((lo, hi))
    return A

def safe_arcs(v, h=F(1, 14)):
    """{ tau : ||v tau|| >= h } = complement of danger_arcs within [0,1)."""
    # complement of the union of danger arcs
    d = sorted([(a, b) for (a, b) in danger_arcs(v, h) if b > a])
    # merge dangers
    merged = []
    cl = ch = None
    for a, b in d:
        if ch is None:
            cl, ch = a, b
        elif a <= ch:
            if b > ch: ch = b
        else:
            merged.append((cl, ch)); cl, ch = a, b
    if ch is not None:
        merged.append((cl, ch))
    # complement
    comp = []
    prev = F(0)
    for a, b in merged:
        if a > prev:
            comp.append((prev, a))
        prev = max(prev, b)
    if prev < 1:
        comp.append((prev, F(1)))
    return comp

def intersect_arc_lists(A, B):
    """Intersection of two arc lists (each a union of sub-arcs of [0,1))."""
    A = sorted(A); B = sorted(B); out = []
    i = j = 0
    while i < len(A) and j < len(B):
        lo = max(A[i][0], B[j][0]); hi = min(A[i][1], B[j][1])
        if lo < hi:
            out.append((lo, hi))
        if A[i][1] < B[j][1]:
            i += 1
        else:
            j += 1
    return out

def G_set_arcs(C, h=F(1, 14)):
    """G_C = intersection over v in C of safe_arcs(v,h). Returns arc list."""
    res = [(F(0), F(1))]
    for v in C:
        res = intersect_arc_lists(res, safe_arcs(v, h))
        if not res:
            return []
    return res

def lcm(a, b):
    return a // gcd(a, b) * b

def lcm_list(xs):
    return reduce(lcm, xs, 1)

# ----------------------------------------------------------------------------
# CORES
# ----------------------------------------------------------------------------
AP_FULL = list(range(1, 14))                       # {1..13} (the tight AP, n=13)
# The 12-speed EASY CORE used by codex (AP-drop-6): {1..13} \ {6}
CORE_AP_DROP6 = [1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13]
# Hardest known cover from task: {1..11,13} (i.e. AP-drop-12) is the 12-core,
# union {84} is the m=1 covering 13-set.
CORE_1to11_13 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13]

def header(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)

# ============================================================================
# PART 1: THE OPTIMUM GEOMETRY -- where does parked w sit at tau*?
# ============================================================================
def part1():
    header("PART 1: Parked-runner position at the OPTIMUM tau*  (binding vs centered)")
    print("S = C u {w}, w == 0 mod 14.  Report M(S), tau*, frac(w tau*), ||w tau*||,")
    print("and whether w is BINDING (||w tau*|| == M) i.e. hugging the EDGE of section 0.")
    print()

    def run(C, w, label):
        S = sorted(set(C + [w]))
        m, t = M(S)
        wf = fracpart(w * t)         # signed position in [0,1)
        wn = nrm(w * t)              # ||w tau*||
        binding = (wn == m)
        # which runners are binding (at distance exactly m)?
        binders = [v for v in S if nrm(v * t) == m]
        # position of w relative to NEAREST grid point of its OWN 0-band:
        # w==0 mod14 => the grid a/14 maps w to integer; safe band center is
        # tau halfway between grid pts. "edge" means ||w t|| small, "middle" ~1/2.
        print(f"  {label:34s} w={w:4d}  M={str(m):9s} ~{float(m):.5f}  tau*={str(t):9s}")
        print(f"      ||w tau*|| = {str(wn):8s} ~{float(wn):.5f}   "
              f"frac(w tau*)={str(wf):9s} ~{float(wf):.5f}")
        print(f"      w BINDING? {binding}   binders={binders}   (M-clear vs 1/14: "
              f"{'>=' if m>=F(1,14) else '<'} ; slack={float(m-F(1,14)):+.6f})")
        return m, t, wn, binding

    print("--- Family w = 84m on core {1..11,13} (the canonical hard cover) ---")
    fam = []
    for mm in range(1, 9):
        w = 84 * mm
        res = run(CORE_1to11_13, w, f"core1..11,13 + 84*{mm}")
        fam.append((mm, w) + res)
    # closed form check: M = 7m/(84m+5)
    print("\n  Closed-form check  M ?= 7m/(84m+5):")
    for mm, w, m, t, wn, binding in fam:
        cf = F(7 * mm, 84 * mm + 5)
        print(f"    m={mm}: M={str(m):9s}  7m/(84m+5)={str(cf):9s}  match={m==cf}")

    print("\n--- w in {14,28,42,...,168} on core {1..11,13} ---")
    for k in range(1, 13):
        w = 14 * k
        # only meaningful if w not already in core; all 14k are not in core
        run(CORE_1to11_13, w, f"core1..11,13 + 14*{k}")

    print("\n--- Same w=84m family on the AP-drop-6 core {1..13}\\{6} ---")
    for mm in range(1, 6):
        run(CORE_AP_DROP6, 84 * mm, f"core(AP\\6) + 84*{mm}")

    print("\n--- Non-84 multiples of 14 (e.g. 14*lcm-ish) on core {1..11,13} ---")
    for w in [70, 98, 126, 154, 182, 210, 280, 420]:
        run(CORE_1to11_13, w, f"core1..11,13 + {w}")

# ============================================================================
# PART 2: CONSTRUCTIVE CENTERING -- meas(T_w cap G_C) and the SLACK
# ============================================================================
def part2():
    header("PART 2: Constructive centering.  meas(T_w(1/14) cap G_C) and centering SLACK")
    print("T_w(h) = {||w tau|| >= h} (safe band of w).  G_C = gap-1/14 lonely set of core C.")
    print("If T_w(1/14) cap G_C is nonempty, ANY tau there is a survivor => M(S)>=1/14.")
    print("CENTERED max = max_{tau in G_C} ||w tau||  (how 'middle' we can push w while")
    print("staying lonely for the core). SLACK = CENTERED max - 1/14.")
    print()

    def run(C, w, label, show_inter=False):
        Gc = G_set_arcs(C, F(1, 14))
        measGc = union_measure(Gc)
        Tw = safe_arcs(w, F(1, 14))
        inter = intersect_arc_lists(Gc, Tw)
        measInter = union_measure(inter)
        # CENTERED max: maximize ||w tau|| over tau in G_C.
        # ||w tau|| is piecewise linear; its max over a finite union of intervals
        # is attained either at an interval endpoint or at a band-center
        # tau=(2j+1)/(2w) lying inside G_C. Enumerate candidates exactly.
        cands = set()
        for a, b in Gc:
            cands.add(a); cands.add(b)
            # band centers of w in [a,b]
            # center_j = (2j+1)/(2w);  a <= (2j+1)/(2w) <= b
            jlo = (2 * a * w - 1)  # >= : (2j+1) >= 2aw  => 2j >= 2aw-1
            # iterate j over plausible range
            import math
            j0 = int(F(2 * a * w - 1, 2))
            j1 = int(F(2 * b * w - 1, 2)) + 2
            for j in range(max(0, j0 - 2), j1 + 2):
                c = F(2 * j + 1, 2 * w)
                if a <= c <= b:
                    cands.add(c)
        best = F(-1); bestt = None
        for c in cands:
            val = nrm(w * c)
            if val > best:
                best = val; bestt = c
        # at endpoints of Gc, the lonely constraint is >= so closed; safe to use.
        slack = best - F(1, 14)
        ok = measInter > 0 or best >= F(1, 14)
        print(f"  {label:30s} w={w:4d}")
        print(f"      meas(G_C)        = {str(measGc):14s} ~{float(measGc):.8f}")
        print(f"      meas(Tw cap G_C) = {str(measInter):14s} ~{float(measInter):.8f}"
              f"   ({'NONEMPTY' if measInter>0 else 'measure-0'})")
        print(f"      CENTERED max ||w tau||_{{G_C}} = {str(best):10s} ~{float(best):.6f}"
              f"  at tau={str(bestt):10s}")
        print(f"      SLACK over 1/14  = {str(slack):14s} ~{float(slack):+.6f}"
              f"   clears 1/14? {best>=F(1,14)}")
        if show_inter:
            print(f"      Tw cap G_C arcs: {[(str(a),str(b)) for a,b in inter][:6]}"
                  f"{' ...' if len(inter)>6 else ''}")
        return measGc, measInter, best, slack

    print("--- core {1..11,13}, w = 84m ---")
    for mm in range(1, 6):
        run(CORE_1to11_13, 84 * mm, f"core1..11,13 +84*{mm}", show_inter=(mm == 1))

    print("\n--- core {1..11,13}, w = 14k ---")
    for k in [1, 2, 3, 5, 6, 7, 12]:
        run(CORE_1to11_13, 14 * k, f"core1..11,13 +14*{k}")

    print("\n--- AP-drop-6 core {1..13}\\{6}, w = 84m (codex's meas(G_C)=7/858 case) ---")
    measGc6 = union_measure(G_set_arcs(CORE_AP_DROP6, F(1, 14)))
    print(f"  meas(G_C) for AP-drop-6 core = {measGc6} ~ {float(measGc6):.8f}"
          f"   (codex sharpened bound 7/858 = {float(F(7,858)):.8f})")
    for mm in [1, 2, 3]:
        run(CORE_AP_DROP6, 84 * mm, f"core(AP\\6) +84*{mm}")

# ============================================================================
# PART 3: THE CANONICAL CENTERING TAU -- (2k+1)/(2w) forced into G_C?
# ============================================================================
def part3():
    header("PART 3: Canonical centering tau=(2j+1)/(2w) -- is it FORCED into G_C?")
    print("Perfect middle of w-band: tau=(2j+1)/(2w) gives ||w tau||=1/2 (max possible).")
    print("Question: does the easy core's structure FORCE some such center into G_C?")
    print("Also relate to dihedral clock events (14k +/- 1)/(14 v).")
    print()

    def run(C, w, label):
        Gc = G_set_arcs(C, F(1, 14))
        # all band-centers of w in (0,1):
        centers = [F(2 * j + 1, 2 * w) for j in range(w)]
        good = []
        for c in centers:
            if c >= 1:
                continue
            inG = any(a <= c <= b for a, b in Gc)
            if inG:
                # verify it's a survivor with gap == core gap (>=1/14)
                gC = min(nrm(v * c) for v in C)
                good.append((c, gC))
        print(f"  {label:30s} w={w:4d}: #band-centers in (0,1)={w}, "
              f"#inside G_C={len(good)}")
        if good:
            # report the center that maximizes the CORE gap (most robust)
            good.sort(key=lambda x: x[1], reverse=True)
            c, gC = good[0]
            S = sorted(set(C + [w]))
            mS = min(nrm(v * c) for v in S)
            print(f"      best center tau={str(c):10s} ~{float(c):.6f}  "
                  f"core-gap={str(gC):8s} ~{float(gC):.5f}  full-gap(min over S)="
                  f"{str(mS):8s} ~{float(mS):.5f}")
            # dihedral form: is c == (14k +/- 1)/(14 v) for small v,k?
            forms = []
            for v in range(1, 14):
                for sgn in (1, -1):
                    # c = (14k + sgn)/(14 v)  => 14 v c = 14k + sgn => k = (14 v c - sgn)/14
                    num = 14 * v * c - sgn
                    if num.denominator == 1 and int(num) % 14 == 0:
                        kk = int(num) // 14
                        forms.append(f"(14*{kk}{'+' if sgn>0 else '-'}1)/(14*{v})")
            print(f"      dihedral-clock forms (14k+/-1)/(14 v): "
                  f"{forms[:4] if forms else 'none small'}")
        else:
            print(f"      NO band-center of w lies in G_C (centering tau not core-forced).")
        return good

    print("--- core {1..11,13} ---")
    for w in [84, 168, 14, 28, 42]:
        run(CORE_1to11_13, w, f"core1..11,13 + w={w}")
    print("\n--- AP-drop-6 core ---")
    for w in [84, 168, 14]:
        run(CORE_AP_DROP6, w, f"core(AP\\6) + w={w}")

    # SPECIAL: the m=1 optimum tau*=37/89. Where is THAT relative to w=84 centers?
    header("PART 3b: The m=1 optimum tau*=37/89 vs w=84 band geometry")
    w = 84
    t = F(37, 89)
    print(f"  tau* = 37/89 ~ {float(t):.6f}")
    print(f"  ||84 tau*|| = {nrm(84*t)} ~ {float(nrm(84*t)):.6f}   "
          f"frac(84 tau*) = {fracpart(84*t)} ~ {float(fracpart(84*t)):.6f}")
    print(f"  nearest 84-band CENTER below/above tau*:")
    j = int(F(2 * t * 84 - 1, 2))
    for jj in [j, j + 1]:
        c = F(2 * jj + 1, 2 * 84)
        print(f"    center j={jj}: tau={c} ~{float(c):.6f}, distance to tau* = "
              f"{abs(c - t)} ~{float(abs(c-t)):.6f}")
    print(f"  => at the optimum, w=84 is at frac {float(fracpart(84*t)):.4f}, i.e. ||.||="
          f"{float(nrm(84*t)):.4f}=M, hugging the EDGE of section 0 (NOT the 1/2 middle).")

if __name__ == "__main__":
    part1()
    part2()
    part3()
