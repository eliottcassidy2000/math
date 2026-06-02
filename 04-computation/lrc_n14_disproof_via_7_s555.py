#!/usr/bin/env python3
"""
DISPROOF HUNT for LRC@14 steered by the n=7 / even-fold structure.
opus-2026-06-01-S555 (remote-control).

Goal: find a 13-speed set with max-collar M = max_t min_i ||v_i t|| < 1/14
(a genuine LRC(14) counterexample).  HONEST PRIOR: LRC is believed true and no
counterexample is known at ANY n; a real one would be historic.  We therefore
VERIFY every candidate with exact arithmetic and report truthfully.

Why steer by 7: the even-fold gives M14(S) <= M(fold(S)) with fold = halved
evens.  For e=|fold| <= 6, LRC(7) forces M(fold) >= 1/7 and the preimage
construction produced witnesses 127/127 (S554) -> a counterexample must live in
the UNPROTECTED e >= 7 regime, where M(fold) is an unproven LRC(8+) quantity.
The sharpest probe: make the even part a DOUBLED tight AP so M(fold) is pushed
to 1/(e+1) (just above 1/14 for e up to 12) and stack the odd coupling.

Targets:
  A. AP with an arbitrary subset of its 6 even runners DOUBLED (generalises V*).
  B. doubled-AP_k even part {2,4,...,2k} plus (13-k) odd speeds (e=k>=7).
  C. mod-7 / 7-resonant random configs.
  D. measure-minimising hill-climb forced to e>=7, 7-seeded, speeds <= 120.
All decisions exact (Fraction); float only screens the hill-climb.
"""

from fractions import Fraction
from math import gcd
from itertools import combinations
import random

THR = Fraction(1, 14)


def M_exact(V):
    c = set()
    for v in V:
        for k in range(2 * v):
            c.add(Fraction(2 * k + 1, 2 * v) % 1)
    for i in range(len(V)):
        for j in range(len(V)):
            for s in (1, -1):
                d = V[i] + s * V[j]
                if d:
                    for k in range(abs(d) + 1):
                        c.add(Fraction(k, d) % 1)
    best = Fraction(0)
    bt = None
    for t in c:
        mn = min(min((v * t) % 1, 1 - (v * t) % 1) for v in V)
        if mn > best:
            best = mn
            bt = t
    return best, bt


def primitive(S):
    g = 0
    for v in S:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in S)) if g else tuple(sorted(S))


def float_measure(V):
    pts = []
    for v in V:
        for k in range(v + 1):
            b = k / v
            pts.append((b + (1 / 14) / v) % 1.0)
            pts.append((b - (1 / 14) / v) % 1.0)
    pts.sort()
    L = len(pts)
    meas = 0.0
    for i in range(L):
        a = pts[i]
        b = pts[i + 1] if i + 1 < L else pts[0] + 1.0
        mid = (a + (b - a) / 2.0) % 1.0
        ok = True
        for v in V:
            x = (v * mid) % 1.0
            if (x if x < 1 - x else 1 - x) < 1 / 14 - 1e-12:
                ok = False
                break
        if ok:
            meas += (b - a)
    return meas


def report(tag, V):
    Vp = primitive(V)
    M, t = M_exact(Vp)
    flag = "  <<< COUNTEREXAMPLE M<1/14 !!!" if M < THR else (
        "  (tight)" if M == THR else "")
    return (Vp, M, t, M < THR, f"{tag}: {Vp}  M={M}={float(M):.5f}  (1/14={float(THR):.5f}){flag}")


def partA():
    print("== A. AP with a subset of its even runners DOUBLED (generalises V*) ==")
    AP = list(range(1, 14))
    evens = [v for v in AP if v % 2 == 0]      # 2,4,6,8,10,12
    best = (Fraction(1), None)
    ce = []
    n_tight = 0
    for r in range(1, len(evens) + 1):
        for sub in combinations(evens, r):
            V = [2 * v if v in sub else v for v in AP]
            if len(set(V)) != 13:
                continue
            Vp, M, t, isce, line = report(f"double{sub}", V)
            if M < best[0]:
                best = (M, Vp)
            if M == THR:
                n_tight += 1
            if isce:
                ce.append(line)
    print(f"   subsets tested; min M found = {best[0]} = {float(best[0]):.6f} at {best[1]}")
    print(f"   how many stayed exactly tight (M=1/14): {n_tight}")
    print(f"   counterexamples: {len(ce)}")
    for l in ce[:10]:
        print("   " + l)
    print()
    return ce


def partB():
    print("== B. doubled-AP_k even part {2,...,2k} + (13-k) odd speeds (e=k>=7) ==")
    ce = []
    glob = (Fraction(1), None)
    for k in range(7, 13):                       # 7..12 even speeds
        evens = [2 * i for i in range(1, k + 1)]  # {2,4,...,2k}, fold={1..k}
        need = 13 - k
        # choose the (13-k) odd speeds to MINIMISE M (search small odds + some big)
        odd_pool = [x for x in range(1, 60) if x % 2 == 1]
        bestk = (Fraction(1), None)
        rng = random.Random(100 + k)
        # exhaustive small if need small, else random
        if need <= 2:
            combos = combinations(odd_pool[:18], need)
        else:
            combos = (tuple(rng.sample(odd_pool, need)) for _ in range(3000))
        for oc in combos:
            V = evens + list(oc)
            if len(set(V)) != 13:
                continue
            Vp = primitive(V)
            if float_measure(Vp) > 1e-7:         # clearly loose, skip exact
                continue
            M, t = M_exact(Vp)
            if M < bestk[0]:
                bestk = (M, Vp)
            if M < THR:
                ce.append(f"k={k}: {Vp} M={M}<1/14 !!!")
        print(f"   k={k} (fold={{1..{k}}}, M(fold)=1/{k+1}): min M over odd choices = "
              f"{bestk[0]}={float(bestk[0]):.6f} at {bestk[1]}")
        if bestk[0] < glob[0]:
            glob = bestk
    print(f"   global min over B = {glob[0]}={float(glob[0]):.6f} at {glob[1]}")
    print(f"   counterexamples: {len(ce)}")
    for l in ce[:10]:
        print("   " + l)
    print()
    return ce


def partC(trials=4000, seed=7):
    print("== C. mod-7 / 7-resonant random configs ==")
    rng = random.Random(seed)
    ce = []
    glob = (Fraction(1), None)
    for _ in range(trials):
        # bias toward multiples of 7 and small residues mod 7
        V = set()
        while len(V) < 13:
            if rng.random() < 0.4:
                V.add(7 * rng.randint(1, 8))
            else:
                V.add(rng.randint(1, 90))
        Vp = primitive(tuple(V))
        if len(Vp) != 13:
            continue
        if float_measure(Vp) > 1e-7:
            continue
        M, t = M_exact(Vp)
        if M < glob[0]:
            glob = (M, Vp)
        if M < THR:
            ce.append(f"{Vp} M={M}<1/14 !!!")
    print(f"   min M = {glob[0]}={float(glob[0]):.6f} at {glob[1]}")
    print(f"   counterexamples: {len(ce)}")
    print()
    return ce


def partD(restarts=300, steps=300, cap=120, seed=5):
    print(f"== D. measure-minimising hill-climb forced e>=7, 7-seeded (cap {cap}) ==")
    rng = random.Random(seed)
    ce = []
    glob = (Fraction(1), None)
    seeds = [tuple(range(1, 14)),
             (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24),
             tuple([2 * i for i in range(1, 13)] + [1]),     # doubled AP_12 + 1
             tuple([2 * i for i in range(1, 8)] + [1, 3, 5, 9, 11, 13])]
    for r in range(restarts):
        if r < len(seeds):
            V = list(seeds[r])
        else:
            # random with forced e>=7 even speeds
            ev = rng.sample([x for x in range(2, cap + 1, 2)], 7)
            od = rng.sample([x for x in range(1, cap + 1, 2)], 6)
            V = ev + od
        V = list(dict.fromkeys(V))
        while len(V) < 13:
            V.append(rng.randint(1, cap))
            V = list(dict.fromkeys(V))
        cur = float_measure(tuple(V))
        for _ in range(steps):
            i = rng.randrange(13)
            old = V[i]
            V[i] = rng.randint(1, cap)
            if len(set(V)) != 13:
                V[i] = old
                continue
            nm = float_measure(tuple(V))
            if nm <= cur:
                cur = nm
            else:
                V[i] = old
        Vp = primitive(tuple(V))
        M, t = M_exact(Vp)
        if M < glob[0]:
            glob = (M, Vp)
        if M < THR:
            ce.append(f"{Vp} M={M}<1/14 !!!")
    print(f"   global min M = {glob[0]}={float(glob[0]):.6f} at {glob[1]}  "
          f"(1/14={float(THR):.6f})")
    print(f"   counterexamples: {len(ce)}")
    for l in ce[:10]:
        print("   " + l)
    print()
    return ce


if __name__ == "__main__":
    allce = []
    allce += partA()
    allce += partB()
    allce += partC()
    allce += partD()
    print("================ VERDICT ================")
    print(f"  total verified counterexamples (M<1/14): {len(allce)}")
    if not allce:
        print("  => NO disproof: every 7-structured / e>=7 / doubled-AP / hill-climb")
        print("     config has M >= 1/14. The minimum collar found is 1/14 (the tight")
        print("     family). Consistent with LRC(14) TRUE; the 7-fold lower-bounds the")
        print("     protected e<=6 half, and the e>=7 hunt is empty.")
    else:
        for l in allce[:20]:
            print("   " + l)
