#!/usr/bin/env python3
"""
lrc14_Bk_BC_only_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)
ONLY the STRETCH (B) and LARGE-SPREAD ADVERSARY (C) tests for k=4..13.
Confirms no same-k large-spread shape beats the cap-14 bounded-spread mu_min(k). Exact Fractions.
"""
import sys, random
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
random.seed(424242)
TWO7 = F(2, 7)
def merge(iv):
    iv = sorted(iv); out = []
    for a, b in iv:
        if out and a <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], b))
        else: out.append((a, b))
    return out
def meas(arcs): return sum((b - a for a, b in arcs), F(0))
def good_set_exact(E):
    E = sorted(set(E)); k = len(E)
    if k == 1: return [(F(0), F(1))]
    diffs = set()
    for a in range(k):
        for b in range(a + 1, k): diffs.add(E[b] - E[a])
    bps = {F(0), F(1)}
    for d in diffs:
        for m in range(0, d + 1): bps.add(F(m, d))
    bps = sorted(x for x in bps if 0 <= x <= 1)
    good = []
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        pts = sorted(((E[t] * xm) % 1, E[t]) for t in range(k))
        order = [e for _, e in pts]; floors = [int((e * xm) // 1) for e in order]
        for idx in range(k):
            e_cur = order[idx]; f_cur = floors[idx]
            if idx < k - 1: e_nx = order[idx + 1]; f_nx = floors[idx + 1]; wrap = F(0)
            else: e_nx = order[0]; f_nx = floors[0]; wrap = F(1)
            A = F(e_nx - e_cur); Cc = F(f_cur - f_nx) + wrap
            if A == 0:
                if Cc > TWO7: good.append((x0, x1))
                continue
            xb = (TWO7 - Cc) / A
            if A > 0: lo = max(x0, xb); hi = x1
            else: lo = x0; hi = min(x1, xb)
            if lo < hi: good.append((lo, hi))
    return merge(good)
def mu(E): return meas(good_set_exact(E))
def is_primitive(E):
    g = 0
    for e in E: g = gcd(g, e)
    return g == 1
def primitive(E):
    g = 0
    for e in E: g = gcd(g, e)
    return tuple(e // g for e in E) if g > 1 else tuple(E)
REF = {
    4: (F(19,21),  (0,1,2,3)), 5: (F(9,14),(0,1,2,3,4)), 6: (F(4,7),(0,1,2,3,4,5)),
    7: (F(13,35),(0,2,3,4,5,6,8)), 8: (F(71,220),(0,2,3,4,5,6,8,11)),
    9: (F(164,735),(0,2,4,5,6,7,8,9,12)), 10:(F(468,2695),(0,1,3,5,6,8,10,11,12,14)),
    11:(F(409,2548),(0,1,2,3,4,6,8,9,11,13,14)), 12:(F(5367,35035),(0,1,2,3,4,5,6,7,9,12,13,14)),
    13:(F(5367,35035),(0,1,2,3,4,5,6,7,8,9,12,13,14)),
}
print("="*90)
print("(B) STRETCH: same-k DILATIONS of cap-14 minimizer never beat mu_min(k); L1 scalar-eq sanity.")
print("="*90)
for k in range(4, 14):
    ref, Estar = REF[k]
    eqfail = any(mu(tuple(c*e for e in Estar)) != ref for c in (2,3,5,7,11))
    worst = (F(2), None); tried = 0
    for num,den in [(3,2),(2,1),(5,2),(3,1),(7,2),(4,1),(9,2),(5,1),(11,2),(6,1),(7,1),(8,1),(10,1),(13,1),(17,1)]:
        lam = F(num,den)
        img = sorted(set(int((e*lam).__round__()) for e in Estar))
        if len(img) != k: continue
        img = primitive(tuple(img))
        if not is_primitive(img): continue
        m = mu(img); tried += 1
        if m < worst[0]: worst = (m, img)
    beat = worst[0] < ref
    print(f"  k={k:2d}: mu_min={str(ref):>14s}={float(ref):.6f}  L1-eq={'OK' if not eqfail else 'FAIL!'}  "
          f"dilations={tried}  min-stretched={float(worst[0]):.6f} "
          f"{'BEATS!!' if beat else 'none below mu_min (OK)'}")
print()
print("="*90)
print("(C) LARGE-SPREAD ADVERSARY: exact random, spread in [3k,140]; can it beat mu_min(k)?")
print("="*90)
overall_beat = False
for k in range(4, 14):
    ref = REF[k][0]; best = (F(2), None)
    ntr = 9000 if k <= 9 else 4500
    for _ in range(ntr):
        sp = random.randint(3*k, 140)
        if sp < k-1: continue
        rest = tuple(sorted(random.sample(range(1, sp+1), k-1)))
        E = primitive((0,)+rest)
        if len(set(E)) != k: continue
        m = mu(E)
        if m < best[0]: best = (m, E)
    beat = best[0] < ref; overall_beat = overall_beat or beat
    print(f"  k={k:2d}: mu_min={float(ref):.6f}  large-spread min={float(best[0]):.6f}  spread={max(best[1])}  "
          f"-> {'BEATS!!' if beat else 'OK (>= mu_min)'}")
print(f"\n  ==> {'large-spread BEAT the bounded min (reduction in DOUBT)' if overall_beat else 'large-spread NEVER beat bounded mu_min (reduction corroborated)'}.")
print("\nDONE.")
