#!/usr/bin/env python3
"""
lrc14_Bk_Fk_and_trend_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)
Final two numbers for the verdict:
 (A) F(k) iid ceiling exact for k=4..13 (the positive limit mu approaches from below at large spread).
 (B) moderate-spread trend of lowest exact mu(13): spread caps 14,18,22,26 (kept moderate so exact
     arithmetic stays fast). Does mu_min(13) keep dropping below 5367/35035 toward F(13)?
"""
import sys, random
from fractions import Fraction as F
from math import gcd, comb
sys.stdout.reconfigure(line_buffering=True)
random.seed(999)
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
def Fk(k):
    tot = F(0); thr = F(2,7)
    for j in range(1, k+1):
        base = 1 - j*thr
        if base <= 0: break
        tot += F((-1)**(j+1)) * comb(k, j) * (base**(k-1))
    return tot
print("="*90)
print("(A) F(k) iid/equidistribution CEILING (exact). mu(E) -> F(k) from below as spread->inf.")
print("    This is the POSITIVE limit; if mu_min(k) descends to F(k), the floor is >= ... let's see.")
print("="*90)
for k in range(4, 14):
    fk = Fk(k)
    print(f"   F({k:2d}) = {str(fk):>22s} = {float(fk):.6f}")
print()
print("="*90)
print("(B) moderate-spread trend lowest mu(13), caps 14,18,22,26 (fast). Drop below 5367/35035?")
print("="*90)
REF13 = (0,1,2,3,4,5,6,7,8,9,12,13,14)
DIL = (0,2,3,4,6,8,9,10,12,14,18,20,21)  # the lambda=3/2 beating shape (spread 21)
prev = None
for S in (14, 18, 22, 26):
    best = (F(2), None)
    # deterministic seeds
    for seed in [tuple(range(13)), REF13, DIL]:
        E = primitive(tuple(seed))
        if len(set(E))==13 and is_primitive(E) and max(E)<=S:
            m = mu(E)
            if m<best[0]: best=(m,E)
    for _ in range(2500):
        sp = random.randint(12, S)
        rest = tuple(sorted(random.sample(range(1, sp+1), 12)))
        E = primitive((0,)+rest)
        if len(set(E))!=13 or max(E)>S: continue
        m = mu(E)
        if m<best[0]: best=(m,E)
    tag = "" if prev is None else ("  DROP" if best[0]<prev else "  =/up")
    print(f"   spread<={S:3d}: lowest mu(13) = {best[0]} = {float(best[0]):.6f}{tag}  spread={max(best[1])}  E={best[1]}")
    prev = best[0]
print(f"\n   F(13) = {Fk(13)} = {float(Fk(13)):.6f}  (the iid ceiling -- mu approaches this from below at large spread).")
print("\nDONE.")
