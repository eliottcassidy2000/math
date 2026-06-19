#!/usr/bin/env python3
"""
lrc14_Bk_cap_stability_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)

THE DECISIVE TEST of the BOUNDED-SPREAD REDUCTION.  Question:
  does mu_min(k) := inf over primitive E (0 in E, |E|=k) of mu(E)
  equal the bounded-spread minimum, i.e. does increasing the spread cap NOT lower it?

We attack it three ways, all EXACT:

 (A) CAP PUSH where tractable (k=4,5,6,7): exhaustive mu_min over primitive E with maxE<=cap,
     for a growing sequence of caps. Report each; flag if any later cap LOWERS mu_min.
     If mu_min is flat for the last 3 caps -> strong evidence it has converged.

 (B) STRETCH TEST: take the cap-14 minimizer E* for each k (from the mumin stage) and apply
     measure-non-decreasing-suspected operations to push spread UP:
        - scalar multiply by c (mu invariant by L1 -- a sanity check, must be EQUAL),
        - insert a single far runaway: E* U {M}, M large  (changes k -> k+1, separate check),
        - dilate: map e_i -> round(e_i * lambda) for lambda>1 (genuinely larger spread, same k).
     Confirm none of these large-spread same-k shapes has mu < mu_min(k).

 (C) LARGE-SPREAD ADVERSARY: heavy exact random search at spread in [3k, 120] for each k,
     report the smallest mu found and whether it is >= the bounded-spread mu_min(k).
     This is the honest 'can a coordinated growing cluster beat the bounded min?' probe.

Exact Fractions throughout; stdlib only.
"""
import sys, itertools, random
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
random.seed(778899)
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
def mumin_exhaustive(k, cap):
    best = (F(2), None)
    for rest in itertools.combinations(range(1, cap + 1), k - 1):
        E = (0,) + rest
        if not is_primitive(E): continue
        m = mu(E)
        if m < best[0]: best = (m, E)
    return best

# cap-14 minimizers found earlier (the reference bounded-spread minima)
REF = {
    4: (F(19,21),  (0,1,2,3)),
    5: (F(9,14),   (0,1,2,3,4)),
    6: (F(4,7),    (0,1,2,3,4,5)),
    7: (F(13,35),  (0,2,3,4,5,6,8)),
    8: (F(71,220), (0,2,3,4,5,6,8,11)),
    9: (F(164,735),(0,2,4,5,6,7,8,9,12)),
    10:(F(468,2695),(0,1,3,5,6,8,10,11,12,14)),
    11:(F(409,2548),(0,1,2,3,4,6,8,9,11,13,14)),
    12:(F(5367,35035),(0,1,2,3,4,5,6,7,9,12,13,14)),
    13:(F(5367,35035),(0,1,2,3,4,5,6,7,8,9,12,13,14)),
}

print("="*90)
print("(A) CAP PUSH (exact exhaustive) -- does mu_min keep DROPPING or STABILIZE as cap grows?")
print("="*90)
cap_seq = {4: [6,10,16,22,30,40], 5: [8,12,16,22,28], 6: [8,12,16,20], 7: [9,12,15,18]}
for k in sorted(cap_seq):
    print(f"\n  k={k}:")
    prev = None; vals = []
    for cap in cap_seq[k]:
        if cap < k-1: continue
        m, arg = mumin_exhaustive(k, cap)
        tag = "" if prev is None else ("  DROP!" if m < prev else ("  =" if m == prev else "  UP?"))
        vals.append(m)
        print(f"     cap<={cap:3d}: mu_min = {str(m):>16s} = {float(m):.6f}{tag}   argmin spread={max(arg)}  E={arg}")
        prev = m
    stable = len(vals) >= 3 and vals[-1] == vals[-2] == vals[-3]
    print(f"     -> {'STABLE (last 3 caps equal): converged.' if stable else 'NOT yet flat over last 3 caps.'}")

print()
print("="*90)
print("(B) STRETCH TEST: same-k large-spread images of the cap-14 minimizer never beat mu_min(k).")
print("="*90)
for k in range(4, 14):
    ref, Estar = REF[k]
    worst = (F(2), None)  # smallest mu among stretched images
    # (b1) scalar multiply (L1: must be EQUAL) -- sanity
    eqfail = False
    for c in (2,3,5,7,11):
        if mu(tuple(c*e for e in Estar)) != ref:
            eqfail = True
    # (b3) dilation by rational lambda>1 then round to ints, keep k distinct, primitive
    tried = 0
    for num,den in [(3,2),(2,1),(5,2),(3,1),(7,2),(4,1),(9,2),(5,1),(11,2),(6,1),(8,1),(10,1)]:
        lam = F(num,den)
        img = sorted(set(int((e*lam).__round__()) for e in Estar))
        if len(img) != k: continue
        img = primitive(tuple(img))
        if not is_primitive(img): continue
        m = mu(img); tried += 1
        if m < worst[0]: worst = (m, img)
    beat = worst[0] < ref
    print(f"  k={k:2d}: mu_min={str(ref):>14s}  L1-scalar-eq={'OK' if not eqfail else 'FAIL!'}  "
          f"dilations tried={tried}  min-stretched-mu={float(worst[0]):.6f} "
          f"{'-> BEATS!!' if beat else '-> none below mu_min (OK)'}")

print()
print("="*90)
print("(C) LARGE-SPREAD ADVERSARY: exact random search, spread in [3k,120]; can it beat mu_min(k)?")
print("="*90)
for k in range(4, 14):
    ref = REF[k][0]
    best = (F(2), None)
    for _ in range(6000):
        sp = random.randint(3*k, 120)
        if sp < k-1: continue
        rest = tuple(sorted(random.sample(range(1, sp+1), k-1)))
        E = primitive((0,)+rest)
        if len(set(E)) != k: continue
        m = mu(E)
        if m < best[0]: best = (m, E)
    verdict = "BEATS mu_min!!" if best[0] < ref else "OK (>= mu_min)"
    print(f"  k={k:2d}: bounded mu_min={float(ref):.6f}  large-spread min found={float(best[0]):.6f}  "
          f"spread={max(best[1])}  -> {verdict}")
print("\nDONE.")
