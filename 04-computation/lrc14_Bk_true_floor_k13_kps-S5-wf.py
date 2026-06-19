#!/usr/bin/env python3
"""
lrc14_Bk_true_floor_k13_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)

The stretch test showed: the cap-14 bounded-spread minimum 5367/35035 is NOT the
true mu_min for k=12,13 -- a DILATION of it gives mu ~ 0.130992 < 0.153190. So the
bounded-spread-14 reduction is REFUTED at large k; the true infimum is lower and lives
at larger spread. We:

 (1) Find the EXACT mu of the beating dilation for k=12,13 (reproduce 0.130992 as a rational).
 (2) Heavy EXACT search for the lowest mu at k=13 over a wide spread range (the true mu_min(13)
     candidate): report the exact rational and shape, and the spread where it lives.
 (3) Test whether mu_min(13) keeps dropping as we allow ever-larger spread, OR stabilizes at
     some moderate spread (the bounded-spread reduction may still hold, just with a LARGER bound
     than 14). Report the trend: lowest mu at spread-cap S for S = 14, 20, 30, 50, 80, 140.
     If it stabilizes -> reduction holds with bigger B(13). If it keeps dropping -> the floor is
     the equidistribution ceiling F(13) and the reduction is to a genuinely asymptotic statement.
 (4) Report F(13) (the iid ceiling) for reference: P(13 iid uniform pts have maxgap>2/7).
Exact Fractions throughout for mu; F(13) exact via the inclusion-exclusion formula.
"""
import sys, random, itertools
from fractions import Fraction as F
from math import gcd, comb
sys.stdout.reconfigure(line_buffering=True)
random.seed(131313)
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

print("="*90)
print("(1) EXACT mu of the beating DILATION of the cap-14 minimizers for k=12,13.")
print("="*90)
REF12 = (0,1,2,3,4,5,6,7,9,12,13,14)
REF13 = (0,1,2,3,4,5,6,7,8,9,12,13,14)
for name, Estar in [("k=12", REF12), ("k=13", REF13)]:
    base = mu(Estar)
    print(f"\n  {name}: cap-14 minimizer E*={Estar}  mu(E*)={base}={float(base):.6f}")
    best = (base, Estar, "identity")
    for num,den in [(3,2),(2,1),(5,2),(3,1),(7,2),(4,1),(9,2),(5,1),(11,2),(6,1),(7,1),(8,1),(10,1),(13,1),(17,1),(19,1),(23,1)]:
        lam = F(num,den)
        img = sorted(set(int((e*lam).__round__()) for e in Estar))
        if len(img) != len(Estar): continue
        img = primitive(tuple(img))
        if not is_primitive(img): continue
        m = mu(img)
        if m < best[0]: best = (m, img, f"lambda={num}/{den}")
    print(f"     lowest dilation mu = {best[0]} = {float(best[0]):.6f}  at E={best[1]} ({best[2]}, spread {max(best[1])})")

print()
print("="*90)
print("(2)-(3) TRUE mu_min(13) trend: lowest EXACT mu over primitive E (0 in E,|E|=13) with")
print("    spread <= S, for S = 14,20,30,50,80,140.  Does it stabilize or keep dropping?")
print("="*90)
def lowest_mu_under_cap(k, cap, ntrial):
    best = (F(2), None)
    # include the structured dilations and consecutive-ish shapes deterministically
    cons = tuple(range(k))
    for seed_shape in [cons, REF13, REF12+(15,) if k==13 else REF12]:
        if len(set(seed_shape))==k:
            E = primitive(tuple(seed_shape))
            if is_primitive(E) and max(E)<=cap:
                m = mu(E)
                if m < best[0]: best=(m,E)
    for _ in range(ntrial):
        sp = random.randint(k-1, cap)
        rest = tuple(sorted(random.sample(range(1, sp+1), k-1)))
        E = primitive((0,)+rest)
        if len(set(E)) != k or max(E)>cap: continue
        m = mu(E)
        if m < best[0]: best = (m, E)
    return best
prev = None
for S in (14, 20, 30, 50, 80, 140):
    b = lowest_mu_under_cap(13, S, ntrial=4000)
    tag = "" if prev is None else ("  DROP" if b[0] < prev else "  =/up")
    print(f"   spread<={S:3d}: lowest mu = {b[0]} = {float(b[0]):.6f}{tag}  spread={max(b[1])}  E={b[1]}")
    prev = b[0]

print()
print("="*90)
print("(4) F(13) iid ceiling = P(13 iid uniform points have circular max-gap > 2/7).")
print("    F(k) = sum_{j>=1} (-1)^{j+1} C(k,j) max(0, 1-2j/7)^{k-1}.")
print("="*90)
def Fk(k):
    tot = F(0); thr = F(2,7)
    for j in range(1, k+1):
        base = 1 - j*thr
        if base <= 0: break
        tot += F((-1)**(j+1)) * comb(k, j) * (base**(k-1))
    return tot
for k in (7, 10, 12, 13):
    fk = Fk(k)
    print(f"   F({k}) = {fk} = {float(fk):.6f}")
print("\n   NOTE: F(13) is the LARGE-spread/equidistribution value; the TRUE mu_min(13) is the")
print("   infimum and lies AT or BELOW the moderate-spread minima found in (2)-(3).  The gap")
print("   between the bounded-spread-14 min (0.1532) and the true infimum is the residual.")
print("\nDONE.")
