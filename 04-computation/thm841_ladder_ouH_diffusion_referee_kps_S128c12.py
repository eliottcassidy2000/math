#!/usr/bin/env python3
"""
thm841_ladder_ouH_diffusion_referee_kps_S128c12.py
==================================================
kind-pasteur-2026-07-15-S128 (cont.12).
 (A) THM-841 referee: support lemma + ladder measures m_r + factorial moments S_rho on Farey gaps
     vs direct computation; k <= 10, rational lambda grid, exact.
 (B) OU on the H axis: is E[Delta H | T] (uniform arc flip) a linear function of H? of log H?
     Exact per-tournament data n=4,5; report best affine law + residuals, or honest negative.
 (C) The diffusion invariant D(T) = E[(dc3)^2 | T]: exact per class; test affine models in
     (c3, c3^2, V) -- is the OU noise class-determined by degree data alone? (It should NOT be:
     Sum_{arcs}(d_u-d_v)^2 is not degree-determined -- quantify the spread within (n, c3) fibers.)
"""
import sys
from fractions import Fraction as F
from math import comb, gcd
from collections import defaultdict
sys.stdout.reconfigure(line_buffering=True)

# ---------- (A) THE VIOLATION LADDER ----------
print("=" * 96)
print("(A) THM-841 violation ladder referee")
print("=" * 96)

def farey(k):
    return sorted(set(F(a, i) for i in range(1, k + 1) for a in range(i)))

def W_direct(k, lam, t):
    return sum(1 for j in range(1, k + 1) if abs(j * t - round(j * t)) < lam)

def ladder_direct(k, lam, r):
    """measure of {t : W(t) < r} by exact interval arithmetic: collect all arcs, count coverage."""
    events = []
    for j in range(1, k + 1):
        rad = lam / j
        for a in range(j):
            c = F(a, j)
            for lo, hi in [(c - rad, c + rad)]:
                events.append((lo, 1)); events.append((hi, -1))
                if lo < 0:
                    events.append((lo + 1, 1)); events.append((hi + 1, -1))
                if hi > 1:
                    events.append((lo - 1, 1)); events.append((hi - 1, -1))
    events.sort()
    tot = F(0); depth = 0; prev = F(0)
    for x, d in events:
        if x > 1: break
        if x > prev:
            if depth < r:
                tot += min(x, F(1)) - prev
            prev = x
        depth += d
    if prev < 1 and depth < r:
        tot += 1 - prev
    return tot

def ladder_gap(k, lam, r):
    """THM-841 formula: per gap, W = two-sided divisor-chain staircase; integrate {W < r} exactly."""
    fr = farey(k)
    tot = F(0)
    for x, y in zip(fr, fr[1:] + [F(1)]):
        i = x.denominator; j = (y.denominator if y != 1 else 1)
        g = y - x
        # thresholds from left (multiples of i) and right (multiples of j)
        L = sorted([lam / (m * i) for m in range(1, k // i + 1)], reverse=True)   # s < L[m]: violation count
        R = sorted([lam / (m * j) for m in range(1, k // j + 1)], reverse=True)
        # W(s) = #{l in L : s < l} + #{rr in R : g - s < rr}; integrate over s in [0, g] where W < r
        cuts = sorted(set([F(0), g] + [l for l in L if 0 < l < g] + [g - rr for rr in R if 0 < g - rr < g]))
        for a, b in zip(cuts, cuts[1:]):
            mid = (a + b) / 2
            W = sum(1 for l in L if mid < l) + sum(1 for rr in R if g - mid < rr)
            if W < r:
                tot += b - a
    return tot

allA = True
for k in [4, 6, 8, 10]:
    lams = [F(1, k + 2), F(1, 2 * k), F(1, 3 * k), F(2, 3 * (k + 1))]
    for lam in lams:
        if lam >= F(1, k + 1):
            continue
        for r in [1, 2, 3, 4]:
            d = ladder_direct(k, lam, r)
            gp = ladder_gap(k, lam, r)
            if d != gp:
                allA = False
                print("  MISMATCH k=%d lam=%s r=%d: direct=%s gap=%s" % (k, lam, r, d, gp))
    print("k=%2d: ladder m_r (r<=4) exact on %d lambdas: %s" % (k, len(lams), allA))
# S_rho factorial moments: S_rho = integral C(W,rho) dt -- direct vs gap
def Srho_direct(k, lam, rho):
    # integrate C(W,rho): use ladder measures: C(W,rho) = sum_{r} [W >= ...]... do directly by depth sweep
    events = []
    for j in range(1, k + 1):
        rad = lam / j
        for a in range(j):
            c = F(a, j)
            events.append((c - rad, 1)); events.append((c + rad, -1))
            if c - rad < 0: events.append((c - rad + 1, 1)); events.append((c + rad + 1, -1))
            if c + rad > 1: events.append((c - rad - 1, 1)); events.append((c + rad - 1, -1))
    events.sort()
    tot = F(0); depth = 0; prev = F(0)
    for x, d in events:
        if x > 1: break
        if x > prev:
            tot += comb(depth, rho) * (min(x, F(1)) - prev)
            prev = x
        depth += d
    if prev < 1:
        tot += comb(depth, rho) * (1 - prev)
    return tot

def Srho_gap(k, lam, rho):
    fr = farey(k)
    tot = F(0)
    for x, y in zip(fr, fr[1:] + [F(1)]):
        i = x.denominator; j = (y.denominator if y != 1 else 1)
        g = y - x
        L = [lam / (m * i) for m in range(1, k // i + 1)]
        R = [lam / (m * j) for m in range(1, k // j + 1)]
        cuts = sorted(set([F(0), g] + [l for l in L if 0 < l < g] + [g - rr for rr in R if 0 < g - rr < g]))
        for a, b in zip(cuts, cuts[1:]):
            mid = (a + b) / 2
            W = sum(1 for l in L if mid < l) + sum(1 for rr in R if g - mid < rr)
            tot += comb(W, rho) * (b - a)
    return tot

for k in [6, 10]:
    lam = F(1, k + 2)
    ok = all(Srho_direct(k, lam, rho) == Srho_gap(k, lam, rho) for rho in [1, 2, 3])
    print("k=%2d lam=1/%d: S_rho (rho<=3) direct==gap: %s  (S_1=%s S_2=%s)"
          % (k, k + 2, ok, Srho_gap(k, lam, 1), Srho_gap(k, lam, 2)))
    allA &= ok
print("(A) VERDICT: %s" % ("ALL EXACT" if allA else "FAILURES"))

# ---------- (B) OU on the H axis ----------
print()
print("=" * 96)
print("(B) OU on H: E[Delta H | T] vs H (n=4,5) -- law or negative")
print("=" * 96)
def ham(B, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1 << n):
        for v in range(n):
            c = dp[S][v]
            if not c or not (S >> v) & 1: continue
            for u in range(n):
                if not (S >> u) & 1 and B[v][u]:
                    dp[S | 1 << u][u] += c
    return sum(dp[(1 << n) - 1][v] for v in range(n))

for n in [4, 5]:
    m2 = comb(n, 2)
    tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]
    tidx = {t: i for i, t in enumerate(tiles)}
    mm = len(tiles)
    pts = defaultdict(list)
    for t in range(1 << mm):
        tv = [(t >> i) & 1 for i in range(mm)]
        B = [[False] * n for _ in range(n)]
        for k2 in range(2, n + 1):
            B[k2 - 1][k2 - 2] = True
        for (x, y), i in tidx.items():
            B[x - 1][y - 1] = tv[i] == 1
            B[y - 1][x - 1] = tv[i] == 0
        H0 = ham(B, n)
        s = F(0)
        for u in range(n):
            for v in range(n):
                if u != v and B[u][v]:
                    B[u][v], B[v][u] = False, True
                    s += ham(B, n) - H0
                    B[u][v], B[v][u] = True, False
        pts[H0].append(s / m2)
    print("n=%d: per-H drift E[dH|T] (H: mean [min..max] count):" % n)
    for H0 in sorted(pts):
        vals = pts[H0]
        mean = sum(vals) / len(vals)
        print("   H=%2d: %8.4f  [%s..%s]  x%d" % (H0, float(mean), min(vals), max(vals), len(vals)))

# ---------- (C) the diffusion invariant ----------
print()
print("=" * 96)
print("(C) D(T) = E[(dc3)^2 | T]: is the OU noise degree-determined? spread within (n, c3) fibers")
print("=" * 96)
for n in [5, 6]:
    tiles = [(x, y) for y in range(1, n - 1) for x in range(n, y + 1, -1) if x - y >= 2]
    tidx = {t: i for i, t in enumerate(tiles)}
    mm = len(tiles)
    m2 = comb(n, 2)
    fib = defaultdict(set)
    for t in range(1 << mm):
        tv = [(t >> i) & 1 for i in range(mm)]
        d = [0] * n
        for k2 in range(2, n + 1):
            d[k2 - 1] += 1
        for (x, y), i in tidx.items():
            d[(x - 1) if tv[i] else (y - 1)] += 1
        c3 = comb(n, 3) - sum(comb(di, 2) for di in d)
        B = [[False] * n for _ in range(n)]
        for k2 in range(2, n + 1):
            B[k2 - 1][k2 - 2] = True
        for (x, y), i in tidx.items():
            B[x - 1][y - 1] = tv[i] == 1
            B[y - 1][x - 1] = tv[i] == 0
        D = sum(F((d[u] - d[v] - 1) ** 2, m2) for u in range(n) for v in range(n) if u != v and B[u][v])
        fib[c3].add(D)
    multi = {c3: sorted(ds) for c3, ds in fib.items() if len(ds) > 1}
    print("n=%d: c3-fibers with MULTIPLE D values: %d of %d  (noise NOT c3-determined: %s)"
          % (n, len(multi), len(fib), bool(multi)))
    for c3 in sorted(multi)[:2]:
        print("   c3=%d: D in %s" % (c3, [str(x) for x in multi[c3]]))
print("\nDONE")
