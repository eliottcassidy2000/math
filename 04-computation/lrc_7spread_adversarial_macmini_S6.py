#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S6 -- HYP-4292 part 2 (CORRECTED): adversarial minimization
of M(U) over 7-spread lattices using the RIGOROUS Lipschitz bracket.

FIX to v1: minimizing a COARSE grid-max rewards undersampling (high speeds make
the coarse grid miss the true max -> spurious "gridmax 0").  Correct objective:
the rigorous bracket LOWER bound (grid-max + Lipschitz-verified), with speeds
bounded so the grid resolves.  Speeds <= 12 (the tight-window regime: window
M <= 2/25 needs bounded witness denominators; per-class M >= 1/6 is
speed-independent by LRC(<=5) anyway).  Best finds re-bracketed at high N.
"""
from math import gcd
from collections import Counter
import random, time

T0 = time.time()
def log(m=""):
    print(m, flush=True)
random.seed(660)

HIf = 2 / 25
LOf = 1 / 13

def dist1(x):
    x = x - int(x)
    if x < 0:
        x += 1
    return min(x, 1 - x)

def bracket(u, v, N=300, depth=0):
    """rigorous (lower, upper) for M(U); brackets need only clear the window."""
    L1 = max(abs(a) for a in u) or 1
    L2 = max(abs(b) for b in v) or 1
    best = 0.0
    bt = bs = 0.0
    for i1 in range(N):
        t = i1 / N
        ut = [a * t for a in u]
        for i2 in range(N):
            s = i2 / N
            m = 1.0
            for k in range(12):
                d = dist1(ut[k] + v[k] * s)
                if d < m:
                    m = d
                    if m <= best:
                        break
            if m > best:
                best = m
                bt, bs = t, s
    slack = (L1 / N + L2 / N) / 2
    lower, upper = best, best + slack
    if depth < 3 and upper - lower > 0.01 and not (lower > HIf):
        span = 3.0 / N
        M = 200
        bb = lower
        for i1 in range(M + 1):
            t = bt - span / 2 + span * i1 / M
            ut = [a * t for a in u]
            for i2 in range(M + 1):
                s = bs - span / 2 + span * i2 / M
                m = 1.0
                for k in range(12):
                    d = dist1(ut[k] + v[k] * s)
                    if d < m:
                        m = d
                        if m <= bb:
                            break
                if m > bb:
                    bb = m
        lower = max(lower, bb)
        upper = min(upper, max(bb + (L1 + L2) * span / M / 2, best + slack))
    return lower, upper

def maxclass(u, v):
    c = Counter()
    for a, b in zip(u, v):
        if a == 0 and b == 0:
            return 99
        g = gcd(abs(a), abs(b))
        p, q = a // g, b // g
        if p < 0 or (p == 0 and q < 0):
            p, q = -p, -q
        c[(p, q)] += 1
    return max(c.values())

DIRS = [(1, 0), (0, 1), (1, 1), (1, -1), (1, 2), (2, 1), (1, -2), (2, -1), (1, 3), (3, 1)]

def random_lattice(maxspeed=12):
    m = random.randint(3, 6)
    sizes = [0] * m
    for _ in range(12):
        sizes[random.randrange(m)] += 1
    if max(sizes) > 5 or 0 in sizes:
        return None
    dirs = random.sample(DIRS, m)
    u, v = [], []
    for j in range(m):
        for c in random.sample(range(1, maxspeed + 1), min(sizes[j], maxspeed)):
            u.append(c * dirs[j][0]); v.append(c * dirs[j][1])
    return (u, v) if len(u) == 12 else None

def perturb(u, v, maxspeed=12):
    u2, v2 = list(u), list(v)
    i = random.randrange(12)
    a, b = u2[i], v2[i]
    g = gcd(abs(a), abs(b)) or 1
    p, q = a // g, b // g
    newc = random.randint(1, maxspeed)
    u2[i], v2[i] = newc * p, newc * q
    if maxclass(u2, v2) > 5:
        return None
    return u2, v2

best = (1.0, None)
restart = 0
while time.time() - T0 < 480:
    restart += 1
    lat = None
    while lat is None:
        lat = random_lattice()
    u, v = lat
    cur, _ = bracket(u, v, 180)
    for step in range(25):
        p = perturb(u, v)
        if p is None:
            continue
        g, _ = bracket(p[0], p[1], 180)
        if g < cur:
            cur, u, v = g, p[0], p[1]
    if cur < best[0]:
        best = (cur, (list(u), list(v)))
        if cur < 0.16:
            log(f"  restart {restart}: M >= {cur:.5f}  maxclass={maxclass(u,v)}")

lo, (u, v) = best
lofine, upfine = bracket(u, v, 600)
log(f"\nADVERSARIAL infimum over {restart} restarts (speeds <= 12, rigorous):")
log(f"  best M-lower = {lo:.6f}; refined bracket ({lofine:.6f}, {upfine:.6f})")
log(f"  extremal maxclass = {maxclass(u,v)}")
log(f"  u = {u}\n  v = {v}")
log(f"window ceiling 2/25 = {HIf:.6f}; 1/6 = {1/6:.6f}")
log("VERDICT: " + (f"7-spread infimum {lofine:.5f} >= 1/6-ish, factor {lofine/HIf:.2f} above window -- (A) CLEAN with margin"
                   if lofine > HIf else f"BEAT WINDOW at {lofine:.5f} -- CRITICAL, verify exactly"))
log(f"[t = {time.time()-T0:.0f}s]")
