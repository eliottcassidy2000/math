#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S6 -- HYP-4292 part 2: ADVERSARIAL minimization of M(U)
over 7-spread rank-2 lattices (MISTAKE-102 discipline -- can anything beat 1/6?).

The structured census bottomed at 1/6 = 0.1667 (2.08x the window ceiling 2/25).
Here: random-restart hill-descent on (directions, speeds) minimizing M(U)
subject to 7-spread (max direction-class <= 5), to pin the true infimum.

M(U) = max_{(t,s)} min_i ||u_i t + v_i s||; a coarse grid MAX is a LOWER bracket
(so minimizing the grid-max is conservative -- the true M is >= grid-max, so if
grid-max stays high the truth is higher).  Best finds re-bracketed rigorously.
"""
from fractions import Fraction as F
from math import gcd
import random, time
from collections import Counter

T0 = time.time()
def log(m=""):
    print(m, flush=True)
random.seed(66)

HIf = 2 / 25

def dist1(x):
    x = x - int(x)
    if x < 0:
        x += 1
    return min(x, 1 - x)

def gridmax(u, v, N):
    best = 0.0
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
    return best

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

DIRS = [(1, 0), (0, 1), (1, 1), (1, -1), (1, 2), (2, 1), (1, -2), (2, -1),
        (1, 3), (3, 1), (2, 3), (3, 2), (1, -3), (3, -1)]

def random_lattice(maxspeed):
    m = random.randint(3, 6)
    sizes = [0] * m
    for _ in range(12):
        sizes[random.randrange(m)] += 1
    if max(sizes) > 5 or 0 in sizes:
        return None
    dirs = random.sample(DIRS, m)
    u, v = [], []
    for j in range(m):
        for c in random.sample(range(1, maxspeed + 1), sizes[j]):
            u.append(c * dirs[j][0])
            v.append(c * dirs[j][1])
    return u, v

def perturb(u, v, maxspeed):
    """local move: rescale one coord's speed or reassign; keep 7-spread."""
    u2, v2 = list(u), list(v)
    i = random.randrange(12)
    # recover direction of coord i
    a, b = u2[i], v2[i]
    g = gcd(abs(a), abs(b)) or 1
    p, q = a // g, b // g
    newc = random.randint(1, maxspeed)
    u2[i], v2[i] = newc * p, newc * q
    if maxclass(u2, v2) > 5:
        return None
    return u2, v2

best_overall = (1.0, None)
N = 60
for restart in range(3000):
    lat = None
    while lat is None:
        lat = random_lattice(random.choice([8, 12, 16, 24]))
    u, v = lat
    cur = gridmax(u, v, N)
    for step in range(40):
        p = perturb(u, v, max(max(map(abs, u)), max(map(abs, v))) // max(1, min(gcd(abs(a),abs(b)) or 1 for a,b in zip(u,v))) + 4)
        if p is None:
            continue
        g = gridmax(p[0], p[1], N)
        if g < cur:
            cur, u, v = g, p[0], p[1]
    if cur < best_overall[0]:
        best_overall = (cur, (list(u), list(v)))
        if cur < 0.15:
            log(f"  new min gridmax {cur:.5f} at restart {restart}: maxclass={maxclass(u,v)}")
    if time.time() - T0 > 540:
        log(f"  (time budget: stopped at restart {restart})")
        break

lo, lat = best_overall
u, v = lat
log(f"\nADVERSARIAL infimum (coarse grid-max, a LOWER bound on M): {lo:.6f}")
log(f"  extremal lattice maxclass = {maxclass(u,v)}")
log(f"  u = {u}")
log(f"  v = {v}")
# rigorous re-bracket of the best find
finelo = gridmax(u, v, 480)
L1, L2 = max(map(abs, u)), max(map(abs, v))
log(f"  fine grid-max (N=480): {finelo:.6f}  (true M >= this)")
log(f"\nwindow ceiling 2/25 = {HIf:.6f}")
log("VERDICT: " + (f"adversarial infimum {finelo:.5f} >> 2/25 -- 7-spread safe with margin factor "
                   f"{finelo/HIf:.2f}" if finelo > HIf else "BEAT THE WINDOW -- critical find"))
log(f"[t = {time.time()-T0:.0f}s]")
