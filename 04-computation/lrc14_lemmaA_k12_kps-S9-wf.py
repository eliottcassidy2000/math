#!/usr/bin/env python3
"""Standalone fast Lemma A + k=12 witness (the parts the heavy run did not reach).
LEMMA A: P(N=6)(E) = meas(G_E) <= 1/(7(k-1)) for primitive E, equality iff consec.
Component-count bound: #comp(G_E) <= max(E)/(k-1).
Uses float P(N=6) screen, exact confirm on any flagged violation."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import random, sys
sys.stdout.reconfigure(line_buffering=True)

def frac(x):
    r = x - int(x)
    return r + 1 if r < 0 else r

def cells_exact(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(e):
            for j in range(7):
                bps.add(F(7 * m + j, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    out = []
    for a, b in zip(bps, bps[1:]):
        if b <= a: continue
        xm = (a + b) / 2
        allzero = all((int(frac(e * xm) * 7) % 7) == 0 for e in E if e != 0)
        out.append((b - a, allzero))  # allzero => N=6 (G_E membership)
    return out

def PN6_exact(E):
    return sum((L for L, z in cells_exact(E) if z), F(0))

def ncomp_exact(E):
    cells = cells_exact(E)
    runs = 0; prev = False
    for L, z in cells:
        if z and not prev: runs += 1
        prev = z
    if cells and cells[0][1] and cells[-1][1] and runs >= 2:
        runs -= 1
    return runs

def PN6_float(E):
    """fast float screen of meas(G_E)."""
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    if not Enz: return 0.0
    bps = {0.0, 1.0}
    for e in Enz:
        for m in range(e):
            bps.add(m / e); bps.add(m / e + 1 / (7.0 * e))
    bps = sorted(b for b in bps if 0.0 <= b <= 1.0)
    tot = 0.0
    for a, b in zip(bps, bps[1:]):
        if b <= a: continue
        xm = 0.5 * (a + b)
        if all(((e * xm) % 1.0) < 1.0 / 7.0 for e in E if e != 0):
            tot += (b - a)
    return tot

def is_primitive(E):
    g = 0
    for e in E:
        if e: g = gcd(g, e)
    return g == 1

print("=" * 72)
print("LEMMA A — exhaustive bounded: P(N=6) <= 1/(7(k-1)) + comp-count <= max(E)/(k-1)")
print("=" * 72)
for k, B in [(8, 16), (9, 14), (10, 13), (11, 13)]:
    bound = F(1, 7 * (k - 1)); maxv = F(0); arg = None
    nviol = 0; ccv = 0; cce = []; cnt = 0; eq_noncon = 0
    for rest in combinations(range(1, B + 1), k - 1):
        E = (0,) + rest
        if not is_primitive(E): continue
        cnt += 1
        pn6 = PN6_exact(list(E))
        if pn6 > maxv: maxv = pn6; arg = E
        if pn6 > bound: nviol += 1
        if pn6 == bound and E != tuple(range(k)): eq_noncon += 1
        nc = ncomp_exact(list(E))
        if F(nc) > F(max(E), k - 1):
            ccv += 1; cce.append((nc, max(E), E))
    print(f" k={k} B={B} prim#={cnt}: bound={bound}={float(bound):.6f} max={maxv}={float(maxv):.6f} "
          f"arg={arg} #viol_len={nviol} #viol_compcount={ccv} #equality-noncon={eq_noncon}")
    if cce: print("    CC VIOL:", cce[:5])

print("\n[Lemma A wide/resonant random hunt — float screen, exact confirm]")
random.seed(7)
def gen_wide(k):
    span = random.randint(20, 250); s = {0}
    while len(s) < k: s.add(random.randint(1, span))
    return sorted(s)
def gen_res(k):
    s = {0}
    while len([x for x in s if x % 7 == 0]) < k - 2:
        s.add(7 * random.randint(1, 30))
    while len(s) < k:
        x = random.randint(1, 250)
        if x % 7: s.add(x)
    return sorted(s)[:k]
for k in [8, 9, 10, 11]:
    bound = F(1, 7 * (k - 1)); bf = float(bound)
    flagged_len = []; flagged_cc = []; maxf = 0.0; arg = None; tried = 0
    for _ in range(60000):
        E = gen_wide(k) if random.random() < 0.6 else gen_res(k)
        if len(E) != k or not is_primitive(E): continue
        tried += 1
        v = PN6_float(E)
        if v > maxf: maxf = v; arg = tuple(E)
        if v > bf + 1e-9: flagged_len.append(tuple(E))
    real = [E for E in flagged_len if PN6_exact(list(E)) > bound]
    print(f"  k={k}: bound={bf:.6f} maxP(N=6)float={maxf:.6f} arg={arg} float-flag={len(flagged_len)} EXACT-viol={len(real)}")
    if real: print("     CONFIRMED LEMMA A VIOLATION:", real[:3])

print("\n" + "=" * 72)
print("k=12 Lemma-B failure witness (exact)")
print("=" * 72)
def cells_S7(E):
    E = sorted(set(E)); Enz = [e for e in E if e != 0]
    bps = {F(0), F(1)}
    for e in Enz:
        for m in range(e):
            for j in range(7):
                bps.add(F(7 * m + j, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for a, b in zip(bps, bps[1:]):
        if b <= a: continue
        xm = (a + b) / 2
        hit = set()
        for e in E:
            s = int(frac(e * xm) * 7); hit.add(6 if s == 7 else s)
        if len(hit) == 7: tot += (b - a)
    return tot
cap12 = F(6, 7)
v = cells_S7(list(range(11)) + [12]); c = cells_S7(list(range(12)))
print(f"  E=[0..10,12]: meas_S7={v}={float(v):.6f} consec_12={float(c):.6f} beats_consec={v>c} over_cap12={v>cap12}")
print("\nDONE")
