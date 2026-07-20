#!/usr/bin/env python3
"""
lrc14_x49_gate_boxeph_S137.py  (HYP-8010)

THE x49 GATE LAW, measured at p=43 (Q = 49*43 = 2107) — the first direct test of
"the cage = iterated x7 layers" (HYP-7990's elimination).

DERIVATION (stated; controls verify):  improper mod 49p at the 1/14 threshold =
no b with all dist_{49p}(v*b) > 49p/14 = 3.5p; for odd p the danger band is
+-(7p-1)/2 (fraction ~1/7, scale-free again).  CRT (mod 49, mod p):
  * pure-p scales (b = 49b''):  level-1 mod p (dist_p <= floor(p/14)) — as at
    every gate;
  * pure-49 scales (b = p b'):  window +-3 spread at modulus 49 (unit b'), and
    the SOME-7-MULTIPLE law at gcd(b',49)=7 scales — the x7 divisibility law
    inherited;
  * mixed unit scales: 7-of-49 SLAVING: the covering runner's mod-49 class must
    lie in a 7-term AP (difference p mod 49) determined by its mod-p class —
    and those 7 values are exactly the SEVEN LIFTS of the runner's mod-7 class
    forced at the x7 gate (x = r + j*p, j over 7 consecutive integers; x mod 7
    sweeps all of Z/7 iff p != 0 mod 7).
So the tower semantics: the x7 gate forces each covering runner's mod-7 class
(1-of-7, zero freedom); the x49 gate asks whether that forced skeleton admits
a mod-49 REFINEMENT (7 choices per runner) covering all 1053 scales.

MEASUREMENT at p=43: (i) x7-gate admit from level-1 (greedy, one-sided) — never
measured before; (ii) for x7 witnesses (c43, sign, c7 fixed), the x49 refinement:
MC random-refinement survival + greedy admit (one-sided).  Compare against the
x14 numbers (S135/S136) — the tower-height localization death-star's |S| needs.

Controls (boundary tuples, no strict witness at any modulus): AP13 and GW must
be improper mod 301 AND mod 2107; their identity refinements give admit.

boxeph-2026-07-19-S137.  Pure Python, exact integers.
"""

import random
from math import gcd

P = 43
Q7, H7 = 301, 150          # 7p, semi-scales
Q49, H49 = 2107, 1053      # 49p, semi-scales
DK = 3                     # floor(43/14)
BAND7 = 21                 # floor((7p-1)/14) = floor(p/2)
BAND49 = 150               # (7p-1)/2 = (301-1)/2

def dist(x, m):
    r = x % m
    return r if r <= m - r else m - r

def improper(vals, Q, band):
    return all(any(dist(v * b, Q) <= band for v in vals) for b in range(1, Q // 2 + 1))

def level1(res):
    return all(any(dist(v * b, P) <= DK for v in res) for b in range(1, 22))

# CRT helpers
I7 = pow(P % 7, -1, 7)          # for mod 7p
I49 = pow(P % 49, -1, 49)       # for mod 49p

def crt7p(rp, c7):
    return (rp + P * (((c7 - rp) * I7) % 7)) % Q7

def crt49p(rp, c49):
    return (rp + P * (((c49 - rp) * I49) % 49)) % Q49

# ---------------- controls ------------------------------------------------------
AP13 = list(range(1, 14)); GW = list(range(1, 12)) + [13, 24]
assert improper(AP13, Q7, BAND7) and improper(GW, Q7, BAND7)
assert improper(AP13, Q49, BAND49) and improper(GW, Q49, BAND49)
print("controls: AP13, GW improper mod 301 AND mod 2107 (boundary tuples)  OK")

# ---------------- level-1 bases (same seed family for comparability) -------------
def random_level1(n_want, seed):
    rng = random.Random(seed)
    hitters = {b: [r for r in range(1, 22) if dist(r * b, P) <= DK] for b in range(1, 22)}
    out = []
    while len(out) < n_want:
        chosen, unc = [], set(range(1, 22))
        while unc and len(chosen) < 13:
            b = min(unc, key=lambda x: len(hitters[x]))
            c = hitters[b][:]; rng.shuffle(c)
            chosen.append(c[0])
            unc = {x for x in unc if dist(c[0] * x, P) > DK}
        if not unc:
            while len(chosen) < 13:
                chosen.append(rng.randrange(1, 22))
            res = sorted(chosen)
            if level1(res):
                out.append(tuple(res))
    return out

bases = random_level1(20, seed=1443)

# ---------------- stage 1: x7-gate witnesses from level-1 ------------------------
FULL7 = (1 << H7) - 1
def x7_witness(base, restarts, seed):
    slots = []
    for r in base:
        opts = []
        for eps in (1, -1):
            rp = (eps * r) % P
            for c7 in range(7):
                v = crt7p(rp, c7)
                m = 0
                for b in range(1, H7 + 1):
                    if dist(v * b, Q7) <= BAND7:
                        m |= 1 << (b - 1)
                opts.append((m, v))
        slots.append(opts)
    rng = random.Random(seed)
    for _ in range(restarts):
        covered, used, free = 0, [None] * 13, list(range(13))
        ok = True
        while covered != FULL7:
            bestsc, cands = -1, []
            for i in free:
                for o in range(14):
                    g = bin(slots[i][o][0] & ~covered).count('1')
                    if g > bestsc: bestsc, cands = g, [(i, o)]
                    elif g == bestsc: cands.append((i, o))
            if not free or bestsc <= 0: ok = False; break
            i, o = cands[rng.randrange(len(cands))]
            used[i] = o; free.remove(i); covered |= slots[i][o][0]
        if ok:
            sol = [slots[i][used[i]][1] if used[i] is not None
                   else crt7p(base[i] % P, 1) for i in range(13)]
            if improper(sol, Q7, BAND7) and any(v % 7 for v in sol):
                return sol
    return None

wit7 = []
for i, base in enumerate(bases):
    w = x7_witness(base, 300, 900 + i)
    if w: wit7.append(w)
print("x7-gate admit from level-1 (greedy, one-sided >=): %d/20 witnesses" % len(wit7))

# ---------------- stage 2: the x49 refinement -------------------------------------
FULL49 = (1 << H49) - 1
rng = random.Random(4907)
mc_ok = mc_try = 0
admit = 0
for w in wit7:
    # fixed: rp = v mod 43, c7 = v mod 7; free: c49 in the 7 lifts of c7
    slots = []
    for v in w:
        rp, c7 = v % P, v % 7
        opts = []
        for t in range(7):
            c49 = c7 + 7 * t
            u = crt49p(rp, c49)
            assert u % 7 == c7 and u % P == rp
            m = 0
            for b in range(1, H49 + 1):
                if dist(u * b, Q49) <= BAND49:
                    m |= 1 << (b - 1)
            opts.append((m, u))
        slots.append(opts)
    # (a) MC random refinement
    for _ in range(300):
        vals = [slots[i][rng.randrange(7)][1] for i in range(13)]
        mc_try += 1
        if improper(vals, Q49, BAND49):
            mc_ok += 1
    # (b) greedy admit
    got = False
    for _ in range(250):
        covered, used, free, ok = 0, [None] * 13, list(range(13)), True
        while covered != FULL49:
            bestsc, cands = -1, []
            for i in free:
                for o in range(7):
                    g = bin(slots[i][o][0] & ~covered).count('1')
                    if g > bestsc: bestsc, cands = g, [(i, o)]
                    elif g == bestsc: cands.append((i, o))
            if not free or bestsc <= 0: ok = False; break
            i, o = cands[rng.randrange(len(cands))]
            used[i] = o; free.remove(i); covered |= slots[i][o][0]
        if ok:
            vals = [slots[i][used[i]][1] if used[i] is not None
                    else slots[i][0][1] for i in range(13)]
            if improper(vals, Q49, BAND49):
                got = True; break
    admit += got
print("x49 refinement of x7 witnesses:")
print("  (a) MC random-refinement survival mod 2107: %d/%d" % (mc_ok, mc_try))
print("  (b) greedy admit (one-sided >=): %d/%d" % (admit, len(wit7)))
print("\nCOMPARISON LADDER (p=43 throughout, all one-sided where greedy):")
print("  x14 gate from level-1:    admit 40/40, random ~0/30000        (S135)")
print("  2nd prime (43&61) at x14: admit 20/20, random 464/8000=5.8%%   (S136)")
print("  x7  gate from level-1:    admit %d/20                          (this run)" % len(wit7))
print("  x49 refinement of x7:     admit %d/%d, random %d/%d=%.2f%%" % (
    admit, len(wit7), mc_ok, mc_try, 100 * mc_ok / max(1, mc_try)))
print("\nREADING: if the x49 refinement admit is high, the second 7-rung is still")
print("soft at p=43 and the cage needs either higher rungs (343p...), larger p, or")
print("the conjunction ACROSS primes of the 7-power layers; if it collapses, the")
print("cage is localized at x49 and death-star's |S| should be measured there first.")
print("DONE.")
