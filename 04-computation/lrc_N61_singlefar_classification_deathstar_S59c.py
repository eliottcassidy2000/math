#!/usr/bin/env python3
"""
death-star-2026-07-19-S59c -- HYP-7900 part C: the COMPLETE single-far
classification at N=61 + budgeted census.

Window W_61 = (1/62, 2/123).  F_4(61) = {1..59,61,240} = 4/247 CONFIRMED inline.
This decides the WHOLE single-far stratum {1..61}\\{i} u {x} (all i, all x):
per-defect sweep-first certificates (safe-interval computation of the base at
theta = 2/123 directly certifies M(B) > theta when an interior interval exists
-- no full M(B) needed, no open-LRC citation), absorption threshold X0, finite
check to max(X0, 265)+15 (range covers x = 240).  Expected: the classification
finds (i=60, x=240) -> 4/247 and nothing else, UNLESS the gate family has
siblings (e.g. other defects) -- that is what this decides.

Plus a small census (bordered/two-defect/targeted/repair species) at N=61.
"""
from fractions import Fraction as F
from math import ceil, gcd
from functools import reduce
from itertools import combinations
import random, sys, time
sys.path.insert(0, '04-computation')
from lrc_singlefar_absorption_atlas_deathstar_S59 import (
    M_exact, M_exact_wit, safe_intervals)

random.seed(6161)
log = lambda s="": print(s, flush=True)
N = 61
TH = F(2, 123); LO = F(1, 62)

def primitive(S):
    g = reduce(gcd, S)
    return tuple(sorted(v//g for v in S))

def part_classify():
    log(f"== N=61 single-far classification (theta = 2/123) ==")
    members = []
    for i in range(1, N+1):
        t0 = time.time()
        B = [v for v in range(1, N+1) if v != i]
        iv, l = safe_intervals(B, TH)
        if l == 0:
            log(f"   i={i:>2}: NO interior safe interval -- absorption inapplicable; "
                f"stratum undecided (and M(B) <= theta: LRC(62) alarm)")
            continue
        X0 = ceil(2*TH / l)
        hi_check = max(X0, 265) + 15
        found = 0
        for x in range(N+1, hi_check+1):
            m = M_exact(B + [x], stop_above=TH)
            if LO < m < TH:
                Mw, q, a = M_exact_wit(B + [x])
                members.append((i, x, Mw, q, a))
                found += 1
                log(f"   i={i:>2}: x={x} MEMBER M={Mw} (q={q}, a={a})")
        log(f"   i={i:>2}: l={float(l):.6f} X0={X0:>4} checked<= {hi_check} "
            f"members={found} ({time.time()-t0:.0f}s)")
    log(f"\n   COMPLETE N=61 single-far members: "
        f"{[(i,x,str(M)) for i,x,M,q,a in members]}")
    return members

def check_fam(S, found, sp):
    Sp = primitive(S)
    if len(set(Sp)) < N: return 0
    m = M_exact(list(Sp), stop_above=TH)
    if LO < m < TH and m not in found:
        found[m] = (Sp, sp)
        log(f"   !! N=61 census MEMBER {m} via {sp}: {Sp}")
    return 1

def part_census():
    log(f"\n== N=61 budgeted census (non-single-far species) ==")
    found = {}; counts = {}
    # B bordered dilated APs
    t0 = time.time(); n = 0
    for d in range(2, 8):
        for a0 in range(1, 40):
            for m in range(45, N):
                spine = [a0 + d*k for k in range(m)]
                nb = N - m
                if nb <= 0: continue
                borders = set()
                for sv in spine:
                    for e in (1, 2, 3):
                        for sg in (1, -1):
                            b = sv + sg*e
                            if b > 0 and b not in spine: borders.add(b)
                borders = sorted(borders)
                if len(borders) < nb: continue
                c = 0
                for combo in combinations(borders, nb):
                    n += check_fam(tuple(sorted(spine+list(combo))), found, "B")
                    c += 1
                    if c >= 1500 or time.time()-t0 > 120: break
                if time.time()-t0 > 120: break
            if time.time()-t0 > 120: break
        if time.time()-t0 > 120: break
    counts["B"] = n
    # C two-defect two-far (random pairs)
    t0 = time.time(); n = 0
    base = list(range(1, N+1))
    pairs = list(combinations(range(1, N+1), 2)); random.shuffle(pairs)
    for i, j in pairs:
        Bb = [v for v in base if v not in (i, j)]
        for w1 in range(N+1, 6*N):
            for w2 in range(w1+1, 6*N):
                n += check_fam(tuple(Bb+[w1, w2]), found, "C")
                if time.time()-t0 > 120: break
            if time.time()-t0 > 120: break
        if time.time()-t0 > 120: break
    counts["C"] = n
    # T targeted structured multiples (incl. the D-graded shapes at 2 defects)
    t0 = time.time(); n = 0
    for i, j in combinations(range(1, N+1), 2):
        Bb = [v for v in base if v not in (i, j)]
        cands = {3*i, 4*i, 3*j, 4*j, 4*i+1, 4*i-1, 4*j+1, 4*j-1, 5*i, 5*j}
        for w1, w2 in combinations(sorted(c for c in cands if c > N), 2):
            n += check_fam(tuple(Bb+[w1, w2]), found, "T")
            if time.time()-t0 > 60: break
        if time.time()-t0 > 60: break
    counts["T"] = n
    # E needle repair at the rung bands (4,247) and (5,308/309)
    t0 = time.time(); n = 0
    vq = [(4, 247), (5, 308), (5, 309), (6, 370)]
    while time.time()-t0 < 90:
        val, q = random.choice(vq)
        if gcd(val, q) != 1: continue
        band = list(range(val, q-val+1))
        res = {val, q-val}
        while len(res) < N: res.add(random.choice(band))
        S = sorted(res); curM = M_exact(S)
        for _ in range(100):
            if LO < curM < TH: break
            idx = random.randrange(len(S)); new = random.choice(band)
            T2 = sorted(set(S[:idx]+S[idx+1:]+[new]))
            if len(T2) < N: continue
            m2 = M_exact(T2, stop_above=TH)
            if m2 < curM: S, curM = T2, m2
        n += check_fam(tuple(primitive(S)), found, "E")
    counts["E"] = n
    log(f"   census counts: {counts}; members beyond classification: "
        f"{sorted(found) if found else 'NONE'}")
    return found, counts

if __name__ == "__main__":
    ms = part_classify()
    found, counts = part_census()
    log("\n== N=61 VERDICT ==")
    log(f"  single-far members (COMPLETE): {[(i,x,str(M)) for i,x,M,q,a in ms]}")
    log(f"  census extras: {sorted(found) if found else 'none'} {counts}")
