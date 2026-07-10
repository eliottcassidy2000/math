#!/usr/bin/env python3
"""
Connected-cascade closure forensics (klein-2026-07-09-S233, companion to
lrc14_connected_cascade_klein_S233.py).

(1) Is M2(u1*u2) = A2(u1/u2) an identity (all 78 pairs, q=139/353)?  This is
    the t=2 agreement that let THM-684(I)'s product-count misidentification
    pass S232's checks.
(2) EXACT full-layer closure: LM(q)/Q = (b/Q)^13 + L2sig + L3sig + R_{>=4},
    LM by 13-way bitmask AND; R_{>=4} = the exact total mass of layers 4..13.
    Decay of R_{>=4} in q = does the cascade truncation converge?
(3) Top-|dev3| triples at q=5003: are GEN's exceptional connected triples its
    ADDITIVE relations (Schur a+b=c / AP 2b=a+c), i.e. additive->multiplicative
    leakage into the cumulant layer?  DIL's should be ratio-coherent.
"""
from itertools import combinations
from math import gcd

GEN = [12, 33, 46, 47, 68, 73, 79, 81, 85, 87, 91, 112, 120]
DIL = [20, 41, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260]


def band(q):
    return (q + 13) // 14, (13 * q) // 14


def masks(S, q):
    lo, hi = band(q)
    out = []
    for v in [x % q for x in S]:
        m = 0
        for c in range(1, q):
            if lo <= c * v % q <= hi:
                m |= 1 << c
        out.append(m)
    return out, hi - lo + 1


def part1(q):
    lo, hi = band(q)
    B = list(range(lo, hi + 1))
    Bs = set(B)
    v = [x % q for x in GEN]
    bm, b = masks(GEN, q)
    neq = 0
    for i, j in combinations(range(13), 2):
        a2 = (bm[i] & bm[j]).bit_count()
        T = v[i] * v[j] % q
        m2 = sum(1 for y in B if T * pow(y, -1, q) % q in Bs)
        if a2 != m2:
            neq += 1
    print(f"(1) q={q}: M2(u1*u2) vs A2(u1/u2) over 78 pairs: "
          f"{78 - neq}/78 equal" + ("  -> IDENTITY (t=2 masks the error)"
                                    if neq == 0 else f"  ({neq} differ)"))


def closure(S, name, primes):
    print(f"(2) {name}: exact layer closure  LM/Q = bud + L2sig + L3sig + R>=4")
    for q in primes:
        bm, b = masks(S, q)
        Q = q - 1
        r = b / Q
        A2 = {}
        for i, j in combinations(range(13), 2):
            A2[i, j] = (bm[i] & bm[j]).bit_count()
        l2 = r ** 11 / Q * sum((Q * a - b * b) / Q for a in A2.values())
        l3 = 0.0
        for i, j, k in combinations(range(13), 3):
            a3 = (bm[i] & bm[j] & bm[k]).bit_count()
            p3 = Q * Q * a3 - b * Q * (A2[i, j] + A2[i, k] + A2[j, k]) + 2 * b ** 3
            l3 += p3 / Q ** 3
        l3 *= r ** 10
        lm = bm[0]
        for m in bm[1:]:
            lm &= m
        lmQ = lm.bit_count() / Q
        rest = lmQ - r ** 13 - l2 - l3
        print(f"    q={q:5d}: LM/Q={lmQ:.4f} bud={r ** 13:.4f} "
              f"L2sig={l2:+.4f} L3sig={l3:+.4f} R>=4={rest:+.4f}")


def classify(a, c, d):
    if a + c == d:
        return "SCHUR"
    if 2 * c == a + d:
        return "AP"
    ok = True
    for x, y in ((a, c), (a, d), (c, d)):
        g = gcd(x, y)
        if y // g > 6 or x // g > 6:
            ok = False
    return "RATIO" if ok else "-"


def part3(S, name, q=5003):
    bm, b = masks(S, q)
    Q = q - 1
    A2 = {}
    for i, j in combinations(range(13), 2):
        A2[i, j] = (bm[i] & bm[j]).bit_count()
    devs = []
    for i, j, k in combinations(range(13), 3):
        a3 = (bm[i] & bm[j] & bm[k]).bit_count()
        p3 = Q * Q * a3 - b * Q * (A2[i, j] + A2[i, k] + A2[j, k]) + 2 * b ** 3
        devs.append((abs(p3 / Q / Q), p3 / Q / Q, (i, j, k)))
    devs.sort(reverse=True)
    print(f"(3) {name} q={q}: top-6 |dev3| triples (speeds, class):")
    for _, d, (i, j, k) in devs[:6]:
        a, c, e = S[i], S[j], S[k]
        print(f"    ({a},{c},{e}) dev3={d:+9.2f}  {classify(a, c, e)}")
    nrel = [(d, t) for _, d, t in devs
            if classify(S[t[0]], S[t[1]], S[t[2]]) in ("SCHUR", "AP", "RATIO")]
    top20 = devs[:20]
    hits = sum(1 for _, d, t in top20
               if classify(S[t[0]], S[t[1]], S[t[2]]) != "-")
    print(f"    relation-carrying triples in top-20: {hits}/20 "
          f"(instance has {len(nrel)} relation triples of 286)")


part1(139)
part1(353)
closure(GEN, "GEN", [139, 353, 1009, 2003, 5003])
closure(DIL, "DIL", [139, 353, 1009, 2003, 5003])
part3(GEN, "GEN")
part3(DIL, "DIL")
