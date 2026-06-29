#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LRC14 floor CV criterion -- UNIFORMITY sweep (mac-mini-2026-06-29-S2).

THM-579: covering floor holds if CV(N_R)^2 < m_Q/(1-m_Q).  Since R,Q are
independent parts, the per-r requirement is
    max_{|R|=13-r}  CV(N_R)^2   <   cap_r/(1-cap_r),
cap_r = min_{|Q|=r} meas(lonely(Q))  (THM-576 pairwise-avoidance cap).
This sweeps many 14-free R (consec / shifted / AP / random / large-speed) to
estimate max CV(N_R)^2 and compare to the threshold.  Robust margin => the
clean criterion is plausibly uniform (closing the floor on the C-S branch).
"""
from __future__ import annotations
import functools, math, random
print = functools.partial(print, flush=True)
from fractions import Fraction as F
from itertools import combinations

W = F(1, 14)


def danger_intervals(p):
    ivs = []
    for k in range(p + 1):
        lo = max(F(k, p) - W / p, F(0)); hi = min(F(k, p) + W / p, F(1))
        if hi > lo:
            ivs.append((lo, hi))
    return ivs


def union_intervals(lists):
    ivs = sorted(iv for lst in lists for iv in lst)
    if not ivs:
        return []
    out = [list(ivs[0])]
    for lo, hi in ivs[1:]:
        if lo > out[-1][1]:
            out.append([lo, hi])
        else:
            out[-1][1] = max(out[-1][1], hi)
    return [(a, b) for a, b in out]


def complement(intervals):
    out = []; cur = F(0)
    for lo, hi in intervals:
        if lo > cur:
            out.append((cur, lo))
        cur = max(cur, hi)
    if cur < 1:
        out.append((cur, F(1)))
    return out


def lonely_set(P):
    return complement(union_intervals([danger_intervals(p) for p in P]))


def measure(intervals):
    return sum((hi - lo for lo, hi in intervals), F(0))


def intersect(a, b):
    out = []
    for (lo1, hi1) in a:
        for (lo2, hi2) in b:
            lo = max(lo1, lo2); hi = min(hi1, hi2)
            if hi > lo:
                out.append((lo, hi))
    return out


def shift(intervals, s):
    res = []
    for lo, hi in intervals:
        length = hi - lo
        a = lo + s
        fa = a - F(math.floor(a))
        end = fa + length
        if end <= 1:
            res.append((fa, end))
        else:
            res.append((fa, F(1))); res.append((F(0), end - 1))
    return union_intervals([res])


def CV2_of_R(R):
    Lr = lonely_set(tuple(sorted(R)))
    m_R = measure(Lr)
    if m_R == 0:
        return None
    EN = 14 * m_R
    EN2 = F(0)
    for d in range(-13, 14):
        EN2 += (14 - abs(d)) * measure(intersect(Lr, shift(Lr, F(d, 14))))
    return (EN2 - EN * EN) / (EN * EN)


def cap_r(r):
    """min_{|Q|=r} meas(lonely(Q)) over Q subset {1..16} (THM-576 minimizers)."""
    best = None
    for Q in combinations(range(1, 17), r):
        m = measure(lonely_set(Q))
        if best is None or m < best:
            best = m
    return best


def main():
    print("=" * 78)
    print("LRC14 floor CV uniformity sweep (mac-mini-S2)")
    print("=" * 78)
    rng = random.Random(12345)
    for r in range(2, 7):
        size = 13 - r
        cap = cap_r(r)
        thr = cap / (1 - cap)
        # diverse 14-free R of given size
        Rs = set()
        Rs.add(tuple(range(1, size + 1)))                       # consec
        Rs.add(tuple(range(2, size + 2)))                       # shifted
        for d in (2, 3):                                        # APs (14-free)
            ap = tuple(x for x in range(1, d * size + 1, d) if x % 14 != 0)[:size]
            if len(ap) == size:
                Rs.add(ap)
        Rs.add(tuple(sorted([1] + list(range(13 - size + 2, 14))))[:size])  # {1}+top
        # large-speed variants: replace some elements by larger 14-free speeds
        base = list(range(1, size + 1))
        for _ in range(400):
            cand = rng.sample([x for x in range(1, 70) if x % 14 != 0], size)
            Rs.add(tuple(sorted(cand)))
        maxcv = None; argR = None
        for R in Rs:
            cv2 = CV2_of_R(R)
            if cv2 is None:
                continue
            if maxcv is None or cv2 > maxcv:
                maxcv = cv2; argR = R
        ok = float(maxcv) < float(thr)
        print(f"\nr={r} (|R|={size}): cap_r={cap} ({float(cap):.4f}) "
              f"threshold cap/(1-cap)={float(thr):.4f}")
        print(f"   max CV(N_R)^2 over {len(Rs)} configs = {float(maxcv):.4f}  at R={argR}")
        print(f"   criterion max CV^2 < threshold : {'HOLDS' if ok else 'FAILS'} "
              f"(margin {float(thr)-float(maxcv):+.4f})")
    print("\n" + "=" * 78)
    print("If 'HOLDS' with margin for all r, the C-S floor criterion is uniform over")
    print("the sampled (bounded-speed) family -> closes the floor on the C-S branch")
    print("modulo (a) unbounded-R discharge and (b) the covering constraint on R.")
    print("=" * 78)


if __name__ == "__main__":
    main()
