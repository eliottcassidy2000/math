#!/usr/bin/env python3
"""
lrc14_c13_lane_scan_boxeph_S22.py — the c=13 dilated-core multi-killer lane (HYP-6295)

boxeph-2026-07-12-S22.  The seam between mac-mini-S68's cases (a) and (b):
families  S_b = 13*{1..12} u {14b},  13 ∤ b.
 - covering for EVERY b (each d in 2..12 divides 13d; 13-duty on the core; 14-duty on 14b);
 - primitive (gcd(13-core, 14b) = 1 when 13 ∤ b);
 - NO 182-multiple  =>  genuinely multi-killer, outside the single-killer censuses;
 - Vmax = max(156, 14b) > 182 for b >= 14  =>  outside klein's ILP for b >= 14.

The adversary tunes b:  14b ≡ -1 (mod 169) blocks the 13^2-base (b ≡ 12 mod 169);
14b ≡ 0 (mod 181) blocks the {1..12}-landscape escape at 181 (b ≡ 181);
b = 170 makes the killer 2380 = 13*183+1 resonant with the 1/13-window.

QUESTION: min_b M(S_b) — does the c=13 lane dip below the covering-min 14/183?
"""

from fractions import Fraction as F
from math import gcd
from functools import reduce
import sys

FLOOR = F(14, 183)

def exact_M(speeds):
    n = len(speeds)
    assert len(set(speeds)) == n and all(v > 0 for v in speeds)
    if all(v % 2 for v in speeds):
        return F(1, 2), 2, 1
    best_n, best_d = 0, 1
    best_wit = (0, 0)
    qs = sorted({speeds[i] + speeds[j] for i in range(n)
                 for j in range(i + 1, n)}, reverse=True)
    for q in qs:
        rs = sorted({u % q for u in speeds})
        half = q // 2
        for m in range(1, q):
            d = q
            ok = True
            for r in rs:
                x = (r * m) % q
                if x > half:
                    x = q - x
                if x < d:
                    d = x
                    if d * best_d <= best_n * q:
                        ok = False
                        break
            if ok and d * best_d > best_n * q:
                best_n, best_d = d, q
                best_wit = (q, m)
    return F(best_n, best_d), best_wit[0], best_wit[1]

def is_covering(s):
    return all(any(v % d == 0 for v in s) for d in range(2, 15))

def main():
    bmax = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    core = [13 * j for j in range(1, 13)]
    print(f"c=13 LANE SCAN: S_b = 13*{{1..12}} u {{14b}}, b = 1..{bmax}, 13∤b")
    print(f"floor 14/183 = {float(FLOOR):.6f}; ILP boundary Vmax<=182 <=> b<=13\n")
    rows = []
    for b in range(1, bmax + 1):
        if b % 13 == 0:
            continue
        w = 14 * b
        if w in core:
            continue
        fam = sorted(core + [w])
        assert is_covering(fam) and reduce(gcd, fam) == 1 \
            and not any(v % 182 == 0 for v in fam)
        M, q, m = exact_M(fam)
        rows.append((M, b, q, m))
    rows_by_M = sorted(rows)
    print("15 LOWEST M over the lane:")
    for (M, b, q, m) in rows_by_M[:15]:
        flag = "BELOW-FLOOR!!" if M < FLOOR else ("=FLOOR" if M == FLOOR else "")
        ilp = "ILP" if 14 * b <= 182 else "OPEN"
        print(f"  b={b:4d} (14b={14*b:5d}, {ilp:4s})  M = {str(M):>10s} = "
              f"{float(M):.6f}  wit q={q} m={m}  {flag}")
    below = [r for r in rows_by_M if r[0] < FLOOR]
    print(f"\nfamilies below 14/183: {len(below)}   (of {len(rows)})")
    if below:
        print("*** THE c=13 LANE DIPS BELOW THE DEEP WELL — list: ***")
        for (M, b, q, m) in below:
            print(f"  b={b}: M = {M} = {float(M):.6f} wit q={q} m={m}")
    mn = rows_by_M[0]
    print(f"\nlane min: M = {mn[0]} at b = {mn[1]} (wit q={mn[2]}, m={mn[3]});"
          f" margin over floor {float(mn[0]-FLOOR):+.6f}")
    # spotlight designed b values
    print("\ndesigned-b spotlight:")
    for b in (12, 170, 181, 350, 519):
        if b % 13 == 0 or b > bmax:
            continue
        r = next((r for r in rows if r[1] == b), None)
        if r:
            print(f"  b={b:4d}: M = {str(r[0]):>10s} = {float(r[0]):.6f}  "
                  f"wit q={r[2]} m={r[3]}")

if __name__ == '__main__':
    main()
