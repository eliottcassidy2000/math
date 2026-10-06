#!/usr/bin/env python3
"""Check Theorem 2.4: first time tau(a,b) at which both runners are >= 1/3 from the start."""
from fractions import Fraction as Fr
from math import gcd

def tau(speeds, n_plus_1):
    d = Fr(1, n_plus_1)
    # candidate first-entry times: left endpoints (m + d)/v of safe intervals, sorted
    cands = sorted({(m + d) / v for v in speeds for m in range(0, v)})
    for t in cands:
        if t <= 0 or t >= 1:
            continue
        ok = True
        for v in speeds:
            x = (t * v) % 1
            if not (d <= x <= 1 - d):
                ok = False; break
        if ok:
            return t
    return None

top = []
viol = 0
for a in range(1, 121):
    for b in range(a + 1, 121):
        if gcd(a, b) != 1:
            continue
        t = tau((a, b), 3)
        if a == 1:
            pred = Fr(1, 3) if b % 3 else Fr(1, 3) + Fr(1, 9 * (b // 3))
            viol += t != pred
        elif b < 2 * a:
            viol += t != Fr(1, 3 * a)
        else:
            viol += t > Fr(2, 3 * a)
        top.append((t, a, b))
top.sort(reverse=True)
print("coprime pairs a<b<=120 checked; violations of Theorem 2.4:", viol)
print("top 8 tau:", [(str(t), a, b) for t, a, b in top[:8]])
print("max over a>=2:", max((t, a, b) for t, a, b in top if a >= 2))
