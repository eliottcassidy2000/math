#!/usr/bin/env python3
"""klein-2026-07-02-S106 -- HYP-4011: F3-SHARP numeric verification.
CLAIM: meas(L_C) >= int_{L_B} D_c(t) dt - eps(N), eps(N) <= 2 c_B (j+1) / N -- no sqrt.
MECHANISM: per component I of L_B, substitute u = Nt; the far-runner joint-safe event in
u sweeps N|I| cycles of the (drifting-offset) tooth pattern; int_I D_c(t) dt integrates
the TRUE drifting density, so per-cycle the two agree EXACTLY (rate_core exactness);
the error is only the two partial cells per component + tooth-boundary cells:
<= 2(j+1)/N per component. Test: deep-cluster cores C = B u {N + c_i} across N,
exact interval arithmetic, fit eps(N) vs the bound and vs sqrt."""
from fractions import Fraction as F
import itertools

R = F(1, 14)

def merged_danger(S):
    ivs = []
    for v in S:
        rv = R / v
        for a in range(v + 1):
            lo, hi = F(a, v) - rv, F(a, v) + rv
            lo, hi = max(lo, F(0)), min(hi, F(1))
            if hi > lo: ivs.append((lo, hi))
    ivs.sort()
    out = []
    for lo, hi in ivs:
        if out and lo <= out[-1][1]: out[-1] = (out[-1][0], max(out[-1][1], hi))
        else: out.append((lo, hi))
    return out

def lonely_measure(S):
    cov = merged_danger(S)
    tot, prev = F(0), F(0)
    for lo, hi in cov:
        if lo > prev: tot += lo - prev
        prev = max(prev, hi)
    if prev < 1: tot += 1 - prev
    return tot

def lonely_components(S):
    cov = merged_danger(S); out = []; prev = F(0)
    for lo, hi in cov:
        if lo > prev: out.append((prev, lo))
        prev = max(prev, hi)
    if prev < 1: out.append((prev, F(1)))
    return out

def D_integral(B, offs, N, sub=6):
    """int_{L_B} D_c(t) dt with D_c(t) = joint safe density of the offset pattern at
    slow time t: D_c(t) = meas_u{all ||u + c_i t|| >= r} (exact per t; integrate by
    fine Simpson over each component -- D_c piecewise linear in t so midpoint suffices
    on small subdivisions)."""
    def D_at(t):
        ivs = []
        for c in offs:
            ctr = (c * t) % 1
            ivs.append(((ctr - R) % 1, (ctr + R) % 1))
        # measure of union on circle of j arcs of width 2r centered ctr
        pts = []
        for a, b in ivs:
            if a <= b: pts.append((a, b))
            else: pts.extend([(F(0), b), (a, F(1))])
        pts.sort(); tot = F(0); cl = ch = None
        for lo, hi in pts:
            if ch is None: cl, ch = lo, hi
            elif lo <= ch: ch = max(ch, hi)
            else: tot += ch - cl; cl, ch = lo, hi
        if ch is not None: tot += ch - cl
        return 1 - tot
    tot = F(0)
    for lo, hi in lonely_components(B):
        n = sub
        h = (hi - lo) / n
        for i in range(n):
            tot += h * D_at(lo + h * (2 * i + 1) / 2)
    return tot

B = [1, 2, 3, 4]
offs = [0, 1, 2]   # far cluster N, N+1, N+2 ; j = 3
j = len(offs)
cB = len(lonely_components(B))
print(f"B = {B} (c_B = {cB} components), offsets {offs} (j = {j}); bound 2 c_B (j+1)/N")
print(f"{'N':>6} {'meas(L_C)':>12} {'int D_c':>12} {'eps':>12} {'2cB(j+1)/N':>12} {'4(j+1)sqrt(D/N)':>16} ok")
for N in (40, 80, 120, 240, 480):
    C = B + [N + c for c in offs]
    m = lonely_measure(C)
    dint = D_integral(B, offs, N, sub=8)
    eps = abs(float(dint - m))
    bnd = 2 * cB * (j + 1) / N
    old = 4 * (j + 1) * (max(offs) / N) ** 0.5 if max(offs) else 0
    print(f"{N:>6} {float(m):>12.6f} {float(dint):>12.6f} {eps:>12.2e} {bnd:>12.2e} {old:>16.2e} {eps <= bnd}")
print("PASS iff eps <= the linear bound at every N (and visibly O(1/N) scaling).")
