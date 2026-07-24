#!/usr/bin/env python3
"""
lrc14_margin_vs_bindingspeed_opus_S4.py   opus-2026-07-23-S4

DECIDES whether the certified-Fejer concentration reduction is PRACTICAL or ASTRONOMICAL.
My Fejer law: certification degree N* ~ (max binding speed at tau*) / delta,  delta = gap - 1/14.
So a uniform bulk degree exists iff  small delta  FORCES  small max binding speed.
Open after my S4 self-correction. Here: scan 13-speed configs, compute EXACT gap, delta, and the
binding set {v : ||v tau*|| = gap}; ask whether large binding speeds coexist with small delta.
"""
from fractions import Fraction as Fr
from math import floor

def ndist(x):
    f = x - floor(x); return min(f, 1 - f)
def gval(V, t): return min(ndist(v * t) for v in V)
def breakpoints(V):
    pts = {Fr(0), Fr(1)}
    for v in V:
        for j in range(0, 2 * v + 1): pts.add(Fr(j, 2 * v))
    for i in range(len(V)):
        for j in range(i + 1, len(V)):
            for d in {abs(V[i] - V[j]), V[i] + V[j]}:
                if d == 0: continue
                for k in range(0, d + 1): pts.add(Fr(k, d))
    return sorted(p for p in pts if 0 <= p <= 1)
def analyze(V):
    bp = breakpoints(V); best = None
    for p in bp:
        gv = gval(V, p)
        if best is None or gv > best[0]: best = (gv, p)
    gap, ts = best
    binding = [v for v in V if ndist(v * ts) == gap]
    return gap, ts, binding

FL = Fr(1, 14)
rows = []
fams = []
# single-far from AP: drop j, add r
for j in [6, 11, 12, 13]:
    for r in range(14, 61):
        V = sorted(set(range(1, 14)) - {j} | {r})
        if len(V) == 13: fams.append((f"AP\{{{j}}}u{{{r}}}", V))
# {1..12, r} and {1..11,13, r}
for r in range(14, 61):
    fams.append((f"{{1..12,{r}}}", list(range(1, 13)) + [r]))
    fams.append((f"{{1..11,13,{r}}}", list(range(1, 12)) + [13, r]))

for nm, V in fams:
    gap, ts, binding = analyze(V)
    d = gap - FL
    if d <= 0: continue                       # tight or below (shouldn't happen except AP/GW)
    rows.append((float(d), d, gap, nm, max(binding), binding, ts))

rows.sort(key=lambda x: x[0])
print("13-speed configs sorted by margin delta = gap - 1/14  (does small delta force small binding speed?)")
print("=" * 104)
print(f"  {'delta':>10} {'gap':>8}  {'config':22s} {'maxBind':>7}  {'binding set':22s} {'tau*':>10}")
print("-" * 104)
for (df, d, gap, nm, mb, binding, ts) in rows[:28]:
    print(f"  {str(d):>10} {str(gap):>8}  {nm:22s} {mb:7d}  {str(binding):22s} {str(ts):>10}")
print("-" * 104)
# the decisive statistic: max binding speed among the SMALLEST-margin configs
import statistics
small = [r for r in rows if r[0] < 0.004]      # gap < 1/14 + 0.004  (i.e. < ~0.0754)
big   = [r for r in rows if r[0] >= 0.004]
print(f"configs with delta < 0.004 : n={len(small)}, max binding speed = "
      f"{max((r[4] for r in small), default=None)}, median = {statistics.median([r[4] for r in small]) if small else None}")
print(f"configs with delta >= 0.004: n={len(big)}, max binding speed = "
      f"{max((r[4] for r in big), default=None)}")
# is there ANY config with small delta AND large binding speed?
hard = [r for r in rows if r[0] < 0.004 and r[4] > 20]
print(f"\nHARD configs (delta<0.004 AND maxBind>20): {len(hard)}")
for r in hard[:12]:
    print(f"   delta={str(r[1]):>10} gap={str(r[2]):>7} {r[3]:22s} maxBind={r[4]} binding={r[5]}")
print("\nIf HARD is EMPTY: small margin forces small binding speed => uniform bulk degree ~13/eps0 (PRACTICAL).")
print("If HARD is nonempty: N* ~ maxBind/delta is unbounded => reduction stays astronomical.")
