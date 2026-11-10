#!/usr/bin/env python3
"""AUDIT check (c): the class-count function N_n(h) from the auditor's own histogram (d_window.c output).
usage: c_hist_analyze.py ORDER HISTFILE [HISTFILE ...]   (several files are summed: res/mod parts)
Compares with the packet's out/hist_ORDER.txt when present."""
import sys, os, math
from collections import defaultdict
order = int(sys.argv[1]); N = defaultdict(int)
for fn in sys.argv[2:]:
    for l in open(fn):
        h, c = map(int, l.split()); N[h] += c
tot = sum(N.values()); mean = sum(h * c for h, c in N.items()) / tot
sd = math.sqrt(sum(h * h * c for h, c in N.items()) / tot - mean ** 2)
print(f"order {order}: classes {tot}; distinct H values {len(N)}; max H {max(N)}; all H odd: {all(h % 2 for h in N)}")
print(f"   class-mean {mean:.3f}, class-sd {sd:.3f}  (labelled mean n!/2^(n-1) = {math.factorial(order) / 2 ** (order - 1):.3f})")
missing = [h for h in range(1, max(N) + 1, 2) if h not in N]
first_missing_other = next(h for h in missing if h not in (7, 21))
print(f"   odd values not occurring at this order, below the first other gap: {[h for h in missing if h < first_missing_other]}; first missing odd value other than 7, 21: {first_missing_other}"
      f"  => every odd h <= {first_missing_other - 2} except 7, 21 occurs")
lo, hi = mean - sd, mean + sd
win = [h for h in range(1, max(N) + 1, 2) if lo <= h <= hi]
least = min(N.get(h, 0) for h in win); arg = [h for h in win if N.get(h, 0) == least]
print(f"   window mean +- sd = [{lo:.2f}, {hi:.2f}] -> odd h from {win[0]} to {win[-1]} ({len(win)} values); least class count {least} at h = {arg}")
zs = []
for h in win:
    nb = (N.get(h - 2, 0) + N.get(h + 2, 0)) / 2
    if nb > 0: zs.append((N.get(h, 0) - nb) / math.sqrt(1.5 * nb))
rms = math.sqrt(sum(z * z for z in zs) / len(zs))
print(f"   rms normalised second difference {rms:.3f}  (variance ratio to independent Poisson {rms * rms:.2f})")
worst = (0, None, None)
for q in (3, 4, 5, 7, 8):
    s = defaultdict(float); k = defaultdict(int); num = defaultdict(float); den = defaultdict(float)
    for h in win:
        nb = (N.get(h - 2, 0) + N.get(h + 2, 0) + N.get(h - 4, 0) + N.get(h + 4, 0)) / 4
        if nb > 0:
            s[h % q] += N.get(h, 0) / nb; k[h % q] += 1; num[h % q] += N.get(h, 0); den[h % q] += nb
    rel = {r: s[r] / k[r] for r in sorted(s)}; rel2 = {r: num[r] / den[r] for r in sorted(s)}
    dev = max(abs(v - 1) for v in rel.values()); r_ = max(rel, key=lambda r: abs(rel[r] - 1))
    if dev > worst[0]: worst = (dev, q, r_)
    print(f"   h mod {q}: mean of local ratios {({r: round(v, 3) for r, v in rel.items()})}; ratio of sums {({r: round(v, 3) for r, v in rel2.items()})}")
print(f"   largest residue deviation {100 * worst[0]:.1f}% at h = {worst[2]} mod {worst[1]}")
here = os.path.dirname(os.path.abspath(__file__)); pk = os.path.join(here, "..", "out", f"hist_{order}.txt")
if os.path.exists(pk):
    P = defaultdict(int)
    for l in open(pk):
        h, c = map(int, l.split()); P[h] += c
    print(f"   identical to the packet's out/hist_{order}.txt: {dict(P) == dict(N)}")
