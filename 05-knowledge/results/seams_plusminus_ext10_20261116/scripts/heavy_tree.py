#!/usr/bin/env python3
"""heavy_tree.py (S-10): the heaviest-parent tree — each class of order n attached to its distinct deletion class of
largest count (ties: least gentourng index, which is the least canonical string order used before).  Uses the tooth
packet's parents_n.txt (distinct parent indices per class) and counts_z.txt.  Reports branching (max, childless, mean)
and the share of first creators whose tree parent is a first creator, orders 4-9, and compares with the least-string tree."""
import os
from collections import defaultdict, Counter
T = '/private/tmp/claude-501/-Users-e-Documents-GitHub-ephrepos-math/d63fb198-a11f-4364-a0ae-e072089f7b73/scratchpad/work/tooth/'
H = defaultdict(list)
for line in open(T + 'counts_z.txt'):
    n, h = map(int, line.split()); H[n].append(h)
seen = set(); first = {}
for n in range(1, 10):
    new = set(H[n]) - seen; first[n] = [h in new for h in H[n]]; seen |= set(H[n])
mono_bad = 0; mono_tot = 0
for n in range(4, 10):
    parents = [list(map(int, l.split())) for l in open(T + f'parents_{n}.txt')]
    for i, ps in enumerate(parents):
        mono_tot += 1
        if max(H[n-1][p] for p in ps) > H[n][i]: mono_bad += 1
print(f"deletion monotonicity H(T - v) <= H(T): {mono_tot} classes of orders 4-9, {mono_bad} violations")
print("order: heaviest-parent tree (max branching, childless parents, first creators with a first-creator parent) | least-string tree (same)")
for n in range(4, 10):
    parents = [list(map(int, l.split())) for l in open(T + f'parents_{n}.txt')]
    cpar = [int(x) for x in open(T + f'cparent_{n}.txt')]
    heavy = [max(ps, key=lambda p: (H[n-1][p], -p)) for ps in parents]
    def stats(par):
        br = Counter(par); mx = max(br.values()); childless = len(H[n-1]) - len(br)
        fc = [i for i in range(len(H[n])) if first[n][i]]; bad = [i for i in fc if not first[n-1][par[i]]]
        return mx, childless, f"{len(fc) - len(bad)}/{len(fc)}" + (f" (exceptions: counts {sorted(set(H[n][i] for i in bad))})" if bad else "")
    a = stats(heavy); b = stats(cpar)
    print(f"  {n}: ({a[0]}, {a[1]}, {a[2]}) | ({b[0]}, {b[1]}, {b[2]})")
