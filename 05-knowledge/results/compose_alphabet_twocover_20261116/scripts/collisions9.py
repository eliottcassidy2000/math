#!/usr/bin/env python3
"""collisions9.py (C-8): (h, P(2)) collision groups among the classes of order 9 and how many are split by P(3)
(covers9.txt from covers9.c)."""
from collections import defaultdict
g = defaultdict(set); cnt = defaultdict(int)
for line in open('covers9.txt'):
    h, p2, p3 = map(int, line.split()); g[(h, p2)].add(p3); cnt[(h, p2)] += 1
coll = [k for k, c in cnt.items() if c >= 2]; split = [k for k in coll if len(g[k]) >= 2]
print(f"(C-8) order 9: classes {sum(cnt.values())}, distinct (h, P(2)) {len(cnt)}, collision groups {len(coll)}, split by P(3): {len(split)}")
print("      split groups (h, P(2)): P(3) values:", {k: sorted(g[k]) for k in sorted(split)[:12]})
diffs = sorted(set(abs(a - b) for k in split for a in g[k] for b in g[k] if a < b))
print("      P(3) differences within split groups:", diffs[:20], "; all multiples of 40:", all(d % 40 == 0 for d in diffs))
