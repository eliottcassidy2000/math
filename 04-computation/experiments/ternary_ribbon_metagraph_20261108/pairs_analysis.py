#!/usr/bin/env python3
"""pairs_analysis.py n -- statistics of (H closed end, H open end) over the lines of order n from lines_n{n}_p*.txt."""
import sys, glob
from collections import Counter
n = int(sys.argv[1])
pairs = Counter(); hc = Counter(); ho = Counter(); ksec = Counter(); meta = []
for path in sorted(glob.glob(f"lines_n{n}_p*.txt")):
    sec = None
    for line in open(path):
        line = line.strip()
        if line.startswith('#pairs'): sec = 'pairs'; continue
        if line.startswith('#Hclosed'): sec = 'hc'; continue
        if line.startswith('#Hopen'): sec = 'ho'; continue
        if line.startswith('#k'): sec = 'k'; continue
        if line.startswith('#'): sec = None; meta.append(line) if line.startswith('# n=') else None; continue
        if not line: continue
        p = line.split()
        if sec == 'pairs': pairs[(int(p[0]), int(p[1]), int(p[2]), int(p[3]))] += int(p[4])
        elif sec == 'hc': hc[int(p[0])] += int(p[1])
        elif sec == 'ho': ho[int(p[0])] += int(p[1])
        elif sec == 'k': ksec[int(p[0])] += int(p[1])
total = sum(pairs.values())
print(f"n={n}: lines {total}; k-distribution of the open end {dict(sorted(ksec.items()))}")
same = sum(c for (a, b, bl, st), c in pairs.items() if a == b)
lt = sum(c for (a, b, bl, st), c in pairs.items() if a < b); gt = sum(c for (a, b, bl, st), c in pairs.items() if a > b)
print(f"  same-value lines {same} ({same/total:.4f}); closed end smaller {lt} ({lt/total:.4f}), larger {gt} ({gt/total:.4f})")
import math
EHc = sum(a * c for (a, b, bl, st), c in pairs.items()) / total; EHo = sum(b * c for (a, b, bl, st), c in pairs.items()) / total
print(f"  mean H closed {EHc:.2f}, mean H open {EHo:.2f}; max H closed {max(hc)}, max H open {max(ho)}")
cov = sum((a - EHc) * (b - EHo) * c for (a, b, bl, st), c in pairs.items()) / total
va = sum((a - EHc) ** 2 * c for (a, b, bl, st), c in pairs.items()) / total; vb = sum((b - EHo) ** 2 * c for (a, b, bl, st), c in pairs.items()) / total
print(f"  correlation of the two end values {cov / math.sqrt(va * vb):.4f}")
# log-ratio distribution
ratios = Counter()
for (a, b, bl, st), c in pairs.items(): ratios[round(math.log(b / a), 1)] += c
print("  distribution of log(H_open/H_closed) rounded to 0.1 (top 12):", sorted(ratios.items(), key=lambda kv: -kv[1])[:12])
# extremes
mn = min((min(a, b), a, b) for (a, b, bl, st) in pairs); mx = max((max(a, b), a, b) for (a, b, bl, st) in pairs)
print(f"  smallest value on any line {mn}; largest {mx}")
print("  lines whose open end has H in {1,3,5,7,9}:", {h: ho.get(h, 0) for h in (1, 3, 5, 7, 9)})
print("  lines whose closed end has H in {3,5,7,9,11,15}:", {h: hc.get(h, 0) for h in (3, 5, 7, 9, 11, 15)})
# blue lines: value pairs
bl_pairs = Counter()
for (a, b, bl, st), c in pairs.items():
    if bl: bl_pairs[(a, b)] += c
print(f"  blue lines: {sum(bl_pairs.values())}; distinct value pairs {len(bl_pairs)}; same-value blue lines {sum(c for (a,b),c in bl_pairs.items() if a==b)}")
