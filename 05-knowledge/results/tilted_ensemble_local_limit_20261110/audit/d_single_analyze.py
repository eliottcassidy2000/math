#!/usr/bin/env python3
"""AUDIT: best coverage of the window [3h,5h) by the extensions of ONE base (d_single_realiser.c output).
usage: d_single_analyze.py DIR   (files DIR/single_M.txt for M = 2..9 and DIR/single_10_r.txt)"""
import sys, glob, os
D = sys.argv[1] if len(sys.argv) > 1 else "/tmp/claude-audit-ensemble"
best = {}; ncls = {}
for fn in sorted(glob.glob(os.path.join(D, "single_*.txt"))):
    base = os.path.basename(fn)[len("single_"):-4]; order = int(base.split("_")[0])
    for l in open(fn):
        h, w, b, c = map(int, l.split())
        ncls[h] = ncls.get(h, 0) + c
        if h not in best or b > best[h][0]: best[h] = (b, order)
orders = sorted({int(os.path.basename(f)[len("single_"):-4].split("_")[0]) for f in glob.glob(os.path.join(D, "single_*.txt"))})
print("base orders read:", orders, "; realised odd h < 1024:", len(best), "; bases with H < 1024:", sum(ncls.values()))
full = sorted(h for h, (b, m) in best.items() if b == h)
print("h for which a SINGLE base fills all h odd values of [3h,5h):", full)
top = sorted(((b / h, h, b, m) for h, (b, m) in best.items() if h >= 9), reverse=True)[:8]
print("largest single-base coverage, h >= 9 (h, covered, base order, fraction):", [(h, b, m, round(f, 3)) for f, h, b, m in top])
big = sorted(((b / h, h, b, m) for h, (b, m) in best.items() if h >= 237), reverse=True)[:6]
print("largest single-base coverage, h >= 237:", [(h, b, m, round(f, 3)) for f, h, b, m in big])
print("=> for every realised odd h in [237, 1023] the best single base covers at most %.1f%% of [3h,5h)" % (100 * big[0][0]))
