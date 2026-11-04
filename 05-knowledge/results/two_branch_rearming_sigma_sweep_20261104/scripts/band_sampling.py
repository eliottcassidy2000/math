#!/usr/bin/env python3
"""band_sampling.py (klein-2026-11-04) -- depth-3 evidence for the full-digit selector with order dropping.
data/samp_11_*.txt and data/samp_12_*.txt hold ready pairs "H sigma" of random strong tournaments of orders 11 and 12
with vertex 1 a sole exit into vertex 0 (scripts/ready_search.c mode `sample`, seeds 1 and 2, 3,000,000 draws each),
kept only when 12000 <= H <= 14000, the band of the lineage's depth-2 leaves (12869..12877, which live at order 15).
Each pair realises 3H + 2 sigma as the Hamilton count of a strong C2 child of order 13 or 14; the script reports how
many odd numbers of the depth-3 window [38001, 39999] are reached, and whether the eight tree values
38609, 38611, 38615, 38617, 38627, 38629, 38633, 38635 are among them."""
import glob, os
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, '..', 'data')
pairs = {}; tot = 0
for f in sorted(glob.glob(os.path.join(DATA, 'samp_1[12]_*.txt'))):
    for l in open(f):
        p = l.split()
        if len(p) != 2: continue
        h, s = int(p[0]), int(p[1]); tot += 1; pairs.setdefault(h, set()).add(s)
sw = set(3 * h + 2 * s for h in pairs for s in pairs[h] if s >= 1)
win = list(range(38001, 40000, 2)); got = sum(1 for v in win if v in sw)
print('band ready pairs (orders 11-12, 12000 <= H <= 14000): %d pairs, %d distinct H of the 1001 odd values in the band, mean digits per H %.1f' % (tot, len(pairs), tot / max(1, len(pairs))))
print('C2 children cover %d of %d odd numbers in [38001, 39999] (%.3f); first uncovered: %s' % (got, len(win), got / len(win), [v for v in win if v not in sw][:12]))
print('tree level-3 values reached:', {v: (v in sw) for v in (38609, 38611, 38615, 38617, 38627, 38629, 38633, 38635)})
