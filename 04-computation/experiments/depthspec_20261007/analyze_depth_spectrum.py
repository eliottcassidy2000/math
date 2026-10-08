#!/usr/bin/env python3
import math, statistics
rows = []
for line in open('depth_spectrum_1001_3999.txt'):
    K, d1, d3, desc = line.split()
    K = int(K); desc = int(desc)
    d1 = None if d1 == 'None' else int(d1); d3 = None if d3 == 'None' else int(d3)
    best = min([d for d in (d1, d3) if d is not None], default=None)
    rows.append((K, d1, d3, best, desc))
N = len(rows)
print(f"odd K in [1001, 3999]: {N}")
print(f"  D=1 merges: {sum(r[1] is not None for r in rows)/N:.3f}; D=3 merges: {sum(r[2] is not None for r in rows)/N:.3f}; either: {sum(r[3] is not None for r in rows)/N:.3f}")
bests = [r[3] for r in rows if r[3] is not None]
print(f"  median best deletion depth {statistics.median(bests)}, mean {statistics.mean(bests):.0f}; median descent depth {statistics.median(r[4] for r in rows)} (= {statistics.median(r[4]/r[0] for r in rows):.2f} K)")
print(f"  deletion certificate shorter than descent: {sum(1 for r in rows if r[3] is not None and r[3] < r[4])/N:.3f}")
print("  survival P(best deletion depth > s), counting 'no merge' as infinite:")
for s in (10, 30, 100, 300, 1000, 3000, 10000):
    p = sum(1 for r in rows if r[3] is None or r[3] > s) / N
    print(f"    s={s:6d}: {p:.4f}   p*sqrt(s) = {p*math.sqrt(s):.2f}")
# depth relative to K for merged ones
rel = sorted(r[3] / r[0] for r in rows if r[3] is not None)
print("  quantiles of best depth / K:", [round(rel[int(q*len(rel))], 3) for q in (0.1, 0.25, 0.5, 0.75, 0.9, 0.99)])
