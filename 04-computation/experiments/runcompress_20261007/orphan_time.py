#!/usr/bin/env python3
"""Available-time test of the orphan law: for odd K in [1000, 6400], compare the post-run Terras length (sigma_T(M_K) - K)/K of
deletion orphans vs non-orphans, and the orphan rate by tercile of that normalized length."""
import statistics
data = {}
for line in open('mersenne_sigma_6400.txt'):
    K, o, t = map(int, line.split()); data[K] = (o, t)
by_odd = {}; rows = []
for K in sorted(data):
    o, t = data[K]
    has = any(t - data[Kp][1] == K - Kp for Kp in by_odd.get(o, []))
    if K % 2 == 1 and 1000 <= K <= 6400: rows.append(((t - K) / K, not has))
    by_odd.setdefault(o, []).append(K)
orf = [x for x, o in rows if o]; nof = [x for x, o in rows if not o]
print(f"normalized post-run length (sigma_T - K)/K: orphans mean {statistics.mean(orf):.3f} (n={len(orf)}), non-orphans {statistics.mean(nof):.3f} (n={len(nof)})")
xs = sorted(x for x, _ in rows); q1, q2 = xs[len(xs)//3], xs[2*len(xs)//3]
for lo, hi, name in ((0, q1, 'short third'), (q1, q2, 'middle third'), (q2, 99, 'long third')):
    sel = [o for x, o in rows if lo <= x < hi]
    print(f"  {name:12s} ({lo:.2f}..{hi:.2f}): orphan rate {sum(sel)/len(sel):.4f}  (n={len(sel)})")
