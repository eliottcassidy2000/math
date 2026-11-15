#!/usr/bin/env python3
"""Wave 38: fibre widths beyond exhaustive orders, by sampling.

Input files s<n>.txt ("H P2" per uniformly random labelled tournament, from w38_sample_pairs.c) and the exact class
file pairs10.txt.  For each order the script takes the samples whose H lies within +-REL of the median H (a narrow
relative window, so the ratio R = P2/H removes the leading scaling) and reports the spread of R: standard deviation,
central quantile ranges, and the full sampled range.  At order 10 the same statistics are computed on the exact
class list and compared with the exact median fibre width / H (3.025, w38_fibre_widths.py), which calibrates how a
sampled window range relates to a true fibre.
Usage: python3 w38_sampled_fibres.py SAMPLE_DIR PAIRS10_FILE
"""
import sys
from statistics import median, pstdev

D, P10 = sys.argv[1], sys.argv[2]
REL = 0.01


def window_stats(pairs, label):
    hs = sorted(h for h, _ in pairs)
    hm = hs[len(hs) // 2]
    R = sorted(p / h for h, p in pairs if abs(h - hm) <= REL * hm)
    k = len(R)

    def q(t):
        return R[min(k - 1, max(0, int(t * k)))]
    print(f"{label}: median H = {hm}, window +-{REL:.0%}: {k} tournaments; R = P2/H mean {sum(R) / k:.3f}, "
          f"sd {pstdev(R):.3f}; q05-q95 {q(0.95) - q(0.05):.3f}; q01-q99 {q(0.99) - q(0.01):.3f}; "
          f"full range {R[-1] - R[0]:.3f}; sd/(2(n-1)) = {pstdev(R) / (2 * (n - 1)):.4f}")
    return pstdev(R)


n = 10
exact = [tuple(map(int, l.split())) for l in open(P10)]
window_stats(exact, "order 10 exact classes")
for n in range(10, 17):
    try:
        pairs = [tuple(map(int, l.split())) for l in open(f"{D}/s{n}.txt")]
    except FileNotFoundError:
        continue
    window_stats(pairs, f"order {n} sampled ({len(pairs)})")
