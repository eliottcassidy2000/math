#!/usr/bin/env python3
"""reader_cde: arithmetic of the cell sum in THM-5141 (D).

For n-bit words, R = ceil(log2 n) + 1, t = ceil(n/4), the cells are (t', s) with t' >= t, s >= 0, w = n - t' - s >= 1
and |s - t'| <= R + 3 (the draft), or t' - R - 3 <= s <= t' + R (the carries bounds as proved).
Per cell the draft uses   B(t', s) = 4 * 4^w + 4 * 2^max(t'-1, s) * Cw * w * 2^(3w/2)   (Cw = 1 here).
We report log2 of  sum_cells B / 2^n  against  the draft's uniform bound  n (2R+7) [4^(R+4) + n (32n)^(5/2)]
and against the candidates n^(5/2) and n^2 (w-free S_w bound), to see which polylog power the cell sum really has.
All in exact integer / float-log arithmetic.
"""
import math, sys

def run(n, drop_w=False, asym=True):
    R = math.ceil(math.log2(n)) + 1
    t = -(-n // 4)
    tot1 = 0.0  # sum of 4*4^w / 2^n
    tot2 = 0.0  # sum of 4 * 2^max * w * 2^(3w/2) / 2^n
    cells = 0
    for tp in range(t, n):
        lo = max(0, tp - R - 3)
        hi = tp + R if asym else tp + R + 3
        for s in range(lo, hi + 1):
            w = n - tp - s
            if w < 1:
                continue
            cells += 1
            tot1 += 2.0 ** (2 + 2 * w - n)
            e = 2 + max(tp - 1, s) + 1.5 * w - n
            tot2 += (1 if drop_w else w) * 2.0 ** e
    return R, cells, tot1, tot2

print("n, R, cells, sum1/2^n, sum1/(n^2), sum2/2^n, sum2/n^(5/2) [with factor w], sum2'/n^(3/2) [w dropped], draft uniform bound / n^(9/2) ln n")
for n in [8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 65536]:
    R, cells, t1, t2 = run(n)
    _, _, _, t2d = run(n, drop_w=True)
    draft = n * (2 * R + 7) * (4.0 ** (R + 4) + n * (32.0 * n) ** 2.5)
    print(f"n={n:6d} R={R:2d} cells={cells:7d}  sum1={t1:.4e} sum1/n^2={t1/n**2:.4f}  sum2={t2:.4e} sum2/n^2.5={t2/n**2.5:.4f}  sum2'/n^1.5={t2d/n**1.5:.4f}  draft/(n^4.5 ln n)={draft/(n**4.5*math.log(n)):.1f}")
