#!/usr/bin/env python3
"""HYP-9214 table: equal-odd-step-time merge of a random first-reset-2 source n (B bits) with (n-1)/2, before 1
(mac-mini, 2026-10-06).  Sources n = 2^(r+1) t - 1, r >= 1, first reset exponent 2 (b = 1 + v_2(3^r t - 1) >= 3).
Prints, per size B: merge rate, the fraction of sources merging within 10 steps beyond the run, the median lag, and
the merge rate split by the parity of the debt height h = b - 2.  Seeds 2026 (main table) and 99 (large sizes).
Run: python3 sevens_20261006_debt_resolution_trend.py   (about 20 s)
"""
import random, math, statistics
def U(x):
    y = 3 * x + 1
    return y >> ((y & -y).bit_length() - 1)
def v2(x): return (x & -x).bit_length() - 1
def meet(n, m):
    j = 0
    while n != 1 and m != 1:
        n, m, j = U(n), U(m), j + 1
        if n == m:
            return j if n != 1 else None
    return None
def run(Bs, seed, N, N_big=None):
    rng = random.Random(seed)
    print(f"seed {seed}:")
    print("     B      N   merge-rate   SE    within-10 (of sources)   median lag   rate (h even / h odd)")
    for B in Bs:
        n_s = N_big if (N_big and B >= 2048) else N
        tot = hit = w10 = 0; lags = []; par = {0: [0, 0], 1: [0, 0]}
        while tot < n_s:
            n = rng.getrandbits(B) | (1 << (B - 1)) | 1
            r = v2(n + 1) - 1
            if r < 1:
                continue
            t = (n + 1) >> (r + 1)
            b = 1 + v2(3 ** r * t - 1)
            if b < 3:
                continue
            tot += 1
            j = meet(n, (n - 1) // 2)
            h = b - 2
            par[h % 2][0] += 1
            if j is not None:
                hit += 1; lags.append(j - r); par[h % 2][1] += 1
                w10 += (j - r) <= 10
        p = hit / tot
        pe = par[0][1] / max(1, par[0][0]); po = par[1][1] / max(1, par[1][0])
        print(f"{B:6d} {tot:6d}     {p:.3f}     {math.sqrt(p * (1 - p) / tot):.3f}          {w10 / tot:.3f}              {statistics.median(lags) if lags else '-':>6}       {pe:.3f} / {po:.3f}", flush=True)
run((16, 32, 64, 128, 256, 512, 1024, 2048), 2026, 1000, 600)
run((2048, 4096, 8192), 99, 400)
