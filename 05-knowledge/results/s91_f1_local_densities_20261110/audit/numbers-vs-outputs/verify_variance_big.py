#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: relative variance var(r)/mean(r)^2 of the F1 representation count on each sheet for
the pre-registered window of odd N in [2^24, 2^27) (direct enumeration, own kernel of verify_misc.py)."""
import numpy as np
from verify_misc import build_q, window

LO, K = 24, 27
A, B = 1 << LO, 1 << K
q = build_q(B // 2 + 2)
nodd = (B - A) // 2
r0 = np.zeros(nodd, np.int32)
r1 = np.zeros(nodd, np.int32)
t0, t1 = window(q, A, B, K, r0, r1)
del q
print(f"[2^{LO}, 2^{K}): triples per odd N: sheet 0 {t0 / nodd:.4f}, sheet 1 {t1 / nodd:.4f}, total {(t0 + t1) / nodd:.4f}; triples {t0 + t1}")
for name, a in (("sheet 0", r0), ("sheet 1", r1)):
    m = a.mean(dtype=np.float64)
    v = (a.astype(np.float64) ** 2).mean() - m * m
    print(f"   {name}: mean {m:.4f}, variance {v:.4f}, var/mean^2 {v / m ** 2:.4f}, var/mean {v / m:.3f}")
r = r0 + r1
m = r.mean(dtype=np.float64)
v = (r.astype(np.float64) ** 2).mean() - m * m
print(f"   both   : mean {m:.4f}, variance {v:.4f}, var/mean^2 {v / m ** 2:.4f}, var/mean {v / m:.3f}; misses {int((r == 0).sum())} ({100.0 * (r == 0).sum() / nodd:.4f} %)")
