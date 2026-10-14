"""Audit E (2026-10-14) of HYP-9263 'Evidence': independent recomputation of 20,000 partial quotients of theta = log2 3
with mpmath at a different precision (22,000 digits), plus a self-consistency check that the precision suffices
(recompute at 23,000 digits and compare all quotients). Reports the Khinchin geometric mean, the Levy constant proxy,
the largest quotients, and max_{n > 1000} log a_(n+1) / log q_n (claimed < 1e-3 in HYP-9263) with all exceedances."""
import time
from mpmath import mp, log, floor, pi
from math import log as mlog, exp as mexp

def cf(N, dps):
    mp.dps = dps
    x = log(3) / log(2); a = []
    for _ in range(N):
        ai = int(floor(x)); a.append(ai); x = 1 / (x - ai)
    return a

N = 20000
t = time.time(); a = cf(N, 22000); t1 = time.time() - t
b = cf(N, 23000)
print("20000 partial quotients at 22000 and 23000 digits agree:", a == b, "(%.1fs)" % t1)
print("first 60:", a[:60])
geo = sum(mlog(x) for x in a[1:N]) / (N - 1)
qprev, qcur = 0, 1; mx = (0, -1); best = (0.0, -1); exceed = []
for i in range(1, N):
    qprev, qcur = qcur, a[i] * qcur + qprev
    if a[i] > mx[0]: mx = (a[i], i)
    if i > 1000 and i + 1 < N:
        r = mlog(a[i + 1]) / mlog(qcur)
        if r > best[0]: best = (r, i)
        if r >= 1e-3: exceed.append((i, a[i + 1], round(r, 5)))
K = 2.685452001065306; L = float(pi ** 2 / (12 * log(2)))
print("geometric mean a_1..a_%d = %.6f (Khinchin %.6f)" % (N - 1, mexp(geo), K))
print("log q_n / n at n = %d: %.6f (Levy %.6f)" % (N - 1, mlog(qcur) / (N - 1), L))
print("largest partial quotient a_%d = %d" % (mx[1], mx[0]))
print("max_{1000 < n < %d} log a_(n+1)/log q_n = %.6f at n = %d (a_(n+1) = %d)" % (N - 1, best[0], best[1], a[best[1] + 1]))
print("number of n in (1000, %d) with log a_(n+1)/log q_n >= 1e-3: %d; first ten (n, a_(n+1), ratio): %s" % (N - 1, len(exceed), exceed[:10]))
print("largest n with ratio >= 1e-3:", exceed[-1] if exceed else None)
top = sorted(range(1, N), key=lambda i: -a[i])[:10]
print("ten largest (index, value):", [(i, a[i]) for i in top])
print("max of a_1..a_40:", max(a[1:41]), "; max of a_0..a_39:", max(a[:40]))
