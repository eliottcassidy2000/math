#!/usr/bin/env python3
"""The owner's question as a table: S19's coupled pair y (odd, 2-adic Haar, here a random 8000-bit odd integer) and
x = 3*2^v*y + 1 (v >= 2, P(v=k) = 2^-(k-1)).  Exponent streams a_s = v2(3 y_s + 1), b_s = v2(3 x_s + 1) in odd-step
time.  For windows of s we report Corr(b_s, a_(s+d)) for d = -2..3, plus the share of pairs whose x-orbit has merged
with the y-orbit (x_s = y_(s+d0) for some fixed lag d0 from then on).  Merges are detected exactly (x_s == y_(s+d))."""
import random, math, sys
def v2(x): return (x & -x).bit_length() - 1
def U(x):
    y = 3 * x + 1; return y >> v2(y)
rng = random.Random(20261007); P = 1500; S = 1200; BITS = 8000
A = []; B = []; MERGE = []
for p in range(P):
    y = rng.getrandbits(BITS) | 1 | (1 << (BITS - 1))
    v = 2
    while rng.random() < 0.5: v += 1
    x = 3 * (1 << v) * y + 1
    ys = [y]; xs = [x]
    for s in range(S + 5):
        ys.append(U(ys[-1])); xs.append(U(xs[-1]))
    a = [v2(3 * t + 1) for t in ys]; b = [v2(3 * t + 1) for t in xs]
    # first s with x_s == y_(s+d) for some d in [-6, 6]
    m = None
    yi = {}
    for i, t in enumerate(ys): yi.setdefault(t, i)
    for s, t in enumerate(xs[:S]):
        if t in yi and abs(yi[t] - s) <= 6: m = (s, yi[t] - s); break
    A.append(a); B.append(b); MERGE.append(m)
lags = [l for _, l in (mm for mm in MERGE if mm)]
from collections import Counter
print(f"pairs {P}, odd steps {S}; merged within {S}: {sum(1 for m in MERGE if m)}/{P}; merge lag d0 counts: {dict(Counter(lags))}")
def corr(xs, ys):
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    sxy = sum((u-mx)*(w-my) for u, w in zip(xs, ys)); sxx = sum((u-mx)**2 for u in xs); syy = sum((w-my)**2 for w in ys)
    return sxy / math.sqrt(sxx * syy)
print("window of s      merged-share   Corr(b_s, a_(s+d)) for d = -2, -1, 0, 1, 2, 3")
for lo, hi in [(5, 15), (30, 50), (100, 150), (300, 400), (800, 1000), (1100, 1195)]:
    share = sum(1 for m in MERGE if m and m[0] <= lo) / P
    row = []
    for d in (-2, -1, 0, 1, 2, 3):
        xs_ = []; ys_ = []
        for p in range(P):
            for s in range(lo, hi):
                if 0 <= s + d < len(A[p]): xs_.append(B[p][s]); ys_.append(A[p][s + d])
        row.append(corr(xs_, ys_))
    print(f"s in [{lo:4d},{hi:4d}):  {share:.3f}        " + "  ".join(f"{c:+.3f}" for c in row), flush=True)
