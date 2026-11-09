"""S19 audit: HYP-9655 (C) one-apex share table.
One-apex values H(W) = 1 + 2a(W), a(0) = 0, a(2W+1) = a(W), a(2W) = a(W) + ceil(W/2).
Exactness of a table computed from W < 2^B: H(2W+1) = H(W), and for even W' = 2W >= 2, a(W') >= ceil(W/2), so
H(W') >= 1 + W'/2. Hence every value v is attained by W = 0 or an even W' <= 2(v - 1); values < 2^(B-1) are complete.
Also recomputes the share of odd targets in [1001, 200000], and checks against the core-{u} closed form of THM-5085
(h = 1 + N - Q(N) - [N odd], the tail values of the one-vertex core).
"""
import sys
import numpy as np

B = int(sys.argv[1]) if len(sys.argv) > 1 else 22
M = 1 << B
a = np.zeros(M, dtype=np.int64)
for W in range(1, M):
    a[W] = a[W >> 1] if W & 1 else a[W >> 1] + ((W >> 1) + 1) // 2
H = 1 + 2 * a
vals = set(H.tolist())
complete_below = 1 << (B - 1)
print(f"W < 2^{B}: values complete below 2^{B-1} = {complete_below} (proof in the docstring)")
rows = []
for k in [6, 8, 10, 12, 14, 16, 17, 18, 19]:
    lo, hi = 1 << k, 1 << (k + 1)
    assert hi <= complete_below
    odd = hi // 2 - lo // 2
    got = sum(1 for v in range(lo + 1, hi, 2) if v in vals)
    rows.append((k, got, odd, got / odd))
for k, got, odd, s in rows:
    print(f"k={k:2d}: {got}/{odd} = {s:.4f}")
t = [v for v in range(1001, 200001, 2)]
got = sum(1 for v in t if v in vals)
print(f"[1001, 200000]: {got}/{len(t)}")
if B >= 22:
    lo, hi = 1 << 20, 1 << 21
    c20 = sum(1 for v in range(lo + 1, hi, 2) if v in vals)
    print(f"k=20: {c20}/{(hi - lo) // 2} = {c20 / ((hi - lo) // 2):.4f} (HYP-9480 records 0.4974 for the root image)")
    c6 = sum(1 for v in vals if v <= 999999)
    print(f"one-apex (root) odd values <= 999,999: {c6} (Session205 finite-oz-tail-bank and Session197 ONE_CORE no-core part: 249,077)")
# cross-check with the closed form of the core {u}: (h0,R0,D0,hu) = (1,1,-1,1), h(N) = 1 + N - Q(N) - [N odd]
Q = np.zeros(M, dtype=np.int64)
for m in range(1, M):
    k = m >> 1
    Q[m] = Q[k] + k + (k & 1) if m & 1 else Q[k]
N = np.arange(M, dtype=np.int64)
h = 1 + N - Q - (N & 1)
vals2 = set(int(v) for v in h if v < complete_below)
print("core-{u} closed-form value set equals the one-apex set below the completeness bound:",
      vals2 == set(v for v in vals if v < complete_below))
