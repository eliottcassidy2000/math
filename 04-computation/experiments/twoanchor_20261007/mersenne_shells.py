#!/usr/bin/env python3
"""Mersenne residual branch by shells m = v2(K-1) (mac-mini-2026-10-07-twoanchor).
Odd K: n = 2^K - 1 has run end x = 2*3^(K-1) - 1 with j = v2(x-1) = 3 + m, J = 1 + floor(m/2), r = 1 + (m mod 2).
For m >= 4 the K -> K-3 chain (child 2^(K-3) - 1) is in the universal state (3, 1-27) at the end of the two-run;
its post-run absorption law should not depend on m (THM-4601 (v)). Also lists the exponents certified by the two
explicit doubly-uniform head pairs (post-run depth 11 for r = 1, 9 for r = 2)."""
import sys
from collections import defaultdict
def v2(x): return (x & -x).bit_length() - 1
def T(x): return (3*x + 1) >> 1 if x & 1 else x >> 1
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12000
S = 200
shell = defaultdict(lambda: [0, 0, 0])   # m -> [count, certified within S, depth<=11 explicit]
explicit = []
for K in range(17, KMAX, 2):
    m = v2(K - 1)
    if m < 4: continue
    J = 1 + m//2; r = 1 + (m % 2)
    x = 2*3**(K-1) - 1; y = (x + 1)//27 - 1
    u, v, k = x, y, 3
    for s in range(2*J):
        k += (u & 1) - (v & 1); u, v = T(u), T(v)
    assert k == 3 and u - 1 == 27*(v - 1)
    dep = None
    for s in range(S + 1):
        if u == v and k == 0: dep = s; break
        k += (u & 1) - (v & 1); u, v = T(u), T(v)
    shell[m][0] += 1
    if dep is not None:
        shell[m][1] += 1
        if dep <= 11: shell[m][2] += 1; explicit.append((K, m, dep))
print(f"odd K < {KMAX}, shells m = v2(K-1) >= 4; K -> K-3 certified within post-run depth {S}:")
for m in sorted(shell):
    c, a, e = shell[m]
    print(f"  m={m:2d} (r={1+(m%2)}): {a}/{c} = {a/c:.3f}   (post-run depth <= 11: {e})")
tot = [sum(shell[m][i] for m in shell if (m % 2) == par) for par in (0, 1) for i in (0, 1)]
print(f"  pooled r=1 (m even): {tot[1]}/{tot[0]} = {tot[1]/tot[0]:.3f};  r=2 (m odd): {tot[3]}/{tot[2]} = {tot[3]/tot[2]:.3f}")
print("  shortest certificates (K, m, post-run depth):", sorted(explicit, key=lambda z: z[2])[:10])
