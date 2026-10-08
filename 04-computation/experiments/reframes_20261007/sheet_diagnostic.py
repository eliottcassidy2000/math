#!/usr/bin/env python3
"""P1 sheet diagnostic for the orphan law (HYP-9242).  Under 3x-1 (T(x) = x/2 or (3x-1)/2) the conjugate of the Mersenne
line is P_K = 2^K + 1 (x -> -x maps 3x+1 on -(2^K+1) to 3x-1 on 2^K+1); deletion children P_(K-D), run ends
x_K = 2*3^(K-1) + 1 with x_K = 3^D x_(K-D) + 1 - 3^D.  Orbits end in one of the positive 3x-1 cycles {1}, {5,7,10},
{17,...,34}.  Partner: P_K and P_K' (K' < K) enter the same cycle at the same point, at the same shifted time
(Terras time to entry minus K equal) with equal odd counts -- equivalent to an equal-time merge before the cycle.
Orphan: odd K with no partner K' < K.  Compare with the 3x+1 Mersenne line (mersenne_sigma_12800.txt)."""
import sys, math
KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
cyc = {}
for start in (1, 5, 17):
    x = start
    while True:
        cyc[x] = start
        x = x // 2 if x % 2 == 0 else (3 * x - 1) // 2
        if x == start: break
def sig(n):
    t = o = 0; x = n
    while x not in cyc:
        if x & 1: x = (3 * x - 1) >> 1; o += 1
        else: x >>= 1
        t += 1
        if t > 10**6: raise RuntimeError(f"no cycle within 1e6 steps for n = 2^{n.bit_length()-1}+1")
    return (cyc[x], x, o, t)
data = {}
for K in range(2, KMAX + 1):
    c, z, o, t = sig((1 << K) + 1)
    data[K] = (c, z, o, t - K)
# residual class (first post-run letter 2): for 3x-1 on 2^K+1 the run end -x_K = -(2*3^(K-1)+1) is 9 mod 16 iff K is EVEN
# (for 3x+1 on 2^K-1 it is K ODD); the other parity always merges with K-1 three steps after the run.
seen = {}; orphans = []; orph_all = []; res = 0; cyc_count = {}
for K in range(2, KMAX + 1):
    key = data[K]
    cyc_count[key[0]] = cyc_count.get(key[0], 0) + 1
    if key not in seen:
        orph_all.append(K)
        if K % 2 == 0: orphans.append(K)
    if K % 2 == 0: res += 1
    seen.setdefault(key, K)
print(f"3x-1 on 2^K+1, K <= {KMAX}: cycle shares {cyc_count}; residual (even) K {res}; orphans among them {len(orphans)}; orphans among all K {len(orph_all)} (odd: {sum(1 for K in orph_all if K % 2)})")
# 3x+1 Mersenne comparison from the exact table
m = {}
for line in open('../runcompress_20261007/mersenne_sigma_12800.txt'):
    K, o, t = map(int, line.split()); m[K] = (o, t - K)
seenm = set(); orph_m = []; orph_m_all = []
for K in range(2, KMAX + 1):
    if m[K] not in seenm:
        orph_m_all.append(K)
        if K % 2 == 1: orph_m.append(K)
    seenm.add(m[K])
print(f"3x+1 on 2^K-1: orphans among all K {len(orph_m_all)} (even: {sum(1 for K in orph_m_all if K % 2 == 0)})")
print(f"3x+1 on 2^K-1, K <= {KMAX}: orphans {len(orph_m)}")
for lo, hi in ((100, 400), (400, 1600), (1600, 6400), (6400, 12801)):
    if lo >= KMAX: break
    hi2 = min(hi, KMAX + 1)
    n_odd = sum(1 for K in range(lo, hi2) if K % 2 == 1)   # residual-class counts are equal (even vs odd K)
    a = sum(1 for K in orphans if lo <= K < hi2); b = sum(1 for K in orph_m if lo <= K < hi2)
    mid = math.sqrt((lo + hi2) / 2)
    print(f"   K in [{lo},{hi2}): residual-class orphan fraction 3x-1 {a/n_odd:.4f} (x sqrt K {a/n_odd*mid:.2f})   3x+1 {b/n_odd:.4f} (x sqrt K {b/n_odd*mid:.2f})")
