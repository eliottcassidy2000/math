#!/usr/bin/env python3
"""audit_C (independent): odd-step count and Terras length to the first 1 of M_K = 2^K - 1, for K in [K0, K1].

Independent implementation (not the session's mersenne_sigma*.py):
  * the first K Terras steps are done in closed form: T^j(2^K - 1) = 3^j 2^(K-j) - 1 (j <= K), so after K odd steps the
    orbit is at 3^K - 1;
  * afterwards 16 Terras steps are applied at once with a residue table:
        T^16(2^16 q + r) = 3^(o(r)) q + T^16(r),   o(r) = #odd steps among the first 16 steps of r,
    valid while x >= 2^17 (then 1 cannot be reached inside the block, since each Terras step at most halves);
  * the tail below 2^17 is done one Terras step at a time.
Output lines: K odd_steps terras_steps (same format as mersenne_sigma_12800.txt).
Usage: python3 c4_mersenne_orbits.py K0 K1 outfile
"""
import sys

B = 16
MASK = (1 << B) - 1
ODD = [0] * (1 << B)
VAL = [0] * (1 << B)
for r in range(1 << B):
    x, o = r, 0
    for _ in range(B):
        if x & 1:
            x = (3 * x + 1) >> 1
            o += 1
        else:
            x >>= 1
    ODD[r] = o
    VAL[r] = x
P3 = [3 ** i for i in range(B + 1)]
LIM = 1 << (B + 1)


def mersenne_orbit(K):
    x = 3 ** K - 1          # T^K(2^K - 1)
    odd = K
    ter = K
    while x >= LIM:
        r = x & MASK
        o = ODD[r]
        x = P3[o] * (x >> B) + VAL[r]
        odd += o
        ter += B
    while x != 1:
        if x & 1:
            x = (3 * x + 1) >> 1
            odd += 1
        else:
            x >>= 1
        ter += 1
    return odd, ter


def plain(K):
    """straightforward reference (used for a self-check on small K)"""
    n = (1 << K) - 1
    odd = ter = 0
    while n != 1:
        if n & 1:
            n = (3 * n + 1) >> 1
            odd += 1
        else:
            n >>= 1
        ter += 1
    return odd, ter


if __name__ == "__main__":
    for K in range(2, 300):
        assert mersenne_orbit(K) == plain(K), K
    K0, K1, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    with open(out, "w") as f:
        for K in range(K0, K1 + 1):
            o, t = mersenne_orbit(K)
            f.write(f"{K} {o} {t}\n")
            if K % 200 == 0:
                f.flush()
