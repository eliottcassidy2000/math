#!/usr/bin/env python3
"""Independent check of sampled witnesses: rebuild tournament #index of tiltsample(n, q, seed) by replaying the
xorshift128+ stream in Python, then recount H and P2 by a separate subset DP.
usage: verify_witnesses.py n q seed index [index ...]   (prints H P2 for each)"""
import sys
M64 = (1 << 64) - 1
def stream(seed):
    s0 = (seed * 0x9E3779B97F4A7C15 + 1) & M64; s1 = seed ^ 0xD1B54A32D192ED03
    st = [s0, s1]
    def rnd():
        a = st[0]; b = st[1]; st[0] = b
        a = (a ^ (a << 23)) & M64
        st[1] = a ^ b ^ (a >> 17) ^ (b >> 26)
        return (st[1] + b) & M64
    for _ in range(20): rnd()
    return rnd
def counts(out, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]; H = [0] * (1 << n)
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        row = dp[S]; t = 0
        for v in range(n):
            x = row[v]
            if not x: continue
            t += x; nb = out[v] & ~S
            while nb:
                u = (nb & -nb).bit_length() - 1; nb &= nb - 1; dp[S | 1 << u][u] += x
        H[S] = t
    return H[full], sum(H[S] * H[full ^ S] for S in range(1, full))
if __name__ == '__main__':
    n, q, seed = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
    want = sorted(int(x) for x in sys.argv[4:])
    thr = int(q * 18446744073709551615.0)
    rnd = stream(seed); it = 0
    for target in want:
        while True:
            out = [0] * n
            for i in range(n):
                for j in range(i + 1, n):
                    if rnd() < thr: out[i] |= 1 << j
                    else: out[j] |= 1 << i
            if it == target:
                h, p2 = counts(out, n); print(target, h, p2); it += 1; break
            it += 1
