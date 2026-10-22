# The 3-adic Sierpinski tournament on an initial interval: Sigma_N = S_infinity restricted to {0, ..., N-1}
# (x -> y iff the unit part of y - x is 1 mod 3).  The residue classes mod 3 inside {0..N-1} are copies of
# Sigma_(N_b) (N_b = #{x < N : x = b mod 3}), class b => class b+1, so
#     Sigma_N = C_3[Sigma_(N_0), Sigma_(N_1), Sigma_(N_2)],
# and the C_3 substitution law (THM-3121; W = (e1+e2+3e3)/(1-e3)) gives the full ordered path-cover profile
# recursively.  Output: H(Sigma_N) for N <= NMAX, checked by Held-Karp for N <= 14.  mac-mini-2026-10-22-sierp.
import sys, math
import numpy as np
from functools import lru_cache

def step(P, M):
    MA, MB, MC = M
    R = np.zeros_like(P)
    R[1:, :, :] += P[:-1, :, :]; R[:, 1:, :] += P[:, :-1, :]; R[:, :, 1:] += P[:, :, :-1]
    R[1:, 1:, :] += P[:-1, :-1, :]; R[:, 1:, 1:] += P[:, :-1, :-1]; R[1:, :, 1:] += P[:-1, :, :-1]
    R[1:, 1:, 1:] += 3 * P[:-1, :-1, :-1]
    Q = R; s = 1; L = min(M)
    while s <= L:
        Q2 = Q.copy(); Q2[s:, s:, s:] += Q[:MA + 1 - s, :MB + 1 - s, :MC + 1 - s]; Q = Q2; s *= 2
    return Q

def covers3(bA, bB, bC, J):
    M = (len(bA) - 1, len(bB) - 1, len(bC) - 1)
    vA = np.array([0] + list(bA[1:]), dtype=object); vB = np.array([0] + list(bB[1:]), dtype=object); vC = np.array([0] + list(bC[1:]), dtype=object)
    P = np.zeros((M[0] + 1, M[1] + 1, M[2] + 1), dtype=object); P[0, 0, 0] = 1
    c = [0] * (J + 1)
    for j in range(1, J + 1):
        P = step(P, M)
        t = np.tensordot(P, vC, axes=([2], [0])); t = np.tensordot(t, vB, axes=([1], [0]))
        c[j] = int(np.dot(t, vA))
    return c

@lru_cache(maxsize=None)
def profile(N):
    """ordered path-cover numbers c_0..c_N of Sigma_N (c_0 unused = 0 except N = 0)"""
    if N == 1: return (0, 1)
    if N == 2: return (0, 1, 2)          # arc 0 -> 1: c_1 = 1, c_2 = 2
    Ns = [len(range(b, N, 3)) for b in range(3)]
    return tuple(covers3(profile(Ns[0]), profile(Ns[1]), profile(Ns[2]), N))

def held_karp(N):
    def arc(x, y):
        z = y - x
        while z % 3 == 0: z //= 3
        return z % 3 == 1
    out = [sum(1 << y for y in range(N) if y != x and arc(x, y)) for x in range(N)]
    dp = [[0] * N for _ in range(1 << N)]
    for v in range(N): dp[1 << v][v] = 1
    for m in range(1, 1 << N):
        for v in range(N):
            c = dp[m][v]
            if not c: continue
            nxt = out[v] & ~m
            while nxt:
                w = (nxt & -nxt).bit_length() - 1; nxt &= nxt - 1
                dp[m | (1 << w)][w] += c
    return sum(dp[(1 << N) - 1])

if __name__ == '__main__':
    NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    for N in range(1, NMAX + 1):
        H = profile(N)[1]
        chk = ''
        if N <= 14: chk = ' HK ok' if held_karp(N) == H else ' HK MISMATCH'
        kap = H * 2 ** (N - 1) / math.factorial(N)
        v3 = 0; t = H
        while t % 3 == 0 and t: t //= 3; v3 += 1
        print(f"N={N:3d} H={H} v3={v3} kappa={kap:.6f}{chk}", flush=True)
