#!/usr/bin/env python3
"""
alpha_from_cc_bin.py -- Load precomputed cc binary (from fast_cycle_cc.c) and run SSC.

Usage:
  python3 alpha_from_cc_bin.py <n> [circulant|paley]

Loads cc_n<n>_<type>.bin from 05-knowledge/results/, runs the numpy SSC
with 2-prime CRT to recover exact alpha_k values.

This separates the slow cycle_cc step (now done in C, ~10s for n=25) from
the numpy SSC step (done in Python, ~2000s for n=25).
"""

import numpy as np
import time, sys
from math import factorial

P1 = (1 << 31) - 1    # 2147483647 (Mersenne prime)
P2 = (1 << 31) - 19   # 2147483629

# ─── CRT ─────────────────────────────────────────────────────────────────────

def _ext_gcd(a, b):
    if b == 0: return a, 1, 0
    g, x, y = _ext_gcd(b, a % b)
    return g, y, x - (a // b) * y

def crt2(r1, r2, m1, m2):
    _, u, _ = _ext_gcd(m1, m2)
    return (r1 + m1 * ((r2 - r1) * u % m2)) % (m1 * m2)

# ─── Fast SOS using reshape trick ────────────────────────────────────────────

def sos_forward_inplace(A_r, n, p):
    N = 1 << n
    for i in range(n):
        bit = 1 << i
        A_r_4d = A_r.reshape(n+1, N // (2 * bit), 2, bit)
        A_r_4d[:, :, 1, :] = (A_r_4d[:, :, 1, :] + A_r_4d[:, :, 0, :]) % p

def sos_inverse_inplace(A_r, n, p):
    N = 1 << n
    for i in range(n):
        bit = 1 << i
        A_r_4d = A_r.reshape(n+1, N // (2 * bit), 2, bit)
        A_r_4d[:, :, 1, :] = (A_r_4d[:, :, 1, :] - A_r_4d[:, :, 0, :] + p) % p

def make_popcount(n):
    N = 1 << n
    pc = np.zeros(N, dtype=np.int32)
    for i in range(n):
        pc += (np.arange(N, dtype=np.int32) >> i) & 1
    return pc

# ─── SSC mod p ───────────────────────────────────────────────────────────────

def ssc_mod(cc_mod, G_mod, n, p, pc):
    N = 1 << n
    nranks = n + 1
    cc_ranks = list(range(3, n+1, 2))

    A_r = np.zeros((nranks, N), dtype=np.int64)
    B_r = np.zeros((nranks, N), dtype=np.int64)
    for k in range(nranks):
        mk = (pc == k)
        A_r[k][mk] = cc_mod[mk]
        B_r[k][mk] = G_mod[mk]

    sos_forward_inplace(A_r, n, p)
    sos_forward_inplace(B_r, n, p)

    C_r = np.zeros((nranks, N), dtype=np.int64)
    for k in range(nranks):
        for j in cc_ranks:
            if j > k:
                break
            kj = k - j
            prod = A_r[j] * B_r[kj] % p
            C_r[k] = (C_r[k] + prod) % p

    sos_inverse_inplace(C_r, n, p)

    result = np.zeros(N, dtype=np.int64)
    for k in range(nranks):
        mk = (pc == k)
        result[mk] = C_r[k][mk]

    return result

def run_mod(n, cc_np, p, pc):
    full = (1 << n) - 1
    kmax = n // 3
    N = 1 << n
    nranks = n + 1

    cc = cc_np % p

    # G_1 = SOS(cc) mod p
    G = cc.copy().reshape(1, N)
    for i in range(n):
        bit = 1 << i
        G_4d = G.reshape(1, N // (2 * bit), 2, bit)
        G_4d[0, :, 1, :] = (G_4d[0, :, 1, :] + G_4d[0, :, 0, :]) % p
    G = G.reshape(N)

    Gfull = {1: int(G[full])}
    print(f"    k=1: G_1[full] ≡ {Gfull[1]}")

    for k in range(2, kmax+1):
        t0 = time.time()
        G = ssc_mod(cc, G, n, p, pc)
        Gfull[k] = int(G[full])
        print(f"    k={k}: G_{k}[full] ≡ {Gfull[k]}  [{time.time()-t0:.1f}s]")
        sys.stdout.flush()

    return Gfull

# ─── Main ─────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 23
    type_name = sys.argv[2] if len(sys.argv) > 2 else 'circulant'
    kmax = n // 3

    cc_file = f"05-knowledge/results/cc_n{n}_{type_name}.bin"
    print(f"n={n}, type={type_name}, kmax={kmax}")
    print(f"Loading cc from {cc_file}...")

    t0 = time.time()
    cc_raw = np.fromfile(cc_file, dtype=np.int64)
    N = 1 << n
    if len(cc_raw) != N:
        print(f"ERROR: expected {N} entries, got {len(cc_raw)}")
        sys.exit(1)
    print(f"  loaded: {time.time()-t0:.1f}s")

    print("Precomputing popcounts...")
    t0 = time.time()
    pc = make_popcount(n)
    print(f"  done: {time.time()-t0:.1f}s")
    sys.stdout.flush()

    print(f"\n=== Modular run 1 (p={P1}) ===")
    t1 = time.time()
    r1 = run_mod(n, cc_raw, P1, pc)
    print(f"  Total run 1: {time.time()-t1:.1f}s")
    sys.stdout.flush()

    print(f"\n=== Modular run 2 (p={P2}) ===")
    t2 = time.time()
    r2 = run_mod(n, cc_raw, P2, pc)
    print(f"  Total run 2: {time.time()-t2:.1f}s")
    sys.stdout.flush()

    print(f"\n=== CRT recovery ===")
    alphas = {}
    for k in range(1, kmax+1):
        Gk = crt2(r1[k], r2[k], P1, P2)
        fk = factorial(k)
        alphas[k] = Gk // fk
        print(f"  α_{k} = {alphas[k]:,}")

    H = 1 + sum(2**k * v for k, v in alphas.items())
    print(f"\n=== FINAL RESULTS n={n} type={type_name} ===")
    print(f"H = {H:,}")

    print(f"\nTerm breakdown:")
    for k in range(1, kmax+1):
        v = alphas.get(k, 0)
        pct = 100 * 2**k * v / (H - 1) if H > 1 else 0
        print(f"  2^{k} × α_{k} = {2**k * v:>35,}   ({pct:.1f}%)   α_{k}={v:,}")

    print(f"\nRatios:")
    a1 = alphas[1]; a2 = alphas[2]; a3 = alphas.get(3, 0)
    print(f"  α₁/(2α₂) = {a1/(2*a2):.6f}")
    if a3:
        print(f"  α₃/α₂    = {a3/a2:.6f}")
        print(f"  8α₃/2α₁  = {8*a3/(2*a1):.6f}  (>1 means 8α₃ > 2α₁)")
