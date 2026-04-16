#!/usr/bin/env python3
"""
alpha_full_ssc_fast.py -- Compute ALL alpha_k via numpy SSC, reshape-SOS trick.

Key optimization over alpha_full_ssc_numpy.py:
  Instead of scatter-gather SOS (A_r[:, sel] += A_r[:, src]),
  reshape A_r to (n+1, 2^{n-i-1}, 2, 2^i) and do:
    A_r_4d[:, :, 1, :] += A_r_4d[:, :, 0, :]
  This accesses CONTIGUOUS memory blocks -- 10-20x faster than scatter-gather.

Speedup estimate:
  n=21: ~3-8s per SSC call (vs ~40s with scatter-gather)
  n=23: ~30-80s per SSC call (feasible: total ~20 min vs ~80 min)
  n=25: ~3-8 min per SSC call (borderline)

Uses two Mersenne-type primes + CRT for exact integer recovery.
P1 = 2^31 - 1 = 2147483647 (Mersenne prime)
P2 = 2^31 - 19 = 2147483629
Product P1*P2 ~ 4.6e18 > expected max G_k[full] for n<=25.

Overflow safety: (P1-1)^2 ≈ 4.6e18 < 9.2e18 = int64 max. Products fit before mod.
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
    """Smallest non-negative x ≡ r1 (mod m1) ≡ r2 (mod m2)."""
    _, u, _ = _ext_gcd(m1, m2)
    return (r1 + m1 * ((r2 - r1) * u % m2)) % (m1 * m2)

# ─── Tournament / cycle counting ─────────────────────────────────────────────

def circulant(n, S):
    adj = [0]*n
    for v in range(n):
        for k in S:
            adj[v] |= 1 << ((v+k) % n)
    return adj

def cycle_cc(adj, n):
    """Count directed odd cycles by vertex-set bitmask."""
    full = (1 << n)-1
    cc = [0] * (1 << n)
    for s in range(n):
        s_bit = 1 << s
        hi_mask = full & ~(s_bit - 1)
        queue = {(s_bit, s): 1}
        while queue:
            nq = {}
            for (mask, v), cnt in queue.items():
                L = bin(mask).count('1')
                if L >= 3 and L % 2 == 1:
                    if (adj[v] >> s) & 1:
                        cc[mask] += cnt
                cands = adj[v] & ~mask & hi_mask
                while cands:
                    ub = cands & -cands; u = ub.bit_length()-1
                    nq[(mask|ub, u)] = nq.get((mask|ub, u), 0) + cnt
                    cands ^= ub
            queue = nq
    return cc

def make_popcount(n):
    N = 1 << n
    pc = np.zeros(N, dtype=np.int32)
    for i in range(n):
        pc += (np.arange(N, dtype=np.int32) >> i) & 1
    return pc

# ─── Fast SOS using reshape trick ────────────────────────────────────────────

def sos_forward_inplace(A_r, n, p):
    """
    Forward SOS mod p on ALL ranks simultaneously using reshape trick.
    A_r: shape (n+1, 2^n), dtype int64, values in [0, p).
    Modifies A_r in-place.

    For dimension i: reshape to (n+1, 2^{n-i-1}, 2, 2^i),
    then A_r_4d[:, :, 1, :] = (A_r_4d[:, :, 1, :] + A_r_4d[:, :, 0, :]) % p.
    This accesses contiguous blocks of 2^i elements at a time.
    """
    N = 1 << n
    for i in range(n):
        bit = 1 << i
        A_r_4d = A_r.reshape(n+1, N // (2 * bit), 2, bit)
        A_r_4d[:, :, 1, :] = (A_r_4d[:, :, 1, :] + A_r_4d[:, :, 0, :]) % p

def sos_inverse_inplace(A_r, n, p):
    """
    Inverse SOS (Möbius) mod p on ALL ranks simultaneously.
    Modifies A_r in-place.
    """
    N = 1 << n
    for i in range(n):
        bit = 1 << i
        A_r_4d = A_r.reshape(n+1, N // (2 * bit), 2, bit)
        A_r_4d[:, :, 1, :] = (A_r_4d[:, :, 1, :] - A_r_4d[:, :, 0, :] + p) % p

# ─── Fast SSC mod p ──────────────────────────────────────────────────────────

def ssc_mod(cc_mod, G_mod, n, p, pc):
    """
    Compute (cc ★ G)[T] mod p for all T using reshape-SOS trick.

    cc_mod, G_mod: numpy int64 arrays of size 2^n, values in [0, p).
    Returns numpy int64 array, values in [0, p).
    """
    N = 1 << n
    nranks = n + 1

    # Non-zero ranks of cc: only ODD >= 3 (only odd-length cycles)
    cc_ranks = list(range(3, n+1, 2))

    # Rank-stratify into (n+1, N) arrays
    A_r = np.zeros((nranks, N), dtype=np.int64)
    B_r = np.zeros((nranks, N), dtype=np.int64)
    for k in range(nranks):
        mk = (pc == k)
        A_r[k][mk] = cc_mod[mk]
        B_r[k][mk] = G_mod[mk]

    # Forward SOS on both A and B
    sos_forward_inplace(A_r, n, p)
    sos_forward_inplace(B_r, n, p)

    # Pointwise convolution: C_r[k] = sum_j A_r[j] * B_r[k-j]
    # Only iterate j over nonzero cc ranks (odd >= 3)
    C_r = np.zeros((nranks, N), dtype=np.int64)
    for k in range(nranks):
        for j in cc_ranks:
            if j > k:
                break
            kj = k - j
            prod = A_r[j] * B_r[kj] % p
            C_r[k] = (C_r[k] + prod) % p

    # Inverse SOS
    sos_inverse_inplace(C_r, n, p)

    # Assemble: result[T] = C_r[popcount(T)][T]
    result = np.zeros(N, dtype=np.int64)
    for k in range(nranks):
        mk = (pc == k)
        result[mk] = C_r[k][mk]

    return result

# ─── Modular run ─────────────────────────────────────────────────────────────

def run_mod(n, cc_list, p, pc):
    """Run full computation mod p. Returns {k: G_k[full] mod p}."""
    full = (1 << n) - 1
    kmax = n // 3
    N = 1 << n
    nranks = n + 1

    # Convert cc to numpy mod p
    cc = np.array(cc_list, dtype=np.int64) % p

    # G_1 = SOS(cc) mod p (using reshape trick on 1D array)
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

    return Gfull

# ─── Main ─────────────────────────────────────────────────────────────────────

def compute_all_alphas(n, adj):
    """Compute exact alpha_k for k=1..n//3 using two-prime CRT."""
    kmax = n // 3

    print(f"Computing cycle_cc for n={n}...")
    t0 = time.time()
    cc_list = cycle_cc(adj, n)
    print(f"  cycle_cc: {time.time()-t0:.1f}s")

    print(f"Precomputing popcounts...")
    t0 = time.time()
    pc = make_popcount(n)
    print(f"  done: {time.time()-t0:.1f}s")

    print(f"\n=== Modular run 1 (p={P1}) ===")
    t1 = time.time()
    r1 = run_mod(n, cc_list, P1, pc)
    print(f"  Total run 1: {time.time()-t1:.1f}s")

    print(f"\n=== Modular run 2 (p={P2}) ===")
    t2 = time.time()
    r2 = run_mod(n, cc_list, P2, pc)
    print(f"  Total run 2: {time.time()-t2:.1f}s")

    print(f"\n=== CRT recovery ===")
    alphas = {}
    for k in range(1, kmax+1):
        Gk = crt2(r1[k], r2[k], P1, P2)
        fk = factorial(k)
        if Gk % fk != 0:
            print(f"  ERROR: G_{k}={Gk} not divisible by {k}!={fk}")
            alphas[k] = Gk // fk
        else:
            alphas[k] = Gk // fk
        print(f"  α_{k} = {alphas[k]:,}")

    return alphas


if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 21

    if n == 17:
        S = [1,2,3,4,5,6,7,8]
        H_ref = 13_689_269_499
    elif n == 19:
        S = list(range(1, (n+1)//2))
        H_ref = 1_184_212_824_763
    elif n == 21:
        S = list(range(1, (n+1)//2))
        H_ref = 125_547_534_942_879
    elif n == 23:
        S = list(range(1, (n+1)//2))
        H_ref = None
    elif n == 25:
        S = list(range(1, (n+1)//2))
        H_ref = None
    else:
        S = list(range(1, (n+1)//2))
        H_ref = None

    print(f"n={n}, S={S}, kmax={n//3}")
    t_total = time.time()
    adj = circulant(n, S)
    alphas = compute_all_alphas(n, adj)
    elapsed = time.time() - t_total

    H = 1 + sum(2**k * v for k, v in alphas.items())
    print(f"\n=== FINAL RESULTS n={n} (total: {elapsed:.1f}s) ===")
    print(f"H = {H:,}")
    if H_ref:
        print(f"H_ref = {H_ref:,}")
        print(f"Match: {H == H_ref}")

    print(f"\nTerm breakdown:")
    for k in range(1, n//3+1):
        v = alphas.get(k, 0)
        pct = 100 * 2**k * v / (H - 1) if H > 1 else 0
        print(f"  2^{k} × α_{k} = {2**k * v:>30,}   ({pct:.1f}%)   α_{k}={v:,}")

    print(f"\nRatios:")
    a1 = alphas[1]; a2 = alphas[2]; a3 = alphas.get(3, 0)
    print(f"  α₁/(2α₂) = {a1/(2*a2):.6f}")
    if a3:
        print(f"  α₃/α₂    = {a3/a2:.6f}")
        print(f"  8α₃/2α₁  = {8*a3/(2*a1):.6f}  (>1 means 8α₃ > 2α₁)")
