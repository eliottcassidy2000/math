#!/usr/bin/env python3
"""
alpha_full_ssc_numpy.py -- Compute ALL alpha_k via iterated SSC, numpy + CRT.

Uses two Mersenne-type primes + CRT to avoid Python big-integer arithmetic.
All hot loops are vectorized numpy on int64 arrays.

For n=21: each SSC call ~5-10s (vs ~hours in pure Python).
Total n=21: ~5-10 minutes (cc: 65s + 6 SSC calls × 2 primes × ~10s).

Algorithm:
  G_1[T] = SOS(cc)[T]
  G_k[T] = (cc ★ G_{k-1})[T]  (Bjorklund SSC)
  alpha_k = G_k[full] / k!     (recovered via CRT of two mod-p runs)

Modular safety:
  P1 = 2^31 - 1 = 2147483647 (Mersenne prime)
  P2 = 2^31 - 19 = 2147483629 (prime)
  P1 * P2 ~ 4.6e18 > expected max(G_k[full]) ~ 24T for n=21
  Overflow check: (P1-1)^2 = 4.6e18 < 9.2e18 = int64 max  (product fits)
"""

import numpy as np
import time, sys
from math import factorial

P1 = (1 << 31) - 1    # 2147483647  (Mersenne prime)
P2 = (1 << 31) - 19   # 2147483629  (prime)

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
    """Count directed odd cycles by vertex-set bitmask. O(n * 2^n) states."""
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

# ─── Numpy helpers ────────────────────────────────────────────────────────────

def make_sos_pairs(n):
    """Precompute (sel, src) index arrays for each SOS dimension i."""
    N = 1 << n
    idx = np.arange(N)
    pairs = []
    for i in range(n):
        bit = 1 << i
        sel = np.where(idx & bit != 0)[0]
        src = sel ^ bit
        pairs.append((sel, src))
    return pairs

def make_popcount(n):
    """Precompute popcount for each index in [0, 2^n)."""
    N = 1 << n
    pc = np.zeros(N, dtype=np.int32)
    for i in range(n):
        pc += (np.arange(N, dtype=np.int32) >> i) & 1
    return pc

# ─── Core SSC mod p ──────────────────────────────────────────────────────────

def ssc_mod(cc_mod, G_mod, n, p, pairs, pc):
    """
    Compute (cc ★ G)[T] mod p for all T.

    cc_mod, G_mod: numpy int64 arrays of size 2^n, values in [0, p).
    Returns numpy int64 array, values in [0, p).

    Key correctness invariant: (P1-1)^2 < int64 max, so product A*B
    never overflows int64 before the % p reduction.
    """
    N = 1 << n

    # Non-zero ranks of cc: only ODD >= 3 (only odd-length cycles exist)
    cc_ranks = list(range(3, n+1, 2))  # 3, 5, 7, ..., up to n or n-1

    # Rank-stratify A (=cc) and B (=G) into (n+1, N) 2D arrays
    A_r = np.zeros((n+1, N), dtype=np.int64)
    B_r = np.zeros((n+1, N), dtype=np.int64)
    for k in range(n+1):
        mk = (pc == k)
        A_r[k][mk] = cc_mod[mk]
        B_r[k][mk] = G_mod[mk]

    # Forward SOS mod p on all ranks simultaneously (vectorized over rank axis)
    for sel, src in pairs:
        A_r[:, sel] = (A_r[:, sel] + A_r[:, src]) % p
        B_r[:, sel] = (B_r[:, sel] + B_r[:, src]) % p

    # Pointwise convolution mod p:  C_r[k][T] = sum_j A_r[j][T] * B_r[k-j][T]
    # Only iterate j over nonzero cc ranks (odd >= 3) for speed
    C_r = np.zeros((n+1, N), dtype=np.int64)
    for k in range(n+1):
        for j in cc_ranks:
            if j > k:
                break
            kj = k - j
            # Product fits in int64: both A_r[j], B_r[kj] < p < 2^31, product < 2^62
            prod = A_r[j] * B_r[kj] % p
            C_r[k] = (C_r[k] + prod) % p

    # Inverse SOS (Möbius) mod p
    for sel, src in pairs:
        C_r[:, sel] = (C_r[:, sel] - C_r[:, src] + p) % p

    # Assemble: result[T] = C_r[popcount(T)][T]
    result = np.zeros(N, dtype=np.int64)
    for k in range(n+1):
        mk = (pc == k)
        result[mk] = C_r[k][mk]

    return result

# ─── Modular run ─────────────────────────────────────────────────────────────

def run_mod(n, cc_list, p, pairs, pc):
    """Run full alpha computation mod p. Returns {k: G_k[full] mod p}."""
    full = (1 << n) - 1
    kmax = n // 3

    # Convert cc to numpy mod p
    cc = np.array(cc_list, dtype=np.int64) % p

    # G_1 = SOS(cc) mod p
    G = cc.copy()
    for sel, src in pairs:
        G[sel] = (G[sel] + G[src]) % p

    Gfull = {1: int(G[full])}
    print(f"    k=1: G_1[full] ≡ {Gfull[1]}")

    for k in range(2, kmax+1):
        t0 = time.time()
        G = ssc_mod(cc, G, n, p, pairs, pc)
        Gfull[k] = int(G[full])
        print(f"    k={k}: G_{k}[full] ≡ {Gfull[k]}  [{time.time()-t0:.1f}s]")

    return Gfull

# ─── Main ─────────────────────────────────────────────────────────────────────

def compute_all_alphas(n, adj):
    """Compute exact alpha_k for k=1..n//3 using two-prime CRT."""
    full = (1 << n) - 1
    kmax = n // 3

    print(f"Computing cycle_cc for n={n}...")
    t0 = time.time()
    cc_list = cycle_cc(adj, n)
    print(f"  cycle_cc: {time.time()-t0:.1f}s")

    print(f"Precomputing index structures...")
    t0 = time.time()
    pc = make_popcount(n)
    pairs = make_sos_pairs(n)
    print(f"  done: {time.time()-t0:.1f}s")

    print(f"\n=== Modular run 1 (p={P1}) ===")
    r1 = run_mod(n, cc_list, P1, pairs, pc)

    print(f"\n=== Modular run 2 (p={P2}) ===")
    r2 = run_mod(n, cc_list, P2, pairs, pc)

    print(f"\n=== CRT recovery (M = P1*P2 ≈ {P1*P2:.3e}) ===")
    alphas = {}
    for k in range(1, kmax+1):
        Gk = crt2(r1[k], r2[k], P1, P2)
        fk = factorial(k)
        if Gk % fk != 0:
            print(f"  ERROR: G_{k}={Gk} not divisible by {k}!={fk}")
            alphas[k] = Gk // fk  # floor division, flag the error
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
        print(f"  2^{k} × α_{k} = {2**k * v:>30,}   (α_{k} = {v:,})")

    print(f"\nRatios:")
    a1 = alphas[1]; a2 = alphas[2]; a3 = alphas.get(3, 0)
    print(f"  α₁/(2α₂) = {a1/(2*a2):.6f}")
    if a3:
        print(f"  α₃/α₂    = {a3/a2:.6f}")
        print(f"  8α₃/2α₁  = {8*a3/(2*a1):.6f}  (>1 means 8α₃ dominates 2α₁)")
