#!/usr/bin/env python3
"""
alpha_full_ssc.py -- Compute ALL alpha_k via iterated Subset Sum Convolution.

G_1[T] = f[T] = SOS(cc)[T]
G_k[T] = (cc ★ G_{k-1})[T] = sum_{S ⊆ T} cc[S] * G_{k-1}[T\S]
alpha_k = G_k[full] / k!

Each SSC call is O(n^2 * 2^n). For n=19 kmax=6: ~6 × 16s = ~96s.

Iterative SSC avoids the need to enumerate combinatorial types like
[3,3,5,5] or [3,3,3,7] separately — handles all in one formula.
"""

import time, sys

def circulant(n, S):
    adj = [0]*n
    for v in range(n):
        for k in S:
            adj[v] |= 1 << ((v+k) % n)
    return adj

def cycle_cc(adj, n):
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
                    key = (mask|ub, u)
                    nq[key] = nq.get(key, 0) + cnt
                    cands ^= ub
            queue = nq
    return cc

def sos(f, n):
    f = list(f)
    for i in range(n):
        for mask in range(1 << n):
            if (mask >> i) & 1:
                f[mask] += f[mask ^ (1 << i)]
    return f

def ssc(A, B, n):
    """
    Compute (A ★ B)[T] = sum_{S ⊆ T} A[S] * B[T\S] for all T.
    Björklund subset sum convolution. O(n^2 * 2^n).
    """
    N = 1 << n

    # Rank-stratify A and B
    A_r = [[0]*N for _ in range(n+1)]
    B_r = [[0]*N for _ in range(n+1)]
    pc = [bin(m).count('1') for m in range(N)]  # precomputed popcounts
    for m in range(N):
        if A[m]: A_r[pc[m]][m] = A[m]
        if B[m]: B_r[pc[m]][m] = B[m]

    # Forward SOS on each rank
    for k in range(n+1):
        ar = A_r[k]; br = B_r[k]
        any_a = any(ar); any_b = any(br)
        for i in range(n):
            bit = 1 << i
            if any_a:
                for mask in range(bit, N):
                    if mask & bit: ar[mask] += ar[mask ^ bit]
            if any_b:
                for mask in range(bit, N):
                    if mask & bit: br[mask] += br[mask ^ bit]

    # Pointwise convolution: C_r[k][T] = sum_{j=0}^{k} A_r[j][T] * B_r[k-j][T]
    C_r = [[0]*N for _ in range(n+1)]
    for k in range(n+1):
        cr = C_r[k]
        for j in range(k+1):
            aj = A_r[j]; bkj = B_r[k-j]
            # Skip if either row is all-zero
            if not any(aj) or not any(bkj): continue
            for T in range(N):
                if aj[T] and bkj[T]:
                    cr[T] += aj[T] * bkj[T]

    # Inverse SOS (Möbius) on each rank
    for k in range(n+1):
        cr = C_r[k]
        if not any(cr): continue
        for i in range(n):
            bit = 1 << i
            for mask in range(bit, N):
                if mask & bit: cr[mask] -= cr[mask ^ bit]

    # Assemble: (A ★ B)[T] = C_r[pc[T]][T]
    result = [C_r[pc[T]][T] for T in range(N)]
    return result

def compute_all_alphas(n, adj):
    """Compute alpha_1, ..., alpha_kmax using iterated SSC."""
    from math import factorial
    full = (1 << n) - 1
    kmax = n // 3

    print(f"Computing cc and f...")
    t0 = time.time()
    cc = cycle_cc(adj, n)
    f = sos(cc, n)
    print(f"  cc+f: {time.time()-t0:.1f}s")

    alpha1 = sum(cc)
    alpha2 = sum(cc[m] * f[(~m) & full] for m in range(1 << n)) // 2
    print(f"  α₁={alpha1:,}  α₂={alpha2:,}")

    # Iterative SSC: G_k = cc ★ G_{k-1}, starting from G_1 = f
    alphas = {1: alpha1, 2: alpha2}
    G = f  # G_1

    for k in range(2, kmax+1):
        print(f"Computing G_{k} = cc ★ G_{k-1}...")
        t0 = time.time()
        G = ssc(cc, G, n)
        alpha_k = G[full] // factorial(k)
        alphas[k] = alpha_k
        print(f"  α_{k} = {alpha_k:,}  ({time.time()-t0:.1f}s)")

    # Verify H
    H_ocf = 1 + sum(2**k * v for k, v in alphas.items())
    return alphas, H_ocf

if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 17

    # Tournament selection
    if n == 17:
        S = [1,2,3,4,5,6,7,8]
        H_ref = 13_689_269_499
    elif n == 19:
        S = list(range(1, (n+1)//2))
        H_ref = 1_184_212_824_763
    elif n == 11:
        S = [1,3,4,5,9]  # Paley T_11
        H_ref = 95_095
    elif n == 13:
        S = [1,3,5,7,9,11]
        H_ref = 3_711_175
    elif n == 15:
        S = [2,4,6,8,10,12,14]
        H_ref = 198_464_295
    else:
        S = list(range(1, (n+1)//2))
        H_ref = None

    print(f"n={n}, S={S}, kmax={n//3}")

    adj = circulant(n, S)
    alphas, H_ocf = compute_all_alphas(n, adj)

    print(f"\n=== RESULTS n={n} ===")
    print(f"H (OCF) = {H_ocf:,}")
    if H_ref:
        print(f"H (ref) = {H_ref:,}")
        print(f"Match: {H_ocf == H_ref}")

    from math import factorial
    print(f"\nTerm breakdown:")
    for k in range(1, n//3+1):
        v = alphas.get(k, 0)
        print(f"  2^{k} * α_{k} = {2**k * v:>25,}  (α_{k}={v:,})")

    print(f"\nRatios:")
    a1, a2 = alphas[1], alphas[2]
    a3 = alphas.get(3, 0)
    print(f"  α₁/(2α₂) = {a1/(2*a2):.5f}")
    if a3:
        print(f"  α₃/α₂ = {a3/a2:.5f}")
        print(f"  α₃/α₁ = {a3/a1:.5f}")
