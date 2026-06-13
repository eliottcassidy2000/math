#!/usr/bin/env python3
"""
ssc_lowmem.py -- Memory-efficient SSC for large n.

Standard SSC requires O(n * 2^n) memory for the stratified arrays,
which is ~20 GB for n=25.

This implementation recomputes A/B stratification rows on-the-fly,
requiring only O(2^n) peak memory (~1 GB for n=25).

Time cost: O(kmax^2 * n * 2^n) -- about 10-60 minutes for n=25.

Usage:
  python3 ssc_lowmem.py <n> [circulant|paley]

Reads cc_n<n>_<type>_p<P1>.i32 and cc_n<n>_<type>_p<P2>.i32.
Outputs alpha_k values via CRT.
"""

import numpy as np
import time, sys
from math import factorial

P1 = 2147483647    # 2^31 - 1
P2 = 2147483629    # 2^31 - 19

def ext_gcd(a, b):
    if b == 0: return a, 1, 0
    g, x, y = ext_gcd(b, a % b)
    return g, y, x - (a // b) * y

def crt2(r1, r2, m1, m2):
    _, u, _ = ext_gcd(m1, m2)
    return (r1 + m1 * ((r2 - r1) * u % m2)) % (m1 * m2)

def make_pc(n):
    N = 1 << n
    pc = np.zeros(N, dtype=np.int32)
    for i in range(n):
        pc += (np.arange(N, dtype=np.int32) >> i) & 1
    return pc

def sos_fwd(a, n, p):
    """Forward SOS transform in-place (a has shape (N,), int64 mod p)."""
    N = 1 << n
    for i in range(n):
        bit = 1 << i
        a_4d = a.reshape(N // (2*bit), 2, bit)
        a_4d[:, 1, :] = (a_4d[:, 1, :] + a_4d[:, 0, :]) % p

def sos_inv(a, n, p):
    """Inverse SOS transform in-place."""
    N = 1 << n
    for i in range(n):
        bit = 1 << i
        a_4d = a.reshape(N // (2*bit), 2, bit)
        a_4d[:, 1, :] = (a_4d[:, 1, :] - a_4d[:, 0, :] + p) % p

def stratify_sos(arr, rank, n, p, pc):
    """Return SOS-forward-transformed version of arr restricted to given rank."""
    N = 1 << n
    out = np.zeros(N, dtype=np.int64)
    mask = (pc == rank)
    out[mask] = arr[mask]
    sos_fwd(out, n, p)
    return out

def modinv(a, p):
    """Modular inverse of a mod p (p prime)."""
    return pow(int(a), int(p)-2, int(p))

def run_ssc_lowmem(n, cc_raw, p, pc, kmax):
    """
    Memory-efficient SSC mod p.
    Returns {k: G_k[full_mask]} for k=1..kmax.

    Peak memory: O(2^n) per array — about 128 MB each at n=25.
    At most 5 arrays in memory simultaneously.
    """
    N = 1 << n
    full = N - 1
    cc_ranks = list(range(3, n+1, 2))  # odd cycle lengths: 3,5,...,n

    cc = (cc_raw % p).astype(np.int64)

    # G_1 = SOS(cc)
    G = cc.copy()
    sos_fwd(G, n, p)
    Gfull = {1: int(G[full])}
    print(f"    k=1: G_1[full] ≡ {Gfull[1]}")
    sys.stdout.flush()

    for step in range(2, kmax+1):
        t0 = time.time()

        # Compute G_step[mask] = C_r[popcount(mask)][mask] / step
        # where C_r[r][mask] = sum_{j in cc_ranks, j ≤ r, r-j ≥ 0}
        #                        A_sos[j][mask] * B_sos[r-j][mask]

        # Active output ranks: r such that some (j, r-j) gives nonzero contribution
        # j ∈ cc_ranks (odd ≥ 3), r-j ≥ 0 → r ≥ 3
        # min r for G_step: 3*step (need step disjoint 3-cycles minimum)
        # max r: n
        active_out_ranks = [r for r in range(3*step, n+1)]

        # For extracting G_step, we accumulate G_new_r for each rank r
        # then combine: G_new[mask] = G_new_r[mask] for mask with popcount==r
        G_new = np.zeros(N, dtype=np.int64)

        for r in active_out_ranks:
            # C_sos_r = sum over j in cc_ranks, j ≤ r, r-j ≥ 0
            C_r = np.zeros(N, dtype=np.int64)
            for j in cc_ranks:
                if j > r:
                    break
                b_rank = r - j
                # Compute A_sos[j] on-the-fly
                A_j = stratify_sos(cc, j, n, p, pc)
                # Compute B_sos[b_rank] on-the-fly from current G
                B_brank = stratify_sos(G, b_rank, n, p, pc)
                # Accumulate
                C_r = (C_r + A_j * B_brank) % p
                del A_j, B_brank

            # Inverse SOS
            sos_inv(C_r, n, p)

            # Extract into G_new at rank r, divided by step
            inv_step = modinv(step, p)
            mask_r = (pc == r)
            G_new[mask_r] = C_r[mask_r] * inv_step % p
            del C_r

        G = G_new
        Gfull[step] = int(G[full])
        elapsed = time.time() - t0
        print(f"    k={step}: G_{step}[full] ≡ {Gfull[step]}  [{elapsed:.1f}s]")
        sys.stdout.flush()

    return Gfull

if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 25
    type_name = sys.argv[2] if len(sys.argv) > 2 else 'circulant'
    kmax = n // 3

    base = f"05-knowledge/results/cc_n{n}_{type_name}"
    f1 = f"{base}_p{P1}.i32"
    f2 = f"{base}_p{P2}.i32"
    print(f"n={n}, type={type_name}, kmax={kmax}")
    print(f"Loading {f1} and {f2} ...")

    t0 = time.time()
    cc1 = np.fromfile(f1, dtype=np.uint32).astype(np.int64)
    cc2 = np.fromfile(f2, dtype=np.uint32).astype(np.int64)
    assert len(cc1) == (1 << n), f"expected {1<<n} entries"
    print(f"  loaded: {time.time()-t0:.1f}s")

    print("Precomputing popcounts...")
    t0 = time.time()
    pc = make_pc(n)
    print(f"  done: {time.time()-t0:.1f}s")
    sys.stdout.flush()

    print(f"\n=== Modular run 1 (p={P1}) ===")
    t1 = time.time()
    r1 = run_ssc_lowmem(n, cc1, P1, pc, kmax)
    print(f"  Total run 1: {time.time()-t1:.1f}s")
    sys.stdout.flush()

    print(f"\n=== Modular run 2 (p={P2}) ===")
    t2 = time.time()
    r2 = run_ssc_lowmem(n, cc2, P2, pc, kmax)
    print(f"  Total run 2: {time.time()-t2:.1f}s")
    sys.stdout.flush()

    print(f"\n=== CRT recovery ===")
    alphas = {}
    for k in range(1, kmax+1):
        Gk = crt2(r1[k], r2[k], P1, P2)
        alphas[k] = Gk // factorial(k)
        print(f"  α_{k} = {alphas[k]:,}")

    H = 1 + sum(2**k * v for k, v in alphas.items())
    print(f"\n=== FINAL RESULTS n={n} type={type_name} ===")
    print(f"H = {H:,}")
    for k in range(1, kmax+1):
        v = alphas.get(k, 0)
        pct = 100 * 2**k * v / (H - 1) if H > 1 else 0
        print(f"  2^{k} × α_{k} = {2**k*v:>35,}  ({pct:.1f}%)  α_{k}={v:,}")

    a1, a2 = alphas[1], alphas[2]
    a3 = alphas.get(3, 0)
    print(f"\nRatios:  α₁/(2α₂) = {a1/(2*a2):.6f}", end="")
    if a3: print(f"   α₃/α₂ = {a3/a2:.6f}", end="")
    print()
