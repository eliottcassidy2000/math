#!/usr/bin/env python3
"""
alpha3_ssc.py -- Compute alpha_3 via Subset Sum Convolution.

Algorithm (Björklund et al. 2007):
  (A ★ B)[T] = sum_{S ⊆ T} A[S] * B[T\S]  for all T, in O(n^2 * 2^n) time.

For alpha_3:
  G_2[T] = (cc ★ f)[T] = sum_{S ⊆ T} cc[S] * f[T\S]    (ordered pairs of disjoint cycles)
  G_3[full] = sum_m cc[m] * G_2[~m & full]               (ordered triples)
  alpha_3 = G_3[full] / 6

Or equivalently via self-convolution:
  h = cc ★ cc: h[T] = sum_{S ⊆ T} cc[S] * cc[T\S]
  G_3[full] = sum_T f[T] * h[full & ~T]
  alpha_3 = G_3[full] / 6

Note: cc ★ f is more direct but requires ranking f.
We use cc ★ cc here since it's one self-convolution.

Memory: (n+1) * 2^n Python integers. For n=19: 20 * 524288 ≈ 10M ints.
Time estimate: O(n^2 * 2^n) ≈ 189M big-int operations. ~60-120s in Python.
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

def popcount(m):
    return bin(m).count('1')

def ssc_self_conv(cc, n):
    """
    Compute h = cc ★ cc where h[T] = sum_{S ⊆ T} cc[S] * cc[T\S].
    Returns h as a list of size 2^n.
    Uses ranked SOS (Björklund subset sum convolution).
    O(n^2 * 2^n) time.
    """
    N = 1 << n
    print(f"  SSC self-conv: n={n}, N={N:,}")

    # Step 1: Rank-stratify cc
    # cc_r[k][m] = cc[m] if popcount(m)=k else 0
    t0 = time.time()
    cc_r = [[0]*N for _ in range(n+1)]
    for m in range(N):
        if cc[m]:
            cc_r[popcount(m)][m] = cc[m]
    print(f"  rank-stratify: {time.time()-t0:.1f}s")

    # Step 2: Forward SOS on each rank
    t0 = time.time()
    for k in range(n+1):
        row = cc_r[k]
        if not any(row): continue  # skip empty ranks
        for i in range(n):
            bit = 1 << i
            for mask in range(bit, N):
                if mask & bit:
                    row[mask] += row[mask ^ bit]
    print(f"  forward SOS: {time.time()-t0:.1f}s")

    # Step 3: Pointwise self-convolution: h_r[k][T] = sum_{j=0}^{k} cc_r[j][T] * cc_r[k-j][T]
    t0 = time.time()
    h_r = [[0]*N for _ in range(n+1)]
    for k in range(n+1):
        hr_k = h_r[k]
        for j in range(k+1):
            cj = cc_r[j]
            ckj = cc_r[k-j]
            if not any(cj) or not any(ckj): continue
            for T in range(N):
                if cj[T] and ckj[T]:
                    hr_k[T] += cj[T] * ckj[T]
    print(f"  pointwise conv: {time.time()-t0:.1f}s")

    # Step 4: Inverse SOS (Möbius) on each rank
    t0 = time.time()
    for k in range(n+1):
        row = h_r[k]
        if not any(row): continue
        for i in range(n):
            bit = 1 << i
            for mask in range(bit, N):
                if mask & bit:
                    row[mask] -= row[mask ^ bit]
    print(f"  inverse SOS: {time.time()-t0:.1f}s")

    # Step 5: Assemble h[T] = h_r[popcount(T)][T]
    t0 = time.time()
    h = [h_r[popcount(T)][T] for T in range(N)]
    print(f"  assemble: {time.time()-t0:.1f}s")

    return h

def compute_alpha3_via_ssc(n, cc, f):
    """
    Compute alpha_3 using:
      h = cc ★ cc
      G_3[full] = sum_T f[T] * h[full & ~T]
      alpha_3 = G_3[full] / 6
    """
    full = (1 << n) - 1

    print("Computing cc ★ cc via SSC...")
    h = ssc_self_conv(cc, n)

    print("Computing G_3[full]...")
    t0 = time.time()
    G3 = sum(f[T] * h[full & ~T] for T in range(1 << n) if f[T] and h[full & ~T])
    print(f"  G_3 sum: {time.time()-t0:.1f}s")

    assert G3 % 6 == 0, f"G3={G3} not divisible by 6!"
    return G3 // 6

# ─── Main ────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 17  # default n=17 for testing

    if n == 17:
        S = [1,2,3,4,5,6,7,8]
        H_ref = 13_689_269_499
        alpha3_ref = 458_011_858
    elif n == 19:
        S = list(range(1, (n+1)//2))  # cyclic interval
        H_ref = 1_184_212_824_763
        alpha3_ref = None
    else:
        S = list(range(1, (n+1)//2))
        H_ref = None; alpha3_ref = None

    print(f"n={n}, S={S}")

    # Build tournament and compute cc, f
    adj = circulant(n, S)
    t0 = time.time()
    cc = cycle_cc(adj, n)
    f = sos(cc, n)
    print(f"cc+f: {time.time()-t0:.1f}s")

    # Verify alpha_1, alpha_2 first
    alpha1 = sum(cc)
    alpha2 = sum(cc[m] * f[(~m) & ((1<<n)-1)] for m in range(1 << n)) // 2
    print(f"α₁={alpha1:,}  α₂={alpha2:,}")

    # Main: compute alpha_3 via SSC
    t0 = time.time()
    alpha3 = compute_alpha3_via_ssc(n, cc, f)
    total_time = time.time() - t0
    print(f"\nα₃ = {alpha3:,}  (total: {total_time:.1f}s)")

    if alpha3_ref is not None:
        print(f"Expected: {alpha3_ref:,}")
        print(f"Match: {alpha3 == alpha3_ref}")

    print(f"\nRatio α₁/(2α₂) = {alpha1/(2*alpha2):.5f}")
    print(f"Ratio α₃/α₂ = {alpha3/alpha2:.5f}")
