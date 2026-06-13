#!/usr/bin/env python3
"""
staircase_allzero_k12_s_monad.py   monad-researcher-2026-06-02-S578

INV-190 / HYP-2095 extension: compute H(all-0 interleaved staircase) at k=12 (n=24)
using memory-efficient Held-Karp DP with array.array('q') (signed 64-bit integers).

For k=12, n=24: dp array has 2^24 * 24 = 402,653,184 entries.
Using array.array('q') instead of a Python list saves ~4x memory
(8 bytes per entry vs 8 bytes pointer + Python object overhead for large ints).

The dp values (path counts for sub-tournaments) fit comfortably in int64:
H(k=12) ≈ 2.19e11 * 34.7 ≈ 7.6e12, well within int64 max (9.2e18).

Known sequence k=2..11:
  5, 29, 233, 2489, 33773, 562685, 11222321, 262755369, 7110764837, 219612027389
"""

import sys
import time
from array import array


def build_staircase_adj(k):
    n = 2 * k
    rank = {}
    for p in range(k):
        rank[2 * p] = p
        rank[2 * p + 1] = k + p

    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            same_pair = (i // 2 == j // 2)
            if same_pair:
                if i % 2 == 1 and j % 2 == 0:
                    A[i][j] = 1
            else:
                if rank[i] < rank[j]:
                    A[i][j] = 1

    outmask = [0] * n
    for v in range(n):
        m = 0
        for u in range(n):
            if A[v][u]:
                m |= (1 << u)
        outmask[v] = m
    return A, n, outmask


def held_karp_H_array(n, outmask):
    """Count Hamiltonian paths using array.array('q') for memory efficiency.
    Uses signed 64-bit integers -- safe because H(k=12) << 2^63.
    """
    full = 1 << n
    size = full * n
    print(f"  Allocating dp array: {size:,} entries = {size * 8 / 1e9:.2f} GB", flush=True)
    t_alloc = time.time()
    dp = array('q', [0]) * size  # compact signed int64 array
    print(f"  Allocation done in {time.time() - t_alloc:.1f}s", flush=True)

    for v in range(n):
        dp[(1 << v) * n + v] = 1

    fullmask = full - 1
    t_last = time.time()
    progress_interval = max(1, full // 20)

    for mask in range(1, full):
        if mask % progress_interval == 0:
            pct = 100 * mask // full
            elapsed = time.time() - t_last
            print(f"  Progress: {pct}% (mask={mask:#010x}), last {progress_interval} masks: {elapsed:.1f}s",
                  flush=True)
            t_last = time.time()

        base = mask * n
        not_mask = ~mask
        m = mask
        while m:
            lsbv = m & (-m)
            v = lsbv.bit_length() - 1
            m ^= lsbv
            val = dp[base + v]
            if val:
                avail = outmask[v] & not_mask
                while avail:
                    lsb = avail & (-avail)
                    avail ^= lsb
                    u = lsb.bit_length() - 1
                    dp[(mask | lsb) * n + u] += val

    fbase = fullmask * n
    return sum(dp[fbase + v] for v in range(n))


def score_seq(A, n):
    return sorted(sum(A[v]) for v in range(n))


def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i + 1, n):
            for r in range(j + 1, n):
                if ((A[i][j] and A[j][r] and A[r][i]) or
                        (A[i][r] and A[r][j] and A[j][i])):
                    c3 += 1
    return c3


def main():
    known = {2: 5, 3: 29, 4: 233, 5: 2489, 6: 33773,
             7: 562685, 8: 11222321, 9: 262755369,
             10: 7110764837, 11: 219612027389}

    kmax = 12
    if len(sys.argv) > 1:
        kmax = int(sys.argv[1])
    kmin = 12
    if len(sys.argv) > 2:
        kmin = int(sys.argv[2])

    print("=" * 64)
    print("All-0 Interleaved Staircase: H via array.array Held-Karp DP")
    print(f"Running k={kmin}..{kmax}")
    print("=" * 64)

    H_vals = []
    for k in range(kmin, kmax + 1):
        n = 2 * k
        print(f"\nk={k}, n={n}", flush=True)
        t0 = time.time()
        A, _, outmask = build_staircase_adj(k)
        H = held_karp_H_array(n, outmask)
        elapsed = time.time() - t0
        sc = score_seq(A, n)
        c3 = count_3cycles(A, n)
        exp = known.get(k, "NEW")
        match = (H == exp) if isinstance(exp, int) else "NEW"
        print(f"\nk={k:2d}, n={n:2d}: H={H}")
        print(f"   expected={exp}  match={match}  total_t={elapsed:.2f}s")
        print(f"   score_seq={sc}")
        print(f"   c3={c3}  k(k-1)={k*(k-1)}  match={c3 == k*(k-1)}")
        H_vals.append((k, H, c3, elapsed))
        sys.stdout.flush()

    print("\n" + "=" * 64)
    print("FULL SEQUENCE k=2..11 (known) + new:")
    all_known = [(k, v) for k, v in sorted(known.items())]
    for k, v in all_known:
        print(f"  k={k:2d}: H={v}")
    for k, H, c3, t in H_vals:
        if k not in known:
            print(f"  k={k:2d}: H={H}  [NEW, t={t:.2f}s]")

    if H_vals:
        print("\nGrowth ratios (all known + new):")
        prev = None
        for k in sorted(list(known.keys()) + [v[0] for v in H_vals]):
            h = known.get(k) or next(v[1] for v in H_vals if v[0] == k)
            if prev is not None:
                print(f"  k={k:2d}: {h / prev:.6f}")
            prev = h


if __name__ == "__main__":
    main()
