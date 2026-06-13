#!/usr/bin/env python3
"""
staircase_allzero_k9_s_monad.py   monad-compute-2026-06-02

INV-190 / HYP-2095 extension: compute H(all-0 interleaved staircase) at
k=9 (n=18) and attempt k=10 (n=20) using an optimized Held-Karp bitmask DP.

Extends staircase_allzero_k7_s577.py (which reached k=8, n=16).
Known sequence k=2..8:  5, 29, 233, 2489, 33773, 562685, 11222321.

The all-0 interleaved staircase at n=2k:
  - Pairs: (0,1), (2,3), ..., (2k-2, 2k-1)
  - Global ranking: rank[2p]=p (dominants), rank[2p+1]=k+p (recessives)
  - Within-pair: odd (recessive) beats even (dominant): 2p+1 -> 2p
  - Between pairs: lower rank beats higher: rank[i]<rank[j] => i->j

Optimization vs the original:
  - flat dp array dp[mask*n + v] (avoids list-of-lists overhead)
  - precomputed out-neighbor bitmask per vertex; iterate available bits
  - masks visited in increasing order (adding a bit always increases mask)
"""

import sys
import time


def build_staircase_adj(k):
    """Return (A, n, outmask) where A[i][j]=1 iff i->j; outmask[v] is the
    bitmask of out-neighbors of v."""
    n = 2 * k
    rank = {}
    for p in range(k):
        rank[2 * p] = p          # dominant vertices
        rank[2 * p + 1] = k + p  # recessive vertices

    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            same_pair = (i // 2 == j // 2)
            if same_pair:
                if i % 2 == 1 and j % 2 == 0:   # recessive beats dominant
                    A[i][j] = 1
            else:
                if rank[i] < rank[j]:           # lower rank beats higher
                    A[i][j] = 1

    outmask = [0] * n
    for v in range(n):
        m = 0
        for u in range(n):
            if A[v][u]:
                m |= (1 << u)
        outmask[v] = m
    return A, n, outmask


def held_karp_H(n, outmask):
    """Count Hamiltonian paths in a tournament via flat-array Held-Karp DP.
    dp[mask*n + v] = #paths ending at v covering exactly the vertices in mask.
    Time O(2^n * n^2), space O(2^n * n)."""
    size = (1 << n) * n
    dp = [0] * size
    for v in range(n):
        dp[(1 << v) * n + v] = 1

    full = 1 << n
    for mask in range(1, full):
        base = mask * n
        not_mask = ~mask
        m = mask
        # iterate set bits of mask
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

    fullmask = full - 1
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
    print("=" * 64)
    print("All-0 Interleaved Staircase: H via optimized Held-Karp DP")
    print("monad-compute extension to k=9 (n=18) [+ k=10 attempt]")
    print("=" * 64)

    known = {2: 5, 3: 29, 4: 233, 5: 2489, 6: 33773, 7: 562685, 8: 11222321}

    # which k to run: default 2..9; pass an int arg for max k (e.g. 10)
    kmax = 9
    if len(sys.argv) > 1:
        kmax = int(sys.argv[1])

    H_vals = []
    for k in range(2, kmax + 1):
        n = 2 * k
        t0 = time.time()
        A, _, outmask = build_staircase_adj(k)
        H = held_karp_H(n, outmask)
        elapsed = time.time() - t0
        sc = score_seq(A, n)
        c3 = count_3cycles(A, n)
        exp = known.get(k, "NEW")
        match = (H == exp) if isinstance(exp, int) else "NEW"
        print(f"\nk={k:2d}, n={n:2d}: H={H}")
        print(f"   expected={exp}  match={match}  t={elapsed:.2f}s")
        print(f"   score_seq={sc}")
        print(f"   c3={c3}  k(k-1)={k*(k-1)}  match={c3 == k*(k-1)}")
        H_vals.append((k, H, c3))
        sys.stdout.flush()

    print("\n" + "=" * 64)
    print("SEQUENCE")
    print("=" * 64)
    print("k    : " + ", ".join(str(k) for k, _, _ in H_vals))
    print("H(k) : " + ", ".join(str(h) for _, h, _ in H_vals))
    print("c3   : " + ", ".join(str(c) for _, _, c in H_vals))

    # growth ratios
    print("\nGrowth ratios H(k)/H(k-1):")
    for i in range(1, len(H_vals)):
        k = H_vals[i][0]
        print(f"  k={k:2d}: {H_vals[i][1] / H_vals[i-1][1]:.6f}")

    # Markov check on consecutive triples
    print("\nMarkov x^2+y^2+z^2 - 3xyz for consecutive triples:")
    hs = [h for _, h, _ in H_vals]
    for i in range(len(hs) - 2):
        x, y, z = hs[i], hs[i+1], hs[i+2]
        ks = (H_vals[i][0], H_vals[i+1][0], H_vals[i+2][0])
        print(f"  k={ks}: {x*x + y*y + z*z - 3*x*y*z}")

    print("\nKEY NEW VALUES (k>=9):")
    for k, h, c3 in H_vals:
        if k >= 9:
            print(f"  H(all-0 staircase, k={k}, n={2*k}) = {h}   c3={c3} (k(k-1)={k*(k-1)})")


if __name__ == "__main__":
    main()
