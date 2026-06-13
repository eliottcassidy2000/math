#!/usr/bin/env python3
"""
staircase_allzero_k9_monad_s8.py   monad-compute-2026-06-04-S8

INV-190 handoff (monad-researcher-S577): compute H(all-0 interleaved
staircase) at k=9 (n=18) via Held-Karp bitmask DP, extending the known
sequence H(k=2..8) = 5, 29, 233, 2489, 33773, 562685, 11222321.

Reuses the staircase construction from staircase_allzero_k7_s577.py and
validates ALL known values k=2..8 before reporting the new k=9 term, so a
build/DP regression can never silently corrupt the new value.

The all-0 interleaved staircase at n=2k:
  - Pairs: (0,1), (2,3), ..., (2k-2, 2k-1)
  - Global ranking: rank[2p]=p (dominants), rank[2p+1]=k+p (recessives)
  - Within-pair: odd (recessive) beats even (dominant): 2p+1 -> 2p
  - Between pairs: lower rank beats higher: rank[i]<rank[j] => i->j
"""

import time


def build_staircase_adj(k):
    """Return adjacency matrix A where A[i][j]=1 iff i->j, plus n."""
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
                if rank[i] < rank[j]:            # lower rank beats higher
                    A[i][j] = 1
    return A, n


def held_karp_H(A, n):
    """
    Count Hamiltonian paths in tournament via Held-Karp bitmask DP.
    dp[mask][v] = #paths ending at v covering exactly the vertices in mask.
    For speed: precompute out-neighbour bitmask per vertex and push forward.
    """
    FULL = (1 << n) - 1
    # out[v] = bitmask of vertices u with v->u
    out = [0] * n
    for v in range(n):
        m = 0
        Av = A[v]
        for u in range(n):
            if Av[u]:
                m |= (1 << u)
        out[v] = m

    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1

    for mask in range(1, 1 << n):
        row = dp[mask]
        rest = FULL ^ mask
        for v in range(n):
            cnt = row[v]
            if cnt == 0:
                continue
            # candidate next vertices: out-neighbours of v not yet used
            cand = out[v] & rest
            while cand:
                lsb = cand & (-cand)
                u = lsb.bit_length() - 1
                dp[mask | lsb][u] += cnt
                cand ^= lsb

    full_row = dp[FULL]
    return sum(full_row[v] for v in range(n))


def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k2 in range(j + 1, n):
                if ((A[i][j] and A[j][k2] and A[k2][i]) or
                        (A[i][k2] and A[k2][j] and A[j][i])):
                    c3 += 1
    return c3


def score_seq(A, n):
    return sorted(sum(A[v]) for v in range(n))


print("=" * 64)
print("INV-190 all-0 interleaved staircase: H(k=9, n=18) via Held-Karp")
print("monad-compute-2026-06-04-S8")
print("=" * 64)

known = {2: 5, 3: 29, 4: 233, 5: 2489, 6: 33773, 7: 562685, 8: 11222321}

results = []
all_ok = True
for k in range(2, 10):   # k=2..9 (n=4..18)
    n = 2 * k
    t0 = time.time()
    A, _ = build_staircase_adj(k)
    H = held_karp_H(A, n)
    dt = time.time() - t0
    c3 = count_3cycles(A, n)
    sc = score_seq(A, n)

    exp = known.get(k)
    if exp is None:
        tag = "NEW"
    elif H == exp:
        tag = "match"
    else:
        tag = f"MISMATCH(expected {exp})"
        all_ok = False

    c3ok = (c3 == k * (k - 1))
    print(f"\nk={k}, n={n}: H={H}  [{tag}]  t={dt:.2f}s")
    print(f"  c3={c3} (k(k-1)={k*(k-1)}, {'OK' if c3ok else 'FAIL'})")
    print(f"  score_seq={sc}")
    results.append((k, n, H, c3))

print("\n" + "=" * 64)
print("VALIDATION")
print("=" * 64)
print(f"Known values k=2..8 reproduced: {'YES' if all_ok else 'NO'}")

H_vals = [r[2] for r in results]
print(f"\nFull sequence H(k=2..9): {H_vals}")

new = [r for r in results if r[0] == 9]
print("\n" + "=" * 64)
print("KEY RESULT")
print("=" * 64)
for k, n, H, c3 in new:
    status = "VALID" if all_ok else "SUSPECT (validation failed)"
    print(f"H(all-0 staircase, k={k}, n={n}) = {H}   [{status}]")
    print(f"  c3 = {c3} (k(k-1) = {k*(k-1)}, {'OK' if c3 == k*(k-1) else 'FAIL'})")

# Growth ratios (asymptotic eyeball)
print("\nGrowth ratios H(k)/H(k-1):")
for i in range(1, len(H_vals)):
    print(f"  k={results[i][0]}: {H_vals[i] / H_vals[i-1]:.6f}")

# Quick no-recurrence sanity: order-3 fit on first terms, test on tail
print("\nOrder-3 linear recurrence check (fit on k=2..7, test k=8,9):")
from fractions import Fraction
H = H_vals
M = [
    [Fraction(H[2]), Fraction(H[1]), Fraction(H[0])],
    [Fraction(H[3]), Fraction(H[2]), Fraction(H[1])],
    [Fraction(H[4]), Fraction(H[3]), Fraction(H[2])],
]
b = [Fraction(H[3]), Fraction(H[4]), Fraction(H[5])]
for col in range(3):
    piv = next((r for r in range(col, 3) if M[r][col] != 0), None)
    M[col], M[piv] = M[piv], M[col]
    b[col], b[piv] = b[piv], b[col]
    for r in range(3):
        if r != col and M[r][col] != 0:
            f = M[r][col] / M[col][col]
            for c2 in range(3):
                M[r][c2] -= f * M[col][c2]
            b[r] -= f * b[col]
sol = [b[r] / M[r][r] for r in range(3)]
print(f"  p={sol[0]}, q={sol[1]}, r={sol[2]}")
ok = all(sol[0]*H[i-1] + sol[1]*H[i-2] + sol[2]*H[i-3] == H[i]
         for i in range(3, len(H)))
print(f"  Holds for all terms incl. k=8,9: {ok}")
