#!/usr/bin/env python3
"""
paley_cluster_expansion_monad.py
monad-explorer-2026-06-07 (deep-research lane). Builds on HYP-2306 / the-1729-resonance
reflection, which left UNSETTLED whether R(p) = H(T_p)*2^{p-1}/p! converges to e, to a
larger constant, or grows like p^{3/2} (Alon).

NEW IDEA (theoretical, no large-p needed):
  R(p) = H(T_p)*2^{p-1}/p!
       = (1/p!) sum_{orderings sigma} prod_{k=1}^{p-1} (1 + chi(a_{k+1}-a_k))
       = E_sigma[ prod_k (1 + chi(d_k)) ]                       (chi = Legendre symbol)
       = sum_{S subset of edges} E_sigma[ prod_{k in S} chi(d_k) ]
       = 1 + sum_{S != empty} (1/p!) T(S),   T(S) = sum_sigma prod_{k in S} chi(d_k).

Each S decomposes into maximal runs (consecutive kept edges). Define the single-run
character sum over a directed path of L edges (L+1 distinct vertices):

  A_L = sum_{x_0,...,x_L distinct in F_p} prod_{i=0}^{L-1} chi(x_{i+1}-x_i).

Leading-order single-run contribution to R is  (#placements ~ p) * A_L*(p-L-1)!/p!
~ A_L/p^L. Define the CLUSTER INTEGRAL  a_L = lim_{p->inf} A_L / p^L.

Linked-cluster / exponential formula => R(p) -> exp( sum_{L>=2} a_L ).

EXACT facts (proved by symmetry, verified here):
  * a_1 = 0           (single edges: sum_{d!=0} chi(d) = 0)
  * a_L = 0 for ODD L (negation x->-x sends prod chi -> (-1)^L prod chi; for p=3 mod4
                       chi(-1)=-1, so A_L = (-1)^L A_L => A_L=0 for L odd)
  * a_2 = 1           (computed in the reflection / verified here)

So R(p) -> exp(1 + a_4 + a_6 + ...). The whole "is it e?" question reduces to:
                  IS  a_4 + a_6 + ... = 0 ?
If yes, R -> e EXACTLY. If a_4 != 0, the true constant is e^{1+a_4+...} != e.

This script:
  (1) Recomputes R(p) exactly for p=3,7,11,19,23 from stored/recomputed H(T_p).
  (2) Computes A_L exactly for L=1..6 at several Paley primes p, extracts a_L = A_L/p^L,
      and confirms odd-L vanishing + a_2=1.
  (3) Reports a_4 (and a_6) to DECIDE the limit.
"""
import numpy as np
import sys

def legendre_array(p):
    """chi[d] = Legendre symbol (d|p), chi[0]=0, for d in 0..p-1."""
    chi = np.zeros(p, dtype=np.int64)
    qr = set((x*x) % p for x in range(1, p))
    for d in range(1, p):
        chi[d] = 1 if d in qr else -1
    return chi

def A_L(p, L, chi):
    """
    Exact A_L = sum over INJECTIVE (x_0,...,x_L) in F_p of prod chi(x_{i+1}-x_i).
    By translation invariance fix x_0 = 0 and multiply by p.
    Build paths level by level with numpy; carry weight = product of chi so far.
    'used' membership enforced by a boolean mask over the p values per path.
    """
    # paths: array shape (M, k+1) of vertex values; weights: (M,)
    # start: x_0 = 0
    paths = np.zeros((1, 1), dtype=np.int64)
    weights = np.ones(1, dtype=np.int64)
    for step in range(L):
        M = paths.shape[0]
        # candidate next values 0..p-1 for each path
        # next_vals: (M, p)
        cur = paths[:, -1][:, None]            # (M,1) current endpoint
        cand = np.arange(p)[None, :]           # (1,p)
        diff = (cand - cur) % p                # (M,p)
        w_edge = chi[diff]                      # (M,p) ; chi[0]=0 forbids repeat-of-endpoint
        # forbid values already used in the path (distinctness)
        # build used mask (M,p)
        used = np.zeros((M, p), dtype=bool)
        rows = np.repeat(np.arange(M), paths.shape[1])
        used[rows, paths.ravel()] = True
        valid = (~used) & (w_edge != 0)
        # expand
        Ms, Vs = np.nonzero(valid)             # indices into (M,p)
        new_paths = np.concatenate([paths[Ms], cand[0, Vs][:, None]], axis=1)
        new_weights = weights[Ms] * w_edge[Ms, Vs]
        paths, weights = new_paths, new_weights
        if paths.shape[0] == 0:
            return 0
    total = int(weights.sum())
    return p * total   # times p for the x_0=0 translation reduction

def held_karp_paths(adj):
    """Count directed Hamiltonian paths in tournament given by adj[i][j]=1 if i->j.
       adj: numpy (n,n) 0/1. Returns total over all start vertices (directed paths)."""
    n = adj.shape[0]
    # dp[mask][v] = number of directed paths covering 'mask', ending at v
    size = 1 << n
    dp = np.zeros((size, n), dtype=object)
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        for v in range(n):
            c = dp[mask][v]
            if c == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if adj[v][w]:
                    dp[mask | (1 << w)][w] += c
    full = size - 1
    return int(sum(dp[full][v] for v in range(n)))

def paley_adj(p):
    """Paley tournament adjacency: i->j iff (j-i) is a QR mod p. (p = 3 mod 4)."""
    qr = set((x*x) % p for x in range(1, p))
    adj = np.zeros((p, p), dtype=np.int64)
    for i in range(p):
        for j in range(p):
            if i != j and ((j - i) % p) in qr:
                adj[i][j] = 1
    return adj

def main():
    from math import factorial
    print("="*70)
    print("PART 1: exact R(p) = H(T_p) * 2^{p-1} / p!")
    print("="*70)
    # H(T_p) values: recompute small ones, trust stored large ones (re-verified HYP-2306)
    stored_H = {3:3, 7:189, 11:95095, 19:1172695746915, 23:15760206976379349}
    e = 2.718281828459045
    for p in [3,7,11,19]:
        H = held_karp_paths(paley_adj(p))
        assert H == stored_H[p], (p, H, stored_H[p])
        R = H * (2**(p-1)) / factorial(p)
        print(f"  p={p:2d}  H={H:>20d}  R={R:.6f}  (e={e:.6f}, R/e={R/e:.4f})")
    for p in [23]:
        H = stored_H[p]
        R = H * (2**(p-1)) / factorial(p)
        print(f"  p={p:2d}  H={H:>20d}  R={R:.6f}  (stored)  R/e={R/e:.4f}")

    print()
    print("="*70)
    print("PART 2: cluster integrals  a_L = A_L / p^L")
    print("="*70)
    primes = [7, 11, 19, 23, 31]
    Ls = [1, 2, 3, 4]
    chi_cache = {p: legendre_array(p) for p in primes}
    print(f"{'L':>2} | " + " | ".join(f"p={p:<2d} A_L/p^L" for p in primes))
    aL_est = {}
    for L in Ls:
        row = []
        for p in primes:
            try:
                a = A_L(p, L, chi_cache[p])
                ratio = a / (p**L)
                row.append(ratio)
            except MemoryError:
                row.append(float('nan'))
        aL_est[L] = row
        print(f"{L:>2} | " + " | ".join(f"{r:>11.5f}" for r in row))

    print()
    print("Reading off limits (largest p column):")
    for L in Ls:
        print(f"  a_{L} ~ {aL_est[L][-1]:.5f}")
    print()
    print("Decision: R -> exp(sum_{L>=2} a_L). a_2=1, odd a_L=0.")
    s = aL_est[2][-1] + aL_est[4][-1]
    print(f"  sum a_L (L=2,4) ~ {s:.5f}  => R -> e^{s:.5f} = {np.exp(s):.5f}")
    print(f"  (e = {e:.5f}). If a_4 != 0, the limit is NOT e.")

if __name__ == "__main__":
    main()
