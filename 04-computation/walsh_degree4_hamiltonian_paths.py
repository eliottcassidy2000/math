#!/usr/bin/env python3
"""
Walsh Degree-4 Structure at n=5: Hamiltonian Paths as Walsh Monomials

At n=5, the degree-4 Walsh component D4 of H(T) has exactly 60 nonzero
coefficients, all +/-1/8. These correspond to the 60 Hamiltonian paths on K_5.

The sign of each coefficient encodes the "orientation parity" of the path:
paths with an even number of "backward" edges (relative to the natural
ordering i < j) get +1/8, and odd get -1/8.

This script verifies the structure exhaustively and connects it to the
"4D hypercube" perspective: tournaments are points on {0,1}^10, and
H(T) is a degree-4 polynomial on this hypercube whose degree-4 part
is supported on exactly the 60 Hamiltonian paths.

The degree-2 Walsh structure (Petersen complement) and degree-4 structure
(Hamiltonian paths) together fully determine H(T), since H has Walsh
degree <= 4.

opus-2026-03-07-S35
"""

import numpy as np
from itertools import combinations, permutations

def bits_to_adj(bits, n):
    A = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> k) & 1: A[i][j] = 1
            else: A[j][i] = 1
            k += 1
    return A

def count_H(A, n):
    dp = {}
    for v in range(n): dp[(1<<v, v)] = 1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask & (1<<v)): continue
            c = dp.get((mask, v), 0)
            if c == 0: continue
            for u in range(n):
                if mask & (1<<u): continue
                if A[v][u]:
                    key = (mask|(1<<u), u)
                    dp[key] = dp.get(key, 0) + c
    return sum(dp.get(((1<<n)-1, v), 0) for v in range(n))

def fast_walsh_hadamard(f, m):
    N = 2**m
    a = f.copy()
    for i in range(m):
        step = 1 << i
        for j in range(0, N, step*2):
            for k in range(j, j + step):
                u, v = a[k], a[k + step]
                a[k] = u + v
                a[k + step] = u - v
    return a / N


def main():
    n = 5
    m = n*(n-1)//2  # 10
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    N = 2**m

    # Compute H for all tournaments
    H_all = np.zeros(N)
    for bits in range(N):
        A = bits_to_adj(bits, n)
        H_all[bits] = count_H(A, n)

    # Walsh transform
    H_hat = fast_walsh_hadamard(H_all.copy(), m)

    # Degree-4 analysis
    print(f"{'='*60}")
    print(f"DEGREE-4 WALSH STRUCTURE AT n={n}")
    print(f"{'='*60}")

    deg4_monomials = []
    for S in range(N):
        if bin(S).count('1') != 4: continue
        if abs(H_hat[S]) < 1e-10: continue
        bits_list = [k for k in range(m) if (S >> k) & 1]
        es = [edges[k] for k in bits_list]
        verts = set()
        for e in es: verts.update(e)
        deg_map = {}
        for v in verts:
            deg_map[v] = sum(1 for e in es if v in set(e))
        deg_seq = tuple(sorted(deg_map.values(), reverse=True))
        is_path = (deg_seq == (2,2,2,1,1) and len(verts) == 5)
        deg4_monomials.append({
            'S': S, 'edges': es, 'verts': verts, 'coeff': H_hat[S],
            'is_path': is_path, 'deg_seq': deg_seq
        })

    print(f"\nNonzero degree-4 Walsh coefficients: {len(deg4_monomials)}")
    print(f"All are Hamiltonian paths: {all(m['is_path'] for m in deg4_monomials)}")

    pos_count = sum(1 for m in deg4_monomials if m['coeff'] > 0)
    neg_count = sum(1 for m in deg4_monomials if m['coeff'] < 0)
    coeff_values = set(round(m['coeff'], 8) for m in deg4_monomials)
    print(f"Coefficient values: {coeff_values}")
    print(f"Positive (+1/8): {pos_count}")
    print(f"Negative (-1/8): {neg_count}")

    # Sign rule for degree-4 coefficients
    print(f"\n{'='*60}")
    print(f"SIGN RULE FOR DEGREE-4 COEFFICIENTS")
    print(f"{'='*60}")

    # PROVED: H_hat[S] = (-1)^{descents(path)} / 8
    # where descents = number of edges traversed from larger to smaller vertex
    # when walking the path starting from its smaller endpoint.
    # This is well-defined because 4 is even, so both traversal directions
    # give the same descent parity.

    all_match = True
    for d4 in deg4_monomials:
        es = d4['edges']
        degree = {}
        for e in es:
            for v in e:
                degree[v] = degree.get(v, 0) + 1
        endpoints = [v for v, d in degree.items() if d == 1]
        path = [min(endpoints)]
        used = set()
        while len(path) < 5:
            curr = path[-1]
            for e in es:
                if e in used: continue
                if curr in e:
                    other = e[0] if e[1] == curr else e[1]
                    path.append(other)
                    used.add(e)
                    break
        descents = sum(1 for i in range(4) if path[i] > path[i+1])
        expected = (-1)**descents / 8
        if abs(d4['coeff'] - expected) > 1e-10:
            all_match = False
            print(f"  MISMATCH: path={path}, coeff={d4['coeff']}, expected={expected}")

    print(f"  Sign = (-1)^descents / 8: {'CONFIRMED (all 60 paths)' if all_match else 'FAILED'}")
    print(f"  Even descents (coeff +1/8): {pos_count}")
    print(f"  Odd descents (coeff -1/8): {neg_count}")

    # Summary of hypercube interpretation
    print(f"\n{'='*60}")
    print(f"HYPERCUBE INTERPRETATION")
    print(f"{'='*60}")
    print(f"""
Tournaments on {n} vertices = points on the Boolean hypercube {{0,1}}^{m}.

H(T) is a multilinear polynomial of degree 4 on this hypercube:
  H(x) = 7.5 + sum_{{|S|=2}} H_hat[S]*chi_S(x) + sum_{{|S|=4}} H_hat[S]*chi_S(x)

where chi_S(x) = prod_{{k in S}} (2*x_k - 1) is the Walsh character.

The degree-2 monomials (30 of them, coeff +/-3/4) form the edge set
of L(K_5), the complement of the Petersen graph. They encode how
adjacent edge pairs (blue lines) correlate with H.

The degree-4 monomials (60 of them, coeff +/-1/8) are the Hamiltonian
paths on K_5. They encode the deepest correlation structure of H —
how 4-edge sets interact to determine the Hamiltonian path count.

The Walsh-Fourier diagonalization (THM-066) shows that:
  - The degree-2 part comes purely from t3 (the 3-cycle count)
  - The degree-4 part comes from the combination 1-t3+2*t5
  - These are the "Fourier coefficients" w_4/16 and w_0 respectively
  - The telescoping identity t5_hat = t3_hat/2 prevents cross-contamination
""")

if __name__ == '__main__':
    main()
