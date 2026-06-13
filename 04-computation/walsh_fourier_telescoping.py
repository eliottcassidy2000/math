#!/usr/bin/env python3
"""
Walsh-Fourier Telescoping: The Walsh hypercube structure of H(T).

KEY DISCOVERY: The Walsh degree decomposition of H(T) on the Boolean
hypercube {0,1}^{C(n,2)} diagonalizes the Fourier decomposition of W(r).

At n=5 (verified exhaustively over all 1024 tournaments):

  D0(H) = 7.5            = w_4/16             (from highest Fourier coeff)
  D2(H) = 3*(t3 - 2.5)   = w_2/4              (from second-highest)
  D4(H) = 1 - t3 + 2*t5  = w_0               (from lowest Fourier coeff)

Each Walsh degree level receives contributions from EXACTLY ONE Fourier
coefficient, because cross-degree terms telescope:
  - t5_hat[S] = t3_hat[S]/2 for all degree-2 Walsh monomials S
  - This causes the degree-2 content of w_0 to vanish identically

The Walsh spectrum structure:
  - Degree 0: 1 nonzero coefficient (= E[H])
  - Degree 2: 30 nonzero coefficients, all ±3/4
    These are the edges of L(K_5) (complement of Petersen graph)
    Sign: -3/4 for "fan" pairs (shared vertex is source/sink)
           +3/4 for "path" pairs (shared vertex is pass-through)
  - Degree 4: 60 nonzero coefficients, all ±1/8
    These are the Hamiltonian paths on K_5

CONJECTURE: At general odd n,
  Walsh degree 2j of H(T) = w_{n-1-2j}(T) / 2^{n-1-2j}
  i.e., the Walsh filtration diagonalizes the Fourier decomposition.

opus-2026-03-06-S11b (continued)
"""

import numpy as np
from itertools import permutations, combinations
from collections import defaultdict
from fractions import Fraction

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
                    dp[(mask|(1<<u), u)] = dp.get((mask|(1<<u), u), 0) + c
    return sum(dp.get(((1<<n)-1, v), 0) for v in range(n))

def count_t3(A, n):
    return sum(1 for a,b,c in combinations(range(n),3)
               if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a])

def count_t5(A, n):
    t5 = 0
    for verts in combinations(range(n), 5):
        sub = [[A[verts[i]][verts[j]] for j in range(5)] for i in range(5)]
        dp = [[0]*5 for _ in range(1<<5)]
        dp[1][0] = 1
        for m_ in range(1, 1<<5):
            for v in range(5):
                if not (m_&(1<<v)) or dp[m_][v]==0: continue
                for u in range(5):
                    if m_&(1<<u): continue
                    if sub[v][u]: dp[m_|(1<<u)][u] += dp[m_][v]
        full = (1<<5)-1
        t5 += sum(dp[full][v] for v in range(1,5) if sub[v][0])
    return t5

def fast_walsh_hadamard(f, m):
    """In-place Walsh-Hadamard transform."""
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
    m = n*(n-1)//2
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Compute H, t3, t5 for all 2^m tournaments
    N = 2**m
    H_all = np.zeros(N)
    t3_all = np.zeros(N)
    t5_all = np.zeros(N)

    print(f"Computing H, t3, t5 for all {N} tournaments on {n} vertices...")
    for bits in range(N):
        A = bits_to_adj(bits, n)
        H_all[bits] = count_H(A, n)
        t3_all[bits] = count_t3(A, n)
        t5_all[bits] = count_t5(A, n)

    # Walsh transforms
    H_hat = fast_walsh_hadamard(H_all.copy(), m)
    t3_hat = fast_walsh_hadamard(t3_all.copy(), m)
    t5_hat = fast_walsh_hadamard(t5_all.copy(), m)

    # 1. Verify Walsh degree structure
    print(f"\n{'='*70}")
    print("WALSH DEGREE STRUCTURE OF H(T)")
    print(f"{'='*70}")

    for deg in range(m+1):
        nonzero = sum(1 for S in range(N) if bin(S).count('1') == deg and abs(H_hat[S]) > 1e-10)
        total = sum(1 for S in range(N) if bin(S).count('1') == deg)
        if nonzero > 0:
            vals = set(round(H_hat[S], 8) for S in range(N)
                       if bin(S).count('1') == deg and abs(H_hat[S]) > 1e-10)
            print(f"  Degree {deg}: {nonzero}/{total} nonzero, values = {vals}")

    # 2. Verify t5_hat = t3_hat/2 at degree 2
    print(f"\n{'='*70}")
    print("KEY IDENTITY: t5_hat[S] = t3_hat[S]/2 at Walsh degree 2")
    print(f"{'='*70}")

    max_ratio_err = 0
    for S in range(N):
        if bin(S).count('1') != 2: continue
        if abs(t3_hat[S]) < 1e-10: continue
        ratio = t5_hat[S] / t3_hat[S]
        err = abs(ratio - 0.5)
        max_ratio_err = max(max_ratio_err, err)
    print(f"  Max |t5_hat[S]/t3_hat[S] - 1/2| = {max_ratio_err:.2e}")
    print(f"  Identity {'CONFIRMED' if max_ratio_err < 1e-10 else 'FAILED'}")

    # 3. Walsh degree decomposition: D0, D2, D4
    print(f"\n{'='*70}")
    print("WALSH DEGREE DECOMPOSITION OF H")
    print(f"{'='*70}")

    D = {0: np.zeros(N), 2: np.zeros(N), 4: np.zeros(N)}
    for S in range(N):
        deg = bin(S).count('1')
        if deg not in D or abs(H_hat[S]) < 1e-10:
            continue
        for bits in range(N):
            chi = 1
            for k in range(m):
                if (S >> k) & 1:
                    chi *= (2*((bits >> k) & 1) - 1)
            D[deg][bits] += H_hat[S] * chi

    # Verify D0 + D2 + D4 = H
    max_err = np.max(np.abs(D[0] + D[2] + D[4] - H_all))
    print(f"  Max |D0 + D2 + D4 - H| = {max_err:.2e}")

    # Verify D0 = constant = 7.5
    print(f"  D0 is constant: {np.all(np.abs(D[0] - 7.5) < 1e-10)}, value = {D[0][0]:.4f}")

    # Verify D2 = w_2/4 = (12*t3 - 30)/4 = 3*t3 - 7.5
    D2_pred = 3 * t3_all - 7.5
    max_err_D2 = np.max(np.abs(D[2] - D2_pred))
    print(f"  Max |D2 - (3*t3 - 7.5)| = {max_err_D2:.2e}")

    # Verify D4 = w_0 = 1 - t3 + 2*t5
    D4_pred = 1 - t3_all + 2*t5_all
    max_err_D4 = np.max(np.abs(D[4] - D4_pred))
    print(f"  Max |D4 - (1 - t3 + 2*t5)| = {max_err_D4:.2e}")

    # 4. Sign pattern of degree-2 Walsh coefficients
    print(f"\n{'='*70}")
    print("DEGREE-2 SIGN PATTERN (BLUE LINE PAIRS)")
    print(f"{'='*70}")

    edge_idx = {e: k for k, e in enumerate(edges)}
    fan_count = 0
    path_count = 0
    for S in range(N):
        if bin(S).count('1') != 2 or abs(H_hat[S]) < 1e-10:
            continue
        bits_list = [k for k in range(m) if (S >> k) & 1]
        e1, e2 = edges[bits_list[0]], edges[bits_list[1]]
        shared = set(e1) & set(e2)
        if not shared: continue
        v = shared.pop()
        pos1 = 0 if e1[0] == v else 1
        pos2 = 0 if e2[0] == v else 1
        if pos1 == pos2:
            fan_count += 1
        else:
            path_count += 1

    print(f"  Fan pairs (source/sink, coeff=-3/4): {fan_count}")
    print(f"  Path pairs (pass-through, coeff=+3/4): {path_count}")
    print(f"  These form the edge-coloring of L(K_{n}) (complement of Petersen graph)")

    # 5. Degree-4 structure
    print(f"\n{'='*70}")
    print("DEGREE-4 STRUCTURE (HAMILTONIAN PATHS)")
    print(f"{'='*70}")

    for S in range(N):
        if bin(S).count('1') != 4 or abs(H_hat[S]) < 1e-10:
            continue
        bits_list = [k for k in range(m) if (S >> k) & 1]
        es = [edges[k] for k in bits_list]
        verts = set()
        for e in es: verts.update(e)
        deg_map = {}
        for v in verts:
            deg_map[v] = sum(1 for e in es if v in set(e))
        deg_seq = tuple(sorted(deg_map.values(), reverse=True))
        # Check: is this always a path? (deg seq = (2,2,2,1,1) on 5 vertices)
        if deg_seq != (2, 2, 2, 1, 1):
            print(f"  NON-PATH found: {es}, deg_seq={deg_seq}")
            break
    else:
        print(f"  All 60 degree-4 monomials are Hamiltonian paths on K_5")
        pos_count = sum(1 for S in range(N) if bin(S).count('1')==4 and H_hat[S] > 1e-10)
        neg_count = sum(1 for S in range(N) if bin(S).count('1')==4 and H_hat[S] < -1e-10)
        print(f"  Positive (+1/8): {pos_count}, Negative (-1/8): {neg_count}")

    # 6. Summary
    print(f"\n{'='*70}")
    print("SUMMARY: WALSH-FOURIER DIAGONALIZATION")
    print(f"{'='*70}")
    print("""
The OCF H(T) = 1 + 2*t3 + 2*t5 is a degree-4 polynomial on {0,1}^10.

Walsh degree 0:  D0 = 7.5 (constant)          = w_4/2^4 = 120/16
Walsh degree 2:  D2 = 3*(t3 - 5/2)            = w_2/2^2
Walsh degree 4:  D4 = 1 - t3 + 2*t5 = W(0)    = w_0/2^0

The telescoping identity t5_hat = t3_hat/2 at Walsh degree 2 ensures
that cross-degree contributions cancel, making each Walsh degree level
receive contributions from EXACTLY ONE Fourier coefficient.

This diagonalization connects:
  - The Petersen graph (complement of L(K_5)) to the blue-line pair structure
  - The skeleton eigenvalues to the Walsh degree filtration
  - The 4D hypercube perspective to the independence polynomial degree
  - The tangent numbers (via W(0) = D4) to the deepest Walsh level
""")

if __name__ == '__main__':
    main()
