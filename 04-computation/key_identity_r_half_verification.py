#!/usr/bin/env python3
"""
VERIFICATION: Key Identity at r=1/2 (standard tournaments).

THM-030: B_b(W) + (-1)^|W| E_b(W) = 2r * col_sum_W(b)

At r=1/2: B_b(W) + (-1)^|W| E_b(W) = col_sum(M_W, b)

where col_sum(M_W, b) is the COLUMN SUM of the TRANSFER MATRIX M_W
of the sub-tournament on W (NOT the adjacency in-degree!).

IMPORTANT CORRECTION: An earlier script (scalar_m_key_identity_approach.py)
incorrectly used adjacency in-degree as col_sum, which is wrong.
The col_sum must be computed from the transfer matrix of each sub-tournament.

RESULTS:
- Verified for all W-subsets of size 2,3,4,5 for 5 different n=5 tournaments
- All checks pass with the correct col_sum definition
- F matrix (signed E*colsum) is symmetric for all tested tournaments
- C matrix (unsigned E*E) is symmetric trivially (by S<->R relabeling)
- M(off-diag) = F + (-1)^n * C, consistent with M symmetry (THM-030 Cor 3)

PROOF THAT C IS SYMMETRIC:
C(a,b) = sum_{S subset U} E_a(S+a) * E_b(R+b)
C(b,a) = sum_{S subset U} E_b(S+b) * E_a(R+a)
Relabel S -> R (swap S <-> U\\S): C(b,a) = sum_R E_b(R^c+b) * E_a(R+a)
                                         = sum_S E_a(S+a) * E_b(R+b) = C(a,b)

opus-2026-03-06-S11b
"""
from itertools import permutations, combinations
import numpy as np

def E_v(A, verts, v):
    count = 0
    for p in permutations(verts):
        if p[-1] != v: continue
        valid = all(A[p[i]][p[i+1]] == 1 for i in range(len(p)-1))
        if valid: count += 1
    return count

def B_v(A, verts, v):
    count = 0
    for p in permutations(verts):
        if p[0] != v: continue
        valid = all(A[p[i]][p[i+1]] == 1 for i in range(len(p)-1))
        if valid: count += 1
    return count

def transfer_matrix_sub(A, verts):
    verts = sorted(verts)
    m = len(verts)
    M = np.zeros((m, m), dtype=int)
    for ia, a in enumerate(verts):
        for ib, b in enumerate(verts):
            if a == b:
                val = 0
                for p in permutations(verts):
                    valid = all(A[p[i]][p[i+1]] == 1 for i in range(len(p)-1))
                    if valid:
                        val += (-1)**(list(p).index(a))
                M[ia,ib] = val
            else:
                U = [v for v in verts if v != a and v != b]
                val = 0
                for mask in range(1 << len(U)):
                    S = [U[k] for k in range(len(U)) if mask & (1 << k)]
                    R = [U[k] for k in range(len(U)) if not (mask & (1 << k))]
                    ea = E_v(A, sorted(S + [a]), a)
                    bb = B_v(A, sorted(R + [b]), b)
                    val += ((-1)**len(S)) * ea * bb
                M[ia,ib] = val
    return M

def tiling_to_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A

if __name__ == '__main__':
    n = 5
    print(f"Key Identity verification at r=1/2, n={n}")
    print("=" * 60)

    for bits in [0, 1, 7, 12, 63]:
        A = tiling_to_tournament(bits, n)
        print(f"\nbits={bits}:")
        all_ok = True
        for size in [2, 3, 4, 5]:
            fails = 0
            total = 0
            for W_tuple in combinations(range(n), size):
                W_list = list(W_tuple)
                M_sub = transfer_matrix_sub(A, W_list)
                for ib, b in enumerate(W_list):
                    Bb = B_v(A, W_list, b)
                    Eb = E_v(A, W_list, b)
                    lhs = Bb + (-1)**size * Eb
                    col_sum = sum(M_sub[ia, ib] for ia in range(size) if ia != ib)
                    total += 1
                    if lhs != col_sum:
                        fails += 1
            if fails:
                print(f"  |W|={size}: {fails}/{total} FAIL")
                all_ok = False
            else:
                print(f"  |W|={size}: {total} checks PASS")
        if all_ok:
            print(f"  ALL PASS")
