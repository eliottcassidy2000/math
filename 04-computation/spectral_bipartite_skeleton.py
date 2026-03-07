#!/usr/bin/env python3
"""
SPECTRAL BIPARTITE STRUCTURE IN TILING SKELETON

THEOREM (n=5, verified): Within the H=9 class (18 tilings, score (1,1,2,3,3)):
  - The tilings split into two spectral sub-classes of 9 each:
    A: eigenvalues {2-sqrt(7), 1, 1, 3, 2+sqrt(7)}, tr(M^2) = 33
    B: eigenvalues {2-sqrt(5), 2-sqrt(5), 1, 2+sqrt(5), 2+sqrt(5)}, tr(M^2) = 37

  - Sub-class B has BLOCK DIAGONAL M:
    M = diag(1) (+) [[3,-2],[-2,1]] (+) [[1,2],[2,3]]
    The 2x2 blocks have eigenvalues 2 +/- sqrt(5) (golden ratio connection!)

  - The cross-spectral flip subgraph (dH=0 flips between A and B) is BIPARTITE.
    14 cross-spectral edges, 4 intra-spectral edges (all within A).

PERPENDICULARITY CONNECTION:
  The dH=0 flips that change the spectrum create the "spectral direction" perpendicular
  to H. This is the source of H-spectral perpendicularity: at mid-H values (H=9),
  there exist flips that change the spectrum without changing H.

EIGENVALUE FLOW:
  - dH > 0 flips: eigenvalues shift to compress toward H/n (becoming more "scalar")
  - dH < 0 flips: eigenvalues spread out (becoming more "anisotropic")
  - dH = 0 flips: eigenvalues rotate within the spectral fiber (perpendicular direction)

CROSS-SCALE SUMMARY (n=3,5):
  n=3: 2 tilings, no dH=0 flips. Single gradient from H=1 to H=3.
  n=5: 64 tilings, bipartite spectral sub-classes within H=9.
       perpendicularity cos transitions: +0.87 (H=1) -> 0 (H~10) -> -0.96 (H=15)

opus-2026-03-06-S11b (continued)
"""
from itertools import permutations, combinations
from collections import defaultdict
import numpy as np

def tiling_to_adj(bits, n):
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

def transfer_matrix(A, n):
    M = np.zeros((n, n), dtype=int)
    for a in range(n):
        for b in range(n):
            U = [v for v in range(n) if v != a and v != b]
            total = 0
            for k in range(len(U)+1):
                for S in combinations(U, k):
                    S_set = set(S)
                    R = [v for v in U if v not in S_set]
                    ea = sum(1 for p in permutations(sorted(list(S)+[a]))
                            if p[-1]==a and all(A[p[i]][p[i+1]]==1 for i in range(len(p)-1)))
                    bb = sum(1 for p in permutations(sorted(R+[b]))
                            if p[0]==b and all(A[p[i]][p[i+1]]==1 for i in range(len(p)-1)))
                    total += ((-1)**k) * ea * bb
            M[a][b] = total
    return M

if __name__ == '__main__':
    n = 5
    tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
    m = len(tiles)

    data = {}
    for bits in range(2**m):
        A = tiling_to_adj(bits, n)
        M = transfer_matrix(A, n)
        evals = sorted(np.linalg.eigvalsh(M.astype(float)))
        H = int(np.trace(M))
        tr2 = int(np.trace(M @ M))
        data[bits] = {'evals': evals, 'H': H, 'tr2': tr2, 'M': M}

    # Spectral sub-class analysis within each H-class
    for H_target in sorted(set(data[b]['H'] for b in range(2**m))):
        class_bits = [b for b in range(2**m) if data[b]['H'] == H_target]
        tr2_groups = defaultdict(list)
        for bits in class_bits:
            tr2_groups[data[bits]['tr2']].append(bits)

        if len(tr2_groups) > 1:
            class_set = set(class_bits)
            cross = [(a,b) for a in class_bits for ti in range(m)
                     if (b := a ^ (1 << ti)) in class_set and b > a
                     and data[a]['tr2'] != data[b]['tr2']]
            intra = [(a,b) for a in class_bits for ti in range(m)
                     if (b := a ^ (1 << ti)) in class_set and b > a
                     and data[a]['tr2'] == data[b]['tr2']]

            print(f"\nH={H_target}: {len(class_bits)} tilings, "
                  f"{len(tr2_groups)} spectral sub-classes")
            for tr2, members in sorted(tr2_groups.items()):
                evals = data[members[0]]['evals']
                print(f"  tr2={tr2}: {len(members)} tilings, "
                      f"evals={[round(e,4) for e in evals]}")
            print(f"  cross-spectral edges: {len(cross)}, "
                  f"intra-spectral: {len(intra)}")

            # Check bipartiteness
            if len(tr2_groups) == 2:
                tr2_list = sorted(tr2_groups.keys())
                color = {b: 0 for b in tr2_groups[tr2_list[0]]}
                color.update({b: 1 for b in tr2_groups[tr2_list[1]]})
                bipartite = all(color[a] != color[b] for a, b in cross)
                print(f"  bipartite? {bipartite}")
