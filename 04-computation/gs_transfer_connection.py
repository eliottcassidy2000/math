#!/usr/bin/env python3
"""
Connection between GS structure and transfer matrix symmetry.
Instance: opus-2026-03-06-S11

Key question: Does M[a,b] restricted to GS tilings automatically give M[a,b]=M[b,a]?

The transfer matrix M[a,b] at a specific tournament value (c=0) counts:
  M[a,b] = sum_S (-1)^|S| * E_a(S+a) * B_b(S+b)
where S ranges over subsets of {v : v != a, v != b}.

At c=0, this becomes a sum over 2-path covers weighted by tournament arc signs.

For a specific tiling (bit pattern b, giving tournament T(b)):
  M[a,b] counts the weighted 2-path covers of T(b).

If b is a GS tiling, then T(b) is self-converse. Does this force M[a,b] = M[b,a]?

We check: for each GS tiling, compute M[a,b] and M[b,a] directly.
"""
from itertools import permutations
from collections import defaultdict

def compute_M_entry(A, n, a, b):
    """Compute M[a,b] for tournament with adjacency matrix A.
    M[a,b] = sum_S (-1)^|S| * H_a(S+a) * H_b(R+b)
    where S union R = {v != a, b}, and
    H_a(S+a) = #Hamiltonian paths in T[S+a] ending at a,
    H_b(R+b) = #Hamiltonian paths in T[R+b] ending at b.

    At c=0 (unit arc weights), this simplifies.
    """
    U = [v for v in range(n) if v != a and v != b]
    total = 0

    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1) ** len(S)

        # H_a(S + a): count Hamiltonian paths in induced subgraph on S+{a} ending at a
        S_a = sorted(S + [a])
        ha = count_paths_ending_at(A, S_a, a)

        # H_b(R + b): count Hamiltonian paths in induced subgraph on R+{b} ending at b
        R_b = sorted(R + [b])
        hb = count_paths_ending_at(A, R_b, b)

        total += sign * ha * hb

    return total

def count_paths_ending_at(A, vertices, target):
    """Count Hamiltonian paths in induced subgraph ending at target."""
    if len(vertices) == 1:
        return 1 if vertices[0] == target else 0

    count = 0
    for p in permutations(vertices):
        if p[-1] != target:
            continue
        valid = True
        for i in range(len(p) - 1):
            if not A[p[i]][p[i+1]]:
                valid = False
                break
        if valid:
            count += 1
    return count

def analyze_gs_transfer(n):
    """Analyze M[a,b] vs M[b,a] for GS tilings."""
    verts = list(range(n, 0, -1))
    tiles = []
    for y in range(1, n-1):
        for x in range(n, y+1, -1):
            tiles.append((x, y))
    m = len(tiles)
    tile_idx = {(x,y): i for i, (x,y) in enumerate(tiles)}

    trans_map = []
    for (x, y) in tiles:
        nx, ny = n - y + 1, n - x + 1
        trans_map.append(tile_idx.get((nx, ny), -1))

    def is_grid_sym(bits):
        for i in range(m):
            j = trans_map[i]
            if j >= 0 and j != i and ((bits >> i) & 1) != ((bits >> j) & 1):
                return False
        return True

    def bits_to_adj(bits):
        A = [[0]*n for _ in range(n)]
        for k in range(n-1):
            A[k][k+1] = 1
        for i, (xL, yL) in enumerate(tiles):
            xi = verts.index(xL)
            yi = verts.index(yL)
            if (bits >> i) & 1 == 0:
                A[xi][yi] = 1
            else:
                A[yi][xi] = 1
        return A

    # Generate GS tilings
    gs_tilings = [b for b in range(1 << m) if is_grid_sym(b)]

    print(f"\nn={n}: {len(gs_tilings)} GS tilings")

    # For each GS tiling, check M[a,b] = M[b,a] for ALL (a,b) pairs
    sym_always = 0
    sym_fails = 0

    for idx, bits in enumerate(gs_tilings[:20]):  # Check first 20
        A = bits_to_adj(bits)
        all_sym = True
        for a in range(n):
            for b in range(a+1, n):
                mab = compute_M_entry(A, n, a, b)
                mba = compute_M_entry(A, n, b, a)
                if mab != mba:
                    all_sym = False
                    break
            if not all_sym:
                break

        if all_sym:
            sym_always += 1
        else:
            sym_fails += 1
            print(f"  GS tiling {bits:0{m}b}: M NOT symmetric!")

    checked = min(20, len(gs_tilings))
    print(f"  Checked {checked} GS tilings: {sym_always} symmetric, {sym_fails} asymmetric")

    # Also check some non-GS tilings for comparison
    non_gs = [b for b in range(1 << m) if not is_grid_sym(b)][:10]
    non_gs_sym = 0
    for bits in non_gs:
        A = bits_to_adj(bits)
        all_sym = True
        for a in range(n):
            for b in range(a+1, n):
                mab = compute_M_entry(A, n, a, b)
                mba = compute_M_entry(A, n, b, a)
                if mab != mba:
                    all_sym = False
                    break
            if not all_sym:
                break
        if all_sym:
            non_gs_sym += 1

    print(f"  Non-GS check: {non_gs_sym}/{len(non_gs)} symmetric")

def main():
    for n in [4, 5]:
        analyze_gs_transfer(n)

if __name__ == '__main__':
    main()
