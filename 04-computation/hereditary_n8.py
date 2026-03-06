#!/usr/bin/env python3
"""
Quick check: hereditary property at n=8.
Since n=8 is even, we expect NO hereditary chain (consistent with n=4,6 pattern).

Find an n=8 SC maximizer (H=661), compute its vertex deletion H-values.
Check if any equal max H(7)=189.

kind-pasteur-2026-03-06-S18f
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import hamiltonian_path_count

def score_seq(T):
    return tuple(sorted(sum(T[i]) for i in range(len(T))))

def generate_sc_one(n, sigma, bits):
    """Generate one SC tournament from bit index."""
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    free = []
    paired_reps = []
    seen = set()
    for i, j in pairs:
        if (i, j) in seen:
            continue
        si, sj = sigma[i], sigma[j]
        img = (min(si, sj), max(si, sj))
        if img == (i, j):
            free.append((i, j))
            seen.add((i, j))
        else:
            paired_reps.append(((i, j), img))
            seen.add((i, j))
            seen.add(img)

    T = [[0]*n for _ in range(n)]
    bit_idx = 0
    for i, j in free:
        val = (bits >> bit_idx) & 1
        T[i][j] = val
        T[j][i] = 1 - val
        bit_idx += 1
    for (i, j), (si, sj) in paired_reps:
        val = (bits >> bit_idx) & 1
        T[i][j] = val
        T[j][i] = 1 - val
        T[sigma[i]][sigma[j]] = 1 - val
        T[sigma[j]][sigma[i]] = val
        bit_idx += 1
    return T

n = 8
sigma = [1, 0, 3, 2, 5, 4, 7, 6]

# Find a maximizer
print("Finding n=8 maximizer (H=661)...")
for bits in range(1 << 16):
    T = generate_sc_one(n, sigma, bits)
    s = score_seq(T)
    if s != (3, 3, 3, 3, 4, 4, 4, 4):
        continue
    h = hamiltonian_path_count(T)
    if h == 661:
        print(f"Found! bits={bits}")

        # Compute deletion H-values
        print(f"\nDeletion H-values (max H(7)=189):")
        for v in range(n):
            verts = [i for i in range(n) if i != v]
            sub = [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
            h_sub = hamiltonian_path_count(sub)
            s_sub = score_seq(sub)
            is_max = h_sub == 189
            print(f"  v={v}: H(T-v)={h_sub}, score={s_sub}, is_max={is_max}")

        # Deletion H-sum
        h_sum = sum(hamiltonian_path_count(
            [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
        ) for v in range(n) for verts in [[i for i in range(n) if i != v]])
        # Wait, that's wrong. Let me redo.
        del_hs = []
        for v in range(n):
            verts = [i for i in range(n) if i != v]
            sub = [[T[verts[i]][verts[j]] for j in range(n-1)] for i in range(n-1)]
            del_hs.append(hamiltonian_path_count(sub))

        print(f"\nDeletion spectrum: {sorted(del_hs, reverse=True)}")
        print(f"Sum H(T-v) = {sum(del_hs)}, H = {h}, ratio = {sum(del_hs)/h:.4f}")
        print(f"All equal? {len(set(del_hs)) == 1}")
        print(f"Any equal to 189? {189 in del_hs}")

        # Also check n -> n-2
        print(f"\nPair deletions to n=6 (max H(6)=45):")
        pair_max_count = 0
        for v1 in range(n):
            for v2 in range(v1+1, n):
                verts = [i for i in range(n) if i != v1 and i != v2]
                sub = [[T[verts[i]][verts[j]] for j in range(n-2)] for i in range(n-2)]
                h_sub = hamiltonian_path_count(sub)
                if h_sub == 45:
                    pair_max_count += 1
        print(f"  {pair_max_count}/{n*(n-1)//2} pair-deletions give max H(6)=45")

        break

print("\nDone.")
