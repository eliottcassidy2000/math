#!/usr/bin/env python3
"""
Test: Can the Paley tournament on 7 vertices be self-flip under ANY base path?

Self-flip means: flipping all non-base-path arcs gives an isomorphic tournament.
This explains why the Paley (regular, high-Aut) tournaments are NOT in the SC+SF kernel.

Instance: opus-2026-03-06-S7
"""

from itertools import permutations

# Build Paley tournament on 7 vertices
# QR(7) = {1, 2, 4}
qr = {1, 2, 4}
n = 7
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j and (j - i) % n in qr:
            A[i][j] = 1

# Count Ham paths
ham_paths = []
for p in permutations(range(n)):
    if all(A[p[k]][p[k+1]] for k in range(n-1)):
        ham_paths.append(p)
H = len(ham_paths)
print(f"H(Paley_7) = {H}")
print(f"|Aut| = 21")

def are_iso(A, B):
    for p in permutations(range(n)):
        if all(A[p[i]][p[j]] == B[i][j] for i in range(n) for j in range(n)):
            return True, p
    return False, None

# For the standard base path 0->1->2->...->6
print("\n=== Standard base path 0->1->...->6 ===")
A_flip = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i == j:
            continue
        if abs(i-j) == 1:
            A_flip[i][j] = A[i][j]
        else:
            A_flip[i][j] = 1 - A[i][j]

iso, perm = are_iso(A, A_flip)
scores_orig = tuple(sorted([sum(A[i]) for i in range(n)], reverse=True))
scores_flip = tuple(sorted([sum(A_flip[i]) for i in range(n)], reverse=True))
H_flip = sum(1 for p in permutations(range(n)) if all(A_flip[p[k]][p[k+1]] for k in range(n-1)))
print(f"Scores orig: {scores_orig}")
print(f"Scores flip: {scores_flip}")
print(f"H(flip) = {H_flip}")
print(f"Isomorphic to original? {iso}")

# Check ALL Ham paths as base paths
print(f"\n=== Checking all {H} Hamiltonian paths as base paths ===")
count_sf = 0
flip_h_values = set()

for pi, path in enumerate(ham_paths):
    path_arcs = set()
    for k in range(n-1):
        path_arcs.add((path[k], path[k+1]))

    A_p_flip = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if (i, j) in path_arcs:
                A_p_flip[i][j] = A[i][j]
            else:
                A_p_flip[i][j] = 1 - A[i][j]

    iso_p, perm_p = are_iso(A, A_p_flip)
    if iso_p:
        count_sf += 1
        print(f"  Path {path}: flip IS isomorphic! via perm {perm_p}")

    if pi < 5 or iso_p:
        h_f = sum(1 for p in permutations(range(n))
                  if all(A_p_flip[p[k]][p[k+1]] for k in range(n-1)))
        flip_h_values.add(h_f)

print(f"\nOut of {H} Hamiltonian paths, {count_sf} give self-flip")
print(f"Sample flip H values: {sorted(flip_h_values)}")

# Now check the OTHER regular tournaments
print("\n\n=== Checking other regular tournaments on 7 vertices ===")
# Find all regular tournaments (score (3,3,3,3,3,3,3))
from collections import defaultdict

# We already know there are 3 iso classes with this score:
# H=189 (Paley), H=175, H=171
# Let me find representatives of the other two

# Generate some regular tournaments by random-ish construction
# Actually, let me just check a known non-Paley regular tournament

# Non-Paley regular: try QR = {1,3,5} (mod 7)
# Check: 1^2=1, 2^2=4, 3^2=2 mod 7. QR={1,2,4}. NQR={3,5,6}.
# {1,3,5} is not QR, not NQR. Let me try a different approach.

# Build from a specific circulant
# Try S = {1, 2, 5} (mod 7)
A2 = [[0]*7 for _ in range(7)]
S = {1, 2, 5}
for i in range(7):
    for j in range(7):
        if i != j and (j-i) % 7 in S:
            A2[i][j] = 1

scores2 = tuple(sorted([sum(A2[i]) for i in range(7)], reverse=True))
print(f"S={{1,2,5}}: scores = {scores2}")
if scores2 == (3,3,3,3,3,3,3):
    H2 = sum(1 for p in permutations(range(7)) if all(A2[p[k]][p[k+1]] for k in range(6)))
    print(f"H = {H2}")
    iso_to_paley, _ = are_iso(A, A2)
    print(f"Isomorphic to Paley? {iso_to_paley}")

# Try S = {1, 3, 4}
A3 = [[0]*7 for _ in range(7)]
S3 = {1, 3, 4}
for i in range(7):
    for j in range(7):
        if i != j and (j-i) % 7 in S3:
            A3[i][j] = 1

scores3 = tuple(sorted([sum(A3[i]) for i in range(7)], reverse=True))
print(f"\nS={{1,3,4}}: scores = {scores3}")
if scores3 == (3,3,3,3,3,3,3):
    H3 = sum(1 for p in permutations(range(7)) if all(A3[p[k]][p[k+1]] for k in range(6)))
    print(f"H = {H3}")
    iso_to_paley, _ = are_iso(A, A3)
    print(f"Isomorphic to Paley? {iso_to_paley}")

    # Check self-flip for this tournament too
    hp3 = [p for p in permutations(range(7)) if all(A3[p[k]][p[k+1]] for k in range(6))]
    sf3 = 0
    for path in hp3[:20]:  # check first 20
        path_arcs = set((path[k], path[k+1]) for k in range(6))
        A3f = [[0]*7 for _ in range(7)]
        for i in range(7):
            for j in range(7):
                if i == j: continue
                if (i,j) in path_arcs:
                    A3f[i][j] = A3[i][j]
                else:
                    A3f[i][j] = 1 - A3[i][j]
        iso3, _ = are_iso(A3, A3f)
        if iso3: sf3 += 1
    print(f"Self-flip count (first 20 paths): {sf3}/20")


if __name__ == '__main__':
    pass
