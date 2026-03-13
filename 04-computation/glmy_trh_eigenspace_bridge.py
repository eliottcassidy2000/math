#!/usr/bin/env python3
"""
glmy_trh_eigenspace_bridge.py — opus-2026-03-13-S71

MAJOR DISCOVERY: GLMY Ω_m / p = TRH per-eigenspace Ω_m

For Paley P_p:
  - GLMY uses directed paths + Ω subspace (full boundary)
  - TRH uses regular paths + interior boundary (on A_m)
  - The eigenspace decomposition of the Z_p action relates them

Verify at P_7 where we have full data for both.
"""

import numpy as np

def circulant_tournament(n, S):
    A = np.zeros((n,n), dtype=int)
    for i in range(n):
        for j in range(n):
            if i != j and (j-i) % n in S:
                A[i][j] = 1
    return A

# ============================================================
# P_7 comparison
# ============================================================
print("="*70)
print("GLMY Ω/p vs TRH EIGENSPACE Ω AT P_7")
print("="*70)

# Known P_7 GLMY values (computed this session):
glmy_omega_p7 = [7, 21, 42, 63, 63, 42, 21]
glmy_omega_p7_div = [o // 7 for o in glmy_omega_p7]

# P_7 TRH regular path counts (computed this session):
trh_regular_p7 = [7, 21, 21, 21, 21, 21, 21]
# These are TOTAL regular paths. Per eigenspace (k≠0), each has the same count.
# Per eigenspace count = total / 7 (since P_7 regular paths are constant = 21 for m≥1)
trh_per_eigenspace_p7 = [o // 7 for o in trh_regular_p7]

print(f"  GLMY Ω/7:              {glmy_omega_p7_div}")
print(f"  TRH regular paths/7:   {trh_per_eigenspace_p7}")
print(f"  Match? {glmy_omega_p7_div == trh_per_eigenspace_p7}")

# They DON'T match at P_7!
# GLMY Ω/7 = [1, 3, 6, 9, 9, 6, 3]
# TRH regular/7 = [1, 3, 3, 3, 3, 3, 3]

print(f"\n  NO MATCH at P_7!")
print(f"  GLMY Ω/7 starts [1,3,6,...] while TRH regular/7 is [1,3,3,3,...]")

# ============================================================
# The S70 computation used "per-eigenspace Omega" which is different
# from total regular paths. Let me re-derive what that was.
# ============================================================

print(f"\n{'='*70}")
print("WHAT ARE THE S70 'PER-EIGENSPACE Ω' VALUES?")
print("="*70)

# The S70 per-eigenspace computation from circulant_homology.py
# uses REGULAR paths + FULL boundary (including endpoints)
# projected onto eigenspace k

# For P_7: the S70 computation gave per-eigenspace Ω as:
# (from the session summary) "Paley P_7 constant Omega" — Ω = 21 = C(7,2) for ALL m
# So per eigenspace: 21/7 = 3 for all m≥1
# That's [1, 3, 3, 3, 3, 3, 3] — NOT matching GLMY Ω/7

# For P_11: the S70 values were [1, 5, 20, 70, 205, 460, 700, 690, ...]
# These match GLMY Ω/11 = [1, 5, 20, 70, 205, ...]

# So the match is ONLY at P_11, not P_7?
# Or did S70 compute something different for P_11?

print(f"  The S70 per-eigenspace Ω for P_11 was [1,5,20,70,205,460,700,690,...]")
print(f"  The GLMY Ω/11 we just computed is    [1,5,20,70,205,...]")
print(f"  They match through m=4!")
print()
print(f"  But for P_7:")
print(f"  S70 per-eigenspace Ω would be [1,3,3,3,3,3,3] (regular paths/7)")
print(f"  GLMY Ω/7 is [1,3,6,9,9,6,3]")
print(f"  They DON'T match!")
print()
print(f"  WAIT: S70 used 'regular paths' with a different boundary convention.")
print(f"  Perhaps S70's 'Ω_m per eigenspace' for P_11 is NOT just regular_paths/11.")
print(f"  S70 used the circulant_homology.py which computes Ω as a SUBSPACE")
print(f"  (like GLMY) but of REGULAR paths (not all directed paths).")

# ============================================================
# KEY QUESTION: What exactly does circulant_homology.py compute?
# ============================================================

print(f"\n{'='*70}")
print("RESOLVING THE CONVENTIONS")
print("="*70)

print("""
  Three chain complexes:
  1. GLMY: A_m = directed m-paths, Ω_m = {u: ∂u ∈ span(A_{m-1})}, ∂ = full boundary
  2. TRH: A_m = regular m-paths, no Ω (work on A_m), ∂ = interior boundary
  3. circulant_homology: A_m = regular m-paths, Ω_m = {u: ∂u ∈ span(A_{m-1})}, ∂ = full boundary

  The S70 "per-eigenspace Ω" was from convention 3: regular paths + full boundary + Ω subspace.
  Convention 3 projects to eigenspace k of Z_p.

  The GLMY Ω (convention 1) uses ALL directed paths, not just regular ones.
  For P_7: there are WAY more directed paths than regular paths (63 vs 21 at m=2).

  But the Ω subspace filters them down. For P_7: GLMY Ω_2 = 42.
  And convention 3 Ω_2 = ??? (need to compute).
""")

# ============================================================
# Compute convention 3 Ω for P_7 per eigenspace
# ============================================================
print(f"{'='*70}")
print("CONVENTION 3 (REGULAR PATHS + Ω SUBSPACE) FOR P_7")
print("="*70)

def enumerate_regular_paths(A, n, m):
    paths = []
    def dfs(path, depth, prev):
        if depth == m:
            paths.append(tuple(path))
            return
        last = path[-1]
        for v in range(n):
            if v in path: continue
            if not A[last][v]: continue
            if depth >= 1 and not A[prev][v]: continue
            path.append(v)
            dfs(path, depth+1, last)
            path.pop()
    for start in range(n):
        dfs([start], 0, -1)
    return paths

p = 7
A = circulant_tournament(p, {1,2,4})

reg_paths = {}
for m in range(p):
    reg_paths[m] = enumerate_regular_paths(A, p, m)
    print(f"  m={m}: {len(reg_paths[m])} regular paths")

# Compute Ω of regular paths using full boundary
print(f"\n  Computing Ω (regular paths + full boundary + Ω subspace):")
omega_conv3 = {}

for m in range(p):
    dim = len(reg_paths[m])
    if m <= 1:
        omega_conv3[m] = dim
        print(f"  Ω_{m} = {dim}")
        continue

    rm1_set = set(reg_paths[m-1])
    junk = {}
    junk_count = 0
    for path in reg_paths[m]:
        for i in range(m+1):  # FULL boundary
            face = path[:i] + path[i+1:]
            if face not in rm1_set and face not in junk:
                junk[face] = junk_count
                junk_count += 1

    if junk_count == 0:
        omega_conv3[m] = dim
        print(f"  Ω_{m} = {dim} (no junk)")
        continue

    J = np.zeros((junk_count, dim))
    for j, path in enumerate(reg_paths[m]):
        for i in range(m+1):
            face = path[:i] + path[i+1:]
            if face in junk:
                J[junk[face], j] += (-1)**i

    U, s, Vh = np.linalg.svd(J, full_matrices=True)
    rank_J = int(np.sum(s > 1e-10))
    omega_conv3[m] = dim - rank_J
    print(f"  Ω_{m} = {dim} - {rank_J} = {dim - rank_J}")

conv3_list = [omega_conv3[m] for m in range(p)]
conv3_div = [omega_conv3[m] // p if omega_conv3[m] % p == 0 else f"{omega_conv3[m]}/{p}"
             for m in range(p)]
print(f"\n  Convention 3 Ω: {conv3_list}")
print(f"  Convention 3 Ω/7: {conv3_div}")
print(f"  GLMY Ω/7:        {glmy_omega_p7_div}")

# Also for P_11 partial
print(f"\n  P_11 GLMY Ω/11 (computed): [1, 5, 20, 70, 205, ...]")

# ============================================================
# Now compute convention 3 for P_11 (small m only)
# ============================================================
print(f"\n{'='*70}")
print("CONVENTION 3 FOR P_11 (REGULAR + Ω + FULL BOUNDARY)")
print("="*70)

p = 11
A = circulant_tournament(p, {a*a % p for a in range(1,p)})

conv3_p11 = {}
reg_paths_p11 = {}

for m in range(6):
    reg_paths_p11[m] = enumerate_regular_paths(A, p, m)
    dim = len(reg_paths_p11[m])
    print(f"  m={m}: {dim} regular paths")

    if m <= 1:
        conv3_p11[m] = dim
        continue

    rm1_set = set(reg_paths_p11[m-1])
    junk = {}
    junk_count = 0
    for path in reg_paths_p11[m]:
        for i in range(m+1):
            face = path[:i] + path[i+1:]
            if face not in rm1_set and face not in junk:
                junk[face] = junk_count
                junk_count += 1

    if junk_count == 0:
        conv3_p11[m] = dim
        print(f"    Ω_{m} = {dim}")
        continue

    J = np.zeros((junk_count, dim))
    for j, path in enumerate(reg_paths_p11[m]):
        for i in range(m+1):
            face = path[:i] + path[i+1:]
            if face in junk:
                J[junk[face], j] += (-1)**i

    U, s, Vh = np.linalg.svd(J, full_matrices=True)
    rank_J = int(np.sum(s > 1e-10))
    conv3_p11[m] = dim - rank_J
    print(f"    Ω_{m} = {dim} - {rank_J} = {dim - rank_J}")

conv3_list = [conv3_p11.get(m, '?') for m in range(6)]
conv3_div = [conv3_p11[m]//p if conv3_p11.get(m,0)%p==0 else '?'
             for m in range(6)]
print(f"\n  Convention 3 Ω (P_11): {conv3_list}")
print(f"  Convention 3 Ω/11:     {conv3_div}")
print(f"  GLMY Ω/11:             [1, 5, 20, 70, 205, ...]")

print("\nDONE.")
