#!/usr/bin/env python3
"""
OCF Arc-Flip Derivative with FULL Odd-Cycle Counting
opus-2026-03-14-S71f

Previous script only counted 3-cycles. Ω(T) has ALL odd directed cycles.
At n=5, this includes 5-cycles (Hamiltonian cycles).
"""
from itertools import permutations, combinations

def make_tournament(bits, n):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_hp(A, n):
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n-1):
            if A[perm[i]][perm[i+1]] != 1:
                ok = False
                break
        if ok: count += 1
    return count

def find_all_odd_cycles(A, n):
    """Find all directed odd cycles as frozensets of vertices.
    Returns list of (vertex_set, length) tuples.
    A directed k-cycle visits k distinct vertices in order."""
    cycles = set()
    
    # 3-cycles
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            cycles.add(frozenset([i,j,k]))
        if A[i][k] and A[k][j] and A[j][i]:
            cycles.add(frozenset([i,j,k]))
    
    # 5-cycles (only at n>=5)
    if n >= 5:
        for verts in combinations(range(n), 5):
            # Check all cyclic orderings
            from itertools import permutations as perms5
            vlist = list(verts)
            # Fix first vertex, try all orderings of rest
            for p in perms5(vlist[1:]):
                order = [vlist[0]] + list(p)
                # Check if this is a directed cycle
                is_cycle = True
                for idx in range(5):
                    if A[order[idx]][order[(idx+1)%5]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.add(frozenset(verts))
                    break  # One direction found is enough for the vertex set
    
    # 7-cycles (only at n>=7)
    if n >= 7:
        for verts in combinations(range(n), 7):
            vlist = list(verts)
            # Too many permutations (720). Sample or use smarter approach.
            # For n=7 only, all vertices are used, so just check permutations
            found = False
            for p in permutations(vlist[1:]):
                if found: break
                order = [vlist[0]] + list(p)
                is_cycle = True
                for idx in range(7):
                    if A[order[idx]][order[(idx+1)%7]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.add(frozenset(verts))
                    found = True
    
    return list(cycles)

def compute_full_alpha(A, n):
    """Compute α₁, α₂ using ALL odd cycles."""
    cycles = find_all_odd_cycles(A, n)
    alpha1 = len(cycles)
    
    # Disjoint pairs
    alpha2 = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if cycles[i].isdisjoint(cycles[j]):
                alpha2 += 1
    
    return alpha1, alpha2

def flip_arc(A, n, i, j):
    B = [row[:] for row in A]
    B[i][j], B[j][i] = B[j][i], B[i][j]
    return B

print("="*70)
print("FULL OCF ARC-FLIP VERIFICATION")
print("="*70)

# n=5 exhaustive
n = 5
total_edges = n*(n-1)//2
num_t = 2**total_edges
print(f"\nn={n}: {num_t} tournaments, {total_edges} edges each")

# First verify OCF
print("\n--- Verifying H = I(Ω,2) = 1 + 2α₁ ---")
print("(α₂=0 at n=5 since two disjoint odd cycles need ≥6 vertices)")
ocf_ok = 0
ocf_fail = 0
for bits in range(num_t):
    A = make_tournament(bits, n)
    h = count_hp(A, n)
    a1, a2 = compute_full_alpha(A, n)
    predicted = 1 + 2*a1 + 4*a2
    if h == predicted:
        ocf_ok += 1
    else:
        ocf_fail += 1
        if ocf_fail <= 5:
            print(f"  FAIL: bits={bits}, H={h}, α₁={a1}, α₂={a2}, predicted={predicted}")

print(f"  OCF verified: {ocf_ok}/{num_t}, failures: {ocf_fail}")

# Now check arc-flip derivatives
print("\n--- Arc-flip derivatives with full odd-cycle α₁ ---")
delta_map = {}
for bits in range(num_t):
    A = make_tournament(bits, n)
    h = count_hp(A, n)
    a1, a2 = compute_full_alpha(A, n)
    
    for i in range(n):
        for j in range(i+1, n):
            B = flip_arc(A, n, i, j)
            hb = count_hp(B, n)
            b1, b2 = compute_full_alpha(B, n)
            dh = hb - h
            da1 = b1 - a1
            da2 = b2 - a2
            key = (dh, da1, da2)
            delta_map[key] = delta_map.get(key, 0) + 1

print("  (ΔH, Δα₁, Δα₂) → count, match?")
mismatches = 0
for k in sorted(delta_map.keys()):
    dh, da1, da2 = k
    predicted = 2*da1 + 4*da2
    match = "✓" if dh == predicted else "✗"
    if dh != predicted: mismatches += 1
    print(f"    ({dh:+3d}, {da1:+3d}, {da2:+3d}): {delta_map[k]:5d} {match}")

print(f"\n  Mismatches: {mismatches}/{len(delta_map)}")
if mismatches == 0:
    print("  ★ ΔH = 2Δα₁ + 4Δα₂ EXACTLY with full odd-cycle counting!")

# Part 2: Breakdown of cycle types
print("\n--- Cycle type breakdown at n=5 ---")
t3_d5_dist = {}
for bits in range(num_t):
    A = make_tournament(bits, n)
    h = count_hp(A, n)
    
    # Count 3-cycles
    t3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]: t3 += 1
        if A[i][k] and A[k][j] and A[j][i]: t3 += 1
    
    # Count 5-cycles (Hamiltonian cycles)
    d5 = 0
    for p in permutations(range(n)):
        if p[0] > p[1]: continue  # Avoid double-counting by direction
        is_cycle = True
        for idx in range(5):
            if A[p[idx]][p[(idx+1)%5]] != 1:
                is_cycle = False
                break
        if is_cycle: d5 += 1
    
    key = (t3, d5)
    if key not in t3_d5_dist:
        t3_d5_dist[key] = {'count': 0, 'H_vals': set()}
    t3_d5_dist[key]['count'] += 1
    t3_d5_dist[key]['H_vals'].add(h)

print("  (t₃, d₅) → count, H values, I(Ω,2)=1+2(t₃+d₅)")
for k in sorted(t3_d5_dist.keys()):
    t3, d5 = k
    info = t3_d5_dist[k]
    predicted = 1 + 2*(t3 + d5)
    h_vals = sorted(info['H_vals'])
    match = "✓" if all(h == predicted for h in h_vals) else "✗"
    print(f"    ({t3:2d}, {d5:2d}): {info['count']:4d} tournaments, H={h_vals}, predicted={predicted} {match}")

print("\n" + "="*70)
print("DONE")
print("="*70)
