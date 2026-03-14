#!/usr/bin/env python3
"""
CORRECT OCF verification with proper odd-cycle counting.
opus-2026-03-14-S71f

Key insight: Ω(T) vertices = ALL directed odd cycles.
A directed k-cycle has k cyclic rotations that represent the same cycle.
So count by fixing start vertex = min(cycle vertices).
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

def count_directed_odd_cycles(A, n):
    """Count ALL directed odd cycles, properly deduplicated.
    A directed k-cycle visits k vertices in a specific order (up to cyclic rotation).
    We canonicalize by starting at the minimum vertex."""
    cycles = set()
    
    for k in range(3, n+1, 2):  # odd lengths: 3, 5, 7, ...
        for verts in combinations(range(n), k):
            # Fix start = min vertex = verts[0]
            # Try all orderings of the remaining k-1 vertices
            for p in permutations(verts[1:]):
                order = [verts[0]] + list(p)
                # Check if this is a directed cycle
                is_cycle = True
                for idx in range(k):
                    if A[order[idx]][order[(idx+1) % k]] != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.add(tuple(order))
    
    return cycles

def compute_ocf(A, n):
    """Compute I(Ω(T), 2) properly."""
    cycles = count_directed_odd_cycles(A, n)
    cycle_list = list(cycles)
    num_cycles = len(cycle_list)
    
    # Build vertex sets for each cycle
    vsets = [frozenset(c) for c in cycle_list]
    
    # Build adjacency (share a vertex)
    adj = [[False]*num_cycles for _ in range(num_cycles)]
    for i in range(num_cycles):
        for j in range(i+1, num_cycles):
            if vsets[i] & vsets[j]:  # share at least one vertex
                adj[i][j] = adj[j][i] = True
    
    # Count independent sets by brute force (small enough at n=5)
    alpha = [0] * (num_cycles + 1)
    alpha[0] = 1
    
    # Enumerate all subsets (2^num_cycles can be large, but at n=5 should be OK)
    if num_cycles <= 20:
        for mask in range(1, 1 << num_cycles):
            bits_set = []
            m = mask
            while m:
                b = m & (-m)
                bits_set.append(b.bit_length() - 1)
                m ^= b
            
            # Check independence
            is_indep = True
            for i_idx in range(len(bits_set)):
                for j_idx in range(i_idx+1, len(bits_set)):
                    if adj[bits_set[i_idx]][bits_set[j_idx]]:
                        is_indep = False
                        break
                if not is_indep:
                    break
            
            if is_indep:
                alpha[len(bits_set)] += 1
    else:
        # Just compute alpha_1
        alpha[1] = num_cycles
    
    # I(Ω, 2)
    result = sum(alpha[k] * (2**k) for k in range(len(alpha)))
    return result, alpha, num_cycles

print("="*70)
print("CORRECT OCF VERIFICATION")
print("="*70)

n = 5
total_edges = n*(n-1)//2
num_t = 2**total_edges

print(f"\nn={n}: {num_t} tournaments")
print("\nVerifying H = I(Ω,2) for all tournaments...")

all_ok = True
cycle_stats = {}
for bits in range(num_t):
    A = make_tournament(bits, n)
    h = count_hp(A, n)
    
    cycles = count_directed_odd_cycles(A, n)
    t3 = sum(1 for c in cycles if len(c) == 3)
    d5 = sum(1 for c in cycles if len(c) == 5)
    
    I_val, alpha, nc = compute_ocf(A, n)
    
    key = (t3, d5)
    if key not in cycle_stats:
        cycle_stats[key] = {'count': 0, 'H': set(), 'I': set(), 'alpha': None}
    cycle_stats[key]['count'] += 1
    cycle_stats[key]['H'].add(h)
    cycle_stats[key]['I'].add(I_val)
    cycle_stats[key]['alpha'] = alpha[:4]
    
    if h != I_val:
        if all_ok:
            print(f"\n  FIRST FAILURE: bits={bits}")
            print(f"    H={h}, I(Ω,2)={I_val}")
            print(f"    t₃={t3}, d₅={d5}, total_cycles={nc}")
            print(f"    α = {alpha[:5]}")
            # Show cycles
            for c in sorted(cycles, key=lambda x: (len(x), x)):
                print(f"    Cycle: {c}")
        all_ok = False

if all_ok:
    print("  ★ H = I(Ω,2) for ALL 1024 tournaments at n=5! ✓")
else:
    # Count failures
    failures = 0
    for bits in range(num_t):
        A = make_tournament(bits, n)
        h = count_hp(A, n)
        I_val = compute_ocf(A, n)[0]
        if h != I_val:
            failures += 1
    print(f"\n  Failures: {failures}/{num_t}")

print("\n--- Cycle type breakdown ---")
for k in sorted(cycle_stats.keys()):
    t3, d5 = k
    info = cycle_stats[k]
    print(f"  (t₃={t3}, d₅={d5}): {info['count']} tournaments, H={sorted(info['H'])}, I={sorted(info['I'])}, α={info['alpha']}")

# At n=5, check if Ω is always complete
print("\n--- Is Ω always complete at n=5? ---")
for bits in [341, 40, 682]:  # Regular tournaments
    A = make_tournament(bits, n)
    h = count_hp(A, n)
    cycles = count_directed_odd_cycles(A, n)
    cycle_list = list(cycles)
    vsets = [frozenset(c) for c in cycle_list]
    
    # Check pairwise adjacency
    non_adjacent = 0
    for i in range(len(vsets)):
        for j in range(i+1, len(vsets)):
            if not (vsets[i] & vsets[j]):
                non_adjacent += 1
    
    print(f"  bits={bits}: H={h}, {len(cycles)} cycles, non-adjacent pairs: {non_adjacent}")

print("\n" + "="*70)
print("DONE")
print("="*70)
