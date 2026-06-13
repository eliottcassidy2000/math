"""
ocf_correct_verification.py -- kind-pasteur-2026-03-13-S62

CRITICAL FIX: Omega(T) vertices are individual DIRECTED odd cycles,
not vertex SETS. A vertex set of size k can support multiple directed
Hamiltonian k-cycles, each a separate vertex in Omega.

Key insight: for vertex set S of odd size k in tournament T:
  c(S) = number of directed Hamiltonian cycles on T[S] = tr(A_S^k)/k
  Each of these c(S) directed cycles is a SEPARATE vertex in Omega.
  All c(S) cycles on S form a clique in Omega (they share all vertices).

GS formula: H = sum over odd-cycle covers of 2^{#cycles}
           = I(Omega, 2) where Omega has directed cycles as vertices.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_directed_ham_cycles_on_subset(A, n, verts):
    """Count directed Hamiltonian cycles on sub-tournament induced by verts."""
    k = len(verts)
    if k < 3:
        return 0
    # Build sub-adjacency matrix
    sub_A = np.zeros((k, k), dtype=int)
    for i in range(k):
        for j in range(k):
            sub_A[i][j] = A[verts[i]][verts[j]]
    # tr(A^k) / k
    Ak = np.linalg.matrix_power(sub_A, k)
    return int(np.trace(Ak)) // k

def build_omega_correct(A, n):
    """Build Omega(T) with vertices = directed odd cycles.

    For each vertex set S of odd size k >= 3, add c(S) vertices
    to Omega (one per directed Hamiltonian cycle on T[S]).
    All vertices from same S form a clique.
    Vertices from different S1, S2 are adjacent iff S1 ∩ S2 != empty.

    Returns: list of (vertex_set, multiplicity) and adjacency matrix.
    """
    # Collect vertex sets with their cycle counts
    cycle_groups = []  # (frozenset, count)

    for size in range(3, n+1, 2):  # odd sizes
        for combo in combinations(range(n), size):
            verts = list(combo)
            c = count_directed_ham_cycles_on_subset(A, n, verts)
            if c > 0:
                cycle_groups.append((frozenset(combo), c))

    # Build expanded vertex list: each group contributes 'count' vertices
    expanded_vertices = []  # index -> (set_index, copy_index)
    set_index_map = {}
    for gi, (vset, count) in enumerate(cycle_groups):
        set_index_map[gi] = len(expanded_vertices)
        for copy in range(count):
            expanded_vertices.append(gi)

    nc = len(expanded_vertices)

    # Build adjacency: two vertices adjacent iff their vertex sets share a vertex
    omega_adj = np.zeros((nc, nc), dtype=int)
    for a in range(nc):
        for b in range(a+1, nc):
            gi_a = expanded_vertices[a]
            gi_b = expanded_vertices[b]
            if cycle_groups[gi_a][0] & cycle_groups[gi_b][0]:
                omega_adj[a][b] = omega_adj[b][a] = 1

    return cycle_groups, expanded_vertices, omega_adj, nc

def independence_poly(adj, nc):
    """Compute I(G, x) coefficients."""
    if nc > 25:
        return None  # Too large
    coeffs = [0] * (nc + 1)
    for mask in range(1 << nc):
        verts = [i for i in range(nc) if mask & (1 << i)]
        independent = True
        for i in range(len(verts)):
            for j in range(i+1, len(verts)):
                if adj[verts[i]][verts[j]]:
                    independent = False
                    break
            if not independent:
                break
        if independent:
            coeffs[len(verts)] += 1
    return coeffs

def eval_poly(coeffs, x):
    return sum(c * x**k for k, c in enumerate(coeffs))

# ============================================================
# TEST 1: n=5 exhaustive
# ============================================================
print("=" * 70)
print("TEST 1: OCF with CORRECT Omega at n=5 (exhaustive)")
print("=" * 70)

n = 5
total_bits = n*(n-1)//2

mismatches = 0
cycle_group_stats = defaultdict(list)

for bits in range(2**total_bits):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    cycle_groups, expanded, omega_adj, nc = build_omega_correct(A, n)

    ip = independence_poly(omega_adj, nc)
    if ip is not None:
        I_at_2 = eval_poly(ip, 2)
    else:
        I_at_2 = -1

    if I_at_2 != H:
        mismatches += 1
        if mismatches <= 5:
            group_info = [(len(g[0]), g[1]) for g in cycle_groups]
            print(f"  MISMATCH: bits={bits}, H={H}, I={I_at_2}, "
                  f"groups={group_info}, nc={nc}")

    # Stats
    for vset, count in cycle_groups:
        cycle_group_stats[len(vset)].append(count)

print(f"\nResult: {mismatches}/{2**total_bits} mismatches")
print(f"=> OCF {'VERIFIED' if mismatches == 0 else 'FAILED'} at n=5")

print(f"\nCycle group statistics (vertex set size -> multiplicity):")
for size in sorted(cycle_group_stats.keys()):
    vals = cycle_group_stats[size]
    print(f"  size {size}: count distribution = {Counter(vals)}")

# ============================================================
# TEST 2: n=6 sampled
# ============================================================
print("\n" + "=" * 70)
print("TEST 2: OCF with CORRECT Omega at n=6 (sampled)")
print("=" * 70)

n = 6
total_bits = n*(n-1)//2
np.random.seed(42)

mismatches = 0
n_samples = 500

for trial in range(n_samples):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    cycle_groups, expanded, omega_adj, nc = build_omega_correct(A, n)

    if nc > 25:
        # Too many cycles; skip
        continue

    ip = independence_poly(omega_adj, nc)
    I_at_2 = eval_poly(ip, 2)

    if I_at_2 != H:
        mismatches += 1
        if mismatches <= 5:
            group_info = [(len(g[0]), g[1]) for g in cycle_groups]
            print(f"  MISMATCH: bits={bits}, H={H}, I={I_at_2}, "
                  f"groups={group_info}, nc={nc}")

print(f"\nResult: {mismatches} mismatches in {n_samples} samples")

# ============================================================
# TEST 3: Understand 5-cycle multiplicity at n=5
# ============================================================
print("\n" + "=" * 70)
print("TEST 3: How many directed 5-cycles per vertex set at n=5?")
print("=" * 70)

n = 5
total_bits = n*(n-1)//2

c5_dist = Counter()

for bits in range(2**total_bits):
    A = bits_to_adj(bits, n)
    c5 = count_directed_ham_cycles_on_subset(A, n, list(range(n)))
    c5_dist[c5] += 1

print(f"  c5 (directed 5-cycles on all 5 vertices) distribution:")
for c, count in sorted(c5_dist.items()):
    print(f"    c5={c}: {count} tournaments")

# For the Paley tournament bits=40 case
print(f"\n  Specific mismatch case (bits=40):")
A = bits_to_adj(40, n)
H = count_ham_paths(A, n)
c5 = count_directed_ham_cycles_on_subset(A, n, list(range(n)))
c3_sets = []
for i in range(n):
    for j in range(i+1, n):
        for k in range(j+1, n):
            c = count_directed_ham_cycles_on_subset(A, n, [i,j,k])
            if c > 0:
                c3_sets.append(([i,j,k], c))

print(f"    H={H}, c5={c5}")
print(f"    3-cycles: {c3_sets}")
print(f"    Total Omega vertices = {sum(c for _, c in c3_sets)} + {c5} = {sum(c for _, c in c3_sets) + c5}")
print(f"    Expected I(Omega, 2) with multiplicity:")
n_c3 = sum(c for _, c in c3_sets)
print(f"      alpha_1 = {n_c3} (3-cycle vertices) + {c5} (5-cycle vertices) = {n_c3 + c5}")

# ============================================================
# TEST 4: Direct GS formula computation
# ============================================================
print("\n" + "=" * 70)
print("TEST 4: Direct GS formula (sum over odd-cycle partitions)")
print("=" * 70)

def gs_formula(A, n):
    """Compute H via the Grinberg-Stanley formula directly.

    H = sum over partitions of V into odd-sized parts where each part
    supports a directed Hamiltonian cycle, weighted by 2^{#parts} * prod(c_k_i).
    """
    total = 0

    def partition_rec(remaining, parts):
        nonlocal total
        if not remaining:
            # Have a complete partition
            # Weight = 2^|parts| * product of cycle counts per part
            weight = 1
            for part in parts:
                c = count_directed_ham_cycles_on_subset(A, n, list(part))
                weight *= c
            if weight > 0:
                total += (2 ** len(parts)) * weight
            return

        # Choose the next part containing the smallest remaining element
        first = min(remaining)
        remaining_list = sorted(remaining)

        # Try all odd-sized subsets containing 'first'
        for size in range(3, len(remaining) + 1, 2):
            others = [v for v in remaining_list if v != first]
            for combo in combinations(others, size - 1):
                part = frozenset([first] + list(combo))
                new_remaining = remaining - part
                partition_rec(new_remaining, parts + [part])

        # Also try: first is a 1-cycle (size 1)
        # But 1-cycles are trivial (not directed cycles in tournament)
        # Actually, does the GS formula count 1-cycles?
        # A fixed point in a permutation is a 1-cycle (odd!).
        # In a tournament, there's no self-loop, so 1-cycles might count as trivially valid?
        # Let me include and see
        part = frozenset([first])
        new_remaining = remaining - part
        partition_rec(new_remaining, parts + [part])

    partition_rec(frozenset(range(n)), [])
    return total

def gs_formula_no_fixed(A, n):
    """GS formula WITHOUT 1-cycles (fixed points)."""
    total = 0

    def partition_rec(remaining, parts):
        nonlocal total
        if not remaining:
            weight = 1
            for part in parts:
                c = count_directed_ham_cycles_on_subset(A, n, list(part))
                weight *= c
            if weight > 0:
                total += (2 ** len(parts)) * weight
            return

        first = min(remaining)
        remaining_list = sorted(remaining)

        for size in range(3, len(remaining) + 1, 2):
            others = [v for v in remaining_list if v != first]
            for combo in combinations(others, size - 1):
                part = frozenset([first] + list(combo))
                new_remaining = remaining - part
                partition_rec(new_remaining, parts + [part])

    partition_rec(frozenset(range(n)), [])
    return total

n = 5
total_bits = n*(n-1)//2

print(f"\n  n={n}, testing GS formula variants:")
print(f"  {'bits':>5} {'H':>5} {'GS+fix':>8} {'GS-fix':>8}")

for bits in [0, 1, 5, 10, 40, 100, 500, 1023]:
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    gs_with_fixed = gs_formula(A, n)
    gs_no_fixed = gs_formula_no_fixed(A, n)

    print(f"  {bits:>5} {H:>5} {gs_with_fixed:>8} {gs_no_fixed:>8}")

# ============================================================
# TEST 5: Count GS formula with 1-cycles
# ============================================================
print("\n" + "=" * 70)
print("TEST 5: Does GS formula match H? (exhaustive n=5)")
print("=" * 70)

n = 5

match_with = 0
match_no = 0

for bits in range(2**total_bits):
    A = bits_to_adj(bits, n)
    H = count_ham_paths(A, n)

    gs_w = gs_formula(A, n)
    gs_n = gs_formula_no_fixed(A, n)

    if gs_w == H:
        match_with += 1
    if gs_n == H:
        match_no += 1

print(f"  GS with 1-cycles = H: {match_with}/{2**total_bits}")
print(f"  GS without 1-cycles = H: {match_no}/{2**total_bits}")

print("\n\nDone.")
